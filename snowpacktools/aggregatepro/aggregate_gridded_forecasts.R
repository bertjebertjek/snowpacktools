#' Worker R script as part of the python module snowpro.aggregatepro.gridded
#' 
#' Call this script from the command line to aggregate gridded forecasts
#' e.g., Rscript ./aggregate_gridded_forecasts.R --config config.ini --mp_csv input/aggregates_groupings/domain-0.csv [--worker_int 1]
#' @AWSOME: it expects to be called from the `~snowpack/gridded-chain/` directory
#' 
#' <fherla>

## ---Setup ---------------------------------------------------------------
## Parse inputs
tryCatch({
  args <- commandArgs()
  configfile <- as.character(args[which(args == "--config") + 1])
  mp_csv  <- as.character(args[which(args == "--mp_csv") + 1])
  worker <- as.character(args[which(args == "--worker_int") + 1])
}, error = function(e) {
  stop("[E] Error when parsing inputs: ", e$message, call. = FALSE)
})

if (length(worker) == 0) worker <- "0"

cat(paste0("[i] (", worker, ") RScript: Working directory set to ", getwd(), "\n"))
# cat(paste0(
# 'For debugging, paste the following into an R session:

# mock_args <- c("--config", "', configfile, '",
#                "--mp_csv", "', mp_csv, '")
# assign("commandArgs", function(...) mock_args, envir = .GlobalEnv)
# source("/opt/awsome/code/snow-cover/postprocessing/snowpacktools/snowpacktools/aggregatepro/aggregate_gridded_forecasts.R", echo = TRUE)
# '
# ))

library(sarp.snowprofile)
library(sarp.snowprofile.alignment)
library(sarp.snowprofile.pyface)
library(configr)
library(stringr)
library(data.table)

## Read config and csv files
config <- configr::read.config(file = configfile)
mp_df <- fread(mp_csv, sep = ",", data.table = FALSE)
vstations <- fread(config$Paths$`_aggregates_vstations_csv_file_processed`, sep = ",", data.table = FALSE)
source(config$Paths$`_aggregates_plotters_path`)

## Parse DTW hyperparameter settings:
config$Advanced$dims = c("gtype", "hardness", "ddate")
config$Advanced$weights = c(as.double(config$Advanced$weights_gtype),
                               as.double(config$Advanced$weights_hardness),
                               as.double(config$Advanced$weights_ddate))
config$Advanced$dims <- config$Advanced$dims[config$Advanced$weights > 0]
config$Advanced$weights <- config$Advanced$weights[config$Advanced$weights > 0]
if (config$Aggregate$saveas_rds) {
  config$Advanced$keepprofiles <- TRUE
} else if (config$Aggregate$saveas_rds == "MINIMAL") {
  config$Advanced$keepprofiles <- FALSE
}

## Initialize progressbar
if (as.logical(config$Aggregate$progressbar) && worker == "0" && requireNamespace("progress", quietly = TRUE)) {
  progressbar = TRUE
  pb <- progress::progress_bar$new(
    format = paste0("(worker 0) [:bar] :percent in :elapsed | eta: :eta"),
    total = nrow(mp_df), clear = FALSE, width= 60)
  cat(paste0("[i] (", worker, ") Showing progressbar only for first worker.\n", 
             "        Information will be more accurate the more different\n",
             "        region--band-asspect combinations this worker will have to carry out.\n",
             "        It might take a while for the progressbar to appear."))
} else {
  progressbar = FALSE
}

## Get filenames of .pro files in _aggregates_snp_pro_dir
file_names <- list.files(path = config$Paths$`_aggregates_snp_pro_dir`, pattern = "\\.pro$", full.names = TRUE)
smet_names <- gsub("\\.pro", ".smet", file_names)
# Extract the id between "VIR" and ".pro"
file_ids <- str_extract(basename(file_names), "(?<=VIR)[^.]+(?=\\.pro)")
## Get filenames of .rds files in _aggregates_output_dir
rds_files <- list.files(path = config$Paths$`_aggregates_output_dir`, pattern = "\\.rds$", full.names = TRUE)
rds_names <- basename(rds_files)

## ---Iterate over combinations of region, band, aspect-----------------------
t0 <- Sys.time()
errorcode <- 0
for (i in seq_len(nrow(mp_df))) {
  if (progressbar) pb$tick()
  iterstatus <- tryCatch({

    ## ---Parse files and query dates-----------------------------------------
    dfrow <- mp_df[i, ]
    vstation_ids <- vstations$vstation[vstations$region_id == mp_df[i, "region_id"] & vstations$band == mp_df[i, "band"] & vstations$aspect == mp_df[i, "aspect"]]
    k_files <- which(file_ids %in% vstation_ids)
    file_names_sub <- file_names[k_files]
    smet_names_sub <- smet_names[k_files]
    file_ids_sub <- file_ids[k_files]
    
    if (length(file_names_sub) >= 2) {
      ##  Extract timezone from SMET Header
      tz  <- sapply(smet_names_sub, function(fn) {
        readSmet(fn, HeaderOnly = TRUE)$tz
      })
      tz_unique <- unique(tz)
      if (length(tz_unique) > 1) stop("Your profiles refer to different time zones! Not supported.")
      tz_string <- ifelse(tz_unique >= 0, paste0("Etc/GMT-", tz_unique), paste0("Etc/GMT+", abs(tz_unique)))
      ##  Generate relevant datetime period to aggregate
      #   based on first file dates and by assuming all .pro files have the same dates
      dtopera <- as.Date(config$General$date_opera)
      dailytime_parts <- strsplit(config$Aggregate$daily_time, ":")[[1]]
      hours <- as.numeric(dailytime_parts[1])
      minutes <- as.numeric(dailytime_parts[2])
      fdatetime <- scanProfileDates(file_names_sub[1], tz = tz_string)
      dt_season_end <- as.POSIXct(format(as.Date(config$General$season_end), paste0("%Y-%m-%d ", hours, ":", minutes)), tz = tz_string)
      dtmax <- min(max(fdatetime), dt_season_end)
      if (config$Aggregate$initialize_from == 'date_opera') {
        dtopera <- as.POSIXct(format(dtopera, paste0("%Y-%m-%d ", hours, ":", minutes)), tz = tz_string)
        if (dtopera >= dt_season_end) {
          stop(paste0("config$initialize_from is set to date_opera, but that's past the season_end.",
         " Either set to season_start or set date_opera to a prior date."))
        } else if (dtopera > dtmax)  {
          print(paste0("[i]  No profiles available yet for date_opera at daily_time. Trying to initialize from",
          " previous day since config$initialize_from is set to date_opera."))
          dtopera <- dtopera - as.difftime(1, units = "days")
          if (dtopera < min(fdatetime)) {
            stop("No profiles available for date_opera@daily_time minus 1 day.")
          }
        }
        dtperiod <- seq(dtopera, dtmax, by = "day")
      } else if (config$Aggregate$initialize_from == 'season_start') {
        dt_season_start = as.POSIXct(format(as.Date(config$General$season_start), paste0("%Y-%m-%d ", hours, ":", minutes)), tz = tz_string)
        if (min(fdatetime) > dt_season_start) {
          dtmin <- dt_season_start + 86400  # adding one day in seconds
        } else {
          dtmin <- dt_season_start
        }
        dtperiod <- seq(dtmin, dtmax, by = "day")
      } else {
        stop("config$Aggregate$initialize_from must be one of ['season_start', 'date_opera']")
      }

      
      ##  ---Read profiles----------------------------------------------------
      if (length(dtperiod) > 1) {
        profileset <- snowprofileSet(do.call("c", lapply(file_names_sub, snowprofilePro, ProfileDate = dtperiod, 
                                                         tz = tz_string, suppressWarnings = TRUE)))
      } else if (length(dtperiod) == 1) {
        profileset <- snowprofileSet(lapply(file_names_sub, snowprofilePro, ProfileDate = dtperiod, 
                                            tz = tz_string, suppressWarnings = TRUE))
      } else {
        cat(paste0("[w] (", worker, ") No profiles at the relevant dates/timestamps for ", mp_df[i, "region_id"], " ", 
                   mp_df[i, "band"], " ", mp_df[i, "aspect"], "\n"))
        quit(save = "no", status = 0)
      }
      wxlist = lapply(smet_names_sub, readSmet)
      ## routine requires unique station names per station:
      sm <- summary(profileset)
      # check when elev changes or when date jumps back into past or when station_id changes
      sm$change <- c(0, diff(sm$elev) != 0 | diff(sm$date) < 0 | !str_detect(sm$station[-1], sm$station[-length(sm$station)]))
      # hack a station_id to satisfy checks in aggregating function
      sm$profile_number <- 1 + cumsum(sm$change)
      sm$station_id <- sm$profile_number
      ## The hack should actually do just fine. Anyway, include check and messaging in case of weird results:
      if (length(unique(sm$station_id)) != length(file_names_sub)) {
        cat(paste0(
          "[w] (", worker, ") No StationName with unique station_id present in .pro file(s)!",
          " This might lead to unexpected errors/bugs. Update your .pro files. \n"
        ))
        cat(paste0(
          "[w] (", worker, ") ", length(file_names_sub), " different profiles available,",
          " but I had to create ", length(unique(sm$station_id)), " different profile_ids \n"
        ))
        ## check whether higher frequency than daily smpling:
        # tmp <- lapply(sm$station_id, function(sid) {
        #   which(duplicated(sm$date[sm$station_id == sid]))
        # })
        # any(unlist(tmp))
      }
      ## write ski_pen into profile meta summary:
      sm$ski_pen = NA
      for (smi in seq(nrow(sm))) {
        sm$ski_pen[smi] = wxlist[[sm$profile_number[smi]]]$data$ski_pen[wxlist[[sm$profile_number[smi]]]$data$timestamp %in% sm$datetime[smi]]
      }

      ## ---Preprocess profiles-----------------------------------------------
      profileset <- computeRTA(profileset)
      profileset <- computePunstable(profileset, ski_pen = sm$ski_pen)
      profileset <- snowprofileSet(lapply(profileset, function(sp) {
        labelPWL(sp, pwl_gtype = c("SH", "DH", "FCxr", "FC"), threshold_gtype = c("FC", "FCxr"), threshold_RTA = 0.8)
      }))

      ## ---Do aggregation----------------------------------------------------
      ## Load existing aggregate profile from file if available
      k_rds <- which(grepl(paste(mp_df[i, "region_id"], mp_df[i, "band"], mp_df[i, "aspect"], sep = "_"), rds_names))
      init <- TRUE
      if (length(k_rds) == 1 & config$Aggregate$initialize_from != 'season_start') {
        avg1 <- readRDS(rds_files[k_rds])
        if ((as.Date(dtopera) > min(avg1$meta$date))
             & (as.Date(dtopera) <= max(avg1$meta$date)+1)) {
          ## delete all past lead-time forecasts
          #  this will essentially re-compute the average profile from DATE_OPERA to the latest available date in .pro files
          init <- FALSE
          tryCatch({
            avg_avgs_dayBefore <- avg1$avgs[[which(avg1$meta$date == (as.Date(dtopera)-1))]]
          }, error = function(e) {
            stop(paste0("Between your average profile on file and your DATE_OPERA seem to lie some days without data. ",
                        "Please set back your DATE_OPERA or delete the relevant .rds file with the average profile ",
                        "if you want to start computations from scratch."))
          })
          
          avg2 <- averageSPalongSeason(profileset, AvgDayBefore = avg_avgs_dayBefore, sm = sm, 
                                       progressbar = FALSE, verbose = FALSE,
                                       keep.profiles = config$Advanced$keepprofiles,
                                       dims = config$Advanced$dims, weights = config$Advanced$weights,
                                       simType = tolower(config$Advanced$sim_type))
          avg <- concat_avgSP_timeseries(avg1, avg2)
        } else {
          init <- TRUE
          cat(paste0(
            "[w] (", worker, ") Looks like the average profile on file is outdated/erroneous. ",
            "I re-initialize the average profile and overwrite the file. \n"
          ))
        }
      }
      if (init) {
        avg <- averageSPalongSeason(profileset, sm = sm,
                                    progressbar = FALSE, verbose = FALSE,
                                    keep.profiles = config$Advanced$keepprofiles,
                                    dims = config$Advanced$dims, weights = config$Advanced$weights,
                                    simType = tolower(config$Advanced$sim_type))
      }

      if (sum(avg$meta$reinitialized) > 1 & sum(avg$meta$reinitialized) > 0.2*nrow(avg$meta)) {
        cat(paste0("[w] (", worker, ") More than 20% of average profiles were re-initialized for ", 
                    mp_df[i, "region_id"], " ", mp_df[i, "band"], " ", mp_df[i, "aspect"], 
                    " --Consider investigating! \n"))
      }
      ## Save to file
      if (config$Aggregate$saveas_rds) {
        saveRDS(avg, paste0(
          config$Paths$`_aggregates_output_dir`,
          "/", mp_df[i, "region_id"], "_", mp_df[i, "band"], "_", mp_df[i, "aspect"], ".rds"
        ))
      }
      


      ## ---Plot figures-------------------------------------------------------

      ## ---hand hardness profile-----------------------------------------------
      ## single hand hardness profile with instability distributions
      if (config$Aggregate$plot_handhardness) {
        for (pdate in avg$meta$date[avg$meta$date >= as.Date(config$General$date_opera) &
                                    avg$meta$date < as.Date(config$General$date_opera) + as.double(config$Aggregate$plot_leadtime_days_handhardness)]) {
          leadtime <- as.numeric(as.Date(pdate) - as.Date(config$General$date_opera)) # (days)
          fname <- paste0(
            config$Paths$`_aggregates_figures_dir`, "/",
            "hhp_", mp_df[i, "region_id"], "_", mp_df[i, "band"], "_", mp_df[i, "aspect"], "_",
            format(as.Date(config$General$date_opera), "%y%m%d"), "+", as.integer(leadtime), "d", ".png"
          )
          png(filename = fname, width = 800, height = 700)
          par(cex.lab = 1.65, cex.axis = 1.8, bg = "white")
          plotTradAvgProfile(avg, pdate)
          dev.off()
        }
      }
      
      ## ---avg timeseries-----------------------------------------------------
      if (config$Aggregate$plot_tsplain) {
        pdate <- max(avg$meta$date)
        leadtime <- as.numeric(pdate - as.Date(config$General$date_opera)) # (days)
        fname <- paste0(
          config$Paths$`_aggregates_figures_dir`, "/",
          "tsplain_", mp_df[i, "region_id"], "_", mp_df[i, "band"], "_", mp_df[i, "aspect"], "_",
          format(as.Date(config$General$date_opera), "%y%m%d"), "+", as.integer(leadtime), "d", ".png"
        )
        png(filename = fname, width = 1200, height = 600)
        par(cex.lab = 1.4, cex.axis = 1.4, bg = "white")
        plotTSplainAvgProfile(avg)
        dev.off()
      }

      ## ---avg timeseries w/ stability overplot-------------------------------
      if (config$Aggregate$plot_tsstability) {
        pdate <- max(avg$meta$date)
        leadtime <- as.numeric(pdate - as.Date(config$General$date_opera)) # (days)
        fname <- paste0(
          config$Paths$`_aggregates_figures_dir`, "/",
          "tsstab_", mp_df[i, "region_id"], "_", mp_df[i, "band"], "_", mp_df[i, "aspect"], "_",
          format(as.Date(config$General$date_opera), "%y%m%d"), "+", as.integer(leadtime), "d", ".png"
        )
        png(filename = fname, width = 1200, height = 600)
        par(cex.lab = 1.4, cex.axis = 1.4, bg = "white")
        plotTSstabilityAvgProfile(avg)
        dev.off()
      }

    } else {
      cat(paste0("[i] (", worker, ") Not aggregating b/c less than two profiles for ", mp_df[i, "region_id"], " ",
                 mp_df[i, "band"], " ", mp_df[i, "aspect"], " \n"))
    }
  }, error = function(e) {
    cat(paste0(e$message, "\n"))
    cat(paste0("\n\n[E] (", worker, ") Error while aggregating profiles for ", mp_df[i, "region_id"], " ",
               mp_df[i, "band"], " ", mp_df[i, "aspect"], " \n\n"))
  })  # END tryCatch
  if (inherits(iterstatus, "error")) {
    errorcode <- 1
    next
  }
}  # END for loop

if (config$Aggregate$debug_mode) cat(paste0("[i] (", worker, ") took ", format(round(Sys.time() - t0, 1)), "\n"))

quit(save = "no", status = errorcode)