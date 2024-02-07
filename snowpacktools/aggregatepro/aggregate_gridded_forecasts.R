#' Call this script from the command line to aggregate gridded forecasts
#' e.g., Rscript ./aggregate_gridded_forecasts.R --config config.ini --mp_csv input/aggregates_mp0.csv [--worker_int 1]
#' AWSOME: it expects to be called from `forecasts` directory (or as specified in forecast_runtime_domain.ini: ['Paths']['_cwd_'])

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

if (length(worker) == 0) worker <- "1"

cat(paste0("[i] (", worker, ") RScript: Working directory set to ", getwd(), "\n"))

library(sarp.snowprofile)
library(sarp.snowprofile.alignment)
library(sarp.snowprofile.pyface)
library(configr)
library(stringr)
library(data.table)

## Read config and csv files
config <- configr::read.config(file = configfile)
mp_df <- fread(mp_csv, sep = ",", data.table = FALSE)
vstations <- fread(config$Paths$`_aggregates_vstations_csv_file`, sep = ",", data.table = FALSE)
source(config$Paths$`_aggregates_plotters_path`)

## Parse DTW weights
config$DTW_weights$dims = c("gtype", "hardness", "ddate")
config$DTW_weights$weights = c(as.double(config$DTW_weights$GTYPE),
                               as.double(config$DTW_weights$HARDNESS),
                               as.double(config$DTW_weights$DDATE))
config$DTW_weights$dims <- config$DTW_weights$dims[config$DTW_weights$weights > 0]
config$DTW_weights$weights <- config$DTW_weights$weights[config$DTW_weights$weights > 0]

## Get filenames of .pro files in _aggregates_snp_pro_dir
file_names <- list.files(path = config$Paths$`_aggregates_snp_pro_dir`, pattern = "\\.pro$", full.names = TRUE)
smet_names <- gsub("\\.pro", ".smet", file_names)
# Extract the id between "VIR" and ".pro"
file_ids <- str_extract(basename(file_names), "(?<=VIR)[^.]+(?=\\.pro)")

## ---Iterate over combinations of region, band, aspect-----------------------
t0 <- Sys.time()
for (i in seq_len(nrow(mp_df))) {
  iterstatus <- tryCatch({

    ## ---Parse files and query dates-----------------------------------------
    dfrow <- mp_df[i, ]
    vstation_ids <- vstations$vstation[vstations$region_id == mp_df[i, "region_id"] & vstations$band == mp_df[i, "band"] & vstations$aspect == mp_df[i, "aspect"]]
    k_files <- which(file_ids %in% vstation_ids)
    file_names_sub <- file_names[k_files]
    smet_names_sub <- smet_names[k_files]
    file_ids_sub <- file_ids[k_files]

    if (length(file_names_sub) >= 3) {
      ##  Extract timezone from SMET Header
      tz  <- sapply(smet_names_sub, function(fn) {
        readSmet(fn, HeaderOnly = TRUE)$tz
      })
      tz_unique <- unique(tz)
      if (length(tz_unique) > 1) stop("Your profiles refer to different time zones! Not supported.")
      tz_string <- ifelse(tz_unique >= 0, paste0("Etc/GMT-", tz_unique), paste0("Etc/GMT+", abs(tz_unique)))
      ##  Generate relevant datetime period to aggregate
      #   based on first file dates and by assuming all .pro files have the same dates
      dtopera <- as.Date(config$Forecast$DATE_OPERA)
      dailytime_parts <- strsplit(config$Aggregate$DAILY_TIME, ":")[[1]]
      hours <- as.numeric(dailytime_parts[1])
      minutes <- as.numeric(dailytime_parts[2])
      fdatetime_max <- max(scanProfileDates(file_names_sub[1], tz = tz_string))
      dtmax <- min(fdatetime_max, as.POSIXct(format(as.Date(config$Forecast$SEASON_END), paste0("%Y-%m-%d ", hours, ":", minutes)), tz = tz_string))
      dtopera <- as.POSIXct(format(dtopera, paste0("%Y-%m-%d ", hours, ":", minutes)), tz = tz_string)
      dtperiod <- seq(dtopera, dtmax, by = "day")

      
      ##  ---Read profiles----------------------------------------------------
      if (length(dtperiod) > 1) {
        profileset <- snowprofileSet(do.call("c", lapply(file_names_sub, snowprofilePro, ProfileDate = dtperiod, 
                                                         tz = config$Forecast$TZONE, suppressWarnings = TRUE)))
      } else if (length(dtperiod) == 1) {
        profileset <- snowprofileSet(lapply(file_names_sub, snowprofilePro, ProfileDate = dtperiod, 
                                            tz = config$Forecast$TZONE, suppressWarnings = TRUE))
      } else {
        cat(paste0("[w] (", worker, ") No profiles at the relevant dates/timestamps for ", mp_df[i, "region_id"], " ", 
                   mp_df[i, "band"], " ", mp_df[i, "aspect"], "\n"))
        quit(save = "no")
      }
      ## routine requires unique station names per station.
      # 1) take from .pro files at StationName
      # 2) if 1) not unique: use smet file
      sm <- summary(profileset)
      if (length(unique(sm$station_id)) == 1) {
        tryCatch({
          sm$station_id <- sapply(smet_names_sub, function(fn) {
            readSmet(fn, HeaderOnly = TRUE)$station_id
          })
        }, error = function(e) {
          if (config$Aggregate$DEBUG_MODE) cat(e$message, "\n")
          stop(paste0("[E] (", worker, ") This error likely occurs when some vstations miss time stamps",
                      " and the station_id is not available from .pro files but retrieved from .smet files. \n"))
        })
      }

      
      ## ---Preprocess profiles-----------------------------------------------
      profileset <- computeRTA(profileset)
      # profileset <- computePunstable(profileset)  # verify unit of ski pen!
      ## Create random ski_pen to test entire framework until ski_pen issue resolved
      warning("Random ski_pen used for testing purposes")
      cat(paste0("[W] (", worker, ") Random ski_pen used for testing purposes \n"))
      profileset <- computePunstable(profileset, ski_pen = rep(0.2, length(profileset)))
      profileset <- snowprofileSet(lapply(profileset, function(sp) {
        labelPWL(sp, pwl_gtype = c("SH", "DH", "FCxr", "FC"), threshold_gtype = c("FC", "FCxr"), threshold_RTA = 0.8)
      }))

      ## ---Do aggregation----------------------------------------------------
      ## Load existing aggregate profile from file if available
      rds_files <- list.files(path = config$Paths$`_aggregates_output_dir`, pattern = "\\.rds$", full.names = TRUE)
      rds_names <- basename(rds_files)
      k_rds <- which(grepl(paste(mp_df[i, "region_id"], mp_df[i, "band"], mp_df[i, "aspect"], sep = "_"), rds_names))
      init <- TRUE
      if (length(k_rds) == 1) {
        avg1 <- readRDS(rds_files[k_rds])
        if ((as.Date(dtopera) > min(avg1$meta$date))
             & (as.Date(dtopera) <= max(avg1$meta$date)+1)) {
          ## delete all past lead-time forecasts
          #  this will essentially re-compute the average profile from DATE_OPERA to the latest available date in .pro files
          init <- FALSE
          avg_avgs_dayBefore <- avg1$avgs[[avg1$meta$date == as.Date(dtopera)-1]]
          avg2 <- averageSPalongSeason(profileset, AvgDayBefore = avg_avgs_dayBefore, sm = sm, 
                                       progressbar = config$Aggregate$DEBUG_MODE, verbose = config$Aggregate$DEBUG_MODE,
                                       dims = config$DTW_weights$dims, weights = config$DTW_weights$weights)
          avg <- concat_avgSP_timeseries(avg1, avg2)
        } else {
          init <- TRUE
          cat(paste(
            "[w] (", worker, ") Looks like the average profile on file is outdated/erroneous.",
            "I re-initialize the average profile and overwrite the file. \n"
          ))
        }
      }
      if (init) {
        avg <- averageSPalongSeason(profileset, sm = sm,
                                    progressbar = config$Aggregate$DEBUG_MODE, verbose = config$Aggregate$DEBUG_MODE,
                                    dims = config$DTW_weights$dims, weights = config$DTW_weights$weights)
      }

      if (sum(avg$meta$reinitialized) > 0.2*nrow(avg$meta)) {
        cat(paste("[w] (", worker, ") More than 20% of average profiles were re-initialized for", 
                    mp_df[i, "region_id"], mp_df[i, "band"], mp_df[i, "aspect"], 
                    "--Consider investigating! \n"))
      }
      ## Save to file
      if (config$Aggregate$SAVEAS_rds) {
        saveRDS(avg, paste0(
          config$Paths$`_aggregates_output_dir`,
          "/", mp_df[i, "region_id"], "_", mp_df[i, "band"], "_", mp_df[i, "aspect"], ".rds"
        ))
      }
      


      ## ---Plot figures-------------------------------------------------------

      ## ---hand hardness profile-----------------------------------------------
      ## single hand hardness profile with instability distributions
      if (config$Aggregate$PLOT_HandHardness) {
        for (pdate in avg$meta$date[avg$meta$date > as.Date(config$Forecast$DATE_OPERA) &
                                    avg$meta$date < as.Date(config$Forecast$DATE_OPERA) + config$Aggregate$PLOT_Leadtime_days_HandHardness]) {
          leadtime <- as.numeric(pdate - as.Date(config$Forecast$DATE_OPERA)) # (days)
          fname <- paste0(
            config$Paths$`_aggregates_figures_dir`, "/",
            "hhp_", mp_df[i, "region_id"], "_", mp_df[i, "band"], "_", mp_df[i, "aspect"], "_",
            format(as.Date(config$Forecast$DATE_OPERA), "%y%m%d"), "+", as.integer(leadtime), "d", ".png"
          )
          png(filename = fname, width = 800, height = 700)
          par(cex.lab = 1.65, cex.axis = 1.8, bg = "white")
          plotTradAvgProfile(avg, pdate)
          dev.off()
        }
      }
      
      ## ---avg timeseries-----------------------------------------------------
      if (config$Aggregate$PLOT_TSplain) {
        pdate <- max(avg$meta$date)
        leadtime <- as.numeric(pdate - as.Date(config$Forecast$DATE_OPERA)) # (days)
        fname <- paste0(
          config$Paths$`_aggregates_figures_dir`, "/",
          "tsplain_", mp_df[i, "region_id"], "_", mp_df[i, "band"], "_", mp_df[i, "aspect"], "_",
          format(as.Date(config$Forecast$DATE_OPERA), "%y%m%d"), "+", as.integer(leadtime), "d", ".png"
        )
        png(filename = fname, width = 1200, height = 600)
        par(cex.lab = 1.4, cex.axis = 1.4, bg = "white")
        plotTSplainAvgProfile(avg)
        dev.off()
      }

      ## ---avg timeseries w/ stability overplot-------------------------------
      if (config$Aggregate$PLOT_TSstability) {
        pdate <- max(avg$meta$date)
        leadtime <- as.numeric(pdate - as.Date(config$Forecast$DATE_OPERA)) # (days)
        fname <- paste0(
          config$Paths$`_aggregates_figures_dir`, "/",
          "tsstab_", mp_df[i, "region_id"], "_", mp_df[i, "band"], "_", mp_df[i, "aspect"], "_",
          format(as.Date(config$Forecast$DATE_OPERA), "%y%m%d"), "+", as.integer(leadtime), "d", ".png"
        )
        png(filename = fname, width = 1200, height = 600)
        par(cex.lab = 1.4, cex.axis = 1.4, bg = "white")
        plotTSstabilityAvgProfile(avg)
        dev.off()
      }

    } else {
      cat(paste0("[i] (", worker, ") Not aggregating b/c less than three profiles for ", mp_df[i, "region_id"], " ",
                 mp_df[i, "band"], " ", mp_df[i, "aspect"], " \n"))
    }
  }, error = function(e) {
    if (config$Aggregate$DEBUG_MODE) cat(paste0(e$message, "\n"))
    cat(paste0("[E] (", worker, ") Error while aggregating profiles for ", mp_df[i, "region_id"], " ",
               mp_df[i, "band"], " ", mp_df[i, "aspect"], " \n"))
  })
  if (inherits(iterstatus, "error")) next
}  # END for loop

if (config$Aggregate$DEBUG_MODE) cat(paste0("[i] (", worker, ") took ", format(round(Sys.time() - a, 1)), "\n"))
