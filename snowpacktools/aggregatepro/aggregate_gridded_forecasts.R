#' Call this script from the command line to aggregate gridded forecasts
#' e.g., Rscript ./aggregate_gridded_forecasts.R --config input/forecast.ini --mp_csv input/aggregates_mp0.csv
#' it expects to be called from `forecasts` directory (or as specified in forecasts.ini: ['Paths']['_cwd_'])

## ---Setup ---------------------------------------------------------------
## Parse inputs
tryCatch({
  args <- commandArgs()
  configfile <- as.character(args[which(args == "--config") + 1])
  mp_csv  <- as.character(args[which(args == "--mp_csv") + 1])
}, error = function(e) {
  stop("[E] Error when parsing inputs: ", e$message, call. = FALSE)
})

library(sarp.snowprofile)
library(sarp.snowprofile.alignment)
library(sarp.snowprofile.pyface)
library(configr)
library(stringr)
library(data.table)

############
# for local testing only (where temp file not avilable):
# setwd("/home/flo/documents/code/awsome/models/SNOWPACK/forecasts")
# config <- configr::read.config(file = "input/forecast.ini")  # for local testing
# config$Forecast$SEASON_START <- as.Date("2022-11-01")
# config$Forecast$SEASON_END <- as.Date("2023-03-01")
# config$Forecast$DATE_OPERA <- as.Date("2023-02-28")
# config$Forecast$TZONE <- "UTC"
# config$Forecast$NTASKS <- 4

# config$Aggregate$DAILY_TIME <- "06:00"
# config$Aggregate$`_vstations_csv_file` <- "./input/vstations-subtirol10000.csv"

# config$Paths$`_snp_ouput_dir` <- "./output/snp-subtirol10000"
# config$Paths$`_aggregates_output_dir`  <- "./output/snp-aggregates-subtirol10000"
# config$Paths$`_aggregates_figures_dir` <- "./output/snp-aggregates-figs-subtirol10000"
# config$Paths$`_aggregates_plotters_path` <- "../profile-aggregation/plotters.R"

# mp_csv <- "./input/aggregates_mp0.csv"
# dtmax <- as.POSIXct(format(as.Date(dtmax), paste0("%Y-%m-%d ", 5, ":", minutes)), tz = config$Forecast$TZONE)

# configfile <- "input/forecast_localHF.ini"
# mp_csv <- "input/aggregates_mp0.csv"
# i  <- 2
#############

## Read config and csv files
config <- configr::read.config(file = configfile)
config_specific <- configr::read.config(file = config$Paths$`_aggregates_ini`)
if (all(!names(config_specific) %in% names(config))) {
  config <- c(config, config_specific)
} else {
  stop("[E] Aggregate ini file contains sections that are already defined in the forecast ini file.")
}
mp_df <- fread(mp_csv, sep = ",", data.table = FALSE)
vstations <- fread(config$Aggregate$`_vstations_csv_file`, sep = ",", data.table = FALSE)
source(config$Paths$`_aggregates_plotters_path`)

## Parse DTW weights
config$Dtw_weights$dims = c("gtype", "hardness", "ddate")
config$Dtw_weights$weights = c(as.double(config$Dtw_weights$GTYPE),
                               as.double(config$Dtw_weights$HARDNESS),
                               as.double(config$Dtw_weights$DDATE))
config$Dtw_weights$dims <- config$Dtw_weights$dims[config$Dtw_weights$weights > 0]
config$Dtw_weights$weights <- config$Dtw_weights$weights[config$Dtw_weights$weights > 0]

## Get filenames of .pro files in _snp_output_dir
file_names <- list.files(path = config$Paths$`_snp_ouput_dir`, pattern = "\\.pro$", full.names = TRUE)
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
    file_ids_sub <- file_ids[k_files]

    if (length(file_names_sub) >= 3) {
      ##  Generate relevant datetime period to aggregate
      #   based on first file dates and by assuming all .pro files have the same dates
      fdatetime_max <- max(scanProfileDates(file_names_sub[1]))
      dtmax <- min(fdatetime_max, as.POSIXct(config$Forecast$SEASON_END, tz = config$Forecast$TZONE))
      dtopera <- as.Date(config$Forecast$DATE_OPERA)
      dailytime_parts <- strsplit(config$Aggregate$DAILY_TIME, ":")[[1]]
      hours <- as.numeric(dailytime_parts[1])
      minutes <- as.numeric(dailytime_parts[2])
      dtopera <- as.POSIXct(format(dtopera, paste0("%Y-%m-%d ", hours, ":", minutes)), tz = config$Forecast$TZONE)
      dtperiod <- seq(dtopera, dtmax, by = "day")

      
      ##  ---Read profiles----------------------------------------------------
      if (length(dtperiod) > 1) {
        profileset <- snowprofileSet(do.call("c", lapply(file_names_sub, snowprofilePro, ProfileDate = dtperiod, suppressWarnings = TRUE)))
      } else if (length(dtperiod) == 1) {
        profileset <- snowprofileSet(lapply(file_names_sub, snowprofilePro, ProfileDate = dtperiod, suppressWarnings = TRUE))
      } else {
        print(paste("[w] No profiles at the relevant dates for", mp_df[i, "region_id"], mp_df[i, "band"], mp_df[i, "aspect"]))
        quit(save = "no")
      }
      ## This hack is necessary until the station_id is written to the
      #  .pro files as StationName:
      tryCatch({
        sm <- summary(profileset)
        # check when elev changes or when date jumps
        sm$change <- c(0, diff(sm$elev) != 0 | diff(sm$date) > 1)
        # hack a station_id to satisfy checks in aggregating function
        sm$station_id <- cumsum(sm$change)
      }, error = function(e){
        stop("Cannot resolve station_id because not all profiles have identical dates available. Need station_id as StationName in .pro files.")
      })
      
      ## ---Preprocess profiles-----------------------------------------------
      profileset <- computeRTA(profileset)
      # profileset <- computePunstable(profileset)  # verify unit of ski pen!
      ## Create random ski_pen to test entire framework until ski_pen issue resolved
      warning("Random ski_pen used for testing purposes")
      print("[W] Random ski_pen used for testing purposes")
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
          ## following steps have been ported to concat_avgSP_timeseries()
          # k_rm <- which(avg1$meta$date > as.Date(dtopera))
          # avg1$avgs[k_rm] <- NULL
          # try({
          #   avg1$sets[k_rm] <- NULL
          # })
          # avg1$meta <- avg1$meta[-k_rm, ]
          avg_avgs_dayBefore <- avg1$avgs[[avg1$meta$date == as.Date(dtopera)-1]]
          avg2 <- averageSPalongSeason(profileset, AvgDayBefore = avg_avgs_dayBefore, sm = sm, 
                                       progressbar = config$Aggregate$debug_mode, verbose = config$Aggregate$debug_mode,
                                       dims = config$Dtw_weights$dims, weights = config$Dtw_weights$weights)
          avg <- concat_avgSP_timeseries(avg1, avg2)
        } else {
          init <- TRUE
          print(paste(
            "[w] Looks like the average profile on file is outdated/erroneous.",
            "I re-initialize the average profile and overwrite the file."
          ))
        }
      }
      if (init) {
        avg <- averageSPalongSeason(profileset, sm = sm,
                                    progressbar = config$Aggregate$debug_mode, verbose = config$Aggregate$debug_mode,
                                    dims = config$Dtw_weights$dims, weights = config$Dtw_weights$weights)
      }

      ## Save to file
      if (sum(avg$meta$reinitialized) > 0.2*nrow(avg$meta)) {
        print(paste("[w] More than 20% of average profiles were re-initialized for", 
                    mp_df[i, "region_id"], mp_df[i, "band"], mp_df[i, "aspect"], 
                    "--Consider investigating!"))
      }
      saveRDS(avg, paste0(
        config$Paths$`_aggregates_output_dir`,
        "/", mp_df[i, "region_id"], "_", mp_df[i, "band"], "_", mp_df[i, "aspect"], ".rds"
      ))


      ## ---Plot figures-------------------------------------------------------
      pdate <- max(avg$meta$date)
      leadtime <- as.numeric(pdate - as.Date(config$Forecast$DATE_OPERA))  # (days)

      ## ---hand hardness profile-----------------------------------------------
      ## single hand hardness profile with instability distributions
      fname <- paste0(
        config$Paths$`_aggregates_figures_dir`, "/",
        "hhp_", mp_df[i, "region_id"], "_", mp_df[i, "band"], "_", mp_df[i, "aspect"], "_",
        format(as.Date(config$Forecast$DATE_OPERA), "%y%m%d"), "+", as.integer(leadtime), "d", ".png"
      )
      png(filename = fname, width = 800, height = 700)
      par(cex.lab = 1.65, cex.axis = 1.8, bg = "white")
      plotTradAvgProfile(avg, pdate)
      dev.off()

      ## ---avg timeseries-----------------------------------------------------
      fname <- paste0(
        config$Paths$`_aggregates_figures_dir`, "/",
        "tsplain_", mp_df[i, "region_id"], "_", mp_df[i, "band"], "_", mp_df[i, "aspect"], "_",
        format(as.Date(config$Forecast$DATE_OPERA), "%y%m%d"), "+", as.integer(leadtime), "d", ".png"
      )
      png(filename = fname, width = 1200, height = 600)
      par(cex.lab = 1.4, cex.axis = 1.4, bg = "white")
      plotTSplainAvgProfile(avg)
      dev.off()

    } else {
      print(paste("[i] Not aggregating b/c less than three profiles for", mp_df[i, "region_id"], mp_df[i, "band"], mp_df[i, "aspect"]))
    }
  }, error = function(e) {
    if (config$Aggregate$debug_mode) print(e$message)
    print(paste("[E] Error while aggregating profiles for", mp_df[i, "region_id"], mp_df[i, "band"], mp_df[i, "aspect"]))
  })
  if (inherits(iterstatus, "error")) next
}  # END for loop

if (config$Aggregate$debug_mode) print(round(Sys.time() - t0, 2))
