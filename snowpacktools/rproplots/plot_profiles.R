
tryCatch({
  args <- commandArgs()
  pro_file <- as.character(args[which(args == "--pro") + 1])
  smet_file  <- as.character(args[which(args == "--smet") + 1])
  outdir <- as.character(args[which(args == "--outdir") + 1])
  dt_nowcast_end <- as.character(args[which(args == "--dt_nowcast_end") + 1])
  dt_forecast_end <- as.character(args[which(args == "--dt_forecast_end") + 1])
  daily_time <- as.character(args[which(args == "--daily_time") + 1])
  worker <- as.character(args[which(args == "--worker") + 1])
}, error = function(e) {
  stop("[E] Error when parsing inputs: ", e$message, call. = FALSE)
})


library(sarp.snowprofile)
library(sarp.snowprofile.pyface)
library(stringr)

pro_name = str_remove(basename(pro_file), '.pro')

if (length(worker) == 0) worker <- "0"
if (length(daily_time) == 0) daily_time <- "23:00"
if (length(smet_file) == 0) smet_file <- paste0(file.path(dirname(pro_file), pro_name), ".smet")
if (length(dt_nowcast_end) == 0) {
    stop("No dt_nowcast_end given")
} else {
    dt_nowcast_end <- as.POSIXct(dt_nowcast_end, format = "%Y-%m-%dT%H:%M", tz = "UTC")
}

dt_start = dt_nowcast_end - as.difftime(10, units = "days")
pro_dates = scanProfileDates(pro_file, tz = "UTC")
if (length(dt_forecast_end) == 0) {
    dt_forecast_end <- max(pro_dates)
    if (dt_forecast_end < dt_nowcast_end) dt_forecast_end = dt_nowcast_end
} else {
    dt_forecast_end <- as.POSIXct(dt_forecast_end, format = "%Y-%m-%dT%H:%M", tz = "UTC")
}
pro_dates = pro_dates[pro_dates >= dt_start & pro_dates <= dt_forecast_end]


sp = snowprofilePro(pro_file, ProfileDate=pro_dates, tz = 'UTC')
meta = summary(sp)
smet = readSmet(smet_file)
smet$data = smet$data[smet$data$timestamp %in% meta$datetime, ]
sp = computePunstable(sp, ski_pen = smet$data$ski_pen)


## TIMESERIES
daterange = seq.POSIXt(dt_start, dt_forecast_end, by='hour')
daterange = daterange[daterange %in% meta$datetime]
if (length(daterange) == 0) stop("No profile data in specified time period.")

fname_base <- paste0(file.path(outdir, pro_name), "_", as.Date(daterange[1]), "_", as.Date(daterange[length(daterange)]), "_timeseries")
main = paste0(pro_name, ", ", format(max(meta$date), "%Y"), ", ", smet$altitude, " m asl, azimuth ", unique(meta$aspect))

png(filename = paste0(fname_base, ".png"), width = 600, height = 600)
par(mar = c(7.1, 4.1, 4.1, 3.1))
hs_max = max(meta$hs[meta$datetime %in% daterange])
plot(sp[meta$datetime %in% daterange], SortMethod = 'time', main = main, box = FALSE, Timeseries_labels='daily 23:00', yaxis=FALSE)
if (hs_max > 20) {
    yaxis_locs = seq(0, hs_max, by = 20) 
} else {
    yaxis_locs = seq(0, hs_max, by = 5)
}
axis(2, at = yaxis_locs, las = 1)
rect(xleft = dt_nowcast_end, xright = dt_forecast_end, ybottom = 0, ytop = hs_max, 
     col = NA, border = "black", lty = "dashed")
text(dt_nowcast_end, hs_max, labels="forecast", adj=c(0, 1))
# abline(v=dt_nowcast_end, lty='--')
# axis(1, at = as.Date(daterange), labels = format(daterange, "%b %d"), las = 2)
# abline(v = as.Date(daterange)[seq(1, length(daterange), by=2)], lty = 3, col = "dark grey")
dev.off()

png(filename = paste0(fname_base, "_punstable.png"), width = 600, height = 600)
par(mar = c(7.1, 4.1, 4.1, 3.1))
hs_max = max(meta$hs[meta$datetime %in% daterange])
plot(sp[meta$datetime %in% daterange], SortMethod = 'time', ColParam = 'p_unstable', main = main, box = FALSE, Timeseries_labels='daily 23:00', yaxis=FALSE)
if (hs_max > 20) {
    yaxis_locs = seq(0, hs_max, by = 20) 
} else {
    yaxis_locs = seq(0, hs_max, by = 5)
}
axis(2, at = yaxis_locs, las = 1)
rect(xleft = dt_nowcast_end, xright = dt_forecast_end, ybottom = 0, ytop = hs_max, 
     col = NA, border = "black", lty = "dashed")
text(dt_nowcast_end, hs_max, labels="forecast", adj=c(0, 1))
# axis(1, at = as.Date(daterange), labels = format(daterange, "%b %d"), las = 2)
# abline(v = as.Date(daterange)[seq(1, length(daterange), by=2)], lty = 3, col = "dark grey")
ClrRamp <- sapply(seq(0, 1, length.out = 100), function(alph) getColoursStability(alph, StabilityIndexThreshold = 0.77, StabilityIndexRange = c(0, 1)))
xleft <- par("usr")[2] - 0.015 * diff(par("usr")[1:2])
xright <- par("usr")[2]
ybottom <- par("usr")[3]
ytop <- par("usr")[4]
colorLevels <- seq(0, 1, length.out = length(ClrRamp))
for (i in 1:(length(colorLevels) - 1)) {
    rect(xleft, ybottom + colorLevels[i] * (diff(par("usr")[3:4]-0.05)),
        xright, ybottom + colorLevels[i + 1] * (diff(par("usr")[3:4]-0.05)),
        col = ClrRamp[i], border = NA
    )
}
axis(4, at = c(0, c(0.2, 0.6, 0.74, 0.87, 0.97)*(diff(par("usr")[3:4])-0.05)), labels = c(0, 0.2, 0.6, 0.77, 0.9, 1))
mtext("P_unstable", side = 4, line = 1.8)
dev.off()

png(filename = paste0(fname_base, "_lwc.png"), width = 600, height = 600)
par(mar = c(7.1, 4.1, 4.1, 3.1))
hs_max = max(meta$hs[meta$datetime %in% daterange])
plot(sp[meta$datetime %in% daterange], SortMethod = 'time', ColParam = 'lwc', main = main, box = FALSE, Timeseries_labels='daily 23:00', yaxis=FALSE)
if (hs_max > 20) {
    yaxis_locs = seq(0, hs_max, by = 20) 
} else {
    yaxis_locs = seq(0, hs_max, by = 5)
}
axis(2, at = yaxis_locs, las = 1)
rect(xleft = dt_nowcast_end, xright = dt_forecast_end, ybottom = 0, ytop = hs_max, 
     col = NA, border = "black", lty = "dashed")
text(dt_nowcast_end, hs_max, labels="forecast", adj=c(0, 1))
# axis(1, at = as.Date(daterange), labels = format(daterange, "%b %d"), las = 2)
# abline(v = as.Date(daterange)[seq(1, length(daterange), by=2)], lty = 3, col = "dark grey")
colorLevels <- seq(0, 6, length.out = 100)
ClrRamp <- getColoursLWC(colorLevels)
xleft <- par("usr")[2] - 0.015 * diff(par("usr")[1:2])
xright <- par("usr")[2]
ybottom <- par("usr")[3] + (diff(par("usr")[3:4])-0.05)*0.035
ytop <- par("usr")[4]
for (i in 1:(length(colorLevels) - 1)) {
    rect(xleft, ybottom + colorLevels[i]/6 * (diff(par("usr")[3:4]-0.05)),
        xright, ybottom + colorLevels[i + 1]/6 * (diff(par("usr")[3:4]-0.05)),
        col = ClrRamp[i], border = NA
    )
}
axis(4, at = c(0, seq(1, 6, by = 1)/6*(diff(par("usr")[3:4])-0.05)), labels = c(0, seq(1, 6, by = 1)))
mtext("LWC (%)", side = 4, line = 1.8, at = ybottom <- par("usr")[3] + (diff(par("usr")[3:4])-0.05)*0.46)
dev.off()








## HAND HARDNESS 

dt_start <- as.POSIXct(
  paste0(format(dt_nowcast_end, "%Y-%m-%d"), " ", daily_time),
  tz = "UTC"
)
if (dt_start > dt_nowcast_end) dt_start <- dt_start - as.difftime(1, units = "days")
daterange = seq.POSIXt(dt_start, dt_forecast_end, by='day')
daterange = daterange[daterange %in% meta$datetime]

for (i in seq_along(daterange)) {
    date = daterange[i]
    fname_base <- paste0(file.path(outdir, pro_name), "_", format(date, "%Y-%m-%dT%H:%M"))
    main = paste0(pro_name, ", ", format(date, "%Y-%m-%d %H:%M"), "\n", smet$altitude, " m asl, azimuth ", unique(meta$aspect))

    png(filename = paste0(fname_base, ".png"), width = 400, height = 600)
    par(mar = c(5.1, 4.1, 10.1, 1.1))
    plot(sp[[which(meta$datetime %in% date)]], main = main, nYTicks = max(2, round(meta$hs[which(meta$datetime %in% date)]/20)))
    mtext("Height (cm)", side = 2, line = 2.7)
    mtext("Hand hardness", side = 1, line = 2.3, at = 3)
    if ("temperature" %in% names(sp[[which(meta$datetime %in% date)]]$layers)) mtext("Temperature (Celsius)", side = 3, line = 2.3, at = 1)
    dev.off()
}
