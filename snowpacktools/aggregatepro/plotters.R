plotTradAvgProfile <- function(avgSP, pdate, btLyrs=NA) {

    ## ---setup plot & preprocess----------------------------------------------
    opar <- par()
    on.exit(suppressWarnings(par(opar)))
    layout(matrix(c(2, 3, 1), 1, 3, byrow = TRUE), c(1.2, 1.2, 1.8))
    # layout(matrix(c(2, 1), 1, 2, byrow = TRUE), c(1.2, 1.8))
    vmar <- 17
    legloc <- -0.28

    k <- which(as.Date(avgSP$meta$date) %in% as.Date(pdate))
    if (length(k) != 1) stop("Either pdate not contained in average profile, or date is ambiguous.")
    avg <- avgSP$avgs[[k]]
    set <- snowprofileSet(avgSP$sets[[k]])
    if (is.na(btLyrs)) {
        btLyrs <- backtrackLayers(avg, profileSet = set)
    }

    ## compute layer distributions
    nlRTASK05 <- (do.call("c", lapply(btLyrs, function(df) {
        sum(df$sk38[df$rta >= 0.8] <= 0.5)
    })))
    nlRTASK95 <- (do.call("c", lapply(btLyrs, function(df) {
        sum(df$sk38[df$rta >= 0.8] <= 0.95)
    })))
    nlRTASK105 <- (do.call("c", lapply(btLyrs, function(df) {
        sum(df$sk38[df$rta >= 0.8] <= 1.05)
    })))
    nlPU077 <- (do.call("c", lapply(btLyrs, function(df) {
        sum(df$p_unstable >= 0.77)
    })))
    
    ## ---plot average profile-------------------------------------------------
    par(mar = c(vmar, 2.1, 4.1, 1.1))
    plot(avg, axes = FALSE, xlab = "")
    axis(1, at = seq(5), labels = c("F", "4F", "1F", "P", "K"))
    axis(2, pretty(c(0, avg$hs)))
    mtext("Hardness", side = 1, line = 5.5, cex = opar$cex.lab)
    mtext(paste0("valid: ", format(avg$date, "%d.%m.%Y")), 
                 side = 1, line = 13.5, cex = opar$cex.lab,
                 col = "gray50")

    ## ---stacked barplot RTA SK38---------------------------------------------
    par(mar = c(vmar, 5.5, 4.1, 0.5))

    barplot(-nlRTASK105 / length(set),
        width = c(avg$layers$height[1], diff(avg$layers$height)), col = "gray70",
        horiz = TRUE, border = NA, space = 0, ylim = c(0, avg$hs), xlim = c(-1, 0), 
        axes = FALSE, col.axis = "transparent", cex.lab = opar$cex.lab
    )
    barplot(-nlRTASK95 / length(set),
        width = c(avg$layers$height[1], diff(avg$layers$height)), col = "gray40",
        horiz = TRUE, border = NA, space = 0, col.axis = "transparent", add = TRUE
    )
    barplot(-nlRTASK05 / length(set),
        width = c(avg$layers$height[1], diff(avg$layers$height)), col = "gray0",
        horiz = TRUE, border = NA, space = 0, col.axis = "transparent", add = TRUE
    )
    mtext("RTA & SK38", side = 3, line = 1.5, cex = opar$cex.lab)
    mtext("Proportion of individual profiles", side = 1, line = 5.5, at = 0,
          cex = opar$cex.lab)
    # mtext("(a)", side = 3, line = 2, cex = 1.7 + cex.offset, at = -1)
    mtext("Height (cm)", side = 2, line = 3.5, cex = opar$cex.lab)
    axis(1, at = -1 * c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
    axis(2, at = pretty(c(0, avg$hs)), las = 1)
    grid()
    legend("bottomleft", legend = c("very poor (<= 0.5)", "poor (<= 0.95)", "fair (<= 1.05)"), 
           fill = c("gray0", "gray40", "gray70"), border = c("gray0", "gray40", "gray70"), 
           bty = "o", box.lwd = 0, cex = opar$cex.lab, xpd = TRUE, inset = c(0, legloc))
    
    ## ---stacked barplot PU---------------------------------------------------
    par(mar = c(vmar, 3.1, 4.1, 0.5))
    barplot(-nlPU077 / length(set),
        width = c(avg$layers$height[1], diff(avg$layers$height)), col = "gray20",
        horiz = TRUE, border = NA, space = 0, ylim = c(0, avg$hs), xlim = c(-1, 0), 
        axes = FALSE, col.axis = "transparent", cex.lab = opar$cex.lab
    )
    mtext("PU", side = 3, line = 1.5, cex = opar$cex.lab)
    # mtext("(e)", side = 3, line = 2, cex = 1.7 + cex.offset, at = -1)
    axis(1, at = -1 * c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c(0, 0.2, 0.4, 0.6, 0.8, 1))
    axis(2, at = pretty(c(0, avg$hs)))
    grid()
    legend("bottomleft", legend = c("poor (>= 0.77)", ""), 
           fill = c("gray20", "transparent"), border = c("gray20", "transparent"), 
           bty = "o", box.lwd = 0, cex = opar$cex.lab, xpd = TRUE, inset = c(0, legloc))
}  # END plotTradAvgProfile



plotTSplainAvgProfile <- function(avgSP) {

    ## ---setup plot & preprocess----------------------------------------------
    opar <- par()
    on.exit(suppressWarnings(par(opar)))
    par(mar = c(6.1, 6.1, 3.1, 0.1))

    ## ---plot average profile-------------------------------------------------
    plot(avgSP$avgs, box = FALSE, xaxs = "i", yaxs = "i", yaxt = "n", ylab = "", yaxis = FALSE)
    lines(avgSP$meta$date, avgSP$meta$hs_median, lwd = 2)
    lines(avgSP$meta$date, avgSP$meta$hs_median-avgSP$meta$thicknessPPDF_median, lwd = 2, lty = "dashed")
    mtext("Height (cm)", side = 2, line = 4, cex = opar$cex.lab)
    # mtext("Aggregated and predominant snowpack conditions", side = 3, line = 1, cex = opar$cex.lab)
    axis(2, at = pretty(c(0, max(avgSP$meta$hs))), las = 1)
    legend("topleft",
        c("<HS>", "<PP/DF>", "SH", "DH", "FC", "FCxr", "RG", "PP", "DF", "MF", "MFcr"),
        lty = c("solid", "dashed", rep(NA, 9)), lwd = 2,
        fill = getColoursGrainType(c(rep(NA, 3), "SH", "DH", "FC", "FCxr", "RG", "PP", "DF", "MF", "MFcr")),
        density = c(rep(0, 3), rep(NA, 9)), border = "transparent",
        horiz = FALSE, bty = "o", box.lwd = 0, cex = opar$cex.lab
    )
}

