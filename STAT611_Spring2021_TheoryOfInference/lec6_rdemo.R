setwd("D:/Google Drive/2021_SPRING/STAT611")
library(fitdistrplus)
data(groundbeef)
serving <- groundbeef$serving

# Plot the empirical distribution
plotdist(serving, histo = TRUE, demp = TRUE)

fit_w  <- fitdist(serving, "weibull")
fit_g  <- fitdist(serving, "gamma")
fit_ln <- fitdist(serving, "lnorm")

rbind(fit_w$estimate, fit_g$estimate, fit_ln$estimate)

# Diagnostics
pdf('diag.pdf', width = 7, height = 6)
par(mfrow=c(2,2))
plot.legend <- c("Weibull", "lognormal", "gamma")
denscomp(list(fit_w, fit_g, fit_ln), legendtext = plot.legend)
cdfcomp (list(fit_w, fit_g, fit_ln), legendtext = plot.legend)
qqcomp  (list(fit_w, fit_g, fit_ln), legendtext = plot.legend)
ppcomp  (list(fit_w, fit_g, fit_ln), legendtext = plot.legend)
dev.off()
