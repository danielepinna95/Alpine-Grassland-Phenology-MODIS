# =============================================================================
# Step 4: Generalised Additive Models for NDVI phenology metrics
# =============================================================================
# One GAM per metric, using all continuous terrain, climate and soil drivers
# as cubic regression spline smooths (bs = "cs", k = 10, REML). The 80/20
# train/test split is stratified on three elevation classes (quantile bins).
#
# Input : data/field_metrics.csv
# Output: outputs/gam_performance.csv  - test-set R2 and RMSE per metric
#         outputs/gam_smooth_terms.csv - edf, F and p-value per smooth term
#         outputs/gam_<metric>.png     - observed vs predicted on the test set
# Usage : Rscript R/04_gam_models.R

suppressPackageStartupMessages({
  library(mgcv)
  library(dplyr)
  library(caret)
  library(ggplot2)
})

data_path <- "data/field_metrics.csv"
out_dir   <- "outputs"
dir.create(out_dir, showWarnings = FALSE)

drivers <- c("net_radiation_j_m2", "mean_temperature_k", "snowfall_total_m_we", "meanElevation",
             "meanSlope", "wind_speed_m_s", "mean_annual_precip_mm",
             "soil_soc", "soil_ph", "soil_wc33", "soil_clay", "soil_silt", "soil_sand")

metrics <- c(mean_ndvi            = "Mean NDVI",
             ndvi_peak            = "NDVI peak",
             normalized_amplitude = "Seasonal variability",
             sos_doy              = "Start of season (DOY)",
             peak_doy             = "Peak day (DOY)",
             eos_doy              = "End of season (DOY)")

r2   <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)
rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))

data <- read.csv(data_path)

performance <- list()
terms       <- list()

for (target in names(metrics)) {
  d <- na.omit(data[, c(drivers, target)]) |>
    mutate(elev_bin = cut(meanElevation, quantile(meanElevation, c(0, 1/3, 2/3, 1)),
                          include.lowest = TRUE, labels = c("Low", "Medium", "High")))

  set.seed(42)
  idx   <- createDataPartition(d$elev_bin, p = 0.8, list = FALSE)
  train <- d[idx, ]
  test  <- d[-idx, ]

  f   <- as.formula(paste(target, "~", paste0("s(", drivers, ", bs = 'cs')", collapse = " + ")))
  gam <- gam(f, data = train, method = "REML")

  test$predicted <- predict(gam, newdata = test)
  performance[[target]] <- data.frame(metric = target, n_train = nrow(train), n_test = nrow(test),
                                      R2 = r2(test[[target]], test$predicted),
                                      RMSE = rmse(test[[target]], test$predicted),
                                      deviance_explained = summary(gam)$dev.expl)

  st <- summary(gam)$s.table
  terms[[target]] <- data.frame(metric = target, term = rownames(st), edf = st[, "edf"],
                                F = st[, "F"], p_value = st[, "p-value"], row.names = NULL)

  p <- ggplot(test, aes(.data[[target]], predicted)) +
    geom_abline(linetype = "dashed", colour = "grey50") +
    geom_point(colour = "darkolivegreen", alpha = 0.8, size = 2.5) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.3, size = 5,
             label = sprintf("R² = %.3f\nRMSE = %.3g", performance[[target]]$R2, performance[[target]]$RMSE)) +
    labs(x = paste("Observed", metrics[[target]]), y = paste("Predicted", metrics[[target]])) +
    theme_classic(base_size = 14)
  ggsave(file.path(out_dir, paste0("gam_", target, ".png")), p, width = 6, height = 5, dpi = 300)
}

performance <- bind_rows(performance)
print(performance, digits = 3)
write.csv(performance, file.path(out_dir, "gam_performance.csv"), row.names = FALSE)
write.csv(bind_rows(terms), file.path(out_dir, "gam_smooth_terms.csv"), row.names = FALSE)
