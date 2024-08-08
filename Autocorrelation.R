##############################################################
# Author: Heng-Chang Chen
##############################################################
# Input: a data frame containing the percentage of GFP-positive cells measured over time (in time series).
# Object: Lag value for time series measurement.
##############################################################

# R functions

autocorrelation <- function(df) {
  df.pacf <- pacf(df$GFP_pos, plot = F)
  df.pacf <- data.frame(df.pacf$acf)
  names(df.pacf) <- c("gfp_autocor")
  df.pacf <- df.pacf %>% dplyr::mutate(lag = c(1:10))
  df_ <- df %>% dplyr::mutate(gfp_pos_lag = lag(GFP_pos), cv_lag = CV_FITC - lag(CV_FITC), sd_lag = SD_FITC - lag(SD_FITC), mean_lag = Mean_FITC - lag(Mean_FITC), median_lag = Median_FITC - lag(Median_FITC)) %>% dplyr::select(gfp_pos_lag, cv_lag, sd_lag, mean_lag, median_lag)
  df_ <- df_[2:11, ]
  df_ <- df_ %>% dplyr::mutate(lag = c(1:10))
  df_fi <- dplyr::inner_join(df_, df.pacf, by = "lag")
  return(df_fi)
}
