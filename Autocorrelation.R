##############################################################
# Author: Heng-Chang Chen
##############################################################
# Input: a data frame containing the percentage of GFP-positive cells measured over time (in time series).
# Object: Lag value for time series measurement.
##############################################################

# R functions

autocorrelation <- function(df) {
  df.pacf <- pacf(df$GFP_pos, plot = F)
  df.pacf <- data.frame(df.pacf$acf) #calculate lag for (%) GFP-positive cells
  names(df.pacf) <- c("gfp_lag")
  df.pacf <- df.pacf %>% dplyr::mutate(lag = c(1:10))
  df_ <- df %>% dplyr::mutate(cv_lag = CV_FITC - lag(CV_FITC), sd_lag = SD_FITC - lag(SD_FITC)) %>% dplyr::select(cv_lag, sd_lag) #calculate lag for CV and SD
  df_ <- df_[2:11, ]
  df_ <- df_ %>% dplyr::mutate(lag = c(1:10))
  df_fi <- dplyr::inner_join(df_, df.pacf, by = "lag")
  return(df_fi)
}
