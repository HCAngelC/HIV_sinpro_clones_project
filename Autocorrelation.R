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
  names(df.pacf) <- c("pacf")
  df.pacf <- df.pacf %>% dplyr::add_row(pacf = 0, .before = 1) %>% dplyr::mutate(Measurement = c(1:11))
  df$Measurement <- as.numeric(df$Measurement)
  df_fi <- dplyr::inner_join(df, df.pacf, by = "Measurement")
  return(df_fi)
}
