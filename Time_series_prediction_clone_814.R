# 1. Create training and testing sets
c814_train <- training(splits_814)
c814_test <- testing(splits_814)

# 2.Modelling Process
#Creating different time series models
## Percentage of GFP cells
### a. Auto ARIMA 
arima_fit_814_gfp <- arima_reg() %>% 
  set_engine("auto_arima") %>% 
  fit(GFP_pos ~ date.tb, data = clone_814)

### b. Boosted ARIMA 
arima_boost_fit_814_gfp <- arima_boost() %>% 
  set_engine("auto_arima_xgboost") %>% 
  fit(GFP_pos ~ date.tb, data = clone_814)

### c. Exponential Smoothing 
ets_fit_814_gfp <- exp_smoothing() %>% 
  set_engine("ets") %>% 
  fit(GFP_pos ~ date.tb, data = clone_814)

### d. Prophet 
prophet_fit_814_gfp <- prophet_reg() %>% 
  set_engine("prophet") %>% 
  fit(GFP_pos ~ date.tb, data = clone_814)

### e. Linear Regression
lm_fit_814_gfp <- linear_reg() %>% 
  set_engine("lm") %>% 
  fit(GFP_pos ~ date.tb, data = clone_814)

# Add all models to a table
models_tbl_814_gfp <- modeltime_table(
  arima_fit_814_gfp,
  arima_boost_fit_814_gfp,
  ets_fit_814_gfp,
  lm_fit_814_gfp,
  prophet_fit_814_gfp
)

# Calibrate models
calibrate_tbl_814_gfp <- models_tbl_814_gfp %>% 
  modeltime_calibrate(new_data = c814_test)

# Make current forecasts
calibrate_tbl_814_gfp %>% 
  modeltime_forecast(
    actual_data = clone_814,
    new_data = c814_test
  ) %>% 
  plot_modeltime_forecast()

# Check results
tbl_814_gfp <- calibrate_tbl_814_gfp %>% 
  modeltime_accuracy()

## Live cells GFP fluo. (Mean)
arima_fit_814_mean.live <- arima_reg() %>% 
  set_engine("auto_arima") %>% 
  fit(Mean_live ~ date.tb, data = clone_814)

### b. Boosted ARIMA 
arima_boost_fit_814_mean.live <- arima_boost() %>% 
  set_engine("auto_arima_xgboost") %>% 
  fit(Mean_live ~ date.tb, data = clone_814)

### c. Exponential Smoothing 
ets_fit_814_mean.live <- exp_smoothing() %>% 
  set_engine("ets") %>% 
  fit(Mean_live ~ date.tb, data = clone_814)

### d. Prophet 
prophet_fit_814_mean.live <- prophet_reg() %>% 
  set_engine("prophet") %>% 
  fit(Mean_live ~ date.tb, data = clone_814)

### e. Linear Regression
lm_fit_814_mean.live <- linear_reg() %>% 
  set_engine("lm") %>% 
  fit(Mean_live ~ date.tb, data = clone_814)

# Add all models to a table
models_tbl_814_mean.live <- modeltime_table(
  arima_fit_814_mean.live,
  arima_boost_fit_814_mean.live,
  ets_fit_814_mean.live,
  lm_fit_814_mean.live,
  prophet_fit_814_mean.live
)

# Calibrate models
calibrate_tbl_814_mean.live <- models_tbl_814_mean.live %>% 
  modeltime_calibrate(new_data = c814_test)

# Make current forecasts
calibrate_tbl_814_mean.live %>% 
  modeltime_forecast(
    actual_data = clone_814,
    new_data = c814_test
  ) %>% 
  plot_modeltime_forecast()

# Check results
tbl_814_mean.live <- calibrate_tbl_814_mean.live %>% 
  modeltime_accuracy()

## gated cell GFP fluo. (Mean)
arima_fit_814_mean.gated <- arima_reg() %>% 
  set_engine("auto_arima") %>% 
  fit(Mean ~ date.tb, data = clone_814)

### b. Boosted ARIMA 
arima_boost_fit_814_mean.gated <- arima_boost() %>% 
  set_engine("auto_arima_xgboost") %>% 
  fit(Mean ~ date.tb, data = clone_814)

### c. Exponential Smoothing 
ets_fit_814_mean.gated <- exp_smoothing() %>% 
  set_engine("ets") %>% 
  fit(Mean ~ date.tb, data = clone_814)

### d. Prophet 
prophet_fit_814_mean.gated <- prophet_reg() %>% 
  set_engine("prophet") %>% 
  fit(Mean ~ date.tb, data = clone_814)

### e. Linear Regression
lm_fit_814_mean.gated <- linear_reg() %>% 
  set_engine("lm") %>% 
  fit(Mean ~ date.tb, data = clone_814)

# Add all models to a table
models_tbl_814_mean.gated <- modeltime_table(
  arima_fit_814_mean.gated,
  arima_boost_fit_814_mean.gated,
  ets_fit_814_mean.gated,
  lm_fit_814_mean.gated,
  prophet_fit_814_mean.gated
)

# Calibrate models
calibrate_tbl_814_mean.gated <- models_tbl_814_mean.gated %>% 
  modeltime_calibrate(new_data = c814_test)

# Make current forecasts
calibrate_tbl_814_mean.gated %>% 
  modeltime_forecast(
    actual_data = clone_814,
    new_data = c814_test
  ) %>% 
  plot_modeltime_forecast()

# Check results
tbl_814_mean.gated <- calibrate_tbl_814_mean.gated %>% 
  modeltime_accuracy()
