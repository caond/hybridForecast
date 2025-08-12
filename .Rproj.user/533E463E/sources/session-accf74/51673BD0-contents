#' hybridForecast: A Hybrid Time Series Forecasting Package
#'
#' This package provides a modular, ensemble-based framework for time series forecasting.
#' It supports trend-seasonal decomposition, lag-based feature engineering, multiple machine learning models
#' (e.g. XGBoost, LightGBM, torch), and ensemble selection based on cross-validated accuracy.
#'
#' @section Main Functions:
#' - \code{\link{hybrid}}: Fit a hybrid forecasting model.
#' - \code{\link{predict_hybrid}}: Forecast future values using a trained hybrid model.
#' - \code{\link{plot.hybrid_train}}: Visualize CV training folds.
#' - \code{\link{simulate_time_series}}: Simulate time series with seasonal and changepoint components.
#'
#' @docType package
#' @name hybridForecast
#' @aliases hybridForecast-package
#' @keywords package
#'
#' @examples
#' # Simulate a time series
#' data <- simulate_time_series(freq = "week", length_out = 156, seed = 42)
#'
#' # Fit hybrid model
#' fit <- hybrid(data)
#'
#' # Plot CV results
#' plot(fit)
#'
#' # Predict future
#' future_df <- make_future_data(fit, horizon = 14)
#' forecast_hybrid <- predict(fit, future_df)
#' plot(forecast_hybrid)
"_PACKAGE"
