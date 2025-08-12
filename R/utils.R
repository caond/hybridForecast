model_used <-c('xgb','lgb','dl')
VERY_LARGE <- 1E34
EPSILON <- 1 / VERY_LARGE

#' Create Future Data Frame for Forecasting
#'
#' Generates a future \code{data.frame} of dates for use with a fitted
#' \code{hybridForecast_model}, with optional exogenous regressors.
#'
#' @param hybrid_model A fitted object from \code{\link{hybrid}} containing
#'   \code{$data} with a \code{ds} column and \code{$seasonal_periods$freq}.
#' @param xreg Optional exogenous regressors for the forecast horizon
#'   (\code{data.frame} or matrix). Must have \code{horizon} rows if provided.
#'   Default \code{NA}.
#' @param horizon Integer, number of future periods to generate (default 7).
#'
#' @return A \code{data.frame} with:
#' \itemize{
#'   \item \code{ds} – future date/time values
#'   \item additional columns – if \code{xreg} is supplied
#' }
#'
#' @details The start date is taken as the maximum \code{ds} in
#'   \code{hybrid_model$data}, incremented according to the model's seasonal frequency.
#'
#' @examples
#' \dontrun{
#' fit <- hybrid(my_data)
#' future_df <- make_future_data(fit, horizon = 14)
#' }
#'
#' @export
make_future_data <-function(hybrid_model,xreg=NA,horizon=7){
  future<-  tail(data.frame(ds=seq(max(hybrid_model$data$ds), by=hybrid_model$seasonal_periods$freq, length.out=horizon+1)),horizon)
  if (!is.na(xreg))
  {
    if (NROW(xreg)!=horizon) stop(paste0('Please make covariates rows= ',horizon))
    future<-cbind(future,xreg)
  }
  future

}

#' @export
simulate_time_series <- function(start_date = "2020-01-01",
                                 length_out = 156,
                                 freq = "week",
                                 trend_slopes = c(0.1, -0.05, 0.2, -0.1, 0.15),
                                 noise_sd = 3,
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  seasonal_periods<-seasonal_periods_fourier(freq)
  # Parse start date
  is_hourly <- freq %in% c("hour", "hours")
  start <- if (is_hourly) as.POSIXct(start_date) else as.Date(start_date)

  # Generate timestamp sequence
  if (is_hourly) {
    ds <- seq.POSIXt(from = start, by = freq, length.out = length_out)
  } else {
    ds <- seq.Date(from = start, by = freq, length.out = length_out)
  }

  # Seasonal components
  seasonal <- matrix(0, nrow = length_out, ncol = length(seasonal_periods$fouriers))
  colnames(seasonal) <- names(seasonal_periods$fouriers)

  for (i in seq_along(seasonal_periods$fouriers)) {
    period <- seasonal_periods$fouriers[[i]]
    seasonal[, i] <- sin(2 * pi * seq_along(ds) / period) * (10 / i)
  }

  seasonal_total <- rowSums(seasonal)

  # Trend component
  n_segments <- length(trend_slopes)
  change_points <- round(seq(1, length_out, length.out = n_segments + 1))
  trend <- numeric(length_out)

  for (i in 1:n_segments) {
    start <- change_points[i]
    end <- change_points[i + 1] - 1
    if (i == 1) {
      trend[start:end] <- trend_slopes[i] * seq(0, end - start)
    } else {
      trend[start:end] <- trend[start - 1] + trend_slopes[i] * seq(1, end - start + 1)
    }
  }

  # Noise
  noise <- rnorm(length_out, mean = 0, sd = noise_sd)

  # Combine
  y <- 50 + seasonal_total + trend + noise
  data.frame(ds = ds, y = y)
}


smape <- function(a, f) {

  a[which(a==Inf)]<- VERY_LARGE
  (1/length(a) * sum(2*abs(f-a) / (abs(a+EPSILON)+abs(f+EPSILON))))
}

auto_max_lag <- function(N, seasonal_periods) {
  stopifnot(N > 10L)
  freq <- seasonal_periods$freq

  # Aggressive caps since Fourier handles seasonality
  lag_cap <- switch(freq,
                    "hour"    = 48L,  # ~2 day
                    "day"     = 30L,
                    "week"    = 8L,
                    "month"   = 12L,
                    "quarter" = 6L,
                    "year"    = 2L,
                    stop("Unknown freq")
  )

  # Small, variance-friendly baseline for PACF
  baseline <- max(3L, floor(6 * log10(N)))
  n_cap    <- max(3L, floor(N / 6))   # stricter than N/4

  # We already use Fourier for seasonality → ignore seasonal lags in PACF
  seasonal_max <- 0L

  as.integer(min(max(baseline, seasonal_max), lag_cap, n_cap))
}



seasonal_periods_fourier<-function(freq)
{

  # Map to seasonal periods
  fouriers <- switch(freq,
                             "hour" = list(
                               day = 24,
                               week = 24 * 7,
                               year = 24 * 365.25
                             ),
                             "day" = list(
                               week = 7,
                               year = 365.25
                             ),
                             "week" = list(
                               year = 52.1775
                             ),
                             "month" = list(
                               year = 12
                             ),
                             "quarter" = list(
                               year = 4
                             ),
                             "year" = list()  # Usually no seasonality
  )

  return(list(freq=freq,fouriers=fouriers))

}

detect_frequency<-function(ds) {
  ds <- sort(unique(ds))
  diff_days <- as.numeric(median(diff(as.Date(ds)), na.rm = TRUE))

  # Define frequency thresholds (in days)
  freq <- case_when(
    diff_days < 1 ~ "hour",
    diff_days < 2 ~ "day",
    diff_days < 8 ~ "week",
    diff_days < 31 ~ "month",
    diff_days < 92 ~ "quarter",
    TRUE ~ "year"
  )
  return(freq)
}

auto_seasonal_periods <- function(ds, min_obs = 3) {
  # ds: vector of Date or POSIXct timestamps

  if (length(ds) < min_obs) stop("Not enough data points to detect frequency")
  freq<-detect_frequency(ds)

  message("Detected frequency: ", freq)

  # Map to seasonal periods
  seasonal_periods <- seasonal_periods_fourier(freq)

  return(seasonal_periods)
}

find_lags<-function(df,max_lag=10)
{
  adf_test<-adf.test(df$y)
  if (adf_test$p.value>0.05)  y <- diff(df$y)
  else y<- df$y
  pacf_vals <- pacf(y, lag.max = max_lag, plot = FALSE)
  threshold <- 2 / sqrt(length(y))
  significant_lags <- which(abs(pacf_vals$acf) > threshold)
  significant_lags
}


create_fourier_terms <- function(t, period, prefix, K) {

  terms <- list()
  for (k in 1:K) {
    terms[[paste0(prefix, "_sin", k)]] <- sin(2 * pi * k * t / period)
    terms[[paste0(prefix, "_cos", k)]] <- cos(2 * pi * k * t / period)
  }
  return(as.data.frame(terms))
}

select_changepoints <- function(y, n_changepoints = NA,freq) {
  n <- length(y)

  if (!is.numeric(y)) stop("Input 'y' must be numeric.")

  if (is.na(n_changepoints)) {

    # sensible caps by frequency
    minseg <- switch(freq,
                     "hour"     = 12,   # 2 days
                     "day"      = 7,   # 1 week
                     "week"     = 4,    # 1 month
                     "month"    = 3,    # 1 quarter
                     "quarter"  = 2,    # 6 months
                     "year"     = 1,    # 1 year
                     stop("Unknown frequency")
    )

    cpt <- changepoint.np::cpt.np(y,method    = "PELT", penalty   = "SIC", minseglen = minseg  )
    cps <- changepoint::cpts(cpt)
  } else {
    n_changepoints <- min(5, n_changepoints)
    t_range <- 1:floor(0.9 * n)
    cps <- quantile(t_range, probs = seq(0.05, 0.95, length.out = n_changepoints))
  }
  cat('\nnum of change points: ', length(cps))
  return(as.numeric(cps))
}

create_trend_features <- function(n, n_changepoints) {
  t <- 1:n

  if (length(n_changepoints) > 0) {
    A <- sapply(n_changepoints, function(sj) pmax(0, t - sj))
    X <- cbind(t, A)
    colnames(X) <- c("t", paste0("cp_", n_changepoints))
  } else {
    X <- matrix(t, ncol = 1)
    colnames(X) <- "t"
  }

  return(X)
}

auto_pick_K_all_components <- function(t, seasonal_periods, max_fraction = 0.05, max_K_cap = 5, min_K = 1) {
  fouriers<-seasonal_periods$fouriers
  freq<-seasonal_periods$freq
  n <- length(t)
  max_total_terms <- floor(n * max_fraction)

  K<-sapply(fouriers, function(period) {
    n_cycles <- n / period
    K <- floor(min(max_total_terms, n_cycles * 2))  # 2 terms per K (sin/cos)
    K <- min(K, max_K_cap)
    K <- max(K, min_K)
    return(K)
  })


  if ( (freq== "month") && ("year" %in% names(fouriers))) {
    K["year"] <- 4L
  }

  K
}

#' Fit trend + seasonal model using ridge regression (glmnet, alpha = 0)
#' @import glmnet
fit_trend_seasonal_model <- function(y, seasonal_periods, n_changepoints) {
  n <- length(y)
  t <- 1:n
  fouriers<-seasonal_periods$fouriers
  freq<-seasonal_periods$freq
  # --- Trend terms ---
  changepoints <- select_changepoints(y, n_changepoints, freq)
  trend_X <- create_trend_features(n, changepoints)  # columns: "t", "cp_<idx>"

  # --- Seasonal terms (Fourier) ---
  K_list <- auto_pick_K_all_components(t, seasonal_periods)
  seasonal_X_list <- lapply(names(fouriers), function(name) {
    create_fourier_terms(t, fouriers[[name]], prefix = name, K = K_list[[name]])
  })
  if (length(seasonal_X_list) > 0) {
    seasonal_X <- do.call(cbind, seasonal_X_list)
  } else {
    seasonal_X <- matrix(, nrow = n, ncol = 0)
  }

  # --- Design matrix ---
  X <- cbind(trend_X, seasonal_X)
  X_mat <- as.matrix(X)

  # --- Lasso regression via cv.glmnet ---

  cvfit <- glmnet::cv.glmnet(
    x = X_mat,
    y = y,
    alpha = 1,
    intercept = FALSE,
    standardize = TRUE,
    nfolds=5
  )

  # Coefficients at lambda.min (you can switch to "lambda.1se" if you prefer stronger shrinkage)
  beta_mat <- coef(cvfit, s = "lambda.min")
  betas <- drop(beta_mat[, 1])                 # numeric vector
  names(betas) <- rownames(beta_mat)
  betas <- betas[setdiff(names(betas), "(Intercept)")]  # safety: drop intercept if present

  # --- Components ---
  # Trend
  b_trend <- betas[colnames(trend_X)]
  b_trend[is.na(b_trend)] <- 0
  trend_vals <- as.vector(as.matrix(trend_X) %*% b_trend)

  # Seasonality per component
  seasonality_outputs <- list()
  if (length(seasonal_X_list) > 0) {
    for (i in seq_along(seasonal_X_list)) {
      mat <- as.matrix(seasonal_X_list[[i]])
      cols <- colnames(mat)
      b_here <- betas[cols]
      b_here[is.na(b_here)] <- 0
      seasonality_outputs[[i]] <- as.vector(mat %*% b_here)
    }
    names(seasonality_outputs) <- names(fouriers)
  }

  return(list(
    model = cvfit,                    # cv.glmnet object
    changepoints = changepoints,
    n_train = n,
    K_list = K_list,                  # chosen Fourier K per component
    coef_names = colnames(X),
    trend = trend_vals,
    fouriers=fouriers,
    seasonality = seasonality_outputs
  ))
}


predict_trend_seasonal_components <- function(model_obj, h) {
  n_future <- model_obj$n_train + h
  t_pred   <- (model_obj$n_train + 1):n_future

  ## --- Trend design (matches training) ---
  if (length(model_obj$changepoints) > 0) {
    trend_X <- sapply(model_obj$changepoints, function(sj) pmax(0, t_pred - sj))
    trend_X <- cbind(t = t_pred, trend_X)
    colnames(trend_X) <- c("t", paste0("cp_", model_obj$changepoints))
  } else {
    trend_X <- matrix(t_pred, ncol = 1)
    colnames(trend_X) <- "t"
  }

  ## --- Seasonal design (matches training) ---
  K_list <- model_obj$K_list
  season_components <- list()
  if (length(model_obj$fouriers) > 0) {
    for (name in names(model_obj$fouriers)) {
      season_components[[name]] <- create_fourier_terms(
        t_pred,
        period = model_obj$fouriers[[name]],
        prefix = name,
        K = K_list[[name]]
      )
    }
  }
  if (length(season_components) > 0) {
    seasonal_X <- do.call(cbind, unname(season_components))
  } else {
    # no seasonal terms
    seasonal_X <- NULL
  }

  ## --- Full design for prediction (same order of columns as training) ---
  full_X <- cbind(trend_X, seasonal_X)
  full_X_mat <- as.matrix(full_X)

  ## --- Predict with ridge (cv.glmnet) at lambda.min ---
  preds <- as.numeric(predict(model_obj$model, newx = full_X_mat, s = "lambda.min"))

  ## --- Components via coefficients ---
  # coef(cvfit, s=...) returns a sparse matrix with rownames including "(Intercept)"
  beta_mat <- coef(model_obj$model, s = "lambda.min")
  betas    <- drop(beta_mat[, 1])
  names(betas) <- rownames(beta_mat)
  betas <- betas[names(betas) != "(Intercept)"]

  ## Trend component
  b_trend <- betas[colnames(trend_X)]
  b_trend[is.na(b_trend)] <- 0
  component_trend <- as.numeric(as.matrix(trend_X) %*% b_trend)

  ## Seasonal components (per seasonal group)
  component_season <- list()
  if (length(season_components) > 0) {
    for (name in names(model_obj$fouriers)) {
      mat <- as.matrix(season_components[[name]])
      cols <- colnames(mat)
      b_here <- betas[cols]
      b_here[is.na(b_here)] <- 0
      component_season[[name]] <- as.numeric(mat %*% b_here)
    }
  }

  list(
    forecast   = preds,
    trend      = component_trend,
    seasonality = component_season,
    t          = t_pred
  )
}


prepare_trend_seasonal_component<-function(train,significant_lags,seasonal_periods,n_changepoints)
{

  ts_model <- fit_trend_seasonal_model(y=train$y,seasonal_periods=seasonal_periods, n_changepoints=n_changepoints)

  train_feat <- cbind(train, trend = ts_model$trend)

  # Add each seasonal component
  for (name in names(ts_model$seasonality)) {
    train_feat[[name]] <- ts_model$seasonality[[name]]
  }

  for (i in significant_lags) {
    train_feat[[paste0("lag", i)]] <- dplyr::lag(train$y, i)
  }
  train_feat <- na.omit(train_feat)
  return(list(train_feat=train_feat,ts_model=ts_model))
}

#' @import torch
#Define the neural network module
Net <- nn_module(
  initialize = function(input_dim, hidden_dim = 64, output_dim = 1) {
    self$fc1 <- nn_linear(input_dim, hidden_dim)
    self$fc2 <- nn_linear(hidden_dim, output_dim)
  },
  forward = function(x) {
    x %>%
      self$fc1() %>%
      nnf_relu() %>%
      self$fc2()
  }
)


# Train the model
train_torch_model <- function(X, y, epochs = 100, lr = 0.01, hidden_dim = 64) {


  # Auto-detect device (GPU if available, otherwise CPU)
  device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")

  # Convert data to tensors
  X_tensor <- torch_tensor(as.matrix(X), dtype = torch_float(), device = device)
  y_tensor <- torch_tensor(as.numeric(y), dtype = torch_float(), device = device)$unsqueeze(2)

  # Initialize model and optimizer
  input_dim <- ncol(X)
  model <- Net(input_dim = input_dim, hidden_dim = hidden_dim)$to(device = device)

  optimizer <- optim_adam(model$parameters, lr = lr)
  loss_fn <- nn_mse_loss()
  # Training loop
  for (epoch in 1:epochs) {
    model$train()
    optimizer$zero_grad()
    output <- model(X_tensor)
    loss <- loss_fn(output, y_tensor)
    loss$backward()
    optimizer$step()
  }

  return(list(model = model, device = device))
}

# Predict using trained model
predict_torch_model <- function(model_obj, X_new) {
  model <- model_obj$model
  device <- model_obj$device

  X_tensor <- torch_tensor(as.matrix(X_new), dtype = torch_float(), device = device)

  model$eval()
  preds <- model(X_tensor)$to(device = "cpu")  # move back to CPU
  as.numeric(preds)
}

