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
  P <- unlist(seasonal_periods$fouriers, use.names = FALSE)

  # generous frequency caps (cover persistence, not seasonality)
  freq_cap <- switch(freq,
                     "hour"    = 168L,  # up to 1 week
                     "day"     = 90L,   # ~3 months
                     "week"    = 26L,   # ~6 months
                     "month"   = 24L,   # 2 years
                     "quarter" = 12L,   # 3 years
                     "year"    = 5L,    # 5 years
                     stop("Unknown freq")
  )

  # size-aware caps
  base_cap <- max(5L, floor(3 * sqrt(N)))
  n_cap    <- max(5L, floor(N / 3))

  # anchor to a few cycles of the shortest Fourier period (if any)
  seasonal_anchor <- if (length(P)) max(6L, floor(3 * min(P))) else 0L

  as.integer(min(max(base_cap, seasonal_anchor), freq_cap, n_cap))
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
  freq <- dplyr::case_when(
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
  adf_test<-tseries::adf.test(df$y)
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
  cat('\nNumber of change points: ', length(cps))
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


prepare_trend_seasonal_ar<-function(train,significant_lags,seasonal_periods,n_changepoints)
{

  ts_model <- fit_trend_seasonal_model(y=train$y,seasonal_periods=seasonal_periods, n_changepoints=n_changepoints)

  train_feat <- cbind(train, trend = ts_model$trend)

  # Add each seasonal component
  for (name in names(ts_model$seasonality)) {
    train_feat[[name]] <- ts_model$seasonality[[name]]
  }
  # Add auto regressive
  for (i in significant_lags) {
    train_feat[[paste0("lag", i)]] <- dplyr::lag(train$y, i)
  }
  train_feat <- na.omit(train_feat)
  return(list(train_feat=train_feat,ts_model=ts_model))
}

#' @import torch
#Define the neural network module
# Net <- nn_module(
#   initialize = function(input_dim, hidden_dim = 64) {
#     self$fc1 <- nn_linear(input_dim, hidden_dim)
#     self$drop <- torch::nn_dropout(p = 0.1)
#     self$fc2 <- nn_linear(hidden_dim, 1)
#   },
#   forward = function(x) {
#     x %>%
#       self$fc1() %>%
#       nnf_relu() %>%
#       self$fc2()
#   }
# )



# Predict using trained model
predict_torch_model <- function(model_obj, X_new) {
  model <- model_obj$model
  device <- model_obj$device

  X_tensor <- torch_tensor(as.matrix(X_new), dtype = torch_float(), device = device)

  model$eval()
  preds <- model(X_tensor)$to(device = "cpu")  # move back to CPU
  as.numeric(preds)
}


# Train the model
train_torch_model <- function(X, y, epochs = 100, lr = 0.001, hidden_dim = 64) {


  # Auto-detect device (GPU if available, otherwise CPU)
  device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")


  # Scale
  #  x_mean <- colMeans(X)
  # x_sd   <- apply(X, 2, sd); x_sd[x_sd == 0] <- 1
  #Xs     <- sweep(sweep(X, 2, x_mean, "-"), 2, x_sd, "/")

  # Convert data to tensors
  X_tensor <- torch_tensor(as.matrix(X), dtype = torch_float(), device = device)
  y_tensor <- torch_tensor(as.numeric(y), dtype = torch_float(), device = device)$unsqueeze(2)

  # Initialize model and optimizer
  input_dim <- ncol(X)
  model <- Net(input_dim, hidden_dim)$to(device = device)

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


Net <- torch::nn_module(
  "Net",
  initialize = function(input_dim, hidden_dim = 64L, p_drop = 0.1) {
    self$fc1  <- torch::nn_linear(input_dim, as.integer(hidden_dim))
    self$drop <- torch::nn_dropout(p = p_drop)
    self$fc2  <- torch::nn_linear(as.integer(hidden_dim), 1L)
  },
  forward = function(x) {
    x %>%
      self$fc1() %>%
      torch::nnf_relu() %>%
      self$drop() %>%          # keep dropout; will auto-disable in eval()
      self$fc2()
  }
)

# ---------- Utilities ----------

make_split <- function(n, val_frac = 0.1) {
  stopifnot(val_frac >= 0, val_frac < 1)
  if (val_frac == 0) return(list(fit_idx = 1:n, val_idx = integer(0)))
  train_end <- floor((1 - val_frac) * n)
  list(fit_idx = 1:train_end, val_idx = (train_end + 1):n)
}

# Deterministic torch device selection (only used by torch trainer)
torch_get_device_deterministic <- function(seed = 1234L) {
  torch::torch_manual_seed(seed)
  if (torch::cuda_is_available()) {
    torch::use_deterministic_algorithms(TRUE)
    torch::torch_backends_cudnn$deterministic <- TRUE
    torch::torch_backends_cudnn$benchmark <- FALSE
    torch::torch_device("cuda")
  } else {
    torch::torch_device("cpu")
  }
}

bayes_optimize <- function(score_fn, bounds,
                           init_points = 5, n_iter = 10,
                           acq = "ucb", kappa = 2.576, seed = 123) {
  set.seed(seed)
  rBayesianOptimization::BayesianOptimization(
    FUN         = function(...) list(Score = score_fn(list(...))),
    bounds      = bounds,
    init_points = init_points,
    n_iter      = n_iter,
    acq         = acq,
    kappa       = kappa,
    verbose     = FALSE
  )$Best_Par
}

# ---------- Generic Bayes-train orchestrator ----------
# Fits with Bayesian Optimization on a holdout slice, then refits on all data.
# - fit_once(params, fit_idx, val_idx) -> list(score=..., model=..., extras=...)
#   During search you should return $score; during refit you can return $model.
# - refit(params) -> list(model=..., extras=...)

bayes_train <- function(X, y, bounds,
                        fit_once, refit,
                        val_frac = 0.1,
                        bo_seed = 123,
                        init_points = 20, n_iter = 10) {
  n <- nrow(X)
  idx <- make_split(n, val_frac)

  # 1) Define scoring fn for BO (maximize)
  score_fn <- function(par_list) {
    out <- fit_once(par_list, idx$fit_idx, idx$val_idx)
    # Expect a scalar numeric to maximize
    as.numeric(out$score)
  }

  # 2) Search
  best <- bayes_optimize(
    score_fn, bounds,
    init_points = init_points, n_iter = n_iter, seed = bo_seed
  )

  # 3) Refit on all data
  final <- refit(as.list(best))

  list(model = final$model, params = best, extras = final$extras)
}

# =====================================================
# ===============  TORCH (nn) TRAINER  ================
# =====================================================

train_torch_model_bayes <- function(X, y,
                                    epochs = 100L,
                                    val_frac = 0.1,
                                    bounds = list(
                                      lr           = c(1e-3, 1e-2),
                                      hidden_dim   = c(32L, 128L),
                                      p_drop       = c(0.05, 0.3),
                                      weight_decay = c(1e-6, 1e-3)
                                    )) {
  set.seed(1234)
  device <- torch_get_device_deterministic(1234L)

  X <- as.matrix(X)
  y <- as.numeric(y)
  input_dim <- ncol(X)
  loss_fn <- torch::nn_mse_loss()

  fit_core <- function(Xt, yt, params) {
    model <- Net(input_dim, as.integer(round(params$hidden_dim)), p_drop = params$p_drop)$to(device = device)
    opt   <- torch::optim_adam(model$parameters, lr = params$lr, weight_decay = params$weight_decay)

    for (ep in seq_len(epochs)) {
      model$train()
      opt$zero_grad()
      pred <- model(Xt)
      loss <- loss_fn(pred, yt)
      loss$backward()
      opt$step()
    }
    model
  }

  fit_once <- function(params, fit_idx, val_idx) {
    X_fit <- torch::torch_tensor(X[fit_idx,, drop = FALSE], dtype = torch::torch_float(), device = device)
    y_fit <- torch::torch_tensor(y[fit_idx], dtype = torch::torch_float(), device = device)$unsqueeze(2)

    model <- fit_core(X_fit, y_fit, params)

    if (length(val_idx)) {
      X_val <- torch::torch_tensor(X[val_idx,, drop = FALSE], dtype = torch::torch_float(), device = device)
      y_val <- torch::torch_tensor(y[val_idx], dtype = torch::torch_float(), device = device)$unsqueeze(2)
      model$eval()
      score <- -as.numeric(loss_fn(model(X_val), y_val)$item())  # maximize negative MSE
      return(list(score = score, model = NULL))
    } else {
      return(list(score = NA_real_, model = model))
    }
  }

  refit <- function(params) {
    X_all <- torch::torch_tensor(X, dtype = torch::torch_float(), device = device)
    y_all <- torch::torch_tensor(y, dtype = torch::torch_float(), device = device)$unsqueeze(2)
    model <- fit_core(X_all, y_all, params)
    list(model = model, extras = list(device = device))
  }

  bayes_train(X, y, bounds, fit_once, refit, val_frac = val_frac)
}


# =====================================================
# ==============  XGBOOST (regression)  ===============
# =====================================================

train_xgboost_model_bayes <- function(X, y,
                                      val_frac = 0.1,
                                      bounds = list(
                                        max_depth        = c(2L, 8L),
                                        min_child_weight = c(1L, 10L),
                                        subsample        = c(0.5, 0.9)
                                      ),
                                      nrounds = 100, eta = 0.1,
                                      early_stopping_rounds = 20) {
  stopifnot(nrow(X) == length(y))
  X <- as.matrix(X); y <- as.numeric(y)

  fit_once <- function(params, fit_idx, val_idx) {
    dtrain <- xgboost::xgb.DMatrix(data = X[fit_idx,, drop = FALSE], label = y[fit_idx])
    params_xgb <- list(
      objective        = "reg:squarederror",
      eval_metric      = "rmse",
      eta              = eta,
      max_depth        = as.integer(params$max_depth),
      min_child_weight = as.integer(params$min_child_weight),
      subsample        = params$subsample,
      colsample_bytree = 0.8
    )

    if (length(val_idx)) {
      dvalid <- xgboost::xgb.DMatrix(data = X[val_idx,, drop = FALSE], label = y[val_idx])
      watch <- list(train = dtrain, valid = dvalid)
      bst <- xgboost::xgb.train(
        params  = params_xgb, data = dtrain, nrounds = nrounds,
        watchlist = watch, early_stopping_rounds = early_stopping_rounds, verbose = 0
      )
      preds <- predict(bst, X[val_idx,, drop = FALSE])
      score <- -smape(y[val_idx], preds)  # maximize
      list(score = score, model = NULL)
    } else {
      bst <- xgboost::xgb.train(params = params_xgb, data = dtrain, nrounds = nrounds, verbose = 0)
      list(score = NA_real_, model = bst)
    }
  }

  refit <- function(params) {
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = y)
    bst <- xgboost::xgb.train(
      params = list(
        objective        = "reg:squarederror",
        eval_metric      = "rmse",
        eta              = eta,
        max_depth        = as.integer(params$max_depth),
        min_child_weight = as.integer(params$min_child_weight),
        subsample        = params$subsample,
        colsample_bytree = 0.8
      ),
      data = dtrain, nrounds = nrounds, verbose = 0
    )
    list(model = bst, extras = list(
      best_iteration = bst$best_iteration %||% xgboost::xgb.attributes(bst)[["best_iteration"]]
    ))
  }

  bayes_train(X, y, bounds, fit_once, refit, val_frac = val_frac)
}

# =====================================================
# ==============  LIGHTGBM (regression)  ==============
# =====================================================

train_lightgbm_model_bayes <- function(X, y,
                                       val_frac = 0.1,
                                       bounds = list(
                                         max_depth        = c(4L, 10L),
                                         num_leaves       = c(8L, 255L),
                                         min_data_in_leaf = c(20L, 600L),
                                         feature_fraction = c(0.6, 1.0),
                                         bagging_fraction = c(0.6, 1.0)
                                       ),
                                       nrounds = 100,
                                       learning_rate = 0.05,
                                       early_stopping_rounds = 20) {
  stopifnot(nrow(X) == length(y))
  X <- as.matrix(X); y <- as.numeric(y)

  fit_once <- function(params, fit_idx, val_idx) {
    dtrain <- lightgbm::lgb.Dataset(data = X[fit_idx,, drop = FALSE], label = y[fit_idx])
    valids <- NULL

    prm <- list(
      objective = "regression",
      metric = "l2",
      learning_rate = learning_rate,
      max_depth = as.integer(params$max_depth),
      num_leaves = as.integer(params$num_leaves),
      min_data_in_leaf = as.integer(params$min_data_in_leaf),
      feature_fraction = params$feature_fraction,
      bagging_fraction = params$bagging_fraction
    )

    if (length(val_idx)) {
      dvalid <- lightgbm::lgb.Dataset(data = X[val_idx,, drop = FALSE], label = y[val_idx])
      valids <- list(train = dtrain, valid = dvalid)
      bst <- lightgbm::lgb.train(
        params = prm, data = dtrain, nrounds = nrounds,
        valids = valids, early_stopping_rounds = early_stopping_rounds, verbose = 0
      )
      preds <- predict(bst, X[val_idx,, drop = FALSE])
      score <- -smape(y[val_idx], preds)  # maximize
      list(score = score, model = NULL)
    } else {
      bst <- lightgbm::lgb.train(params = prm, data = dtrain, nrounds = nrounds, verbose = 0)
      list(score = NA_real_, model = bst)
    }
  }

  refit <- function(params) {
    dtrain <- lightgbm::lgb.Dataset(data = X, label = y)
    bst <- lightgbm::lgb.train(
      params = list(
        objective = "regression",
        metric = "l2",
        learning_rate = learning_rate,
        max_depth = as.integer(params$max_depth),
        num_leaves = as.integer(params$num_leaves),
        min_data_in_leaf = as.integer(params$min_data_in_leaf),
        feature_fraction = params$feature_fraction,
        bagging_fraction = params$bagging_fraction
      ),
      data = dtrain, nrounds = nrounds, verbose = 0
    )
    list(model = bst, extras = list(best_iter = bst$best_iter))
  }

  bayes_train(X, y, bounds, fit_once, refit, val_frac = val_frac)
}

# ---------- Random Search Optimizer ----------
random_optimize <- function(score_fn, bounds, n_iter = 15, seed = 123) {
  set.seed(seed)
  best_score <- -Inf
  best_params <- NULL

  for (i in 1:n_iter) {
    # Sample random parameters from bounds
    params <- lapply(names(bounds), function(param_name) {
      b <- bounds[[param_name]]

      # Check if integer bounds (by checking if both values are integers)
      if (is.integer(b[1]) && is.integer(b[2])) {
        as.integer(round(runif(1, b[1], b[2])))
      } else {
        runif(1, b[1], b[2])
      }
    })
    names(params) <- names(bounds)

    # Evaluate score
    score <- tryCatch({
      score_fn(params)
    }, error = function(e) {
      warning(paste("Error in iteration", i, ":", e$message))
      -Inf
    })

    # Update best
    if (is.finite(score) && score > best_score) {
      best_score <- score
      best_params <- params
    }
  }

  if (is.null(best_params)) {
    stop("Random search failed: no valid parameter combinations found")
  }

  best_params
}

# ---------- Generic Random-Search-train orchestrator ----------
# Fits with Random Search on a holdout slice, then refits on all data.
# - fit_once(params, fit_idx, val_idx) -> list(score=..., model=..., extras=...)
#   During search you should return $score; during refit you can return $model.
# - refit(params) -> list(model=..., extras=...)

random_train <- function(X, y, bounds,
                         fit_once, refit,
                         val_frac = 0.1,
                         search_seed = 123,
                         n_iter = 15) {
  n <- nrow(X)
  idx <- make_split(n, val_frac)

  # 1) Define scoring fn for random search (maximize)
  score_fn <- function(par_list) {
    out <- fit_once(par_list, idx$fit_idx, idx$val_idx)
    # Expect a scalar numeric to maximize
    as.numeric(out$score)
  }

  # 2) Search
  best <- random_optimize(
    score_fn, bounds,
    n_iter = n_iter, seed = search_seed
  )

  # 3) Refit on all data
  final <- refit(as.list(best))

  list(model = final$model, params = best, extras = final$extras)
}

# =====================================================
# ===============  TORCH (nn) TRAINER  ================
# =====================================================

train_torch_model_random <- function(X, y,
                                     epochs = 100L,
                                     search_epochs = 30L,
                                     val_frac = 0.1,
                                     n_iter = 15,
                                     bounds = list(
                                       lr           = c(1e-3, 1e-2),
                                       hidden_dim   = c(32L, 128L),
                                       p_drop       = c(0.05, 0.3),
                                       weight_decay = c(1e-6, 1e-3)
                                     )) {
  set.seed(1234)
  device <- torch_get_device_deterministic(1234L)

  X <- as.matrix(X)
  y <- as.numeric(y)
  input_dim <- ncol(X)
  loss_fn <- torch::nn_mse_loss()

  fit_core <- function(Xt, yt, params, n_epochs) {
    model <- Net(input_dim, as.integer(round(params$hidden_dim)), p_drop = params$p_drop)$to(device = device)
    opt   <- torch::optim_adam(model$parameters, lr = params$lr, weight_decay = params$weight_decay)

    for (ep in seq_len(n_epochs)) {
      model$train()
      opt$zero_grad()
      pred <- model(Xt)
      loss <- loss_fn(pred, yt)
      loss$backward()
      opt$step()
    }
    model
  }

  fit_once <- function(params, fit_idx, val_idx) {
    X_fit <- torch::torch_tensor(X[fit_idx,, drop = FALSE], dtype = torch::torch_float(), device = device)
    y_fit <- torch::torch_tensor(y[fit_idx], dtype = torch::torch_float(), device = device)$unsqueeze(2)

    model <- fit_core(X_fit, y_fit, params, search_epochs)

    if (length(val_idx)) {
      X_val <- torch::torch_tensor(X[val_idx,, drop = FALSE], dtype = torch::torch_float(), device = device)
      y_val <- torch::torch_tensor(y[val_idx], dtype = torch::torch_float(), device = device)$unsqueeze(2)
      model$eval()
      score <- -as.numeric(loss_fn(model(X_val), y_val)$item())  # maximize negative MSE
      return(list(score = score, model = NULL))
    } else {
      return(list(score = NA_real_, model = model))
    }
  }

  refit <- function(params) {
    X_all <- torch::torch_tensor(X, dtype = torch::torch_float(), device = device)
    y_all <- torch::torch_tensor(y, dtype = torch::torch_float(), device = device)$unsqueeze(2)
    model <- fit_core(X_all, y_all, params, epochs)
    list(model = model, extras = list(device = device))
  }

  random_train(X, y, bounds, fit_once, refit, val_frac = val_frac, n_iter = n_iter)
}


# =====================================================
# ==============  XGBOOST (regression)  ===============
# =====================================================

train_xgboost_model_random <- function(X, y,
                                       val_frac = 0.1,
                                       n_iter = 15,
                                       bounds = list(
                                         max_depth        = c(2L, 8L),
                                         min_child_weight = c(1L, 10L),
                                         subsample        = c(0.5, 0.9)
                                       ),
                                       nrounds = 100,
                                       search_nrounds = 30,
                                       eta = 0.1,
                                       early_stopping_rounds = 10) {
  stopifnot(nrow(X) == length(y))
  X <- as.matrix(X); y <- as.numeric(y)

  fit_once <- function(params, fit_idx, val_idx) {
    dtrain <- xgboost::xgb.DMatrix(data = X[fit_idx,, drop = FALSE], label = y[fit_idx])
    params_xgb <- list(
      objective        = "reg:squarederror",
      eval_metric      = "rmse",
      eta              = eta,
      max_depth        = as.integer(params$max_depth),
      min_child_weight = as.integer(params$min_child_weight),
      subsample        = params$subsample,
      colsample_bytree = 0.8
    )

    if (length(val_idx)) {
      dvalid <- xgboost::xgb.DMatrix(data = X[val_idx,, drop = FALSE], label = y[val_idx])
      watch <- list(train = dtrain, valid = dvalid)
      bst <- xgboost::xgb.train(
        params  = params_xgb, data = dtrain, nrounds = search_nrounds,
        watchlist = watch, early_stopping_rounds = early_stopping_rounds, verbose = 0
      )
      preds <- predict(bst, X[val_idx,, drop = FALSE])
      score <- -smape(y[val_idx], preds)  # maximize
      list(score = score, model = NULL)
    } else {
      bst <- xgboost::xgb.train(params = params_xgb, data = dtrain, nrounds = search_nrounds, verbose = 0)
      list(score = NA_real_, model = bst)
    }
  }

  refit <- function(params) {
    dtrain <- xgboost::xgb.DMatrix(data = as.matrix(X), label = y)
    bst <- xgboost::xgb.train(
      params = list(
        objective        = "reg:squarederror",
        eval_metric      = "rmse",
        eta              = eta,
        max_depth        = as.integer(params$max_depth),
        min_child_weight = as.integer(params$min_child_weight),
        subsample        = params$subsample,
        colsample_bytree = 0.8
      ),
      data = dtrain, nrounds = nrounds, verbose = 0
    )
    list(model = bst, extras = list(
      best_iteration = bst$best_iteration %||% xgboost::xgb.attributes(bst)[["best_iteration"]]
    ))
  }

  random_train(X, y, bounds, fit_once, refit, val_frac = val_frac, n_iter = n_iter)
}


# =====================================================
# ==============  LIGHTGBM (regression)  ==============
# =====================================================

train_lightgbm_model_random <- function(X, y,
                                        val_frac = 0.1,
                                        n_iter = 15,
                                        bounds = list(
                                          max_depth        = c(4L, 10L),
                                          num_leaves       = c(8L, 255L),
                                          min_data_in_leaf = c(20L, 600L),
                                          feature_fraction = c(0.6, 1.0),
                                          bagging_fraction = c(0.6, 1.0)
                                        ),
                                        nrounds = 100,
                                        search_nrounds = 30,
                                        learning_rate = 0.05,
                                        early_stopping_rounds = 10) {
  stopifnot(nrow(X) == length(y))
  X <- as.matrix(X); y <- as.numeric(y)

  fit_once <- function(params, fit_idx, val_idx) {
    dtrain <- lightgbm::lgb.Dataset(data = X[fit_idx,, drop = FALSE], label = y[fit_idx])
    valids <- NULL
    prm <- list(
      objective = "regression",
      metric = "l2",
      learning_rate = learning_rate,
      max_depth = as.integer(params$max_depth),
      num_leaves = as.integer(params$num_leaves),
      min_data_in_leaf = as.integer(params$min_data_in_leaf),
      feature_fraction = params$feature_fraction,
      bagging_fraction = params$bagging_fraction
    )

    if (length(val_idx)) {
      dvalid <- lightgbm::lgb.Dataset(data = X[val_idx,, drop = FALSE], label = y[val_idx])
      valids <- list(train = dtrain, valid = dvalid)
      bst <- lightgbm::lgb.train(
        params = prm, data = dtrain, nrounds = search_nrounds,
        valids = valids, early_stopping_rounds = early_stopping_rounds, verbose = 0
      )
      preds <- predict(bst, X[val_idx,, drop = FALSE])
      score <- -smape(y[val_idx], preds)  # maximize
      list(score = score, model = NULL)
    } else {
      bst <- lightgbm::lgb.train(params = prm, data = dtrain, nrounds = search_nrounds, verbose = 0)
      list(score = NA_real_, model = bst)
    }
  }

  refit <- function(params) {
    dtrain <- lightgbm::lgb.Dataset(data = X, label = y)
    bst <- lightgbm::lgb.train(
      params = list(
        objective = "regression",
        metric = "l2",
        learning_rate = learning_rate,
        max_depth = as.integer(params$max_depth),
        num_leaves = as.integer(params$num_leaves),
        min_data_in_leaf = as.integer(params$min_data_in_leaf),
        feature_fraction = params$feature_fraction,
        bagging_fraction = params$bagging_fraction
      ),
      data = dtrain, nrounds = nrounds, verbose = 0
    )
    list(model = bst, extras = list(best_iter = bst$best_iter))
  }

  random_train(X, y, bounds, fit_once, refit, val_frac = val_frac, n_iter = n_iter)
}

