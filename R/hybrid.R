
predict_one_step <- function(hybrid_model,newdata_feat,xreg) {

  n_test <- nrow(newdata_feat)
  preds <- list()
  use_lags <- length(hybrid_model$significant_lags) > 0

  for (model in model_used)
  {

    preds[[model]]<-numeric(n_test)

    if (use_lags) {
      max_lags <- max(hybrid_model$significant_lags)
      lag_names <- paste0("lag", hybrid_model$significant_lags)
      recent_y <- tail(hybrid_model$train_feat$y, max_lags)
    }

    for (i in 1:n_test) {

      # Prepare lag features
      if (use_lags) {
        lag_values <- sapply(hybrid_model$significant_lags, function(lag) {
          recent_y[length(recent_y) - lag + 1]
        })
        lag_data<- as.data.frame(as.list(lag_values))
        names(lag_data) <- lag_names
        input_row <- cbind(newdata_feat[i, c("trend", names(hybrid_model$seasonal_periods$fouriers))], lag_data)
      } else {
        input_row <- newdata_feat[i, c("trend", names(hybrid_model$seasonal_periods$fouriers))]
      }

      if (!is.na(xreg)) input_row<-cbind(newdata_feat[i,xreg, drop=FALSE],input_row)

      if (model=="dl"){
        preds[[model]][i] <- predict_torch_model(hybrid_model[[model]],  as.matrix(input_row))
      } else {
        preds[[model]][i] <- predict(hybrid_model[[model]], as.matrix(input_row))
      }

      # Update recent_y for next iteration
      if (use_lags) {
        recent_y <- tail(c(recent_y, preds[[model]][i]), max_lags)
      }
    }
  }

  return(preds)
}




predict_hybrid_1<-function(hybrid_model,newdata,xreg)
{

  # prepare test data with prophet predict on trend and season

  res <- predict_trend_seasonal_components(hybrid_model$ts_model, h = NROW(newdata))
  newdata_feat <- cbind(newdata, trend = res$trend)

  # Add each seasonal component
  for (name in names(res$seasonality)) {
    newdata_feat[[name]] <- res$seasonality[[name]]
  }



  # then step by step predict the test

  yhat<-predict_one_step(hybrid_model,newdata_feat,xreg)

  # Initialize containers
  preds <-list()
  train <- list()
  residuals <- list()
  sigma <- list()
  yhat_lower <- list()
  yhat_upper <- list()

  # Forecast horizon and confidence interval setup
  h <- NROW(newdata_feat)
  ci_scale <- sqrt(1:h)
  z <- qnorm(0.975)  # 95% CI

  # Models to process

  # Generate predictions and compute residuals
  for (model in model_used) {

    #  train[[model]] <- predict(hybrid_model[[model]], xgb.DMatrix(as.matrix(hybrid_model$train_feat[, hybrid_model$xgb_model$feature_names])))
    if (model=='dl')
    {
      train[[model]] <- predict_torch_model(hybrid_model[[model]], as.matrix(hybrid_model$train_feat[,!names(hybrid_model$train_feat) %in% c('ds','y')]))
    }else{
      train[[model]] <- predict(hybrid_model[[model]], as.matrix(hybrid_model$train_feat[,!names(hybrid_model$train_feat) %in% c('ds','y')]))
    }

    if (hybrid_model$log_transform)
    {
      residuals[[model]] <- exp(hybrid_model$train_feat$y) - exp(train[[model]])
      yhat[[model]]  <- exp(yhat[[model]]) - EPSILON
    }else{
      residuals[[model]] <- hybrid_model$train_feat$y - train[[model]]
    }
    sigma[[model]] <- sd(residuals[[model]], na.rm = TRUE)
    yhat_lower[[model]] <- yhat[[model]] - z * sigma[[model]] * ci_scale
    yhat_upper[[model]] <- yhat[[model]] + z * sigma[[model]] * ci_scale

    preds[[model]]$yhat<-pmax(0,yhat[[model]])
    preds[[model]]$yhat_lower<-pmax(0,yhat_lower[[model]])
    preds[[model]]$yhat_upper<-pmax(0,yhat_upper[[model]])


  }
  return(preds)

}

#' Predict method for hybridForecast_model objects
#'
#' Generates ensemble forecasts from a fitted \code{hybridForecast_model}
#' on new data, combining base model predictions using the model's
#' stored ensemble weights.
#'
#' @method predict hybridForecast_model
#'
#' @param object A fitted \code{hybridForecast_model} returned by \code{\link{hybrid}}.
#' @param newdata A \code{data.frame} of future periods to forecast, including
#'   a \code{ds} column (and any required regressors if applicable).
#' @param xreg Optional exogenous regressors (\code{vector of column names in newdata}) or \code{NA} (default) if none.
#' @param ... Additional arguments (currently unused).
#'
#' @return An object of class \code{"hybridForecast"} containing:
#' \itemize{
#'   \item \code{data} – original training data from the fitted model
#'   \item \code{emsemble} – name of the selected model combination
#'   \item \code{mape} – mean sMAPE from cross-validation
#'   \item \code{forecast} – \code{newdata} with appended columns:
#'     \code{yhat}, \code{yhat_lower}, \code{yhat_upper}
#' }
#'
#'
#' @examples
#' \dontrun{
#' fit <- hybrid(my_data)
#' future_df <- make_future_data(fit, horizon = 14)
#' fc <- predict(fit, newdata = future_df)
#' head(fc$forecast)
#' }
#'
#' @export
predict.hybridForecast_model<-function(object,newdata,xreg=NA, ...)
{
  test_result<-predict_hybrid_1(object,newdata,xreg)
  test_matrix<-do.call(rbind,test_result)
  yhat<-do.call(rbind,test_matrix[,'yhat'])
  yhat_lower<-do.call(rbind,test_matrix[,'yhat_lower'])
  yhat_upper<-do.call(rbind,test_matrix[,'yhat_upper'])
  yhat_ensemble<-object$ensemble_weight %*% yhat
  yhat_ensemble_lower<-object$ensemble_weight %*% yhat_lower
  yhat_ensemble_upper<-object$ensemble_weight %*% yhat_upper

  hybrid_forecast<-list()
  hybrid_forecast$data<-object$data
  hybrid_forecast$emsemble<-object$emsemble
  hybrid_forecast$mape<-object$mape
  newdata$yhat<-t(yhat_ensemble)
  newdata$yhat_lower<-t(yhat_ensemble_lower)
  newdata$yhat_upper<-t(yhat_ensemble_upper)
  hybrid_forecast$forecast<-newdata
  class(hybrid_forecast)<-'hybridForecast'
  return(hybrid_forecast)
}


hybrid_core<-function(train, xreg, seasonal_periods,n_changepoints,max_lag,log_transform)
{
  if(log_transform) train$y<-log(train$y+EPSILON)
  if (is.na(xreg)) train<-train[,c('ds','y')]
  significant_lags<-find_lags(train,max_lag)
  #prepare features

  feat_model<-prepare_trend_seasonal_component(train=train,significant_lags=significant_lags,seasonal_periods=seasonal_periods,n_changepoints=n_changepoints)

  train_feat<-feat_model$train_feat
  features <- setdiff(colnames(train_feat),c("ds","y"))
  #if (!is.na(xreg)) features<-c(features,xreg)

  X_train<-train_feat[, features]
  y_train <- train_feat$y
  ## train xgboost model on train
  dtrain <- xgb.DMatrix(data = as.matrix(X_train), label = y_train)
  dtrain_lgb <- lgb.Dataset(data = as.matrix(X_train), label = y_train)
  set.seed(123)
  xgb_model <- xgb.train(
    params = list(
      objective = "reg:squarederror",
      max_depth = 6,
      eta = 0.05,
      subsample = 0.8,
      colsample_bytree = 0.8,
      min_child_weight = 5
    ),
    data = dtrain,
    nrounds = 100,
  )
  if (FALSE)
{
  formula_str <- paste("y ~", paste(features, collapse = " + "))
  dl_model <- nnet(
    formula =as.formula(formula_str),
    data=train_feat,
    size= max(10, floor(sqrt(NROW(train_feat) / (length(features) + 1)))),
    decay = 0.1,
    MaxNWts = 5000,
    linout = TRUE
  )
  }

  dl_model <- train_torch_model(X_train, y_train)


  #lm_model<-lm(as.formula(formula_str),data=train_feat)
  lgb_model<-lgb.train(
    params =
      list(
        objective = "regression",
        metric = "l2",
        learning_rate = 0.05,
        num_leaves = 31,
        max_depth=-1,
        min_data_in_leaf = 20,
        feature_fraction = 0.9,
        bagging_fraction = 0.8,
        bagging_freq = 1
      ),
    data = dtrain_lgb,
    nrounds = 100,
    verbose = -1
  )

  return(list(xgb=xgb_model,lgb=lgb_model,dl=dl_model,train_feat=train_feat,significant_lags=significant_lags,seasonal_periods=seasonal_periods,ts_model=feat_model$ts_model,log_transform=log_transform))
}

#' Hybrid Forecasting Model
#'
#' Fits multiple base forecasting models and combines them using
#' cross-validated inverse-sMAPE weighting to form an ensemble.
#'
#' @param data Data frame with columns \code{ds} (date/time) and \code{y} (numeric target). Other columns will be used as covariates.
#' @param xreg Optional exogenous regressors (\code{vector of column names in data} or \code{NA}).
#' @param freq Seasonal frequency: hour, day, week, month, quarter or year. If \code{NA}, inferred automatically.
#' @param n_changepoints Number of changepoints for trend. If \code{NA}, inferred automatically.
#' @param max_lag Max lag for autoregressive features. If \code{NA}, chosen automatically.
#' @param horizon Forecast horizon for cross-validation folds. If \code{NA}, 20% of data will be used for testing.
#' @param max_fold Number of CV folds (default 5).
#' @param log_transform Logical (default FALSE); whether to log-transform \code{y} before modelling.
#'
#' @return An object of class \code{"hybridForecast_model"} containing:
#' \itemize{
#'   \item \code{data} – training data
#'   \item \code{cv_data} – fold-level forecasts and errors
#'   \item \code{mape} – mean sMAPE of selected ensemble
#'   \item \code{ensemble_weight} – weights for each base model
#'   \item \code{emsemble} – selected model combination name
#' }
#'
#' @details Uses rolling-origin CV, tests all non-empty combinations of
#' \code{model_used}, selects the one with lowest mean sMAPE, and refits on full data.
#'
#' @export
hybrid <- function(data,
                   xreg=NA,
                   freq = NA,
                   n_changepoints = NA,
                   max_lag = NA,
                   horizon = NA,
                   max_fold = 5,
                   log_transform=FALSE
                   ) {

  data<-data%>%arrange(ds)
  N <- nrow(data)

  if (is.na(freq)) {
    seasonal_periods <- auto_seasonal_periods(data$ds)

  }else{
    seasonal_periods<-seasonal_priods_fourier(freq)
  }


  if (is.na(max_lag)) {
    max_lag <- auto_max_lag(N, seasonal_periods)
    cat('\nnum of lags:', max_lag)
  }

  model_choice <- rep(list(0:1), length(model_used))
  names(model_choice) <- model_used
  combinations <- expand.grid(model_choice) %>% tail(-1) %>% as.matrix()

  folds <- vector("list", max_fold)
  if (is.na(horizon) || (N - horizon<30) ) horizon<- min(80,floor(0.2*N))
  for (fold in seq_len(max_fold)) {
    cat("\nfold", fold)

    train_end <- N - horizon - fold + 1
    test_range <- (train_end + 1):(train_end + horizon)
    train_data <- data[1:train_end, ]
    test_data <- data[test_range, ]

    hybrid_model <- hybrid_core(train_data,xreg, seasonal_periods,n_changepoints, max_lag,log_transform)

    test_result <- predict_hybrid_1(hybrid_model, test_data,xreg)

    y <- test_data$y
    ds <- test_data$ds

    for (model in model_used) {
        test_result[[model]]$mape <- smape(y, test_result[[model]]$yhat)
    }

    # Combine forecasts
    test_matrix <- do.call(rbind, test_result)
    mape <- do.call(rbind, test_matrix[,"mape"])
    yhat <- do.call(rbind, test_matrix[,"yhat"])
    yhat_lower <- do.call(rbind, test_matrix[,"yhat_lower"])
    yhat_upper <- do.call(rbind, test_matrix[,"yhat_upper"])
    weights_matrix <- combinations
    for (i in seq_len(nrow(weights_matrix))) {
      idx <- which(weights_matrix[i, ] > 0)
      inv_mape <- 1 / mape[idx]
      inv_mape[!is.finite(inv_mape)] <- 1
      weights_matrix[i, idx] <- inv_mape
      weights_matrix[i, ] <- weights_matrix[i, ] / sum(weights_matrix[i, ])
      rownames(weights_matrix)[i] <- paste(model_used[idx], collapse = "|")
    }

    yhat_ensemble <- weights_matrix %*% yhat
    yhat_ensemble_lower <- weights_matrix %*% yhat_lower
    yhat_ensemble_upper <- weights_matrix %*% yhat_upper

    mape_by_method<-apply(yhat_ensemble,1,function(fcst){ smape(fcst,y)})

    folds[[fold]] <- list(
      y = y, ds = ds,
      yhat_ensemble = yhat_ensemble,
      yhat_ensemble_lower = yhat_ensemble_lower,
      yhat_ensemble_upper = yhat_ensemble_upper,
      mape_by_method = mape_by_method,
      combinations = weights_matrix
    )
  }

  # Combine folds
  mape_matrix <- do.call(cbind, lapply(folds, `[[`, "mape_by_method"))
  mean_mape <- rowMeans(mape_matrix)
  best_idx <- which.min(mean_mape)
  best_method <- names(mean_mape)[best_idx]
  best_mape <- mean_mape[best_idx]

  all_combinations <- do.call(rbind, lapply(folds, `[[`, "combinations"))
  best_weights <- colMeans(all_combinations[rownames(all_combinations) == best_method, , drop = FALSE])

  yhat_all <- do.call(cbind, lapply(folds, `[[`, "yhat_ensemble"))
  yhat_all_lower <- do.call(cbind, lapply(folds, `[[`, "yhat_ensemble_lower"))
  yhat_all_upper <- do.call(cbind, lapply(folds, `[[`, "yhat_ensemble_upper"))

  ds_all <- do.call(c,lapply(folds, `[[`, "ds"))
  y_all <- unlist(lapply(folds, `[[`, "y"))
  mape_all <- rep(mape_matrix[best_method, ], each = horizon)

  cv_data <- data.frame(
    fold = rep(1:max_fold, each = horizon),
    ds = ds_all,
    y = y_all,
    yhat = yhat_all[rownames(yhat_all) == best_method, ],
    yhat_lower = yhat_all_lower[rownames(yhat_all_lower) == best_method, ],
    yhat_upper = yhat_all_upper[rownames(yhat_all_upper) == best_method, ],
    mape = mape_all
  ) %>% group_by(fold) %>% mutate(fold_label = paste0("Fold ", fold, " (MAPE = ", round(mape[1] * 100, 1), "%)"))

  # Refit on full data

  final_model <- hybrid_core(data, xreg, seasonal_periods,n_changepoints, max_lag,log_transform)

  final_model$data<-data
  final_model$cv_data <- cv_data
  final_model$mape <- best_mape
  final_model$ensemble_weight <- best_weights
  final_model$emsemble <- best_method
  final_model$seasonal_periods <- seasonal_periods

  class(final_model)<-'hybridForecast_model'
  return(final_model)
}
