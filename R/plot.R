#' @method plot hybridForecast_model
#' @export
plot.hybridForecast_model <- function(x, x_lab = NULL, y_lab = NULL, ...) {

    cv_data <- x$cv_data
    freq<-detect_frequency(cv_data$ds)


    # Decide scale and format
    if (freq == "hour") {
      cv_data$ds <- as.POSIXct(cv_data$ds)
      scale_func <- scale_x_datetime
      xbreaks <- list(date_labels = "%H:%M", date_breaks = "1 hour")
      if (is.null(x_lab)) x_lab <- "Hour"
    } else if (freq %in% c("day", "week")) {
      cv_data$ds <- as.Date(cv_data$ds)
      scale_func <- scale_x_date
      xbreaks <- list(date_labels = "%m-%d", date_breaks = "1 day")
      if (is.null(x_lab)) x_lab <- "Date (mm-dd)"
    } else if (freq %in% c("month", "quarterly")) {
      cv_data$ds <- as.Date(cv_data$ds)
      scale_func <- scale_x_date
      xbreaks <- list(date_labels = "%m-%d", date_breaks = "1 month")
      if (is.null(x_lab)) x_lab <- "Date (mm-dd)"
    } else {  # year
      cv_data$ds <- as.Date(cv_data$ds)
      scale_func <- scale_x_date
      xbreaks <- list(date_labels = "%Y", date_breaks = "1 year")
      if (is.null(x_lab)) x_lab <- "Year"
    }

    if (is.null(y_lab)) y_lab <- "Value"

    if (!all(c("yhat_lower", "yhat_upper") %in% names(cv_data))) {
      stop("data must contain 'yhat_lower' and 'yhat_upper'.")
    }

    # Build plot
    p <- ggplot(cv_data, aes(x = ds)) +
      geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper, group = fold_label),
                  fill = "lightblue", alpha = 0.4) +
      geom_line(aes(y = y, group = fold_label), color = "black", size = 1) +
      geom_line(aes(y = yhat, group = fold_label), color = "blue", linetype = "dashed", size = 1) +
      facet_wrap(~fold_label, scales = "free_x", ncol = 2) +
      do.call(scale_func, xbreaks) +
      labs(
        title = paste0("Actual vs Predicted (Ensemble = ", x$emsemble,
                       ", MAPE = ", round(x$mape * 100, 1), "%)"),
        x = x_lab,
        y = y_lab,
        ...
      ) +
      theme_minimal() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    return(p)
}

#' @method plot hybridForecast
#' @export
plot.hybridForecast <- function(x, x_lab = NULL, y_lab = NULL, ...) {

  data_combined <-dplyr::bind_rows( tail(x$data,min(NROW(x$forecast),NROW(x$data))),   x$forecast)

  freq <- detect_frequency(data_combined$ds)
  # Decide scale and format
  if (freq == "hour") {
    data_combined$ds<-as.POSIXct(data_combined$ds)
    scale_func <- scale_x_datetime
    xbreaks <- list(date_labels = "%H:%M", date_breaks = "1 hour")
    if (is.null(x_lab)) x_lab <- "Hour"
  } else if (freq %in% c("day", "week")) {
    data_combined$ds <- as.Date(data_combined$ds)
    scale_func <- scale_x_date
    xbreaks <- list(date_labels = "%m-%d", date_breaks = "1 day")
    if (is.null(x_lab)) x_lab <- "Date (mm-dd)"
  } else if (freq %in% c("month", "quarterly")) {
    data_combined$ds <- as.Date(data_combined$ds)
    scale_func <- scale_x_date
    xbreaks <- list(date_labels = "%m-%d", date_breaks = "1 month")
    if (is.null(x_lab)) x_lab <- "Date (mm-dd)"
  } else {  # year
    data_combined$ds <- as.Date(data_combined$ds)
    scale_func <- scale_x_date
    xbreaks <- list(date_labels = "%Y", date_breaks = "1 year")
    if (is.null(x_lab)) x_lab <- "Year"
  }

  if (is.null(y_lab)) y_lab <- "Value"

  if (!all(c("yhat_lower", "yhat_upper") %in% names(data_combined))) {
    stop("data must contain 'yhat_lower' and 'yhat_upper'.")
  }

  # Build plot
  p <- ggplot(data_combined, aes(x = ds)) +
    geom_ribbon(aes(ymin = yhat_lower, ymax = yhat_upper),
                fill = "lightblue", alpha = 0.4) +
    geom_line(aes(y = y), color = "black", size = 1) +
    geom_line(aes(y = yhat), color = "blue", linetype = "dashed", size = 1) +
    do.call(scale_func, xbreaks) +
    labs(
      title = paste0("Forecasts from Ensemble (Ensemble = ", x$emsemble,
                     ", MAPE = ", round(x$mape * 100, 1), "%)"),
      x = x_lab,
      y = y_lab,
      ...
    ) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  return(p)
}

