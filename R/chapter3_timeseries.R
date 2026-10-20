#' Centred Moving Average of order k (Slides 11-21)
#' Odd k uses equal weights 1/k. Even k uses 0.5 at both extremes!
#' @param x Numeric vector or time series.
#' @param k Order of moving average (e.g. 3, 4, 12).
#' @examples
#' xx <- c(6, 6.2, 5.1, 4.9, 5.2, 5.8, 7)
#' moving_avg(xx, 3)
#' moving_avg(xx, 4)
#' @export
moving_avg <- function(x, k) {
  if (k %% 2 == 1) {
    weights <- rep(1, k) / k
  } else {
    weights <- c(0.5, rep(1, k - 1), 0.5) / k
  }
  as.numeric(stats::filter(x, weights))
}

#' Complete Time Series Decomposition, Deseasonalization & Trend Line (Slides 24-55)
#' Follows exact professor's R code from Chapter 3!
#' @param xx Numeric vector of observations.
#' @param freq Seasonal period (4 for quarterly, 12 for monthly).
#' @param start_yr Start vector c(year, subperiod), e.g., c(1980, 1).
#' @param type "multiplicative" (X = T*E*C*A) or "additive" (X = T+E+C+A).
#' @param predict_t Optional vector of future time indices t to predict X.
#' @examples
#' xx <- c(178.2, 153.2, 185.9, 163.6, 196.3, 156.9, 197.9, 166.4,
#'         197.3, 159.7, 202.6, 175.6, 209.5, 169.5, 202.4, 179.8,
#'         200.0, 168.6, 216.1, 178.5)
#' ts_full_analysis(xx, freq = 4, start_yr = c(1980, 1))
#' @export
ts_full_analysis <- function(xx, freq = 4, start_yr = c(1, 1), type = "multiplicative", predict_t = NULL) {
  uk <- stats::ts(xx, start = start_yr, frequency = freq)
  
  # 1. Trend-Cycle (TC) via moving average of order = freq (Slide 30)
  if (freq %% 2 == 1) {
    w <- rep(1, freq) / freq
  } else {
    w <- c(0.5, rep(1, freq - 1), 0.5) / freq
  }
  tend <- stats::filter(uk, w)
  
  if (type == "multiplicative") {
    # 2. Eliminate Trend-Cycle: E * A = X / TC (Slide 33)
    est_ale <- uk / tend
    # 3. Unnormalized Seasonal Indices (Slide 38)
    est_raw <- colMeans(matrix(est_ale, ncol = freq, byrow = TRUE), na.rm = TRUE)
    # 4. Normalized / Corrected Seasonal Indices (mean == 1) (Slide 40)
    est_norm <- est_raw / mean(est_raw)
    # 5. Random component: A = (E*A) / E (Slide 42)
    ale <- est_ale / est_norm
    # 6. Deseasonalized series: X / E (Slide 46)
    desest <- uk / est_norm
  } else {
    # Additive case: subtract instead of divide (Slide 32)
    est_ale <- uk - tend
    est_raw <- colMeans(matrix(est_ale, ncol = freq, byrow = TRUE), na.rm = TRUE)
    est_norm <- est_raw - mean(est_raw)
    ale <- est_ale - est_norm
    desest <- uk - est_norm
  }
  
  # 7. Linear Regression on Deseasonalized series: True Trend T = a + b*t (Slide 53)
  t_idx <- seq_along(uk)
  model <- stats::lm(as.numeric(desest) ~ t_idx)
  cf <- stats::coef(model) # Intercept a = cf[1], Slope b = cf[2]
  true_trend <- stats::ts(stats::predict.lm(model), start = start_yr, frequency = freq)
  
  # 8. Optional predictions for future t: Trend(t) * Seasonal_Index
  preds <- NULL
  if (!is.null(predict_t)) {
    trend_pred <- cf[1] + cf[2] * predict_t
    season_pos <- ((predict_t - 1) %% freq) + 1
    if (type == "multiplicative") {
      preds <- trend_pred * est_norm[season_pos]
    } else {
      preds <- trend_pred + est_norm[season_pos]
    }
  }
  
  list(
    ts_data = uk,
    trend_cycle_MA = tend,
    seasonal_random_EA = est_ale,
    seasonal_uncorrected = est_raw,
    seasonal_corrected = est_norm,
    random_A = ale,
    deseasonalized = desest,
    regression_a = unname(cf[1]),
    regression_b = unname(cf[2]),
    true_trend_T = true_trend,
    predictions = preds
  )
}

#' Autocorrelation of order n to detect cycles (Slide 61)
#' @param data Numeric vector of time series.
#' @param n Lag / period to check (e.g., n = 3).
#' @examples
#' x <- c(1.2, 1.4, 2.0, 1.23, 1.37, 2.07, 1.25, 1.47, 2.09)
#' autoCorrel(x, 3)
#' @export
autoCorrel <- function(data, n) {
  xx1 <- utils::head(data, -n)
  xx2 <- utils::tail(data, -n)
  cor(xx1, xx2)
}

#' Simple Index Numbers and Deflation via CPI (Slides 63-72)
#' @param x Numeric vector of values over time.
#' @param base_pos Index of the base period (default 1).
#' @examples
#' simple_index(c(1254, 1345, 1408, 1451), base_pos = 1)
#' @export
simple_index <- function(x, base_pos = 1) {
  x / x[base_pos]
}

#' Deflate Nominal Value to Real Value using CPI (Slides 68-72)
#' Real Value at base t0 = Nominal Value at t / (CPI_{t/t0} in per-unit)
#' @param nominal Nominal monetary value at time t.
#' @param cpi_per_unit CPI of time t with respect to base t0 (if 135%, pass 1.35!).
#' @examples
#' # Product cost 1000 in 1998, CPI_{98/90} = 135% (1.35) -> real value in 1990:
#' deflate_val(1000, 1.35)
#' @export
deflate_val <- function(nominal, cpi_per_unit) {
  nominal / cpi_per_unit
}