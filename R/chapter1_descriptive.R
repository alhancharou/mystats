#' Sample Variance dividing by N (s_n^2) and Standard Deviation (s_n)
#' @param x Numeric vector of values (or class marks).
#' @param n_i Numeric vector of frequencies (default is 1 for each value).
#' @examples
#' pop_var(c(12, 10, 9, 9, 10))
#' pop_var(c(0, 1), c(0.3, 0.7))
#' @export
pop_var <- function(x, n_i = rep(1, length(x))) {
  f <- n_i / sum(n_i)
  m <- sum(x * f)
  sum(f * (x - m)^2)
}

#' @rdname pop_var
#' @export
pop_sd <- function(x, n_i = rep(1, length(x))) {
  sqrt(pop_var(x, n_i))
}

#' Weighted / Grouped Mean
#' @param x Numeric vector of values or class marks.
#' @param n_i Numeric vector of frequencies or weights.
#' @examples
#' grouped_mean(c(2.6, 3.7, 5.1, 4.9, 6.4), c(1, 1, 1, 2, 3))
#' @export
grouped_mean <- function(x, n_i = rep(1, length(x))) {
  sum(x * n_i) / sum(n_i)
}

#' Quadratic Mean / Root Mean Square (RMS) (Slide 36)
#' @param x Numeric vector of values.
#' @param n_i Numeric vector of frequencies.
#' @examples
#' rms_val(c(2, 3, 1, 0, 4, 3))
#' @export
rms_val <- function(x, n_i = rep(1, length(x))) {
  sqrt(sum(n_i * (x^2)) / sum(n_i))
}

#' Mean Deviation MD(p) and Mean Squared Error MSE(p) (Slides 49-50)
#' @param x Numeric vector of values.
#' @param p Reference point (e.g., mean or median).
#' @param n_i Numeric vector of frequencies.
#' @examples
#' x <- c(5, 2, 3, 3, 3, 5, 7)
#' mean_dev(x, mean(x))
#' mse_point(x, median(x))
#' @export
mean_dev <- function(x, p, n_i = rep(1, length(x))) {
  sum(abs(x - p) * n_i) / sum(n_i)
}

#' @rdname mean_dev
#' @export
mse_point <- function(x, p, n_i = rep(1, length(x))) {
  sum(((x - p)^2) * n_i) / sum(n_i)
}

#' Pearson's Coefficient of Variation CV = sd / |mean| (Slide 61)
#' @param x Numeric vector of values.
#' @param n_i Numeric vector of frequencies.
#' @examples
#' cv_pearson(c(12, 10, 9, 9, 10))
#' @export
cv_pearson <- function(x, n_i = rep(1, length(x))) {
  m <- grouped_mean(x, n_i)
  s <- pop_sd(x, n_i)
  s / abs(m)
}

#' 1D Moment of order r with respect to point c (Slides 63-68)
#' For ordinary moment m_r: set c = 0. For central moment mu_r: set c = mean.
#' @param x Numeric vector of values.
#' @param r Order of moment (1, 2, 3, 4...).
#' @param c Point of reference (default 0).
#' @param n_i Numeric vector of frequencies.
#' @examples
#' x <- c(0, 0, 2, 2, 2, 2, 3, 3, 4, 4, 7, 7)
#' moment_1d(x, r = 2, c = median(x))
#' @export
moment_1d <- function(x, r = 1, c = 0, n_i = rep(1, length(x))) {
  sum(((x - c)^r) * n_i) / sum(n_i)
}

#' Skewness (Pearson A_p and Fisher g_1) and Fisher Kurtosis g_2 (Slides 70-72)
#' @param x Numeric vector of values.
#' @param n_i Numeric vector of frequencies.
#' @examples
#' skew_fisher(c(2, 5, 3, 4, 7, 0, 11, 2, 3, 8))
#' kurt_fisher(c(2, 5, 3, 4, 7, 0, 11, 2, 3, 8))
#' @export
skew_fisher <- function(x, n_i = rep(1, length(x))) {
  m <- grouped_mean(x, n_i)
  s <- pop_sd(x, n_i)
  mu3 <- moment_1d(x, r = 3, c = m, n_i = n_i)
  mu3 / (s^3)
}

#' @rdname skew_fisher
#' @export
kurt_fisher <- function(x, n_i = rep(1, length(x))) {
  m <- grouped_mean(x, n_i)
  s <- pop_sd(x, n_i)
  mu4 <- moment_1d(x, r = 4, c = m, n_i = n_i)
  (mu4 / (s^4)) - 3
}

#' @rdname skew_fisher
#' @export
skew_pearson <- function(x, mode_val, n_i = rep(1, length(x))) {
  m <- grouped_mean(x, n_i)
  s <- pop_sd(x, n_i)
  (m - mode_val) / s
}

#' Quantile for Continuous Grouped Data in Intervals (Slide 44)
#' Formula: L_{i-1} + ((c*N - N_{i-1}) / n_i) * a_i
#' @param L_low Vector of lower bounds of intervals (L_{i-1}).
#' @param L_high Vector of upper bounds of intervals (L_i).
#' @param n_i Vector of absolute frequencies.
#' @param c_prob Probability in (0, 1), e.g., 0.5 for Median, 0.25 for Q1.
#' @examples
#' quantile_grouped(c(0, 3, 6, 9, 12), c(3, 6, 9, 12, 15), c(40, 30, 10, 5, 5), 0.5)
#' @export
quantile_grouped <- function(L_low, L_high, n_i, c_prob = 0.5) {
  N <- sum(n_i)
  N_cum <- cumsum(n_i)
  idx <- which(N_cum >= c_prob * N)[1]
  N_prev <- ifelse(idx == 1, 0, N_cum[idx - 1])
  a_i <- L_high[idx] - L_low[idx]
  L_low[idx] + ((c_prob * N - N_prev) / n_i[idx]) * a_i
}