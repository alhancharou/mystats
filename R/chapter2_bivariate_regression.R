# --- PART 1: BIVARIATE TABLES (Slides 1-47) ---

#' Mode and Median from probabilities/frequencies
#' @param vals Numeric vector of values.
#' @param probs Numeric vector of probabilities or frequencies.
#' @examples
#' prob_mode(c(0, 1), c(0.3, 0.7))
#' prob_median(c(0, 1), c(0.3, 0.7))
#' @export
prob_mode <- function(vals, probs) {
  p <- probs / sum(probs)
  vals[which.max(p)]
}

#' @rdname prob_mode
#' @export
prob_median <- function(vals, probs) {
  p <- probs / sum(probs)
  vals[which(cumsum(p) >= 0.5)[1]]
}

#' Check independence of two variables in a joint table (Slide 40)
#' @param tab Matrix or table of joint frequencies.
#' @examples
#' tab <- matrix(c(0.03, 0.07, 0.06, 0.14, 0.09, 0.21), nrow = 3, byrow = TRUE)
#' check_indep(tab)
#' @export
check_indep <- function(tab) {
  tab_correct <- tab / sum(tab)
  expected <- outer(rowSums(tab_correct), colSums(tab_correct))
  isTRUE(all.equal(as.matrix(tab_correct), as.matrix(expected), check.attributes = FALSE))
}

#' Covariance and Correlation from a joint frequency table (Slides 29-30, 80)
#' @param tab Matrix or table of joint frequencies.
#' @param row_vals Numeric vector of values for ROW variable.
#' @param col_vals Numeric vector of values for COLUMN variable.
#' @examples
#' tab <- matrix(c(3, 4, 3, 1, 0, 1, 2, 2, 2), nrow = 3, byrow = TRUE)
#' table_cov(tab, row_vals = c(-2, 0, 2), col_vals = c(-1, 0, 1))
#' @export
table_cov <- function(tab, row_vals, col_vals) {
  tab_correct <- tab / sum(tab)
  mean_r <- sum(row_vals * rowSums(tab_correct))
  mean_c <- sum(col_vals * colSums(tab_correct))
  mean_rc <- sum(outer(row_vals, col_vals) * tab_correct)
  mean_rc - (mean_r * mean_c)
}

#' @rdname table_cov
#' @export
table_cor <- function(tab, row_vals, col_vals) {
  cov_val <- table_cov(tab, row_vals, col_vals)
  sd_r <- pop_sd(row_vals, rowSums(tab))
  sd_c <- pop_sd(col_vals, colSums(tab))
  cov_val / (sd_r * sd_c)
}

#' General Bivariate Moment M_{r,s}(a, b) (Slide 25)
#' @param tab Matrix or table of joint frequencies.
#' @param row_vals Numeric vector of ROW values.
#' @param col_vals Numeric vector of COLUMN values.
#' @param r,s Powers for ROW and COL variables.
#' @param a,b Shift points for ROW and COL variables (default 0).
#' @export
table_moment <- function(tab, row_vals, col_vals, r = 1, s = 1, a = 0, b = 0) {
  tab_correct <- tab / sum(tab)
  term_r <- (row_vals - a)^r
  term_c <- (col_vals - b)^s
  sum(outer(term_r, term_c) * tab_correct)
}


# --- PART 2: REGRESSION, R^2, MSE & TRANSFORMATIONS (Slides 56-111) ---

#' Goodness of Fit: SSE, MSE, Variance Explained, and R^2 (Slides 63-66)
#' Works for ANY model! Just pass actual y and predicted y_est.
#' @param y Actual observed values of Y.
#' @param y_est Predicted (estimated) values y_est = f(x).
#' @examples
#' # After fitting any model: fit_metrics(y, y_est)
#' @export
fit_metrics <- function(y, y_est) {
  e <- y - y_est
  sse <- sum(e^2)
  mse <- mean(e^2)        # divides by n!
  var_y <- pop_var(y)     # s_n^2(Y) dividing by n!
  r2 <- 1 - (mse / var_y)
  v_ex <- r2 * var_y
  list(SSE = sse, MSE = mse, RMSE = sqrt(mse), Var_Y = var_y, Var_Explained = v_ex, R2 = r2)
}

#' Both Linear Regression Lines: Y on X (y = a + bx) and X on Y (x = a' + b'y) (Slide 74)
#' Works for both raw vectors (when tab = NULL) and 2D frequency tables!
#' @param x_vals Vector of X values.
#' @param y_vals Vector of Y values.
#' @param tab Optional: joint frequency matrix (assuming X is Rows, Y is Cols).
#' @examples
#' x <- c(0, 1, 2, 3, 6); y <- c(7, 5, 4, 3.5, 3)
#' reg_lines(x, y)
#' @export
reg_lines <- function(x_vals, y_vals, tab = NULL) {
  if (is.null(tab)) {
    mx <- mean(x_vals); my <- mean(y_vals)
    vx <- pop_var(x_vals); vy <- pop_var(y_vals)
    cov_xy <- mean(x_vals * y_vals) - mx * my
  } else {
    f_tab <- tab / sum(tab)
    px <- rowSums(f_tab); py <- colSums(f_tab)
    mx <- sum(x_vals * px); my <- sum(y_vals * py)
    vx <- pop_var(x_vals, px); vy <- pop_var(y_vals, py)
    cov_xy <- table_cov(f_tab, x_vals, y_vals)
  }
  # Line Y/X: y = a + b*x
  b <- cov_xy / vx
  a <- my - b * mx
  # Line X/Y: x = a_prime + b_prime*y
  b_prime <- cov_xy / vy
  a_prime <- mx - b_prime * my
  r <- cov_xy / sqrt(vx * vy)
  mse_yx <- vy * (1 - r^2)
  
  list(
    mean_x = mx, mean_y = my, var_x = vx, var_y = vy, cov = cov_xy,
    r = r, R2 = r^2, MSE_Y_on_X = mse_yx, RMSE_Y_on_X = sqrt(mse_yx),
    Y_on_X_a = a, Y_on_X_b = b,
    X_on_Y_a_prime = a_prime, X_on_Y_b_prime = b_prime
  )
}

#' Nonlinear Curve Fitting via Transformations (Slides 103-110)
#' Supports: "exp" (y = a*e^(bx)), "hyp" (y = 1/(a+bx)),
#' "power" (y = a*x^b), "abx" (y = a*b^x), "quad" (y = a + bx + cx^2)
#' @param x Numeric vector X.
#' @param y Numeric vector Y.
#' @param type Model type: "linear", "exp", "hyp", "power", "abx", "quad".
#' @param new_x Optional vector of X values to predict Y.
#' @examples
#' x <- c(0, 1, 2, 3, 6); y <- c(7, 5, 4, 3.5, 3)
#' fit_curve(x, y, type = "exp")
#' fit_curve(x, y, type = "hyp")
#' @export
fit_curve <- function(x, y, type = "linear", new_x = x) {
  if (type == "linear") {
    mod <- lm(y ~ x)
    cf <- coef(mod) # a = cf[1], b = cf[2]
    y_est <- cf[1] + cf[2] * x
    pred  <- cf[1] + cf[2] * new_x
    params <- c(a = unname(cf[1]), b = unname(cf[2]))
  } else if (type == "exp") {
    # y = a * exp(b*x) => ln(y) = A + b*x, where a = exp(A) (Slide 103)
    mod <- lm(log(y) ~ x)
    cf <- coef(mod)
    a <- exp(unname(cf[1])); b <- unname(cf[2])
    y_est <- a * exp(b * x)
    pred  <- a * exp(b * new_x)
    params <- c(A_ln_a = unname(cf[1]), a = a, b = b)
  } else if (type == "hyp") {
    # y = 1 / (a + b*x) => 1/y = a + b*x (Slide 105)
    mod <- lm(I(1 / y) ~ x)
    cf <- coef(mod)
    a <- unname(cf[1]); b <- unname(cf[2])
    y_est <- 1 / (a + b * x)
    pred  <- 1 / (a + b * new_x)
    params <- c(a = a, b = b)
  } else if (type == "power") {
    # y = a * x^b => ln(y) = ln(a) + b*ln(x) (Slide 109)
    mod <- lm(log(y) ~ log(x))
    cf <- coef(mod)
    a <- exp(unname(cf[1])); b <- unname(cf[2])
    y_est <- a * (x^b)
    pred  <- a * (new_x^b)
    params <- c(a = a, b = b)
  } else if (type == "abx") {
    # y = a * b^x => ln(y) = ln(a) + x*ln(b) (Slide 110)
    mod <- lm(log(y) ~ x)
    cf <- coef(mod)
    a <- exp(unname(cf[1])); b <- exp(unname(cf[2]))
    y_est <- a * (b^x)
    pred  <- a * (b^new_x)
    params <- c(a = a, b = b)
  } else if (type == "quad") {
    # y = a + b*x + c*x^2 (Slide 101, 111)
    mod <- lm(y ~ x + I(x^2))
    cf <- coef(mod)
    y_est <- cf[1] + cf[2] * x + cf[3] * (x^2)
    pred  <- cf[1] + cf[2] * new_x + cf[3] * (new_x^2)
    params <- c(a = unname(cf[1]), b = unname(cf[2]), c = unname(cf[3]))
  }
  metrics <- fit_metrics(y, y_est)
  list(type = type, parameters = params, y_est = y_est, predictions = pred, metrics = metrics)
}

#' Multiple Linear Regression via Normal Equations Matrix: A = (M^t * M)^(-1) * M^t * Y (Slides 94-101)
#' @param M Design matrix (include column of 1s if intercept a_0 is needed!).
#' @param y Response vector Y.
#' @examples
#' # Plane fit z = a + b*x + c*y (Slide 98):
#' x <- c(0, 1, 0, -1, 0, 1, 1, 2)
#' y <- c(0, 1, 1, 1, 3, 2, 2, 0)
#' z <- c(2, 3, 4, 5, 7, 4, 5, 0)
#' M <- cbind(1, x, y)
#' matrix_ols(M, z)
#' @export
matrix_ols <- function(M, y) {
  M <- as.matrix(M)
  coeffs <- as.numeric(solve(t(M) %*% M, t(M) %*% y))
  y_est  <- as.numeric(M %*% coeffs)
  list(
    coefficients = coeffs,
    Mt_M = t(M) %*% M,
    Mt_Y = as.numeric(t(M) %*% y),
    y_est = y_est,
    metrics = fit_metrics(y, y_est)
  )
}

#example of usage

# Fit models and immediately request a prediction for x = 6 (new_x = 6):
m_lin <- fit_curve(x, y, type = "linear", new_x = 6)
m_exp <- fit_curve(x, y, type = "exp",    new_x = 6)
m_hyp <- fit_curve(x, y, type = "hyp",    new_x = 6)
# (if y = a*x^b is requested -> type = "power"; if y = a*b^x -> type = "abx"; for parabola -> "quad")

# Compare R2 and MSE:
m_lin$metrics$R2;  m_lin$metrics$MSE
m_exp$metrics$R2;  m_exp$metrics$MSE
m_hyp$metrics$R2;  m_hyp$metrics$MSE

# Suppose m_exp has the largest R2. Let's look at its parameters (a, b) and the prediction for x = 6:
m_exp$parameters   # Returns coefficients a and b
m_exp$predictions  # Ready-to-use prediction for Y when x = 6