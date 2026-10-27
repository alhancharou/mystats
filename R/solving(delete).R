freq_table <- function(x) {
  n_i_table <- table(x)
  N <- sum(n_i_table)
  
  n_i <- as.numeric(n_i_table)
  f_i <- n_i / N
  
  data.frame(
    x_i = as.numeric(names(n_i_table)), # Уникальные значения
    n_i = n_i,                          # Абсолютная частота
    f_i = f_i,                          # Относительная частота
    N_i = cumsum(n_i),                  # Абсолютная кумулятивная частота
    F_i = cumsum(f_i)                   # Относительная кумулятивная частота
  )
}

# Mean 
calc_group_mean <- function(x, n_i = rep(1, length(x))) {
  sum(x * n_i) / sum(n_i)
}

# Среднее квадратическое (Root Mean Square, RMS)
calc_quad_mean <- function(x, n_i = rep(1, length(x))) {
  sqrt(sum((x^2) * n_i) / sum(n_i))
}

# Нахождение статистической моды (Mo)
calc_mode <- function(x, n_i = rep(1, length(x))) {
  if (all(n_i == 1)) {
    # Если переданы сырые данные, R сам считает частоты
    freq_table <- table(x)
    max_freq <- max(freq_table)
    as.numeric(names(freq_table)[freq_table == max_freq])
  } else {
    # Если переданы уникальные значения (x_i) и их частоты (n_i)
    max_freq <- max(n_i)
    x[n_i == max_freq]
  }
}

# Нахождение медианы (Me) для сырых и сгруппированных (точечных) данных
calc_median <- function(x, n_i = rep(1, length(x))) {
  # rep(x, times = n_i) дублирует каждое значение x_i ровно n_i раз
  expanded_data <- rep(x, times = n_i)
  median(expanded_data)
}

calc_quantiles <- function(x, probs, n_i = rep(1, length(x))) {
  # Разворачиваем частоты в сырой ряд
  expanded_data <- rep(x, times = n_i)
  
  # Считаем нужные квантили
  quantile(expanded_data, probs = probs, type = 2)
}

# Квантили для ИНТЕРВАЛЬНЫХ данных (Continuous Case)
calc_interv_quantile <- function(c, breaks, n_i) {
  N <- sum(n_i)
  cN <- c * N                      # Позиция квантиля
  N_i <- cumsum(n_i)               # Накопленные частоты
  
  # Ищем первый интервал, накопленная частота которого превышает cN
  idx <- which(N_i > cN)[1]
  
  # Извлекаем компоненты для формулы со слайда:
  L_prev <- breaks[idx]            # L_{i-1} (нижняя граница интервала)
  L_i <- breaks[idx + 1]           # Верхняя граница
  a_i <- L_i - L_prev              # Ширина интервала
  n_interval <- n_i[idx]           # Абсолютная частота этого интервала
  
  # N_{i-1} (накопленная частота предыдущего интервала)
  N_prev <- ifelse(idx == 1, 0, N_i[idx - 1])
  
  # Сама формула со слайда:
  L_prev + ((cN - N_prev) / n_interval) * a_i
}

# Меры размаха (Range, Interquartile Range, Intercentile Range)
calc_ranges <- function(x, n_i = rep(1, length(x))) {
  expanded_data <- rep(x, times = n_i)
  
  R <- max(expanded_data) - min(expanded_data)
  
  # Считаем нужные квантили с type = 2
  qs <- quantile(expanded_data, probs = c(0.25, 0.75, 0.01, 0.99), type = 2)
  
  R_Q <- qs["75%"] - qs["25%"]  # Q3 - Q1
  R_P <- qs["99%"] - qs["1%"]   # P99 - P1
  
  # Выводим все три значения
  c(Range = R, Interquartile_Range = R_Q, Intercentile_Range = R_P)
}


# Среднее линейное отклонение (Mean Deviation)
calc_mean_deviation <- function(x, n_i = rep(1, length(x)), p = NULL) {
  N <- sum(n_i)
  
  # Если параметр p не задан вручную, по умолчанию используем среднее арифметическое
  if (is.null(p)) {
    p <- sum(x * n_i) / N
  }
  
  # Формула со слайда: сумма модулей (abs) разностей, умноженных на частоты, деленная на N
  sum(abs(x - p) * n_i) / N
}

# Средняя квадратическая ошибка (Mean Squared Error)
calc_mse <- function(x, n_i = rep(1, length(x)), p = NULL) {
  N <- sum(n_i)
  
  # Если параметр p не задан явно, R посчитает среднее арифметическое
  if (is.null(p)) {
    p <- sum(x * n_i) / N
  }
  
  # Расчет по формуле со слайда
  sum(n_i * (x - p)^2) / N
}

# Дисперсия (Variance, V или генеральная sigma^2)
calc_variance <- function(x, n_i = rep(1, length(x))) {
  N <- sum(n_i)
  mu <- sum(x * n_i) / N # Среднее арифметическое
  
  # Формула со слайда: сумма квадратов отклонений от среднего, деленная на N
  sum(n_i * (x - mu)^2) / N
}

# Стандартное отклонение (Standard Deviation, sigma)
calc_sd <- function(x, n_i = rep(1, length(x))) {
  # Корень из дисперсии
  sqrt(calc_variance(x, n_i))
}

# Выборочная квазидисперсия (Sample quasi-variance, s^2)
calc_quasi_variance <- function(x, n_i = rep(1, length(x))) {
  expanded_data <- rep(x, times = n_i)
  
  # Базовая функция var() автоматически делит на n-1
  var(expanded_data)
}

# Квази-стандартное отклонение (Квадратный корень из s^2)
calc_quasi_sd <- function(x, n_i = rep(1, length(x))) {
  expanded_data <- rep(x, times = n_i)
  
  # Базовая функция sd() автоматически использует n-1
  sd(expanded_data)
}

# Стандартизованная переменная (Z-score)
calc_standardised <- function(x, n_i = rep(1, length(x))) {
  N <- sum(n_i)
  mu <- sum(x * n_i) / N         # Среднее арифметическое (mu)
  sigma <- calc_sd(x, n_i)       # Генеральное стандартное отклонение (sigma)
  
  # Формула со слайда: возвращаем вектор стандартизованных значений
  (x - mu) / sigma
}

# Коэффициент вариации Пирсона (Pearson's Coefficient of Variation)
calc_cv <- function(x, n_i = rep(1, length(x)), as_percentage = FALSE) {
  N <- sum(n_i)
  mu <- sum(x * n_i) / N
  sigma <- calc_sd(x, n_i) # Используем нашу генеральную функцию
  
  # Формула со слайда: sigma / |mu|
  cv <- sigma / abs(mu)
  
  # Возвращаем результат в процентах или в долях
  if (as_percentage) {
    return(cv * 100)
  } else {
    return(cv)
  }
}

# Момент порядка r относительно точки c
calc_moment <- function(x, r, c = 0, n_i = rep(1, length(x))) {
  N <- sum(n_i)
  
  # Формула со слайда: сумма частот, умноженных на отклонения в степени r, деленная на N
  sum(n_i * (x - c)^r) / N
}

# Коэффициент асимметрии Пирсона (Pearson's Skewness Coefficient)
calc_pearson_skewness <- function(x, n_i = rep(1, length(x))) {
  N <- sum(n_i)
  mu <- sum(x * n_i) / N         # Среднее (mu)
  
  # Находим моду и берем первый элемент для страховки от мультимодальности
  mo <- calc_mode(x, n_i)[1]        
  
  sigma <- calc_sd(x, n_i)       # Генеральное стандартное отклонение (sigma)
  
  # Формула со слайда
  (mu - mo) / sigma
}

# Коэффициент асимметрии Фишера (Fisher's Skewness Coefficient)
calc_fisher_skewness <- function(x, n_i = rep(1, length(x))) {
  # 1. Находим среднее (mu)
  N <- sum(n_i)
  mu <- sum(x * n_i) / N
  
  # 2. Считаем 3-й центральный момент (mu_3), используя нашу функцию calc_moment
  mu_3 <- calc_moment(x, r = 3, c = mu, n_i = n_i)
  
  # 3. Считаем генеральное стандартное отклонение (sigma)
  sigma <- calc_sd(x, n_i)
  
  # 4. Применяем формулу со слайда
  mu_3 / (sigma^3)
}

# Коэффициент эксцесса Фишера (Fisher's Kurtosis Coefficient)
calc_fisher_kurtosis <- function(x, n_i = rep(1, length(x))) {
  N <- sum(n_i)
  mu <- sum(x * n_i) / N
  
  # Считаем 4-й центральный момент (mu_4), где c = mu
  mu_4 <- calc_moment(x, r = 4, c = mu, n_i = n_i)
  
  # Считаем генеральное стандартное отклонение (sigma)
  sigma <- calc_sd(x, n_i)
  
  # Формула со слайда
  (mu_4 / sigma^4) - 3
}

#1
x <- c(6, 8, 6, 8, 6, 8, 12, 6, 8, 8, 6, 8, 8, 8, 12,
       12, 8, 8, 12, 6, 8, 6, 6, 8, 12, 6, 6, 6, 6, 6)

freq_table(x)

#2

x_i <- c(13.725, 13.975, 14.225, 14.725)
n_i <- c(7, 14, 18, 6)

calc_group_mean(x_i, n_i)

#3

grades <- c(2.6, 3.7, 5.1, 4.9, 6.4)
weights <- c(1, 1, 1, 2, 3)

calc_group_mean (grades, weights)

#4

calls <- c(2, 3, 1, 0, 4, 3)
calc_quad_mean (calls)

#5

data <- c(7, 11, 11, 8, 12, 9, 6, 6)
data2 <- c(7, 11, 8, 12, 9, 6)
calc_mode(data)
calc_mode(data2)

#6

data3 <- c(11, 11, 16, 17, 25)
calc_median(data3)
data4 <- c(1, 4, 8, 8, 10, 16, 16, 19)
calc_median(data4)

#7

data5 <- c(2, 5, 3, 4, 7, 0, 11, 2, 3, 8)
q_p <- c(0.25, 0.50, 0.75, 0.37, 0.68)
calc_quantiles(data5, q_p)
calc_median(data5)

#8

data6 <- c(10, 20, 35, 44, 55, 64, 75, 81, 87, 99,
           11, 22, 36, 48, 56, 68, 76, 82, 89, 101,
           13, 23, 38, 49, 57, 69, 76, 83, 90, 102,
           15, 23, 41, 50, 60, 70, 78, 83, 94, 105,
           18, 30, 44, 50, 63, 73, 80, 85, 96, 107)

q <- c(0.05, 0.95, 0.25, 0.50, 0.75)
calc_quantiles(data6, q)

#9

data7 <- c(5, 2, 3, 3, 3, 5, 7)
my_mean <- calc_group_mean(data7)
my_median <- calc_median(data7)
md_mean <- calc_mean_deviation(data7, p = my_mean); md_mean
md_median <- calc_mean_deviation(data7, p = my_median); md_median
mse_mean <- calc_mse(data7, p = my_mean); mse_mean
mse_median <- calc_mse(data7, p = my_median); mse_median

#10

x_i <- c(0, 0, 0, 1, 1, 1, 3, 3, 3)
x_i_star <- c(0, 1, 3, 0, 1, 3, 0, 1, 3)
n_i <- c(6, 3, 2, 3, 5, 1, 1, 3, 7)

error <- x_i_star - x_i

#если в задаче фигурирует "ошибка" (error / predicted vs actual), ты всегда принудительно задаешь p = 0
my_mean <- calc_group_mean(error, n_i); my_mean
my_quad_mean <- calc_quad_mean(error, n_i); my_quad_mean
md_error <- calc_mean_deviation(error, n_i, 0); md_error
mse_error <- calc_mse(error, n_i, 0); mse_error

#11

list_a <- c(12, 10, 9, 9, 10)
list_b <- c(5, 10, 16, 15, 4)
x <- calc_sd(list_a); x  
y <- calc_sd(list_b); y

#12

# Студент A: оценка 8.5, среднее 7.9, отклонение 0.8
Z_A <- (8.5 - 7.9) / 0.8; Z_A
# Результат: 0.75

# Студент B: оценка 7.4, среднее 7.0, отклонение 0.5
Z_B <- (7.4 - 7.0) / 0.5; Z_B
# Результат: 0.8

#13

data <- c(0, 0, 2, 2, 2, 2, 3, 3, 4, 4, 7, 7)
my_median <- calc_median(data); my_median
x <- calc_moment(data, 1, my_median); x
y <- calc_moment(data, 2, my_median); y

#ex1

raw <- c(21.8, 40.8, 25.1, 34.6, 39.1, 37.8, 38.4, 36.0, 34.4, 30.1,
         25.1, 37.4, 38.0, 25.1, 37.8, 29.4, 34.3, 27.2, 23.8, 29.9,
         26.8, 36.8, 25.0, 36.1, 24.9, 38.5, 23.8, 36.8, 31.7, 20.6,
         29.8, 24.3, 36.5, 39.0, 40.9, 27.4, 22.2, 28.2, 31.1, 31.5,
         38.9, 24.6, 22.3, 39.1, 25.9, 21.6, 24.1, 29.4, 22.1, 28.9,
         42.4, 38.5, 43.8, 25.9, 26.5, 33.8, 20.8, 37.5, 37.6, 38.4,
         27.4, 28.5, 31.5, 28.5, 25.1, 31.3, 22.3, 25.3, 22.6, 22.3,
         20.4, 31.1, 21.2, 28.7, 24.9, 23.2, 26.4, 20.6, 27.8, 26.9,
         30.2, 20.5, 31.1, 30.0, 26.7, 25.4, 25.1, 30.3, 28.8, 23.6,
         25.3, 26.6, 23.8, 28.5, 27.0, 23.9, 21.9, 31.9, 26.0, 23.5)

mean_raw <- calc_group_mean(raw); mean_raw
med_raw <- calc_median(raw); med_raw
IQR_raw <- calc_ranges(raw); IQR_raw
s_d_raw <- calc_sd(raw); s_d_raw
cv_raw <- calc_cv(raw); cv_raw
skew_raw <- calc_pearson_skewness(raw); skew_raw
skew_raw2 <- calc_fisher_skewness(raw); skew_raw2
fish_raw <- calc_fisher_kurtosis(raw); fish_raw


#ex2
x <- c(5, 2, 4, 9, 5, 7, 4, 5, 6, 5, 7, 7, 5, 5, 2, 10, 5, 6, 5, 4,
       5, 8, 8, 4, 0, 8, 4, 8, 6, 6, 3, 6, 7, 6, 6, 7, 6, 7, 3, 5,
       6, 9, 6, 1, 4, 6, 3, 5, 5, 6, 7)

x_freq <- freq_table(x); x_freq

#ex3

x_raw <- c(3, 4, 8, 10)
n_i <- c(3, 15, 2, 4)
calc_median(x_raw, n_i)

#ex4

x_raw <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
n_i <- c(2, 9, 14, 20, 18, 15, 9, 6, 4, 2, 1)

calc_group_mean(x_raw, n_i)
calc_median(x_raw, n_i)
calc_mode(x_raw, n_i)
calc_quantiles(x_raw, 0.25, n_i)
calc_quantiles(x_raw, 0.75, n_i)

#problem1








# 1. Значения X (при N=2, mu=10, дисперсии=4, отклонение равно sqrt(4) = 2)
x_data <- c(10 - 2, 10 + 2) # Получаем вектор c(8, 12)

# 2. Сравнение разброса (Коэффициент вариации: CV = sd / |mu|)
# Для X используем нашу функцию, для Y считаем напрямую (sd=2, mu=8)
cv_x <- calc_cv(x_data)     # Результат: 0.20
cv_y <- sqrt(4) / 8         # Результат: 0.25 (Разброс Y больше)

# 3. Сравнение эксцесса (Kurtosis)
# Для X просто отдаем сырой вектор нашей функции
kurt_x <- calc_fisher_kurtosis(x_data) # Результат: -2
kurt_x
# Для Y сырых данных нет, считаем mu_4 (центральный момент) через обычные (m2, m3, m4)
m2 <- 4 + 8^2                          # 68
m2
m3 <- 3 * m2 * 8 - 2 * (8^3)           # 608 (т.к. mu_3 = 0 для симметричных)
m3
mu4 <- 5648 - 4*m3*8 + 6*m2*8^2 - 3*8^4  # 16
mu4
kurt_y <- (mu4 / 4^2) - 3              # Результат: -2 (Эксцесс одинаковый)
kurt_y