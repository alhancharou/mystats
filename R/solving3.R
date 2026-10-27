# Универсальная функция для расчета скользящей средней (Trend)
calc_moving_average <- function(x, k) {
  if (k %% 2 != 0) {
    # Для нечетного k (например, 3, 5): обычное среднее
    weights <- rep(1/k, k)
  } else {
    # Для четного k (например, 4, 12): центрированное скользящее среднее 
    # Крайние значения умножаются на 0.5, как требуют лекции
    weights <- c(0.5, rep(1, k - 1), 0.5) / k
  }
  
  # Применяем фильтр и возвращаем простой вектор
  ma <- stats::filter(x, weights, sides = 2)
  return(as.numeric(ma))
}

# Пример использования на экзамене:
# series <- c(2, 6, 1, 5, 3, 7, 2)
# trend_k3 <- calc_moving_average(series, k = 3)


calc_seasonal_indices <- function(X, trend_cycle, k) {
  # Шаг 1: Очищаем от тренда (получаем E * A)
  E_A <- X / trend_cycle
  
  # Шаг 2: Создаем матрицу, где количество колонок = k (сезоны/месяцы)
  # R расставит данные по годам (строкам) и периодам (колонкам)
  matrix_EA <- matrix(E_A, ncol = k, byrow = TRUE)
  
  # Шаг 3: Нескорректированные индексы (усредняем колонки, игнорируя NA)
  uncorrected <- colMeans(matrix_EA, na.rm = TRUE)
  
  # Шаг 4: Скорректированные индексы (нормализация)
  corrected <- uncorrected / mean(uncorrected)
  
  return(corrected)
}


calc_regression_line <- function(deseasonalized_series) {
  # Автоматически создаем вектор времени t от 1 до N (количества элементов в ряду)
  t <- 1:length(deseasonalized_series)
  
  # Строим модель линейной регрессии (зависимость данных от времени t)
  model <- lm(deseasonalized_series ~ t)
  
  # Возвращаем два коэффициента: свободный член (Intercept) и наклон (Slope)
  return(c(intercept = as.numeric(model$coefficients[1]),
           slope = as.numeric(model$coefficients[2])))
}


# Функция для вычисления чистого Тренда (T) на основе регрессии
calc_pure_trend <- function(deseasonalized_series) {
  # 1. Задаем вектор t от 1 до длины ряда
  t <- 1:length(deseasonalized_series)
  
  # 2. Строим модель
  model <- lm(deseasonalized_series ~ t)
  
  # 3. Функция predict автоматически рассчитывает значения a + b*t для каждого t
  predicted_trend <- predict(model)
  
  return(as.numeric(predicted_trend))
}



# 1. Функция прогноза на h шагов вперед
predict_time_series <- function(deseasonalized_series, seasonal_indices, h = 1) {
  n <- length(deseasonalized_series)
  k <- length(seasonal_indices)
  
  # Строим модель регрессии
  t_vec <- 1:n
  model <- lm(deseasonalized_series ~ t_vec)
  a <- model$coefficients[1]
  b <- model$coefficients[2]
  
  # Будущие временные точки t (например, n+1, n+2, ...)
  future_t <- (n + 1):(n + h)
  
  # Будущий тренд T
  future_trend <- a + b * future_t
  
  # Определяем соответствующий сезон для каждой будущей точки
  # При делении по модулю k остаток 0 соответствует последнему сезону k
  season_idx <- ((future_t - 1) %% k) + 1
  future_season <- seasonal_indices[season_idx]
  
  # Мультипликативный прогноз: T * E
  forecast_values <- future_trend * future_season
  names(forecast_values) <- paste0("t=", future_t, "_season_", season_idx)
  
  return(forecast_values)
}


# 2. Функция расчета автокорреляции для сдвига (лага) T
calc_autocorrelation <- function(x, lag) {
  n <- length(x)
  if (lag >= n) stop("Лаг не может быть больше длины ряда")
  
  # Сдвигаем ряд: сопоставляем x_t и x_{t - lag}
  x_current <- x[(lag + 1):n]
  x_lagged  <- x[1:(n - lag)]
  
  # Считаем классический коэффициент корреляции Пирсона
  r <- cor(x_current, x_lagged)
  return(r)
}


# Функция для вычисления простого индекса и процентного изменения
calc_simple_index <- function(current_val, base_val) {
  index <- current_val / base_val
  variation_pct <- (index - 1) * 100
  
  return(c(index = index, 
           variation_pct = variation_pct))
}


# Функция дефлятирования (перевод номинальных значений в реальные)
calc_deflation <- function(nominal_values, cpi) {
  # Используем классический if, чтобы не обрезать вектор
  if (mean(cpi, na.rm = TRUE) > 5) {
    cpi_factor <- cpi / 100
  } else {
    cpi_factor <- cpi
  }
  
  # Расчет реальной стоимости в ценах базового года
  real_values <- nominal_values / cpi_factor
  
  # Расчет изменения покупательской способности (в процентах к первому периоду)
  purchasing_power_change <- ((real_values / real_values[1]) - 1) * 100
  
  return(data.frame(
    Nominal = nominal_values,
    CPI = cpi,
    Real = real_values,
    PP_Change_pct = purchasing_power_change
  ))
}








#1 (3, 14)

X <- c(2.6767, 3.2050, 4.8714, 4.6275, 5.9008, 5.5200, 5.1186,
       4.9507, 4.6612, 6.2572)

trend_k3 <- calc_moving_average(X, 3); trend_k3
trend_k4 <- calc_moving_average(X, 4); trend_k4

#2 (3, 18)

t <- c(1959, 1960, 1961, 1962, 1963, 1964, 1965)
X <- c(6, 6.2, 5.1, 4.9, 5.2, 5.8, 7)

trend_k3 <- calc_moving_average(X, 3); trend_k3
trend_k4 <- calc_moving_average(X, 4); trend_k4

result_table <- data.frame(
  Year = t,
  Original_X = X,
  Trend_Order3 = trend_k3,
  Trend_Order4 = trend_k4
); result_table

#3 (3, 24)

options(digits = 15)

# 1. "Разворачиваем" таблицу в единый вектор X хронологически
X <- c(178.2, 153.2, 185.9, 163.6,  # 1980
       196.3, 156.9, 197.9, 166.4,  # 1981
       197.3, 159.7, 202.6, 175.6,  # 1982
       209.5, 169.5, 202.4, 179.8,  # 1983
       200.0, 168.6, 216.1, 178.5)  # 1984

# 2. Вычисляем фактор Trend-Cycle (порядок 4 для сезонов)
trend_cycle <- calc_moving_average(X, k = 4)
trend_ts <- ts(trend_cycle, start = c(1980, 1), frequency = 4)
# 3. Смотрим результат
print(trend_cycle)

series_ts <- ts(X, start = c(1980, 1), frequency = 4)
plot(series_ts)

lines(trend_ts, col = 2, lwd = 2)

E_A <- X/trend_cycle; E_A

indices <- calc_seasonal_indices(X, trend_cycle, 4); indices

random_comp <- E_A / indices; random_comp

deseasonalised_series <- X / indices; deseasonalised_series

X_ts <- ts(X, start = c(1980, 1), frequency = 4)
deseasonalised_ts <- ts(deseasonalised_series, start = c(1980, 1), frequency = 4)

plot(X_ts)
lines(deseasonalised_ts)

reg_line <- calc_regression_line(deseasonalised_series); reg_line
trend <- calc_pure_trend(deseasonalised_series); trend
trend_ts <- ts(trend, start = c(1980, 1), frequency = 4)

plot(deseasonalised_ts)
lines(trend_ts)

#4 (3, 56)

options (digits = 15)

xx <- c(33742, 35401, 36147, 38067, 34149, 36176, 37485,
        39047, 34783, 37042, 38008, 40132, 36466, 38680, 39976, 42273,
        39131, 40780, 41852, 43684, 38729, 40427, 41576, 43886, 39131,
        40394, 40956, 42959, 38714, 40062, 41152, 43460, 38695, 39780,
        40923, 44093, 40777, 41778, 43160, 45897, 41947, 44061, 44378,
        47237, 43315, 43396, 44843, 46835, 42833, 43548, 44637, 47107,
        42552, 43526, 45039, 47940, 43740, 45007, 46667, 49325, 44878,
        46234, 47055, 50318, 46354, 47260, 48883, 52605, 48527, 50237,
        51592, 55152, 50451, 52294, 54633, 58802, 53990, 55477, 57850,
        61978)

xx_ts <- ts(xx, start = c(1969, 1), frequency = 4); xx_ts
plot(xx_ts)

trend_cycle <- calc_moving_average(xx, 4); trend_cycle
indices <- calc_seasonal_indices(xx, trend_cycle, 4); indices
random_comp <- (xx/trend_cycle)/indices; random_comp
deseasonalized <- xx/indices; deseasonalized
reg_line <- calc_regression_line(deseasonalized); reg_line
pure_trend <- calc_pure_trend(deseasonalized); pure_trend


#5(3, 59)

data <- c(1.2, 1.4, 2.0, 1.23, 1.37, 2.07, 1.25, 1.47, 2.09)

r <- calc_autocorrelation(data, 3); r

#6 (3, 64)

simple_ind <- calc_simple_index(148, 123); simple_ind

#7 (3, 65)

base_period <- 1254
vals <- c(1254, 1345, 1408, 1451)
simple_ind <- calc_simple_index(vals, base_period); simple_ind

#8 (3, 69)

salaries <- c(100, 150)
cpi_vals <- c(100, 200)

ans <- calc_deflation(salaries, cpi_vals); ans

#9 (3, 71)

cpi <- 135/100

ans_1 <- 1000 / cpi; ans_1
ans_2 <- 1000 * cpi; ans_2

#10

cpi_90_69 <- 10
cpi_90_06 <- 0.75

real_90 <- 24000 * cpi_90_06
real_69 <- real_90 / cpi_90_69; real_69
