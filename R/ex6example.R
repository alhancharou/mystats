rm(list = ls())


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






#---------------------- EX 6 -----------------------------

options(digits = 12)
q_s <- c(25, 31, 40, 28,
         27, 34, 43, 30,
         29, 36, 47, 33)

q_s_ts <- ts(q_s, start = (2021), frequency = 4); q_s_ts

ma <- calc_moving_average(q_s_ts, 4); ma

ind <- calc_seasonal_indices(q_s_ts, ma, 4); ind

#random_comp <- (xx/trend_cycle)/indices; random_comp
#deseasonalized <- xx/indices; deseasonalized

deseasonalized <- q_s/ind
calc_pure_trend(deseasonalized)
calc_regression_line(deseasonalized)
#$T_t = 29.2935 + 0.6547 x t$

q12024 <- predict_time_series(deseasonalized, ind, 1); q12024
q32024 <- predict_time_series(deseasonalized, ind, 3); q32024

ac <- calc_autocorrelation(q_s_ts, 4); ac

annual2023 <- c(29, 36, 47, 33)

annual2021 <- c(25, 31, 40, 28)

an2023 <- sum(annual2023)
an2021 <- sum(annual2021)
nom <- an2023/an2021; nom

#real 145 / 1.092 = 132.78. Затем находим рост: (132.78 / 124) - 1)