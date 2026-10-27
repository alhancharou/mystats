#' Полная таблица частот для сырых данных (x_i, n_i, f_i, N_i, F_i)
#' @param x Вектор исходных сырых данных
#' @export

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

#Если в условии УЖЕ дали вектор частот $n_i$ (например, n_i <- c(10, 15, 5))
#и просят найти $f_i$, никакая функция не нужна, просто пиши:

#f_i <- n_i / sum(n_i)

#Если в задаче УЖЕ дали готовую таблицу частот:

#N_i <- cumsum(n_i)       
#F_i <- N_i / sum(n_i)


#' Базовые расчеты для интервального ряда: ширина (a_i) и середина (x_i)
#' @param breaks Вектор всех границ интервалов по порядку
#' @param n_i Вектор абсолютных частот
#' @export
interv_base <- function(breaks, n_i) {
  L_prev <- breaks[-length(breaks)] 
  L_i <- breaks[-1]                 
  
  a_i <- L_i - L_prev               # Ширина (Interval Width)
  x_i <- (L_prev + L_i) / 2         # Середина (Class Mark)
  
  N <- sum(n_i)
  f_i <- n_i / N                    # Относительная частота
  
  # Формулы со слайда (Heights):
  h_i <- n_i / a_i                  # Высота по абсолютной частоте
  h_f <- f_i / a_i                  # Высота по относительной частоте
  
  data.frame(
    L_prev = L_prev,
    L_i = L_i,
    a_i = a_i,
    x_i = x_i,
    n_i = n_i,
    f_i = f_i,
    h_i = h_i,
    h_f = h_f
  )
}

#Mean 
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

# Нахождение Модального интервала (интервала с максимальной плотностью/высотой)
calc_modal_interval <- function(breaks, n_i) {
  # Ширины интервалов
  a_i <- diff(breaks)
  
  # Плотности (высоты)
  h_i <- n_i / a_i
  
  # Находим индекс максимальной плотности
  max_idx <- which.max(h_i)
  
  # Возвращаем границы модального интервала
  c(Lower = breaks[max_idx], Upper = breaks[max_idx + 1])
}


# Нахождение медианы (Me) для сырых и сгруппированных (точечных) данных
calc_median <- function(x, n_i = rep(1, length(x))) {
  # rep(x, times = n_i) дублирует каждое значение x_i ровно n_i раз
  expanded_data <- rep(x, times = n_i)
  median(expanded_data)
}

# Универсальная функция для квантилей (ДИСКРЕТНЫЕ ДАННЫЕ)
calc_quantiles <- function(x, probs, n_i = rep(1, length(x))) {
  # Разворачиваем частоты в сырой ряд
  expanded_data <- rep(x, times = n_i)
  
  # Считаем нужные квантили
  quantile(expanded_data, probs = probs, type = 2)
}

# Квантили для ИНТЕРВАЛЬНЫХ данных (Continuous Case)
calc_continuous_quantile <- function(c, breaks, n_i) {
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

"The interpretation is:
  AP > 0 Right or positive skew
AP = 0 Symmetry
AP < 0 Left or negative skew"



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

"The interpretation is:
  g1 > 0 Right or positive skew
g1 = 0 Symmetric
g1 < 0 Left or negative skew"



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

"The interpretation is:
  g2 < 0 Less peaked than the normal (Platykurtic).
g2 = 0 Same peakedness as the normal (Mesokurtic).
g2 > 0 More peaked than the normal (Leptokurtic)."









# Шаблон для качественных данных: Bar chart и Pie chart

# 1. Вводим данные
# Категории (строки пишутся в кавычках):
#categories <- c("Chemistry", "Mathematics", "Physics", "Biology", "Geology")
# Частоты (числа):
#frequencies <- c(1500, 750, 1000, 500, 250)

# 2. Строим Bar chart (Столбиковая диаграмма)
# Аргумент names.arg привязывает названия категорий к столбикам
#barplot(frequencies, 
 #       names.arg = categories, 
  #      col = "steelblue",               # Цвет столбиков
   #     main = "Enrolled Students by Section", # Заголовок
    #    ylab = "Number of Students")     # Подпись оси Y

# 3. Строим Pie chart (Круговая диаграмма)
# Аргумент labels добавляет подписи к долям
#pie(frequencies, 
 #   labels = categories, 
  #  col = rainbow(length(frequencies)),  # Автоматическая раскраска разными цветами
   # main = "Enrolled Students by Section")

#rel_freq <- frequencies / sum(frequencies) и передаешь в график rel_freq



# Подсчет процента наблюдений, удовлетворяющих логическому условию
# condition_vec - логический вектор (например, x <= 10)
calc_percentage <- function(condition_vec, n_i = rep(1, length(condition_vec))) {
  # Суммируем частоты только для тех x_i, где условие выполняется (TRUE)
  count_satisfying <- sum(n_i[condition_vec])
  total_count <- sum(n_i)
  
  (count_satisfying / total_count) * 100
}


#rm(list = ls())