# Функция для создания простой таблицы частот двумерной переменной
create_bivariate_table <- function(x, y) {
  # 1. table() считает частоты пар, а as.data.frame() делает из этого колонки
  freq_df <- as.data.frame(table(X = x, Y = y))
  
  # 2. Оставляем только те пары, которые реально встречались (n_i > 0)
  freq_df <- freq_df[freq_df$Freq > 0, ]
  
  # 3. Добавляем колонку относительной частоты (f_i)
  freq_df$Rel_Freq <- freq_df$Freq / sum(freq_df$Freq)
  
  # 4. Переименовываем колонки в точности как на слайде
  colnames(freq_df) <- c("Variable_X", "Variable_Y", "Absolute_Freq_n", "Relative_Freq_f")
  
  # Сбрасываем нумерацию строк для красоты
  rownames(freq_df) <- NULL
  
  return(freq_df)
}


# Создание двумерной таблицы абсолютных частот с маргиналами
calc_joint_freq_table <- function(x, y) {
  # Считаем базовые частоты (n_ij) и добавляем итоговые суммы (n_i., n_.j и n)
  joint_table <- addmargins(table(X = x, Y = y))
  
  return(joint_table)
}


# Создание двумерной таблицы относительных частот с маргиналами
calc_joint_rel_freq_table <- function(x, y) {
  # 1. Считаем базовые абсолютные частоты (n_ij)
  base_table <- table(X = x, Y = y)
  
  # 2. Переводим в относительные частоты (f_ij = n_ij / n)
  rel_table <- prop.table(base_table)
  
  # 3. Добавляем итоговые суммы (в правом нижнем углу будет 1)
  joint_rel_table <- addmargins(rel_table)
  
  return(joint_rel_table)
}


# 1. Маргинальное распределение для переменной X (сумма по строкам)
calc_marginal_x <- function(base_table) {
  n_i <- margin.table(base_table, margin = 1) # margin = 1 означает строки
  f_i <- prop.table(n_i)
  
  # Объединяем абсолютные и относительные частоты для удобства
  rbind(Absolute_n = n_i, Relative_f = f_i)
}

# 2. Маргинальное распределение для переменной Y (сумма по столбцам)
calc_marginal_y <- function(base_table) {
  n_j <- margin.table(base_table, margin = 2) # margin = 2 означает столбцы
  f_j <- prop.table(n_j)
  
  rbind(Absolute_n = n_j, Relative_f = f_j)
}


# 1. Условное распределение X при заданном Y (срез по столбцу)
calc_cond_x_given_y <- function(base_matrix, y_val) {
  # Превращаем значение в строку для поиска нужного столбца по имени
  col_name <- as.character(y_val)
  
  # Абсолютные частоты: извлекаем конкретный столбец
  abs_freq <- base_matrix[, col_name]
  
  # Отбрасываем элемент с именем "Sum", если он случайно попал в срез
  if ("Sum" %in% names(abs_freq)) {
    abs_freq <- abs_freq[names(abs_freq) != "Sum"]
  }
  
  # Относительные частоты: делим значения столбца на его фактическую сумму (n.j)
  rel_freq <- prop.table(abs_freq)
  
  return(rbind(Absolute_n = abs_freq, Relative_f = rel_freq))
}

# 2. Условное распределение Y при заданном X (срез по строке)
calc_cond_y_given_x <- function(base_matrix, x_val) {
  row_name <- as.character(x_val)
  
  # Извлекаем строку
  abs_freq <- base_matrix[row_name, ]
  
  # Отбрасываем элемент с именем "Sum", если он есть, чтобы prop.table не сошел с ума
  if ("Sum" %in% names(abs_freq)) {
    abs_freq <- abs_freq[names(abs_freq) != "Sum"]
  }
  
  rel_freq <- prop.table(abs_freq)
  
  return(rbind(Absolute_n = abs_freq, Relative_f = rel_freq))
}


# Функция для расчета двумерного момента порядка (r, s) относительно точки (a, b)
calc_bivariate_moment <- function(freq_matrix, r, s, a = 0, b = 0) {
  # 1. Извлекаем численные значения X и Y из названий строк и столбцов
  x_vals <- as.numeric(rownames(freq_matrix))
  y_vals <- as.numeric(colnames(freq_matrix))
  
  # 2. Рассчитываем матрицу относительных частот f_ij
  f_ij <- prop.table(freq_matrix)
  
  # 3. Считаем отклонения для X и Y согласно формуле
  x_part <- (x_vals - a)^r
  y_part <- (y_vals - b)^s
  
  # 4. Формируем матрицу произведений отклонений и умножаем на частоты f_ij
  moment_matrix <- outer(x_part, y_part) * f_ij
  
  # 5. Двойная сумма - это общая сумма всех элементов получившейся матрицы
  return(sum(moment_matrix))
}


# Функция для прямого расчета ковариации из таблицы частот
calc_covariance <- function(freq_matrix) {
  # 1. Извлекаем значения X и Y
  x_vals <- as.numeric(rownames(freq_matrix))
  y_vals <- as.numeric(colnames(freq_matrix))
  
  # 2. Переводим абсолютные частоты в относительные (f_ij)
  f_ij <- prop.table(freq_matrix)
  
  # 3. Находим маргинальные распределения (относительные частоты)
  f_i <- margin.table(f_ij, margin = 1)
  f_j <- margin.table(f_ij, margin = 2)
  
  # 4. Рассчитываем маргинальные средние (x_bar и y_bar)
  mean_x <- sum(x_vals * f_i)
  mean_y <- sum(y_vals * f_j)
  
  # 5. Считаем обычный момент m_11 (сумма x * y * f_ij)
  m_11 <- sum(outer(x_vals, y_vals) * f_ij)
  
  # 6. Применяем альтернативную формулу: m_11 - (mean_x * mean_y)
  cov_xy <- m_11 - (mean_x * mean_y)
  
  return(cov_xy)
}


# Функция для проверки переменных на независимость
calc_independence <- function(freq_matrix) {
  # 1. Получаем матрицу относительных частот f_ij
  f_ij <- prop.table(freq_matrix)
  
  # 2. Считаем маргинальные частоты
  f_i <- margin.table(f_ij, 1) # по строкам
  f_j <- margin.table(f_ij, 2) # по столбцам
  
  # 3. Создаем "теоретическую" матрицу независимости (перемножение f_i и f_j)
  # outer() перемножит каждый элемент f_i на каждый элемент f_j
  expected_f_ij <- outer(f_i, f_j)
  
  # 4. Проверяем, равны ли матрицы (используем all.equal для обхода проблем с округлением)
  is_independent <- isTRUE(all.equal(f_ij, expected_f_ij, check.attributes = FALSE))
  
  if (is_independent) {
    print("Вывод: Переменные НЕЗАВИСИМЫ (f_ij = f_i * f_j).")
  } else {
    print("Вывод: Переменные ЗАВИСИМЫ (f_ij != f_i * f_j).")
  }
  
  # Возвращаем TRUE или FALSE для дальнейшего использования в коде
  return(is_independent)
}


calc_functional_dependency <- function(base_matrix) {
  # Проверяем условие X = f(Y) (каждый столбец имеет ровно 1 ненулевой элемент)
  # colSums(base_matrix != 0) считает количество не-нулей по столбцам
  x_depends_on_y <- all(colSums(base_matrix != 0) == 1)
  
  # Заодно добавим проверку обратной зависимости Y = f(X)
  # (каждая строка имеет ровно 1 ненулевой элемент)
  y_depends_on_x <- all(rowSums(base_matrix != 0) == 1)
  
  # Возвращаем понятный результат
  return(list(
    "X = f(Y) (Columns check)" = x_depends_on_y,
    "Y = f(X) (Rows check)" = y_depends_on_x
  ))
}



calc_model_fit <- function(y_actual, y_est, var_y) {
  # 1. Вектор остатков e_i = y_i - y_est
  e <- y_actual - y_est
  
  # 2. Среднеквадратичная ошибка (MSE)[cite: 9]
  n <- length(y_actual)
  mse <- sum(e^2) / n
  
  # 3. Коэффициент детерминации (R^2)[cite: 9]
  r_squared <- 1 - (mse / var_y)
  
  # 4. Объясненная дисперсия (V_ex)
  v_ex <- r_squared * var_y
  
  return(list(
    "MSE" = mse,
    "R_squared" = r_squared,
    "Variance_Explained" = v_ex
  ))
}


#y = a + bx
calc_simple_linear_regression <- function(x, y) {
  n <- length(x)
  
  # Центры тяжести (средние значения)
  mean_x <- mean(x)
  mean_y <- mean(y)
  
  # Смещенные (генеральные) дисперсии и ковариация (деление на n)
  cov_xy <- sum((x - mean_x) * (y - mean_y)) / n
  var_x <- sum((x - mean_x)^2) / n
  var_y <- sum((y - mean_y)^2) / n
  
  # 1. Линия Y на X (y = a + bx)
  # Наклон b и свободный член a
  b <- cov_xy / var_x
  a <- mean_y - b * mean_x
  
  # 2. Линия X на Y (x = a' + b'y)
  # Наклон b' и свободный член a'
  b_prime <- cov_xy / var_y
  a_prime <- mean_x - b_prime * mean_y
  
  # 3. Линейный коэффициент корреляции Пирсона (r)
  r <- cov_xy / (sqrt(var_x) * sqrt(var_y))
  
  # 4. Расчет оценок и остатков для модели Y/X
  y_est <- a + b * x
  e <- y - y_est
  
  return(list(
    "Center_of_Gravity" = c(X = mean_x, Y = mean_y),
    "Line_Y_on_X_Coefficients" = c(a = a, b = b),
    "Line_X_on_Y_Coefficients" = c(a_prime = a_prime, b_prime = b_prime),
    "Pearson_Correlation" = r,
    "Estimated_Y" = y_est,
    "Residuals" = e
  ))
}


calc_freq_regression <- function(freq_matrix) {
  n_total <- sum(freq_matrix)
  x_mid <- as.numeric(rownames(freq_matrix))
  y_mid <- as.numeric(colnames(freq_matrix))
  
  n_x <- rowSums(freq_matrix)
  n_y <- colSums(freq_matrix)
  
  mean_x <- sum(x_mid * n_x) / n_total
  mean_y <- sum(y_mid * n_y) / n_total
  
  var_x <- sum(n_x * (x_mid - mean_x)^2) / n_total
  var_y <- sum(n_y * (y_mid - mean_y)^2) / n_total
  
  # Стандартные отклонения для формулы Пирсона
  sd_x <- sqrt(var_x)
  sd_y <- sqrt(var_y)
  
  cov_xy <- 0
  for(i in 1:length(x_mid)) {
    for(j in 1:length(y_mid)) {
      if(freq_matrix[i, j] > 0) {
        cov_xy <- cov_xy + freq_matrix[i, j] * (x_mid[i] - mean_x) * (y_mid[j] - mean_y)
      }
    }
  }
  cov_xy <- cov_xy / n_total
  
  # Линии регрессии
  b <- cov_xy / var_x
  a <- mean_y - b * mean_x
  
  b_prime <- cov_xy / var_y
  a_prime <- mean_x - b_prime * mean_y
  
  # Коэффициент корреляции Пирсона (r)
  r <- cov_xy / (sd_x * sd_y)
  
  # Ошибка линейной модели Y на X (RMSE)
  # Исходя из декомпозиции дисперсии: var_y = MSE + b^2 * var_x
  mse_y_x <- var_y - (b^2 * var_x)
  rmse_y_x <- sqrt(max(0, mse_y_x)) # max(0) спасает от микро-погрешностей округления в R
  
  # RMSE X on Y
  mse_x_y <- var_x - (b_prime^2 * var_y)
  rmse_x_y <- sqrt(max(0, mse_x_y))
  
  return(list(
    "Covariance" = cov_xy,
    "Line_Y_on_X_Coefficients" = c(a = a, b = b),
    "Line_X_on_Y_Coefficients" = c(a_prime = a_prime, b_prime = b_prime),
    "Pearson_Correlation" = r,
    "RMSE_Y_on_X" = rmse_y_x,
    "RMSE_X_on_Y" = rmse_x_y
  ))
}


calc_cond_distribution <- function(base_matrix, margin = "Y|X", cond_vals) {
  # Приводим значения к тексту для безопасного поиска по именам
  cond_str <- as.character(cond_vals)
  
  if (margin == "Y|X") {
    # УСЛОВИЕ НА X (а X у нас теперь в СТРОКАХ). 
    # Вырезаем нужные строки и суммируем их по столбцам, чтобы получить распределение Y.
    subset_mat <- base_matrix[cond_str, , drop = FALSE]
    freqs <- colSums(subset_mat)
    vals <- as.numeric(colnames(base_matrix)) # Значения Y берем из столбцов
    cond_name <- paste0("Y | X in {", paste(cond_str, collapse = ", "), "}")
    
  } else if (margin == "X|Y") {
    # УСЛОВИЕ НА Y (а Y у нас теперь в СТОЛБЦАХ). 
    # Вырезаем нужные столбцы и суммируем их по строкам, чтобы получить распределение X.
    subset_mat <- base_matrix[ , cond_str, drop = FALSE]
    freqs <- rowSums(subset_mat)
    vals <- as.numeric(rownames(base_matrix)) # Значения X берем из строк
    cond_name <- paste0("X | Y in {", paste(cond_str, collapse = ", "), "}")
    
  } else {
    stop("Параметр margin должен быть 'Y|X' или 'X|Y'")
  }
  
  N <- sum(freqs)
  if (N == 0) return("Ошибка: Сумма частот для этого условия равна нулю.")
  
  # Базовые характеристики
  probs <- freqs / N
  mean_val <- sum(vals * freqs) / N
  var_val <- sum(freqs * (vals - mean_val)^2) / N
  
  # Центральные моменты и форма распределения
  m3 <- sum(freqs * (vals - mean_val)^3) / N
  m4 <- sum(freqs * (vals - mean_val)^4) / N
  
  skewness <- ifelse(var_val == 0, NA, m3 / (var_val^(3/2)))
  kurtosis <- ifelse(var_val == 0, NA, m4 / (var_val^2))
  
  # Формируем красивую таблицу распределения
  dist_table <- data.frame(Value = vals, Freq = freqs, Prob = probs)
  
  return(list(
    "Condition" = cond_name,
    "Distribution" = dist_table,
    "Mean" = mean_val,
    "Variance" = var_val,
    "Skewness" = skewness,
    "Kurtosis_Pearson" = kurtosis,
    "Excess_Kurtosis" = kurtosis - 3
  ))
}


calc_pearson_correlation <- function(freq_matrix) {
  n_total <- sum(freq_matrix)
  x_mid <- as.numeric(rownames(freq_matrix))
  y_mid <- as.numeric(colnames(freq_matrix))
  
  n_x <- rowSums(freq_matrix)
  n_y <- colSums(freq_matrix)
  
  mean_x <- sum(x_mid * n_x) / n_total
  mean_y <- sum(y_mid * n_y) / n_total
  
  var_x <- sum(n_x * (x_mid - mean_x)^2) / n_total
  var_y <- sum(n_y * (y_mid - mean_y)^2) / n_total
  
  cov_xy <- 0
  for(i in 1:length(x_mid)) {
    for(j in 1:length(y_mid)) {
      if(freq_matrix[i, j] > 0) {
        cov_xy <- cov_xy + freq_matrix[i, j] * (x_mid[i] - mean_x) * (y_mid[j] - mean_y)
      }
    }
  }
  cov_xy <- cov_xy / n_total
  
  # Формула Пирсона
  r <- cov_xy / sqrt(var_x * var_y)
  return(r)
}


calc_rmse <- function(freq_matrix, margin = "Y|X") {
  n_total <- sum(freq_matrix)
  x_mid <- as.numeric(rownames(freq_matrix))
  y_mid <- as.numeric(colnames(freq_matrix))
  
  n_x <- rowSums(freq_matrix)
  n_y <- colSums(freq_matrix)
  
  mean_x <- sum(x_mid * n_x) / n_total
  mean_y <- sum(y_mid * n_y) / n_total
  
  var_x <- sum(n_x * (x_mid - mean_x)^2) / n_total
  var_y <- sum(n_y * (y_mid - mean_y)^2) / n_total
  
  cov_xy <- 0
  for(i in 1:length(x_mid)) {
    for(j in 1:length(y_mid)) {
      if(freq_matrix[i, j] > 0) {
        cov_xy <- cov_xy + freq_matrix[i, j] * (x_mid[i] - mean_x) * (y_mid[j] - mean_y)
      }
    }
  }
  cov_xy <- cov_xy / n_total
  
  if (margin == "Y|X") {
    b <- cov_xy / var_x
    mse <- var_y - (b^2 * var_x)
  } else if (margin == "X|Y") {
    b_prime <- cov_xy / var_y
    mse <- var_x - (b_prime^2 * var_y)
  } else {
    stop("Параметр margin должен быть 'Y|X' или 'X|Y'")
  }
  
  return(sqrt(mse))
}


#pear_cor <- cor(D, V) #Pearson correlation когда просто 2 вектора 



# Функция принимает коэффициенты уравнений в формате Ax + By = C
# eq1 = c(A, B, C), eq2 = c(A, B, C)
calc_bivariate_regression <- function(eq1, eq2, N, var_x, pt_remove = NULL) {
  
  # 1. Поиск средних значений через систему уравнений (правило Крамера)
  det_main <- eq1[1]*eq2[2] - eq1[2]*eq2[1]
  mean_x <- (eq1[3]*eq2[2] - eq1[2]*eq2[3]) / det_main
  mean_y <- (eq1[1]*eq2[3] - eq1[3]*eq2[1]) / det_main
  
  # 2. Определение, какая линия есть какая (проверяем произведение наклонов)
  b_A <- -eq1[1]/eq1[2]
  d_A <- -eq2[2]/eq2[1]
  r2_A <- b_A * d_A
  
  if (r2_A >= 0 && r2_A <= 1) {
    b <- b_A; d <- d_A; r2 <- r2_A
  } else {
    b <- -eq2[1]/eq2[2]
    d <- -eq1[2]/eq1[1]
    r2 <- b * d
  }
  
  # 3. Вычисление базовых статистик
  r <- sign(b) * sqrt(r2)
  cov_xy <- b * var_x
  var_y <- cov_xy / d
  
  cat("--- ORIGINAL STATISTICS (N =", N, ") ---\n")
  cat("Means:          X =", mean_x, "| Y =", mean_y, "\n")
  cat("Correlation (r):", r, "\n")
  cat("Variance of Y:  ", var_y, "\n")
  cat("Covariance:     ", cov_xy, "\n\n")
  
  # 4. Пересчет при удалении точки
  if (!is.null(pt_remove)) {
    x_pt <- pt_remove[1]
    y_pt <- pt_remove[2]
    
    # Восстанавливаем оригинальные суммы
    sum_x <- mean_x * N
    sum_y <- mean_y * N
    sum_x2 <- (var_x + mean_x^2) * N
    sum_y2 <- (var_y + mean_y^2) * N
    sum_xy <- (cov_xy + mean_x*mean_y) * N
    
    # Вычитаем удаленную точку
    N_new <- N - 1
    sum_x_new <- sum_x - x_pt
    sum_y_new <- sum_y - y_pt
    sum_x2_new <- sum_x2 - x_pt^2
    sum_y2_new <- sum_y2 - y_pt^2
    sum_xy_new <- sum_xy - x_pt*y_pt
    
    # Считаем новые средние, дисперсии и ковариацию
    mx_new <- sum_x_new / N_new
    my_new <- sum_y_new / N_new
    vx_new <- (sum_x2_new / N_new) - mx_new^2
    vy_new <- (sum_y2_new / N_new) - my_new^2
    cov_new <- (sum_xy_new / N_new) - mx_new*my_new
    
    # Новые уравнения
    b_new <- cov_new / vx_new
    a_new <- my_new - b_new * mx_new
    
    d_new <- cov_new / vy_new
    c_new <- mx_new - d_new * my_new
    
    cat("--- NEW REGRESSION LINES (N =", N_new, ") ---\n")
    cat("New Y on X: y =", a_new, "+", b_new, "* x\n")
    cat("New X on Y: x =", c_new, "+", d_new, "* y\n")
  }
  
  return(invisible(list(Mean_X = mean_x, Mean_Y = mean_y, r = r, Var_Y = var_y, Cov_XY = cov_xy)))
}


compare_regression_lines <- function(x = NULL, y = NULL, var_x = NULL, var_y = NULL) {
  # Если переданы векторы данных, считаем дисперсии
  if (!is.null(x) && !is.null(y)) {
    var_x <- mean(x^2) - mean(x)^2
    var_y <- mean(y^2) - mean(y)^2
  }
  
  if (var_y < var_x) {
    verdict <- "R_Y/X is BETTER (lower MSE)"
  } else if (var_x < var_y) {
    verdict <- "R_X/Y is BETTER (lower MSE)"
  } else {
    verdict <- "Both lines have the SAME MSE"
  }
  
  return(list(
    var_X = var_x,
    var_Y = var_y,
    best_line = verdict
  ))
}


# Универсальная функция для восстановления моментов и статистик 
# Обязательные параметры: b (наклон Y on X), d (наклон X on Y)
# Все остальные параметры указывай только если они даны в задаче
calc_bivariate_moments <- function(b, d, 
                                   sd_x = NULL, sd_y = NULL, 
                                   var_x = NULL, var_y = NULL, 
                                   a = NULL, c = NULL,
                                   mean_x = NULL, mean_y = NULL) {
  
  # 1. Корреляция r
  r2 <- b * d
  r <- sign(b) * sqrt(r2)
  
  # 2. Восстанавливаем отклонения из дисперсий (если дали дисперсии)
  if (!is.null(var_x)) sd_x <- sqrt(var_x)
  if (!is.null(var_y)) sd_y <- sqrt(var_y)
  
  # 3. Находим недостающее отклонение через наклоны (b/d)
  if (!is.null(sd_x) && is.null(sd_y)) {
    sd_y <- sd_x * sqrt(b / d)
  } else if (!is.null(sd_y) && is.null(sd_x)) {
    sd_x <- sd_y * sqrt(d / b)
  }
  
  # 4. Вычисляем центральные моменты
  m20 <- sd_x^2       # Дисперсия X (var X)
  m02 <- sd_y^2       # Дисперсия Y (var Y)
  m11 <- b * m20      # Ковариация (можно и через d * m02, они равны)
  
  # 5. Восстанавливаем средние значения, если известны intercept (a или c)
  m_x <- mean_x
  m_y <- mean_y
  
  # Если дали среднее X и свободный член 'a' из уравнения y = a + bx
  if (!is.null(m_x) && is.null(m_y) && !is.null(a)) m_y <- a + b * m_x
  
  # Если дали среднее Y и свободный член 'c' из уравнения x = c + dy
  if (!is.null(m_y) && is.null(m_x) && !is.null(c)) m_x <- c + d * m_y
  
  return(list(
    r = r,
    sd_X = sd_x,
    sd_Y = sd_y,
    m20_var_X = m20,
    m02_var_Y = m02,
    m11_cov = m11,
    mean_X = m_x,
    mean_Y = m_y
  ))
}

#Золотое правило №1: Коэффициент детерминации $r^2 = b x d, где b — наклон линии $Y(X)$,
#а $d$ — наклон линии $X(Y)$. При этом $r^2$ всегда должен быть от 0 до 1




# Универсальная функция линейной регрессии (МНК)
# Y - вектор зависимой переменной
# X - матрица или вектор независимых переменных (предикторов)
# intercept = TRUE: автоматически добавляет столбец единиц (свободный член a0)
# intercept = FALSE: свободный член не добавляется (для моделей вида a*sin + b*cos)
fit_regression_model <- function(Y, X, intercept = TRUE) {
  n <- length(Y)
  X_mat <- as.matrix(X)
  
  # Формируем матрицу плана M
  if (intercept) {
    M <- cbind(Intercept = 1, X_mat)
  } else {
    M <- X_mat
  }
  
  # 1. Нормальные уравнения: (M^t * M) * A = M^t * Y
  MtM <- t(M) %*% M
  MtY <- t(M) %*% Y
  
  # 2. Оценка коэффициентов A
  A <- solve(MtM, MtY)
  
  # 3. Расчет ошибок и качества модели
  Y_est <- as.vector(M %*% A)
  errors <- Y - Y_est
  SSE <- sum(errors^2)
  MSE <- SSE / n
  var_Y <- mean(Y^2) - mean(Y)^2
  R2 <- 1 - (MSE / var_Y)
  
  # Красивый вывод нормальных уравнений в консоль
  cat("--- Normal Equations: (M^t * M) * A = M^t * Y ---\n")
  cat("Left Matrix (M^t * M):\n")
  print(round(MtM, 4))
  cat("\nRight Vector (M^t * Y):\n")
  print(round(MtY, 4))
  cat("-------------------------------------------------\n\n")
  
  coef_vec <- as.vector(A)
  names(coef_vec) <- colnames(M)
  
  return(list(
    Coefficients = coef_vec,
    SSE = SSE,
    MSE = MSE,
    Var_Y = var_Y,
    R_squared = R2,
    Fitted = Y_est
  ))
}



# Экспоненциальная регрессия: y = a * exp(b * x)
fit_exponential_model <- function(x, y) {
  n <- length(x)
  
  # 1. Линеаризация целевой переменной
  Y_ln <- log(y)
  
  # 2. Решение линейной системы (Y_ln = A + bx)
  M <- cbind(1, x)
  MtM <- t(M) %*% M
  MtY <- t(M) %*% Y_ln
  coeffs <- solve(MtM, MtY)
  
  A <- coeffs[1]
  b <- coeffs[2]
  
  # 3. Восстановление коэффициента 'a'
  a <- exp(A)
  
  # 4. Расчет прогнозов и ошибок для ИСХОДНЫХ (нелогарифмированных) данных
  y_est <- a * exp(b * x)
  errors <- y - y_est
  
  SSE <- sum(errors^2)
  MSE <- SSE / n
  var_y <- mean(y^2) - mean(y)^2
  R2 <- 1 - (MSE / var_y)
  
  cat("--- Exponential Fit: y = a * e^(bx) ---\n")
  cat("Linearized Eq: ln(y) =", round(A, 4), "+", round(b, 4), "* x\n")
  cat("Final Eq:         y =", round(a, 4), "* e^(", round(b, 4), "* x)\n\n")
  
  return(list(
    a = a,
    b = b,
    SSE = SSE,
    MSE = MSE,
    R_squared = R2,
    Fitted_y = y_est
  ))
}



# Гиперболическая регрессия: y = 1 / (a + b * x)
fit_hyperbolic_model <- function(x, y) {
  n <- length(x)
  
  # 1. Линеаризация целевой переменной
  Y_inv <- 1 / y
  
  # 2. Решение линейной системы (Y_inv = a + bx)
  M <- cbind(1, x)
  MtM <- t(M) %*% M
  MtY <- t(M) %*% Y_inv
  coeffs <- solve(MtM, MtY)
  
  a <- coeffs[1]
  b <- coeffs[2]
  
  # 3. Расчет прогнозов и ошибок для ИСХОДНЫХ (не перевернутых) данных
  y_est <- 1 / (a + b * x)
  errors <- y - y_est
  
  SSE <- sum(errors^2)
  MSE <- SSE / n
  var_y <- mean(y^2) - mean(y)^2
  R2 <- 1 - (MSE / var_y)
  
  cat("--- Hyperbolic Fit: y = 1 / (a + bx) ---\n")
  cat("Linearized Eq: 1/y =", round(a, 4), "+", round(b, 4), "* x\n")
  cat("Final Eq:        y = 1 / (", round(a, 4), round(b, 4), "* x )\n\n")
  
  return(list(
    a = a,
    b = b,
    SSE = SSE,
    MSE = MSE,
    R_squared = R2,
    Fitted_y = y_est
  ))
}


# Степенная регрессия: y = a * x^b
fit_power_model <- function(x, y) {
  n <- length(x)
  
  # 1. Линеаризация (логарифмируем И x, И y)
  X_ln <- log(x)
  Y_ln <- log(y)
  
  # 2. Решение линейной системы (Y_ln = A + b * X_ln)
  M <- cbind(1, X_ln)
  MtM <- t(M) %*% M
  MtY <- t(M) %*% Y_ln
  coeffs <- solve(MtM, MtY)
  
  A <- coeffs[1]
  b <- coeffs[2]
  
  # 3. Восстановление коэффициента 'a'
  a <- exp(A)
  
  # 4. Расчет прогнозов и ошибок для ИСХОДНЫХ данных
  y_est <- a * (x^b)
  errors <- y - y_est
  
  SSE <- sum(errors^2)
  MSE <- SSE / n
  var_y <- mean(y^2) - mean(y)^2
  R2 <- 1 - (MSE / var_y)
  
  cat("--- Power Fit: y = a * x^b ---\n")
  cat("Linearized Eq: ln(y) =", round(A, 4), "+", round(b, 4), "* ln(x)\n")
  cat("Final Eq:         y =", round(a, 4), "* x ^", round(b, 4), "\n\n")
  
  return(list(
    a = a,
    b = b,
    SSE = SSE,
    MSE = MSE,
    R_squared = R2,
    Fitted_y = y_est
  ))
}