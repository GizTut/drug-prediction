# ==============================================================================
# 03_run_nested_cv.R — 6 model × 7 ilaç × 3 gen sayısı nested CV
# ==============================================================================
# Bu EN UZUN script. Bir kez çalıştır, sonuçları CSV'ye düşür.
#
# Yapacağı:
#   - Her ilaç için, her gen sayısı için (50/100/150):
#     - 6 model (KNN, ENet, NNet, SVM, XGBoost, RF) için nested 5x5 CV
#     - Outer fold metrikleri (RMSE, MAE, R²)
#     - Inner fold'da hiperparametre tuning
#   - Sonuçlar her ilaç klasöründeki results/ altına yazılır
#
# Bağımlı: 00_setup.R + 02_preprocess_per_drug.R çıktıları
# Tahmini süre: ~9-13 saat (7 ilaç × 1.5 saat)
#
# ⚠️  Mac'te uyku modunu kapat: System Settings → Battery → Prevent sleep
# ==============================================================================
setwd(SCRIPTS_PATH)

source("00_setup.R")
cat("\n[03_run_nested_cv] Nested CV başlıyor (UZUN — 9-13 saat)...\n")


# ==============================================================================
# NESTED CV FONKSİYONLARI
# ==============================================================================

# ---- KNN ----
nested_knn <- function(train_file, outer_folds = OUTER_FOLDS,
                       inner_folds = INNER_FOLDS, seed = GLOBAL_SEED) {
  set.seed(seed)
  df <- read.csv(train_file, row.names = 1, check.names = FALSE)
  df$IC50 <- as.numeric(df$IC50)
  outer_idx <- createFolds(df$IC50, k = outer_folds, returnTrain = FALSE)
  tune_grid <- expand.grid(
    kmax = seq(3, 35, by = 2), distance = c(1, 2),
    kernel = c("rectangular", "triangular", "epanechnikov", "gaussian")
  )
  outer_results <- list(); chosen_each <- list()
  for (i in seq_along(outer_idx)) {
    test_ids  <- outer_idx[[i]]
    train_ids <- setdiff(seq_len(nrow(df)), test_ids)
    ot <- df[train_ids, , drop = FALSE]
    oe <- df[test_ids,  , drop = FALSE]
    ctrl <- trainControl(method = "cv", number = inner_folds)
    set.seed(seed + i)
    fit <- caret::train(IC50 ~ ., data = ot, method = "kknn",
                 trControl = ctrl, tuneGrid = tune_grid, metric = "RMSE")
    chosen_each[[i]] <- fit$bestTune
    p <- predict(fit, newdata = oe); y <- oe$IC50
    outer_results[[i]] <- data.frame(
      Fold = i, MAE = mae(y, p), RMSE = rmse(y, p), R2 = cor(y, p)^2
    )
  }
  outer_summary <- bind_rows(outer_results)
  chosen_df <- bind_rows(chosen_each)
  final_params <- chosen_df %>% dplyr::count(kmax, distance, kernel, sort = TRUE) %>%
    dplyr::slice(1) %>% dplyr::select(kmax, distance, kernel)
  list(mean_metrics = c(MAE  = mean(outer_summary$MAE,  na.rm = TRUE),
                        RMSE = mean(outer_summary$RMSE, na.rm = TRUE),
                        R2   = mean(outer_summary$R2,   na.rm = TRUE)),
       outer_fold_metrics = outer_summary,
       final_params = final_params)
}

# ---- Elastic Net ----
nested_enet <- function(train_file, outer_folds = OUTER_FOLDS,
                        inner_folds = INNER_FOLDS, seed = GLOBAL_SEED) {
  set.seed(seed)
  df <- read.csv(train_file, row.names = 1, check.names = FALSE)
  df$IC50 <- as.numeric(df$IC50)
  outer_idx <- createFolds(df$IC50, k = outer_folds, returnTrain = FALSE)
  tune_grid <- expand.grid(lambda = seq(0, 1, by = 0.1),
                           fraction = seq(0.1, 1, by = 0.1))
  outer_results <- list(); chosen_each <- list()
  for (i in seq_along(outer_idx)) {
    test_ids  <- outer_idx[[i]]
    train_ids <- setdiff(seq_len(nrow(df)), test_ids)
    ot <- df[train_ids, , drop = FALSE]
    oe <- df[test_ids,  , drop = FALSE]
    ctrl <- trainControl(method = "cv", number = inner_folds)
    set.seed(seed + i)
    fit <- caret::train(IC50 ~ ., data = ot, method = "enet",
                 trControl = ctrl, tuneGrid = tune_grid, metric = "RMSE")
    chosen_each[[i]] <- fit$bestTune
    p <- predict(fit, newdata = oe); y <- oe$IC50
    r2 <- suppressWarnings(cor(y, p, use = "complete.obs")^2)
    if (is.na(r2) || is.nan(r2)) r2 <- NA_real_
    outer_results[[i]] <- data.frame(
      Fold = i, MAE = mae(y, p), RMSE = rmse(y, p), R2 = r2
    )
  }
  outer_summary <- bind_rows(outer_results)
  chosen_df <- bind_rows(chosen_each)
  final_params <- chosen_df %>% dplyr::count(lambda, fraction, sort = TRUE) %>%
    dplyr::slice(1) %>% dplyr::select(lambda, fraction)
  list(mean_metrics = c(MAE  = mean(outer_summary$MAE,  na.rm = TRUE),
                        RMSE = mean(outer_summary$RMSE, na.rm = TRUE),
                        R2   = mean(outer_summary$R2,   na.rm = TRUE)),
       outer_fold_metrics = outer_summary,
       final_params = final_params)
}

# ---- Neural Net ----
nested_nnet <- function(train_file, outer_folds = OUTER_FOLDS,
                        inner_folds = INNER_FOLDS, seed = GLOBAL_SEED) {
  set.seed(seed)
  df <- read.csv(train_file, row.names = 1, check.names = FALSE)
  df$IC50 <- as.numeric(df$IC50)
  outer_idx <- createFolds(df$IC50, k = outer_folds, returnTrain = FALSE)
  tune_grid <- expand.grid(size = c(1, 3, 5), decay = c(0, 0.001, 0.01))
  outer_results <- list(); chosen_each <- list()
  for (i in seq_along(outer_idx)) {
    test_ids  <- outer_idx[[i]]
    train_ids <- setdiff(seq_len(nrow(df)), test_ids)
    ot <- df[train_ids, , drop = FALSE]
    oe <- df[test_ids,  , drop = FALSE]
    ctrl <- trainControl(method = "cv", number = inner_folds)
    set.seed(seed + i)
    fit <- caret::train(IC50 ~ ., data = ot, method = "nnet",
                 trControl = ctrl, tuneGrid = tune_grid,
                 linout = TRUE, trace = FALSE, maxit = 500, metric = "RMSE")
    chosen_each[[i]] <- fit$bestTune
    p <- predict(fit, newdata = oe); y <- oe$IC50
    r2 <- suppressWarnings(cor(y, p, use = "complete.obs")^2)
    if (is.na(r2) || is.nan(r2)) r2 <- NA_real_
    outer_results[[i]] <- data.frame(
      Fold = i, MAE = mae(y, p), RMSE = rmse(y, p), R2 = r2
    )
  }
  outer_summary <- bind_rows(outer_results)
  chosen_df <- bind_rows(chosen_each)
  final_params <- chosen_df %>% dplyr::count(size, decay, sort = TRUE) %>%
    dplyr::slice(1) %>% dplyr::select(size, decay)
  list(mean_metrics = c(MAE  = mean(outer_summary$MAE,  na.rm = TRUE),
                        RMSE = mean(outer_summary$RMSE, na.rm = TRUE),
                        R2   = mean(outer_summary$R2,   na.rm = TRUE)),
       outer_fold_metrics = outer_summary,
       final_params = final_params)
}

# ---- SVM ----
nested_svm <- function(train_file, outer_folds = OUTER_FOLDS,
                       inner_folds = INNER_FOLDS, seed = GLOBAL_SEED) {
  set.seed(seed)
  df <- read.csv(train_file, row.names = 1, check.names = FALSE)
  df$IC50 <- as.numeric(df$IC50)
  outer_idx <- createFolds(df$IC50, k = outer_folds, returnTrain = FALSE)
  tune_grid <- expand.grid(C = 2^(-2:5), sigma = 2^(-5:2))
  outer_results <- list(); chosen_each <- list()
  for (i in seq_along(outer_idx)) {
    test_ids  <- outer_idx[[i]]
    train_ids <- setdiff(seq_len(nrow(df)), test_ids)
    ot <- df[train_ids, , drop = FALSE]
    oe <- df[test_ids,  , drop = FALSE]
    ctrl <- trainControl(method = "cv", number = inner_folds)
    set.seed(seed + i)
    fit <- caret::train(IC50 ~ ., data = ot, method = "svmRadial",
                 trControl = ctrl, tuneGrid = tune_grid, metric = "RMSE")
    chosen_each[[i]] <- fit$bestTune
    p <- predict(fit, newdata = oe); y <- oe$IC50
    r2 <- suppressWarnings(cor(y, p, use = "complete.obs")^2)
    if (is.na(r2) || is.nan(r2)) r2 <- NA_real_
    outer_results[[i]] <- data.frame(
      Fold = i, MAE = mae(y, p), RMSE = rmse(y, p), R2 = r2
    )
  }
  outer_summary <- bind_rows(outer_results)
  chosen_df <- bind_rows(chosen_each)
  final_params <- chosen_df %>% dplyr::count(C, sigma, sort = TRUE) %>%
    dplyr::slice(1) %>% dplyr::select(C, sigma)
  list(mean_metrics = c(MAE  = mean(outer_summary$MAE,  na.rm = TRUE),
                        RMSE = mean(outer_summary$RMSE, na.rm = TRUE),
                        R2   = mean(outer_summary$R2,   na.rm = TRUE)),
       outer_fold_metrics = outer_summary,
       final_params = final_params)
}

# ---- XGBoost ----
nested_xgb <- function(train_file, outer_folds = OUTER_FOLDS,
                       inner_folds = INNER_FOLDS, seed = GLOBAL_SEED,
                       nrounds_max = 2000, early_stop = 30) {
  set.seed(seed)
  df <- read.csv(train_file, row.names = 1, check.names = FALSE)
  df$IC50 <- as.numeric(df$IC50)
  X_all <- as.matrix(df[, setdiff(colnames(df), "IC50"), drop = FALSE])
  y_all <- df$IC50
  outer_idx <- createFolds(y_all, k = outer_folds, returnTrain = FALSE)
  param_grid <- expand.grid(eta = c(0.05, 0.1, 0.3),
                            max_depth = c(3, 5, 7),
                            subsample = c(0.7, 1.0),
                            colsample_bytree = c(0.7, 1.0),
                            min_child_weight = c(1, 3),
                            gamma = c(0, 0.1))
  outer_results <- list(); chosen_each <- list()
  for (i in seq_along(outer_idx)) {
    test_ids  <- outer_idx[[i]]
    train_ids <- setdiff(seq_len(nrow(df)), test_ids)
    X_tr <- X_all[train_ids, , drop = FALSE]; y_tr <- y_all[train_ids]
    X_te <- X_all[test_ids,  , drop = FALSE]; y_te <- y_all[test_ids]
    dtr <- xgb.DMatrix(data = X_tr, label = y_tr)
    dte <- xgb.DMatrix(data = X_te, label = y_te)
    best_rmse <- Inf; best_params <- NULL; best_nrounds <- NULL
    for (r in seq_len(nrow(param_grid))) {
      params <- list(objective = "reg:squarederror", eval_metric = "rmse",
                     eta = param_grid$eta[r],
                     max_depth = param_grid$max_depth[r],
                     subsample = param_grid$subsample[r],
                     colsample_bytree = param_grid$colsample_bytree[r],
                     min_child_weight = param_grid$min_child_weight[r],
                     gamma = param_grid$gamma[r])
      set.seed(seed + i + r)
      cv <- tryCatch(xgb.cv(params = params, data = dtr, nrounds = nrounds_max,
                            nfold = inner_folds, early_stopping_rounds = early_stop,
                            verbose = 0), error = function(e) NULL)
      if (is.null(cv) || is.null(cv$evaluation_log)) next
      rv <- cv$evaluation_log$test_rmse_mean
      if (length(rv) == 0 || all(is.na(rv))) next
      rmse_min <- min(rv, na.rm = TRUE)
      nb <- which.min(rv)
      if (length(nb) == 0 || nb < 1) next
      if (!is.finite(rmse_min)) next
      if (is.null(best_params) || rmse_min < best_rmse) {
        best_rmse <- rmse_min; best_params <- params; best_nrounds <- nb
      }
    }
    if (is.null(best_params)) {
      stop(sprintf("Outer fold %d: hiçbir XGB parametresi geçerli değil.", i))
    }
    chosen_each[[i]] <- data.frame(
      Fold = i, eta = best_params$eta, max_depth = best_params$max_depth,
      subsample = best_params$subsample,
      colsample_bytree = best_params$colsample_bytree,
      min_child_weight = best_params$min_child_weight,
      gamma = best_params$gamma, best_nrounds = best_nrounds,
      inner_best_rmse = best_rmse
    )
    set.seed(seed + 1000 + i)
    model <- xgb.train(params = best_params, data = dtr,
                       nrounds = best_nrounds, verbose = 0)
    p <- predict(model, dte, iteration_range = c(1, best_nrounds))
    r2 <- suppressWarnings(cor(y_te, p, use = "complete.obs")^2)
    if (is.na(r2) || is.nan(r2)) r2 <- NA_real_
    outer_results[[i]] <- data.frame(
      Fold = i, MAE = mae(y_te, p), RMSE = rmse(y_te, p), R2 = r2
    )
  }
  outer_summary <- bind_rows(outer_results)
  chosen_df <- bind_rows(chosen_each)
  fp_row <- chosen_df %>%
    dplyr::count(eta, max_depth, subsample, colsample_bytree,
                 min_child_weight, gamma, sort = TRUE) %>%
    dplyr::slice(1)
  # Final nrounds için tüm veride CV
  fp_final <- list(objective = "reg:squarederror", eval_metric = "rmse",
                   eta = fp_row$eta, max_depth = fp_row$max_depth,
                   subsample = fp_row$subsample,
                   colsample_bytree = fp_row$colsample_bytree,
                   min_child_weight = fp_row$min_child_weight,
                   gamma = fp_row$gamma)
  dfull <- xgb.DMatrix(data = X_all, label = y_all)
  set.seed(seed + 9999)
  cv_full <- xgb.cv(params = fp_final, data = dfull, nrounds = nrounds_max,
                    nfold = inner_folds, early_stopping_rounds = early_stop,
                    verbose = 0)
  final_nrounds <- which.min(cv_full$evaluation_log$test_rmse_mean)
  final_params <- data.frame(
    eta = fp_row$eta, max_depth = fp_row$max_depth,
    subsample = fp_row$subsample,
    colsample_bytree = fp_row$colsample_bytree,
    min_child_weight = fp_row$min_child_weight,
    gamma = fp_row$gamma, nrounds = final_nrounds
  )
  list(mean_metrics = c(MAE  = mean(outer_summary$MAE,  na.rm = TRUE),
                        RMSE = mean(outer_summary$RMSE, na.rm = TRUE),
                        R2   = mean(outer_summary$R2,   na.rm = TRUE)),
       outer_fold_metrics = outer_summary,
       final_params = final_params)
}

# ---- Random Forest ----
nested_rf <- function(train_file, outer_folds = OUTER_FOLDS,
                      inner_folds = INNER_FOLDS, seed = GLOBAL_SEED) {
  set.seed(seed)
  df <- read.csv(train_file, row.names = 1, check.names = FALSE)
  df$IC50 <- as.numeric(df$IC50)
  outer_idx <- createFolds(df$IC50, k = outer_folds, returnTrain = FALSE)
  n_features <- ncol(df) - 1
  mtry_grid <- unique(round(c(sqrt(n_features), n_features / 3,
                              n_features / 5, n_features / 2)))
  mtry_grid <- mtry_grid[mtry_grid >= 1 & mtry_grid <= n_features]
  tune_grid <- data.frame(mtry = mtry_grid)
  outer_results <- list(); chosen_each <- list()
  for (i in seq_along(outer_idx)) {
    test_ids  <- outer_idx[[i]]
    train_ids <- setdiff(seq_len(nrow(df)), test_ids)
    ot <- df[train_ids, , drop = FALSE]
    oe <- df[test_ids,  , drop = FALSE]
    ctrl <- trainControl(method = "cv", number = inner_folds)
    set.seed(seed + i)
    fit <- caret::train(IC50 ~ ., data = ot, method = "rf",
                 trControl = ctrl, tuneGrid = tune_grid,
                 metric = "RMSE", ntree = 500)
    chosen_each[[i]] <- fit$bestTune
    p <- predict(fit, newdata = oe); y <- oe$IC50
    r2 <- suppressWarnings(cor(y, p, use = "complete.obs")^2)
    if (is.na(r2) || is.nan(r2)) r2 <- NA_real_
    outer_results[[i]] <- data.frame(
      Fold = i, MAE = mae(y, p), RMSE = rmse(y, p), R2 = r2
    )
  }
  outer_summary <- bind_rows(outer_results)
  chosen_df <- bind_rows(chosen_each)
  final_params <- chosen_df %>% dplyr::count(mtry, sort = TRUE) %>%
    dplyr::slice(1) %>% dplyr::select(mtry)
  list(mean_metrics = c(MAE  = mean(outer_summary$MAE,  na.rm = TRUE),
                        RMSE = mean(outer_summary$RMSE, na.rm = TRUE),
                        R2   = mean(outer_summary$R2,   na.rm = TRUE)),
       outer_fold_metrics = outer_summary,
       final_params = final_params)
}


# ==============================================================================
# 7 İLAÇ İÇİN NESTED CV ÇALIŞTIR
# ==============================================================================

run_one_model <- function(model_name, nested_fn, train_files) {
  cat("    -", model_name, "\n")
  details <- list(); summaries <- list()
  for (file in train_files) {
    res <- nested_fn(file)
    n_genes <- gsub("data_train_|\\.csv", "", basename(file))
    details[[paste0(model_name, "_", n_genes)]] <- data.frame(
      Model = model_name, Genes = as.numeric(n_genes),
      Fold = res$outer_fold_metrics$Fold,
      MAE  = res$outer_fold_metrics$MAE,
      RMSE = res$outer_fold_metrics$RMSE,
      R2   = res$outer_fold_metrics$R2
    )
    fp <- res$final_params
    summary_row <- data.frame(File = basename(file), Genes = as.numeric(n_genes),
                              MAE = res$mean_metrics["MAE"],
                              RMSE = res$mean_metrics["RMSE"],
                              R2 = res$mean_metrics["R2"])
    if (is.data.frame(fp)) summary_row <- cbind(summary_row, fp)
    summaries[[basename(file)]] <- summary_row
  }
  list(details = bind_rows(details), summary = bind_rows(summaries))
}

run_nested_cv_for_drug <- function(drug_folder) {
  cat("\n", strrep("=", 60), "\n", sep = "")
  cat("İLAÇ: ", drug_folder, "\n", sep = "")
  cat(strrep("=", 60), "\n", sep = "")
  start_time <- Sys.time()

  drug_path    <- paste0(BASE_PATH, drug_folder, "/")
  results_path <- paste0(drug_path, "results/")
  dir.create(results_path, showWarnings = FALSE, recursive = TRUE)

  train_files <- c(paste0(drug_path, "data_train_50.csv"),
                   paste0(drug_path, "data_train_100.csv"),
                   paste0(drug_path, "data_train_150.csv"))

  if (!all(file.exists(train_files))) {
    cat("⚠️  Train dosyaları yok, bu ilacı atlıyorum.\n"); return(NULL)
  }

  cat("Nested CV (6 model × 3 gen sayısı)...\n")
  knn_run  <- run_one_model("KNN",  nested_knn,  train_files)
  enet_run <- run_one_model("ENET", nested_enet, train_files)
  nnet_run <- run_one_model("NNET", nested_nnet, train_files)
  svm_run  <- run_one_model("SVM",  nested_svm,  train_files)
  xgb_run  <- run_one_model("XGB",  nested_xgb,  train_files)
  rf_run   <- run_one_model("RF",   nested_rf,   train_files)

  fold_details_df <- bind_rows(knn_run$details, enet_run$details,
                               nnet_run$details, svm_run$details,
                               xgb_run$details, rf_run$details)

  write.csv(knn_run$summary,  paste0(results_path, "knn_results.csv"),  row.names = FALSE)
  write.csv(enet_run$summary, paste0(results_path, "enet_results.csv"), row.names = FALSE)
  write.csv(nnet_run$summary, paste0(results_path, "nnet_results.csv"), row.names = FALSE)
  write.csv(svm_run$summary,  paste0(results_path, "svm_results.csv"),  row.names = FALSE)
  write.csv(xgb_run$summary,  paste0(results_path, "xgb_results.csv"),  row.names = FALSE)
  write.csv(rf_run$summary,   paste0(results_path, "rf_results.csv"),   row.names = FALSE)
  write.csv(fold_details_df,
            paste0(results_path, "fold_details_all_models.csv"),
            row.names = FALSE)

  end_time <- Sys.time()
  duration <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 1)
  cat("✅", drug_folder, "tamamlandı (", duration, "dk)\n")
}


# ==============================================================================
# ANA LOOP
# ==============================================================================

total_start <- Sys.time()

for (drug in DRUG_FOLDER_NAMES) {
  tryCatch(run_nested_cv_for_drug(drug),
           error = function(e) {
             cat("⚠️ HATA:", drug, "-", conditionMessage(e), "\n")
           })
}

total_end <- Sys.time()
total_h <- round(as.numeric(difftime(total_end, total_start, units = "hours")), 2)

cat("\n", strrep("=", 60), "\n", sep = "")
cat("✅ 03_run_nested_cv.R TAMAMLANDI\n")
cat("Toplam süre: ", total_h, " saat\n", sep = "")
cat(strrep("=", 60), "\n", sep = "")
cat("\n📌 Sonraki adım: 04_statistical_tests.R\n")
