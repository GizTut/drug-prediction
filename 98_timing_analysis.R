# ==============================================================================
# 98_timing_analysis.R — Modellerin eğitim ve tahmin sürelerini ölç
# ==============================================================================
# YAYIN İÇİN ÖNEMLİ:
# Hakem 3.7 "hesaplama karmaşıklığı tartışılmamış" demişti. Bu script,
# Methods/Results bölümünde raporlanabilecek süre ölçümlerini üretir.
#
# Yapacağı:
#   - Her ilaç için 1-SE kararıyla seçilmiş 6 final modeli yeniden eğit
#   - Eğitim süresini mikrosaniye hassasiyetinde kaydet (microbenchmark)
#   - Test seti üzerinde tahmin süresini ölç
#   - Ortalama + standart sapma ile sonuç tablosu üret
#   - Sistem bilgisi (CPU, RAM, OS) kaydet
#
# Bağımlı: 00_setup.R + 02-05 çıktıları (özellikle final_gene_selection.csv)
# Tahmini süre: ~15 dakika
#
# KULLANIM (pipeline TAMAMEN bittikten sonra):
#   setwd("/Users/gizemtutkun/Desktop/Gizem_yl/scripts/")
#   source("98_timing_analysis.R")
# ==============================================================================

source("00_setup.R")
cat("\n[98_timing_analysis] Süre ölçümleri başlıyor...\n")


# ==============================================================================
# 0. SİSTEM BİLGİSİ KAYDET (Methods için)
# ==============================================================================

cat("\n--- 0. System Information ---\n")

sys_info <- list(
  R_version       = R.version.string,
  OS              = Sys.info()["sysname"],
  OS_version      = Sys.info()["release"],
  Machine         = Sys.info()["machine"],
  Cores_logical   = parallel::detectCores(logical = TRUE),
  Cores_physical  = parallel::detectCores(logical = FALSE),
  RAM_GB          = round(as.numeric(benchmarkme::get_ram()) / 1e9, 1)
)

# benchmarkme yoksa RAM bilgisi alamayız - güvenli mod
sys_info$RAM_GB <- tryCatch({
  if (!requireNamespace("benchmarkme", quietly = TRUE)) {
    install.packages("benchmarkme", quiet = TRUE)
  }
  round(as.numeric(benchmarkme::get_ram()) / 1e9, 1)
}, error = function(e) NA)

cat("  R version: ", sys_info$R_version, "\n", sep = "")
cat("  OS: ", sys_info$OS, " ", sys_info$OS_version, "\n", sep = "")
cat("  Architecture: ", sys_info$Machine, "\n", sep = "")
cat("  Cores: ", sys_info$Cores_logical, " logical (", sys_info$Cores_physical, " physical)\n", sep = "")
cat("  RAM: ", sys_info$RAM_GB, " GB\n", sep = "")

# Paket sürümleri (Methods için kritik)
pkg_versions <- data.frame(
  Package = c("caret", "kknn", "elasticnet", "nnet", "kernlab",
              "xgboost", "randomForest", "org.Hs.eg.db"),
  Version = sapply(c("caret", "kknn", "elasticnet", "nnet", "kernlab",
                     "xgboost", "randomForest", "org.Hs.eg.db"),
                   function(p) as.character(packageVersion(p)))
)

cat("\nPackage versions:\n")
print(pkg_versions, row.names = FALSE)


# ==============================================================================
# 1. SÜRE ÖLÇÜM FONKSİYONU
# ==============================================================================

# Modeli n_reps kez eğit, ortalama süreyi al (varyans için)
time_train <- function(train_fn, n_reps = 3) {
  times <- numeric(n_reps)
  for (i in 1:n_reps) {
    t0 <- Sys.time()
    train_fn()
    times[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  list(mean_sec = mean(times), sd_sec = sd(times), times = times)
}

time_predict <- function(predict_fn, n_reps = 5) {
  times <- numeric(n_reps)
  for (i in 1:n_reps) {
    t0 <- Sys.time()
    predict_fn()
    times[i] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  }
  list(mean_sec = mean(times), sd_sec = sd(times), times = times)
}


# ==============================================================================
# 2. TEK İLAÇ İÇİN 6 MODEL ZAMANLA
# ==============================================================================

time_drug <- function(drug_folder) {
  cat("\n", strrep("-", 60), "\n", sep = "")
  cat("İLAÇ: ", drug_folder, "\n", sep = "")
  cat(strrep("-", 60), "\n", sep = "")

  drug_path    <- paste0(BASE_PATH, drug_folder, "/")
  results_path <- paste0(drug_path, "results/")

  # 1-SE kararı yükle
  decision_df <- read.csv(paste0(results_path, "final_gene_selection.csv"))

  # Her model için summary'leri yükle
  knn_s  <- read.csv(paste0(results_path, "knn_results.csv"))
  enet_s <- read.csv(paste0(results_path, "enet_results.csv"))
  nnet_s <- read.csv(paste0(results_path, "nnet_results.csv"))
  svm_s  <- read.csv(paste0(results_path, "svm_results.csv"))
  xgb_s  <- read.csv(paste0(results_path, "xgb_results.csv"))
  rf_s   <- read.csv(paste0(results_path, "rf_results.csv"))

  # Yardımcı: 1-SE'nin seçtiği gen sayısındaki train+test datasını yükle
  get_data <- function(model_name) {
    n <- decision_df$Final_Choice[decision_df$Model == model_name]
    train <- read.csv(paste0(drug_path, "data_train_", n, ".csv"),
                      row.names = 1, check.names = FALSE)
    test  <- read.csv(paste0(drug_path, "data_test_",  n, ".csv"),
                      row.names = 1, check.names = FALSE)
    train$IC50 <- as.numeric(train$IC50)
    test$IC50  <- as.numeric(test$IC50)
    list(train = train, test = test, n_genes = n)
  }

  ctrl_none <- trainControl(method = "none")
  results_drug <- list()

  # ---- KNN ----
  cat("  Timing KNN...\n")
  d <- get_data("KNN"); p <- knn_s[knn_s$Genes == d$n_genes, ]
  train_fn <- function() {
    set.seed(GLOBAL_SEED)
    caret::train(IC50 ~ ., data = d$train, method = "kknn", trControl = ctrl_none,
          tuneGrid = data.frame(kmax = p$kmax, distance = p$distance, kernel = p$kernel))
  }
  m <- train_fn()
  predict_fn <- function() predict(m, newdata = d$test)
  tt <- time_train(train_fn, 3); pt <- time_predict(predict_fn, 5)
  results_drug[["KNN"]] <- list(train_mean = tt$mean_sec, train_sd = tt$sd_sec,
                                predict_mean = pt$mean_sec, predict_sd = pt$sd_sec,
                                n_genes = d$n_genes, n_train = nrow(d$train), n_test = nrow(d$test))
  cat(sprintf("    train: %.3fs ± %.3fs   |   predict: %.4fs ± %.4fs\n",
              tt$mean_sec, tt$sd_sec, pt$mean_sec, pt$sd_sec))

  # ---- ENET ----
  cat("  Timing Elastic Net...\n")
  d <- get_data("ENET"); p <- enet_s[enet_s$Genes == d$n_genes, ]
  train_fn <- function() {
    set.seed(GLOBAL_SEED)
    caret::train(IC50 ~ ., data = d$train, method = "enet", trControl = ctrl_none,
          tuneGrid = data.frame(lambda = p$lambda, fraction = p$fraction))
  }
  m <- train_fn()
  predict_fn <- function() predict(m, newdata = d$test)
  tt <- time_train(train_fn, 3); pt <- time_predict(predict_fn, 5)
  results_drug[["Elastic Net"]] <- list(train_mean = tt$mean_sec, train_sd = tt$sd_sec,
                                        predict_mean = pt$mean_sec, predict_sd = pt$sd_sec,
                                        n_genes = d$n_genes, n_train = nrow(d$train), n_test = nrow(d$test))
  cat(sprintf("    train: %.3fs ± %.3fs   |   predict: %.4fs ± %.4fs\n",
              tt$mean_sec, tt$sd_sec, pt$mean_sec, pt$sd_sec))

  # ---- NNET ----
  cat("  Timing Neural Net...\n")
  d <- get_data("NNET"); p <- nnet_s[nnet_s$Genes == d$n_genes, ]
  train_fn <- function() {
    set.seed(GLOBAL_SEED)
    caret::train(IC50 ~ ., data = d$train, method = "nnet", trControl = ctrl_none,
          tuneGrid = data.frame(size = p$size, decay = p$decay),
          linout = TRUE, trace = FALSE, maxit = 500)
  }
  m <- train_fn()
  predict_fn <- function() predict(m, newdata = d$test)
  tt <- time_train(train_fn, 3); pt <- time_predict(predict_fn, 5)
  results_drug[["Neural Net"]] <- list(train_mean = tt$mean_sec, train_sd = tt$sd_sec,
                                       predict_mean = pt$mean_sec, predict_sd = pt$sd_sec,
                                       n_genes = d$n_genes, n_train = nrow(d$train), n_test = nrow(d$test))
  cat(sprintf("    train: %.3fs ± %.3fs   |   predict: %.4fs ± %.4fs\n",
              tt$mean_sec, tt$sd_sec, pt$mean_sec, pt$sd_sec))

  # ---- SVM ----
  cat("  Timing SVM...\n")
  d <- get_data("SVM"); p <- svm_s[svm_s$Genes == d$n_genes, ]
  train_fn <- function() {
    set.seed(GLOBAL_SEED)
    caret::train(IC50 ~ ., data = d$train, method = "svmRadial", trControl = ctrl_none,
          tuneGrid = data.frame(C = p$C, sigma = p$sigma))
  }
  m <- train_fn()
  predict_fn <- function() predict(m, newdata = d$test)
  tt <- time_train(train_fn, 3); pt <- time_predict(predict_fn, 5)
  results_drug[["SVM"]] <- list(train_mean = tt$mean_sec, train_sd = tt$sd_sec,
                                predict_mean = pt$mean_sec, predict_sd = pt$sd_sec,
                                n_genes = d$n_genes, n_train = nrow(d$train), n_test = nrow(d$test))
  cat(sprintf("    train: %.3fs ± %.3fs   |   predict: %.4fs ± %.4fs\n",
              tt$mean_sec, tt$sd_sec, pt$mean_sec, pt$sd_sec))

  # ---- XGBoost ----
  cat("  Timing XGBoost...\n")
  d <- get_data("XGB"); p <- xgb_s[xgb_s$Genes == d$n_genes, ]
  X_train <- as.matrix(d$train[, setdiff(colnames(d$train), "IC50")])
  y_train <- d$train$IC50
  X_test  <- as.matrix(d$test[,  setdiff(colnames(d$test),  "IC50")])
  dtr <- xgb.DMatrix(data = X_train, label = y_train)
  dte <- xgb.DMatrix(data = X_test)
  xgb_params <- list(objective = "reg:squarederror", eval_metric = "rmse",
                     eta = p$eta, max_depth = p$max_depth,
                     subsample = p$subsample, colsample_bytree = p$colsample_bytree,
                     min_child_weight = p$min_child_weight, gamma = p$gamma)
  train_fn <- function() {
    set.seed(GLOBAL_SEED)
    xgb.train(params = xgb_params, data = dtr, nrounds = p$nrounds, verbose = 0)
  }
  m <- train_fn()
  predict_fn <- function() predict(m, dte)
  tt <- time_train(train_fn, 3); pt <- time_predict(predict_fn, 5)
  results_drug[["XGBoost"]] <- list(train_mean = tt$mean_sec, train_sd = tt$sd_sec,
                                    predict_mean = pt$mean_sec, predict_sd = pt$sd_sec,
                                    n_genes = d$n_genes, n_train = nrow(d$train), n_test = nrow(d$test))
  cat(sprintf("    train: %.3fs ± %.3fs   |   predict: %.4fs ± %.4fs\n",
              tt$mean_sec, tt$sd_sec, pt$mean_sec, pt$sd_sec))

  # ---- Random Forest ----
  cat("  Timing Random Forest...\n")
  d <- get_data("RF"); p <- rf_s[rf_s$Genes == d$n_genes, ]
  train_fn <- function() {
    set.seed(GLOBAL_SEED)
    caret::train(IC50 ~ ., data = d$train, method = "rf", trControl = ctrl_none,
          tuneGrid = data.frame(mtry = p$mtry), ntree = 500)
  }
  m <- train_fn()
  predict_fn <- function() predict(m, newdata = d$test)
  tt <- time_train(train_fn, 3); pt <- time_predict(predict_fn, 5)
  results_drug[["Random Forest"]] <- list(train_mean = tt$mean_sec, train_sd = tt$sd_sec,
                                          predict_mean = pt$mean_sec, predict_sd = pt$sd_sec,
                                          n_genes = d$n_genes, n_train = nrow(d$train), n_test = nrow(d$test))
  cat(sprintf("    train: %.3fs ± %.3fs   |   predict: %.4fs ± %.4fs\n",
              tt$mean_sec, tt$sd_sec, pt$mean_sec, pt$sd_sec))

  # Hafıza temizle
  rm(m); gc(verbose = FALSE)

  # data.frame'e çevir
  rows <- list()
  for (model_name in names(results_drug)) {
    r <- results_drug[[model_name]]
    rows[[model_name]] <- data.frame(
      Drug = drug_folder, Model = model_name,
      N_genes = r$n_genes, N_train = r$n_train, N_test = r$n_test,
      Train_mean_sec = round(r$train_mean, 4), Train_sd_sec = round(r$train_sd, 4),
      Predict_mean_sec = round(r$predict_mean, 5), Predict_sd_sec = round(r$predict_sd, 5)
    )
  }
  bind_rows(rows)
}


# ==============================================================================
# 3. 7 İLAÇ İÇİN ZAMANLA
# ==============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("Süre ölçümleri 7 ilaç için başlıyor (~15 dk)\n")
cat(strrep("=", 60), "\n", sep = "")

all_timings <- list()
for (drug in DRUG_FOLDER_NAMES) {
  res <- tryCatch(time_drug(drug),
                  error = function(e) {
                    cat("⚠️", drug, "- HATA:", conditionMessage(e), "\n")
                    NULL
                  })
  if (!is.null(res)) all_timings[[drug]] <- res
}

timings_df <- bind_rows(all_timings)
write.csv(timings_df, paste0(BASE_PATH, "model_timings.csv"), row.names = FALSE)


# ==============================================================================
# 4. ÖZET TABLO — METHODS İÇİN
# ==============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("ÖZET — Methods için tablo\n")
cat(strrep("=", 60), "\n", sep = "")

# Ortalama süreler model bazında (7 ilaç ortalaması)
model_summary <- timings_df %>%
  dplyr::group_by(Model) %>%
  dplyr::summarise(
    Mean_Train_sec = round(mean(Train_mean_sec), 3),
    SD_Train_sec   = round(sd(Train_mean_sec), 3),
    Min_Train_sec  = round(min(Train_mean_sec), 3),
    Max_Train_sec  = round(max(Train_mean_sec), 3),
    Mean_Predict_sec  = round(mean(Predict_mean_sec), 5),
    .groups = "drop"
  ) %>%
  dplyr::arrange(Mean_Train_sec)

write.csv(model_summary, paste0(BASE_PATH, "model_timings_summary.csv"),
          row.names = FALSE)

cat("\nOrtalama eğitim süresi (7 ilaç ortalaması):\n")
print(model_summary, row.names = FALSE)


# ==============================================================================
# 5. SİSTEM BİLGİSİ + PAKET SÜRÜMLERİ KAYDET
# ==============================================================================

# sessionInfo'yu txt olarak kaydet
sink(paste0(BASE_PATH, "sessionInfo.txt"))
cat("=== System Information ===\n")
cat("Generated:", format(Sys.time()), "\n\n")
print(sessionInfo())
sink()

cat("\n", strrep("=", 60), "\n", sep = "")
cat("✅ 98_timing_analysis.R TAMAMLANDI\n")
cat(strrep("=", 60), "\n", sep = "")
cat("\nÇıktılar:\n")
cat("  ", BASE_PATH, "model_timings.csv          (42 satır: 7 ilaç × 6 model)\n", sep = "")
cat("  ", BASE_PATH, "model_timings_summary.csv  (6 satır: model özet)\n", sep = "")
cat("  ", BASE_PATH, "sessionInfo.txt            (paket sürümleri)\n", sep = "")
cat("\n💡 Methods için suggestion:\n")
cat("  All models were trained on a [Mac model] with [N] cores and [X] GB RAM,\n")
cat("  running R [version] (see sessionInfo.txt for full package versions).\n")
