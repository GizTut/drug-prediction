# ==============================================================================
# 05_train_final_models.R — 1-SE kararına göre 6 final model eğit
# ==============================================================================
# Bu script:
#   - Her ilaç için 04_statistical_tests'in seçtiği gen sayısı + hiperparametreleri
#     kullanarak final 6 modeli eğitir
#   - Her ilaç klasörünün models/ alt klasörüne kaydeder
#
# Bağımlı: 00_setup.R + 03 + 04 çıktıları
# Tahmini süre: ~10 dakika
# ==============================================================================
setwd(SCRIPTS_PATH)

source("00_setup.R")
cat("\n[05_train_final_models] Final modeller eğitiliyor...\n")


# ==============================================================================
# TEK İLAÇ İÇİN FİNAL MODEL EĞİTİMİ
# ==============================================================================

train_finals_for_drug <- function(drug_folder) {

  cat("\n", strrep("-", 60), "\n", sep = "")
  cat("İLAÇ: ", drug_folder, "\n", sep = "")
  cat(strrep("-", 60), "\n", sep = "")

  drug_path    <- paste0(BASE_PATH, drug_folder, "/")
  results_path <- paste0(drug_path, "results/")
  models_path  <- paste0(drug_path, "models/")
  dir.create(models_path, showWarnings = FALSE, recursive = TRUE)

  # 1-SE kararı ve hiperparametreleri yükle
  decision_df <- read.csv(paste0(results_path, "final_gene_selection.csv"),
                          stringsAsFactors = FALSE)

  knn_summary  <- read.csv(paste0(results_path, "knn_results.csv"))
  enet_summary <- read.csv(paste0(results_path, "enet_results.csv"))
  nnet_summary <- read.csv(paste0(results_path, "nnet_results.csv"))
  svm_summary  <- read.csv(paste0(results_path, "svm_results.csv"))
  xgb_summary  <- read.csv(paste0(results_path, "xgb_results.csv"))
  rf_summary   <- read.csv(paste0(results_path, "rf_results.csv"))

  # Yardımcı: 1-SE seçtiği gen sayısındaki train datasını yükle
  get_train <- function(model_name) {
    n <- decision_df$Final_Choice[decision_df$Model == model_name]
    df <- read.csv(paste0(drug_path, "data_train_", n, ".csv"),
                   row.names = 1, check.names = FALSE)
    df$IC50 <- as.numeric(df$IC50)
    list(df = df, n = n)
  }

  get_params <- function(summary_df, n) {
    summary_df %>% dplyr::filter(Genes == n)
  }

  ctrl_none <- trainControl(method = "none")

  # ---- KNN ----
  knn_data <- get_train("KNN")
  knn_p <- get_params(knn_summary, knn_data$n)
  set.seed(GLOBAL_SEED)
  knn_final <- caret::train(IC50 ~ ., data = knn_data$df, method = "kknn",
                     trControl = ctrl_none,
                     tuneGrid = data.frame(kmax = knn_p$kmax,
                                           distance = knn_p$distance,
                                           kernel = knn_p$kernel))
  save(knn_final, file = paste0(models_path, "model_knn.RData"))
  cat("  ✓ KNN (", knn_data$n, " gen)\n", sep = "")

  # ---- SVM ----
  svm_data <- get_train("SVM")
  svm_p <- get_params(svm_summary, svm_data$n)
  set.seed(GLOBAL_SEED)
  svm_final <- caret::train(IC50 ~ ., data = svm_data$df, method = "svmRadial",
                     trControl = ctrl_none,
                     tuneGrid = data.frame(C = svm_p$C, sigma = svm_p$sigma))
  save(svm_final, file = paste0(models_path, "model_svm.RData"))
  cat("  ✓ SVM (", svm_data$n, " gen)\n", sep = "")

  # ---- ENET ----
  enet_data <- get_train("ENET")
  enet_p <- get_params(enet_summary, enet_data$n)
  set.seed(GLOBAL_SEED)
  enet_final <- caret::train(IC50 ~ ., data = enet_data$df, method = "enet",
                      trControl = ctrl_none,
                      tuneGrid = data.frame(lambda = enet_p$lambda,
                                            fraction = enet_p$fraction))
  save(enet_final, file = paste0(models_path, "model_enet.RData"))
  cat("  ✓ Elastic Net (", enet_data$n, " gen)\n", sep = "")

  # ---- NNET ----
  nnet_data <- get_train("NNET")
  nnet_p <- get_params(nnet_summary, nnet_data$n)
  set.seed(GLOBAL_SEED)
  nnet_final <- caret::train(IC50 ~ ., data = nnet_data$df, method = "nnet",
                      trControl = ctrl_none,
                      tuneGrid = data.frame(size = nnet_p$size,
                                            decay = nnet_p$decay),
                      linout = TRUE, trace = FALSE, maxit = 500)
  save(nnet_final, file = paste0(models_path, "model_nnet.RData"))
  cat("  ✓ Neural Net (", nnet_data$n, " gen)\n", sep = "")

  # ---- XGBoost ----
  xgb_data <- get_train("XGB")
  xgb_p <- get_params(xgb_summary, xgb_data$n)
  X_xgb <- as.matrix(xgb_data$df[, setdiff(colnames(xgb_data$df), "IC50"),
                                 drop = FALSE])
  y_xgb <- xgb_data$df$IC50
  dtrain <- xgb.DMatrix(data = X_xgb, label = y_xgb)
  xgb_params <- list(objective = "reg:squarederror", eval_metric = "rmse",
                     eta = xgb_p$eta, max_depth = xgb_p$max_depth,
                     subsample = xgb_p$subsample,
                     colsample_bytree = xgb_p$colsample_bytree,
                     min_child_weight = xgb_p$min_child_weight,
                     gamma = xgb_p$gamma)
  set.seed(GLOBAL_SEED)
  xgb_final <- xgb.train(params = xgb_params, data = dtrain,
                         nrounds = xgb_p$nrounds, verbose = 0)
  xgb.save(xgb_final, paste0(models_path, "model_xgb.json"))
  cat("  ✓ XGBoost (", xgb_data$n, " gen)\n", sep = "")

  # ---- Random Forest ----
  rf_data <- get_train("RF")
  rf_p <- get_params(rf_summary, rf_data$n)
  set.seed(GLOBAL_SEED)
  rf_final <- caret::train(IC50 ~ ., data = rf_data$df, method = "rf",
                    trControl = ctrl_none,
                    tuneGrid = data.frame(mtry = rf_p$mtry),
                    ntree = 500)
  save(rf_final, file = paste0(models_path, "model_rf.RData"))
  cat("  ✓ Random Forest (", rf_data$n, " gen)\n", sep = "")
}


# ==============================================================================
# 7 İLAÇ İÇİN
# ==============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("Final modeller 7 ilaç için eğitiliyor...\n")
cat(strrep("=", 60), "\n", sep = "")

for (drug in DRUG_FOLDER_NAMES) {
  tryCatch(train_finals_for_drug(drug),
           error = function(e) {
             cat("⚠️ HATA:", drug, "-", conditionMessage(e), "\n")
           })
}

cat("\n", strrep("=", 60), "\n", sep = "")
cat("✅ 05_train_final_models.R TAMAMLANDI\n")
cat(strrep("=", 60), "\n", sep = "")
cat("\n📌 Sonraki adım: 06_test_and_validate.R\n")
