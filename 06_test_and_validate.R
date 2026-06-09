# ==============================================================================
# 06_test_and_validate.R — Test seti değerlendirme + validation tahminleri
# ==============================================================================
# Bu script:
#   - Final modelleri yükler
#   - Test setinde tahmin yapar (MAE, RMSE, R²)
#   - Validation setinde tahmin yapar (gerçek IC50 yok, sadece tahmin)
#   - Her ilaç için results/ klasörüne kaydeder
#
# Bağımlı: 00_setup.R + 05_train_final_models.R çıktıları
# Tahmini süre: ~2 dakika
# ==============================================================================
setwd(SCRIPTS_PATH)

source("00_setup.R")
cat("\n[06_test_and_validate] Test ve validation tahminleri başlıyor...\n")


# ==============================================================================
# TEK İLAÇ İÇİN TEST + VALIDATION
# ==============================================================================

evaluate_drug <- function(drug_folder) {

  drug_path    <- paste0(BASE_PATH, drug_folder, "/")
  results_path <- paste0(drug_path, "results/")
  models_path  <- paste0(drug_path, "models/")

  # 1-SE kararını yükle (hangi gen sayısı?)
  decision_df <- read.csv(paste0(results_path, "final_gene_selection.csv"))

  get_n_genes <- function(model_name) {
    decision_df$Final_Choice[decision_df$Model == model_name]
  }

  # Test verilerini yükle
  test_50  <- read.csv(paste0(drug_path, "data_test_50.csv"),
                       row.names = 1, check.names = FALSE)
  test_100 <- read.csv(paste0(drug_path, "data_test_100.csv"),
                       row.names = 1, check.names = FALSE)
  test_150 <- read.csv(paste0(drug_path, "data_test_150.csv"),
                       row.names = 1, check.names = FALSE)
  test_50$IC50  <- as.numeric(test_50$IC50)
  test_100$IC50 <- as.numeric(test_100$IC50)
  test_150$IC50 <- as.numeric(test_150$IC50)
  pick_test <- function(n) {
    if (n == 50) test_50 else if (n == 100) test_100 else test_150
  }

  # Validation
  val_50  <- read.csv(paste0(drug_path, "data_val_50.csv"),
                      row.names = 1, check.names = FALSE)
  val_100 <- read.csv(paste0(drug_path, "data_val_100.csv"),
                      row.names = 1, check.names = FALSE)
  val_150 <- read.csv(paste0(drug_path, "data_val_150.csv"),
                      row.names = 1, check.names = FALSE)
  pick_val <- function(n) {
    if (n == 50) val_50 else if (n == 100) val_100 else val_150
  }

  # Modelleri yükle
  load(paste0(models_path, "model_knn.RData"))
  load(paste0(models_path, "model_svm.RData"))
  load(paste0(models_path, "model_enet.RData"))
  load(paste0(models_path, "model_nnet.RData"))
  load(paste0(models_path, "model_rf.RData"))
  xgb_final <- xgb.load(paste0(models_path, "model_xgb.json"))

  # ---- Test tahminleri ----
  results_test <- data.frame(
    Model = character(0), Genes = numeric(0),
    MAE = numeric(0), RMSE = numeric(0), R2 = numeric(0)
  )

  predict_one <- function(model_obj, model_label, n) {
    test_df <- pick_test(n)
    if (model_label == "XGBoost") {
      dt <- xgb.DMatrix(data = as.matrix(
        test_df[, setdiff(colnames(test_df), "IC50")]))
      preds <- predict(model_obj, dt)
    } else {
      preds <- predict(model_obj, newdata = test_df)
    }
    data.frame(
      Model = model_label, Genes = n,
      MAE  = round(mae(test_df$IC50, preds), 4),
      RMSE = round(rmse(test_df$IC50, preds), 4),
      R2   = round(cor(test_df$IC50, preds)^2, 4)
    )
  }

  results_test <- rbind(
    predict_one(knn_final,  "KNN",          get_n_genes("KNN")),
    predict_one(svm_final,  "SVM",          get_n_genes("SVM")),
    predict_one(enet_final, "Elastic Net",  get_n_genes("ENET")),
    predict_one(nnet_final, "Neural Net",   get_n_genes("NNET")),
    predict_one(xgb_final,  "XGBoost",      get_n_genes("XGB")),
    predict_one(rf_final,   "Random Forest",get_n_genes("RF"))
  )

  write.csv(results_test, paste0(results_path, "test_performances.csv"),
            row.names = FALSE)

  # ---- Validation tahminleri (gerçek IC50 yok, sadece tahmin) ----
  predict_val <- function(model_obj, model_label, n) {
    val_df <- pick_val(n)
    if (model_label == "XGBoost") {
      dv <- xgb.DMatrix(data = as.matrix(val_df))
      preds <- predict(model_obj, dv)
    } else {
      preds <- predict(model_obj, newdata = val_df)
    }
    round(preds, 3)
  }

  validation_predictions <- data.frame(
    Sample_ID = rownames(val_50),
    IC50_KNN  = predict_val(knn_final,  "KNN",         get_n_genes("KNN")),
    IC50_SVM  = predict_val(svm_final,  "SVM",         get_n_genes("SVM")),
    IC50_ENET = predict_val(enet_final, "Elastic Net", get_n_genes("ENET")),
    IC50_NNET = predict_val(nnet_final, "Neural Net",  get_n_genes("NNET")),
    IC50_XGB  = predict_val(xgb_final,  "XGBoost",     get_n_genes("XGB")),
    IC50_RF   = predict_val(rf_final,   "Random Forest", get_n_genes("RF"))
  )

  write.csv(validation_predictions,
            paste0(results_path, "validation_IC50_predictions.csv"),
            row.names = FALSE)

  cat("✓", drug_folder, "— test (R² =",
      round(mean(results_test$R2), 3), ") + validation tahminleri\n")

  # Hafızayı temizle
  rm(knn_final, svm_final, enet_final, nnet_final, rf_final, xgb_final)
  gc(verbose = FALSE)

  return(results_test)
}


# ==============================================================================
# 7 İLAÇ İÇİN
# ==============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("Test + validation 7 ilaç için çalışıyor...\n")
cat(strrep("=", 60), "\n", sep = "")

for (drug in DRUG_FOLDER_NAMES) {
  tryCatch(evaluate_drug(drug),
           error = function(e) {
             cat("⚠️ HATA:", drug, "-", conditionMessage(e), "\n")
           })
}

cat("\n", strrep("=", 60), "\n", sep = "")
cat("✅ 06_test_and_validate.R TAMAMLANDI\n")
cat(strrep("=", 60), "\n", sep = "")
cat("\n📌 Sonraki adım: 07_visualize.R\n")
