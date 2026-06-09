# ==============================================================================
# 08_feature_importance.R — Permutation feature importance (6 model × 7 ilaç)
# ==============================================================================
# YAYIN İÇİN ÖNEMLİ:
# Hakem 3.6 "SVM için özellik önemi veya SHAP açıklaması" istemişti.
# Bu script TÜM modeller için tek tip (permutation) importance hesaplar:
# - Model-agnostic: aynı yöntem, hepsine adil
# - Yorumlanabilir: "Bu geni karıştırırsak RMSE ne kadar artar?"
# - SHAP'a benzer felsefe ama daha basit ve hızlı
#
# Yapacağı:
#   - Her ilaç için 6 final modeli yükle
#   - Her gen için "permute → predict → RMSE artışı" ölç
#   - 10 tekrar (varyans için)
#   - Top 10 gen listesi (ilaç × model bazında)
#   - 7 ilaç × 6 model × top 10 = 420 satırlık ana tablo
#   - Bonus: "Universal genes" (3+ modelde top 10'da olanlar)
#
# Bağımlı: 00_setup + 05 (final modeller) + 06 (test sets)
# Tahmini süre: ~30-40 dakika (7 ilaç × 6 model × 50-150 gen × 10 permutation)
# ==============================================================================
setwd(SCRIPTS_PATH)

source("00_setup.R")
cat("\n[08_feature_importance] Permutation importance başlıyor...\n")
cat("⏰ Tahmini süre: ~30-40 dakika\n\n")


# ==============================================================================
# 1. PERMUTATION IMPORTANCE FONKSİYONU
# ==============================================================================

# Bir genin permutation importance'ını hesapla
# (n_repeats kez karıştır, ortalama RMSE artışını ver)
permutation_importance <- function(model, X_test, y_test, model_type,
                                    n_repeats = 10, seed = GLOBAL_SEED) {
  # Predict fonksiyonu (model_type'a göre)
  predict_fn <- function(X) {
    if (model_type == "xgb") {
      dt <- xgb.DMatrix(data = as.matrix(X))
      predict(model, dt)
    } else {
      predict(model, newdata = X)
    }
  }

  # Baseline RMSE
  baseline_pred <- predict_fn(X_test)
  baseline_rmse <- rmse(y_test, baseline_pred)

  # Her gen için permute et
  feature_names <- colnames(X_test)
  importance_results <- data.frame(
    Gene = feature_names,
    Importance_mean = NA_real_,
    Importance_sd = NA_real_,
    stringsAsFactors = FALSE
  )

  for (g_idx in seq_along(feature_names)) {
    rmse_increases <- numeric(n_repeats)
    original_col <- X_test[[feature_names[g_idx]]]

    for (rep in 1:n_repeats) {
      set.seed(seed + rep)
      X_perm <- X_test
      X_perm[[feature_names[g_idx]]] <- sample(original_col)
      perm_pred <- predict_fn(X_perm)
      perm_rmse <- rmse(y_test, perm_pred)
      rmse_increases[rep] <- perm_rmse - baseline_rmse
    }

    importance_results$Importance_mean[g_idx] <- mean(rmse_increases)
    importance_results$Importance_sd[g_idx]   <- sd(rmse_increases)
  }

  importance_results <- importance_results[order(-importance_results$Importance_mean), ]
  importance_results$Rank <- seq_len(nrow(importance_results))

  list(baseline_rmse = baseline_rmse,
       importance = importance_results)
}


# ==============================================================================
# 2. TEK İLAÇ İÇİN 6 MODELİN IMPORTANCE'INI HESAPLA
# ==============================================================================

importance_for_drug <- function(drug_folder) {
  cat("\n", strrep("-", 60), "\n", sep = "")
  cat("İLAÇ: ", drug_folder, "\n", sep = "")
  cat(strrep("-", 60), "\n", sep = "")

  drug_path    <- paste0(BASE_PATH, drug_folder, "/")
  results_path <- paste0(drug_path, "results/")
  models_path  <- paste0(drug_path, "models/")

  # 1-SE kararı yükle
  decision_df <- read.csv(paste0(results_path, "final_gene_selection.csv"))

  # Test verilerini yükle
  test_50  <- read.csv(paste0(drug_path, "data_test_50.csv"),  row.names = 1, check.names = FALSE)
  test_100 <- read.csv(paste0(drug_path, "data_test_100.csv"), row.names = 1, check.names = FALSE)
  test_150 <- read.csv(paste0(drug_path, "data_test_150.csv"), row.names = 1, check.names = FALSE)
  test_50$IC50  <- as.numeric(test_50$IC50)
  test_100$IC50 <- as.numeric(test_100$IC50)
  test_150$IC50 <- as.numeric(test_150$IC50)

  pick_test <- function(n) {
    if (n == 50) test_50 else if (n == 100) test_100 else test_150
  }

  # Her modelin 1-SE seçimi
  get_n <- function(model_name) {
    decision_df$Final_Choice[decision_df$Model == model_name]
  }

  # Modelleri yükle
  load(paste0(models_path, "model_knn.RData"))
  load(paste0(models_path, "model_svm.RData"))
  load(paste0(models_path, "model_enet.RData"))
  load(paste0(models_path, "model_nnet.RData"))
  load(paste0(models_path, "model_rf.RData"))
  xgb_final <- xgb.load(paste0(models_path, "model_xgb.json"))

  # 6 model için permutation importance
  all_imp <- list()

  models_config <- list(
    list(name = "KNN",           model = knn_final,  type = "caret", n_key = "KNN"),
    list(name = "Elastic Net",   model = enet_final, type = "caret", n_key = "ENET"),
    list(name = "Neural Net",    model = nnet_final, type = "caret", n_key = "NNET"),
    list(name = "SVM",           model = svm_final,  type = "caret", n_key = "SVM"),
    list(name = "XGBoost",       model = xgb_final,  type = "xgb",   n_key = "XGB"),
    list(name = "Random Forest", model = rf_final,   type = "caret", n_key = "RF")
  )

  for (cfg in models_config) {
    cat("  Computing", cfg$name, "permutation importance...\n")
    n_genes <- get_n(cfg$n_key)
    test_df <- pick_test(n_genes)

    X_test <- test_df[, setdiff(colnames(test_df), "IC50"), drop = FALSE]
    y_test <- test_df$IC50

    t0 <- Sys.time()
    res <- permutation_importance(cfg$model, X_test, y_test, cfg$type,
                                  n_repeats = 10)
    t1 <- Sys.time()
    duration <- round(as.numeric(difftime(t1, t0, units = "secs")), 1)

    cat("    baseline RMSE:", round(res$baseline_rmse, 4),
        "| top gene:", res$importance$Gene[1],
        "(impact +", round(res$importance$Importance_mean[1], 4), ") |",
        duration, "sn\n")

    res$importance$Drug <- drug_folder
    res$importance$Model <- cfg$name
    res$importance$N_genes_used <- n_genes
    all_imp[[cfg$name]] <- res$importance
  }

  # Hafıza temizle
  rm(knn_final, svm_final, enet_final, nnet_final, rf_final, xgb_final)
  gc(verbose = FALSE)

  bind_rows(all_imp)
}


# ==============================================================================
# 3. 7 İLAÇ İÇİN TÜM IMPORTANCE'LARI TOPLA
# ==============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("Permutation importance hesabı 7 ilaç için başlıyor...\n")
cat(strrep("=", 60), "\n", sep = "")

all_drug_imp <- list()
total_start <- Sys.time()

for (drug in DRUG_FOLDER_NAMES) {
  res <- tryCatch(importance_for_drug(drug),
                  error = function(e) {
                    cat("⚠️ HATA:", drug, "-", conditionMessage(e), "\n")
                    NULL
                  })
  if (!is.null(res)) all_drug_imp[[drug]] <- res
}

total_end <- Sys.time()
total_min <- round(as.numeric(difftime(total_end, total_start, units = "mins")), 1)

importance_full <- bind_rows(all_drug_imp)


# ==============================================================================
# 4. KAYDET — TAM TABLO
# ==============================================================================

# Sütun sırasını düzenle
importance_full <- importance_full %>%
  dplyr::select(Drug, Model, N_genes_used, Rank, Gene,
                Importance_mean, Importance_sd)

write.csv(importance_full,
          paste0(BASE_PATH, "permutation_importance_full.csv"),
          row.names = FALSE)


# ==============================================================================
# 5. TOP 10 GENLER (yayın için ana tablo)
# ==============================================================================

top10_per_model <- importance_full %>%
  dplyr::group_by(Drug, Model) %>%
  dplyr::slice_min(order_by = Rank, n = 10) %>%
  dplyr::ungroup()

write.csv(top10_per_model,
          paste0(BASE_PATH, "permutation_top10_per_drug_model.csv"),
          row.names = FALSE)


# ==============================================================================
# 6. UNIVERSAL GENS — Birden fazla modelde top 10'da olanlar
# ==============================================================================

# Her ilaç için: bir gen kaç modelin top 10'unda?
universal_per_drug <- top10_per_model %>%
  dplyr::group_by(Drug, Gene) %>%
  dplyr::summarise(
    N_models = dplyr::n(),
    Models = paste(sort(unique(Model)), collapse = ", "),
    Mean_importance = round(mean(Importance_mean), 4),
    .groups = "drop"
  ) %>%
  dplyr::filter(N_models >= 2) %>%
  dplyr::arrange(Drug, -N_models, -Mean_importance)

write.csv(universal_per_drug,
          paste0(BASE_PATH, "universal_genes_per_drug.csv"),
          row.names = FALSE)

# Tüm ilaçlar arasında: bir gen kaç (ilaç × model) çiftinin top 10'unda?
cross_drug_universal <- top10_per_model %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    Total_appearances = dplyr::n(),
    N_drugs = dplyr::n_distinct(Drug),
    N_models = dplyr::n_distinct(Model),
    Mean_importance = round(mean(Importance_mean), 4),
    .groups = "drop"
  ) %>%
  dplyr::arrange(-Total_appearances) %>%
  dplyr::slice_head(n = 30)

write.csv(cross_drug_universal,
          paste0(BASE_PATH, "cross_drug_universal_genes.csv"),
          row.names = FALSE)


# ==============================================================================
# 7. KEGG/GO İÇİN GİRDİ — Eşsiz genler listesi
# ==============================================================================

# Her ilaç için, tüm modellerin top 10'unun birleşim kümesi
genes_for_enrichment <- top10_per_model %>%
  dplyr::group_by(Drug) %>%
  dplyr::summarise(
    Top_genes_union = paste(sort(unique(Gene)), collapse = ","),
    N_unique_genes = dplyr::n_distinct(Gene),
    .groups = "drop"
  )

write.csv(genes_for_enrichment,
          paste0(BASE_PATH, "genes_for_kegg_go_enrichment.csv"),
          row.names = FALSE)


# ==============================================================================
# 8. ÖZET YAZILARI
# ==============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("✅ 08_feature_importance.R TAMAMLANDI (", total_min, " dakika)\n", sep = "")
cat(strrep("=", 60), "\n", sep = "")

cat("\nÜretilen dosyalar:\n")
cat("  ", BASE_PATH, "permutation_importance_full.csv  (her gen × ilaç × model)\n", sep = "")
cat("  ", BASE_PATH, "permutation_top10_per_drug_model.csv  (yayın ana tablo)\n", sep = "")
cat("  ", BASE_PATH, "universal_genes_per_drug.csv  (her ilaç için 'evrensel' genler)\n", sep = "")
cat("  ", BASE_PATH, "cross_drug_universal_genes.csv  (tüm ilaçlarda en sık görülen 30 gen)\n", sep = "")
cat("  ", BASE_PATH, "genes_for_kegg_go_enrichment.csv  (KEGG/GO girdisi)\n", sep = "")

# Konsolda özet göster
cat("\n=== TÜM İLAÇLARDA EN SIK GÖRÜLEN GENLER (top 30) ===\n")
print(head(cross_drug_universal, 20), row.names = FALSE)

cat("\n=== HER İLAÇ İÇİN EŞSİZ GEN SAYISI (KEGG için) ===\n")
print(genes_for_enrichment %>% dplyr::select(Drug, N_unique_genes),
      row.names = FALSE)

cat("\n📌 Sonraki adım: 09_kegg_go_enrichment.R\n")
