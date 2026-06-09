# ==============================================================================
# 04_statistical_tests.R — Friedman + Wilcoxon + 1-SE kuralı
# ==============================================================================
# Bu script:
#   - Her ilaç için, her model için: Friedman testi (3 gen sayısı arasında)
#   - Wilcoxon signed-rank pairwise (50vs100, 100vs150, 50vs150)
#   - 1-SE kuralı ile son gen sayısı seçimi
#
# Bağımlı: 00_setup.R + 03_run_nested_cv.R çıktıları
# Tahmini süre: ~10 saniye
# ==============================================================================
setwd(SCRIPTS_PATH)

source("00_setup.R")
cat("\n[04_statistical_tests] İstatistiksel testler başlıyor...\n")


# ==============================================================================
# TEK İLAÇ İÇİN STATİSTİKSEL TESTLER
# ==============================================================================

run_stats_for_drug <- function(drug_folder) {

  results_path <- paste0(BASE_PATH, drug_folder, "/results/")
  fold_file <- paste0(results_path, "fold_details_all_models.csv")

  if (!file.exists(fold_file)) {
    cat("⚠️ ", drug_folder, " — fold_details bulunamadı, atlanıyor\n", sep = "")
    return(NULL)
  }

  fold_details_df <- read.csv(fold_file, check.names = FALSE)

  # ---- Friedman ----
  friedman_results <- list()
  for (m in unique(fold_details_df$Model)) {
    sub <- fold_details_df %>% dplyr::filter(Model == m)
    wide_rmse <- sub %>% dplyr::select(Fold, Genes, RMSE) %>%
      pivot_wider(names_from = Genes, values_from = RMSE)
    mat <- as.matrix(wide_rmse[, -1])
    ft <- friedman.test(mat)
    friedman_results[[m]] <- data.frame(
      Model = m, Chi_sq = round(ft$statistic, 3),
      df = ft$parameter, p_value = round(ft$p.value, 4),
      Sonuc = ifelse(ft$p.value < 0.05, "Anlamlı fark var", "Fark anlamsız")
    )
  }
  friedman_df <- bind_rows(friedman_results)

  # ---- Wilcoxon Pairwise ----
  pairs <- list(c(50, 100), c(100, 150), c(50, 150))
  wilcox_results <- list()
  for (m in unique(fold_details_df$Model)) {
    for (p in pairs) {
      g1 <- fold_details_df %>%
        dplyr::filter(Model == m, Genes == p[1]) %>%
        dplyr::arrange(Fold) %>% dplyr::pull(RMSE)
      g2 <- fold_details_df %>%
        dplyr::filter(Model == m, Genes == p[2]) %>%
        dplyr::arrange(Fold) %>% dplyr::pull(RMSE)
      wt <- suppressWarnings(wilcox.test(g1, g2, paired = TRUE))
      wilcox_results[[paste(m, p[1], p[2], sep = "_")]] <- data.frame(
        Model = m, Karsilastirma = paste0(p[1], " vs ", p[2]),
        Mean_RMSE_1 = round(mean(g1), 4), Mean_RMSE_2 = round(mean(g2), 4),
        V_stat = wt$statistic, p_value = round(wt$p.value, 4),
        Anlamli = ifelse(wt$p.value < 0.05, "EVET", "hayır")
      )
    }
  }
  wilcox_df <- bind_rows(wilcox_results)

  # ---- 1-SE Kuralı ----
  decision_results <- list()
  for (m in unique(fold_details_df$Model)) {
    sub <- fold_details_df %>% dplyr::filter(Model == m)
    summary_stats <- sub %>% dplyr::group_by(Genes) %>%
      dplyr::summarise(mean_RMSE = mean(RMSE),
                       se_RMSE = sd(RMSE) / sqrt(n()),
                       .groups = "drop")
    best_idx   <- which.min(summary_stats$mean_RMSE)
    best_rmse  <- summary_stats$mean_RMSE[best_idx]
    best_se    <- summary_stats$se_RMSE[best_idx]
    best_genes <- summary_stats$Genes[best_idx]
    threshold  <- best_rmse + best_se
    candidates <- summary_stats %>% dplyr::filter(mean_RMSE <= threshold)
    selected   <- min(candidates$Genes)
    decision_results[[m]] <- data.frame(
      Model = m, Best_Genes = best_genes,
      Best_RMSE = round(best_rmse, 4), SE = round(best_se, 4),
      Threshold_1SE = round(threshold, 4),
      Final_Choice = selected,
      Note = ifelse(selected == best_genes, "En iyi = en basit",
                    paste0("Daha basit seçildi (", best_genes,
                           " yerine ", selected, ")"))
    )
  }
  decision_df <- bind_rows(decision_results)

  # Kaydet
  write.csv(friedman_df, paste0(results_path, "friedman_test.csv"),
            row.names = FALSE)
  write.csv(wilcox_df, paste0(results_path, "wilcoxon_pairwise.csv"),
            row.names = FALSE)
  write.csv(decision_df, paste0(results_path, "final_gene_selection.csv"),
            row.names = FALSE)

  cat("✓", drug_folder, "— Friedman, Wilcoxon, 1-SE kaydedildi\n")
  return(decision_df)
}


# ==============================================================================
# 7 İLAÇ İÇİN ÇALIŞTIR
# ==============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("İstatistiksel testler 7 ilaç için çalıştırılıyor...\n")
cat(strrep("=", 60), "\n", sep = "")

all_decisions <- list()
for (drug in DRUG_FOLDER_NAMES) {
  res <- run_stats_for_drug(drug)
  if (!is.null(res)) {
    res$Drug <- drug
    all_decisions[[drug]] <- res
  }
}

# Toplu özet
if (length(all_decisions) > 0) {
  all_decisions_df <- bind_rows(all_decisions) %>%
    dplyr::select(Drug, Model, Best_Genes, Best_RMSE, Final_Choice, Note)
  write.csv(all_decisions_df,
            paste0(BASE_PATH, "ALL_DRUGS_GENE_SELECTION.csv"),
            row.names = FALSE)
}

cat("\n", strrep("=", 60), "\n", sep = "")
cat("✅ 04_statistical_tests.R TAMAMLANDI\n")
cat(strrep("=", 60), "\n", sep = "")
cat("\n📌 Sonraki adım: 05_train_final_models.R\n")
