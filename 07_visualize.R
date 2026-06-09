# ==============================================================================
# 07_visualize.R — Yayın için figürler + tablolar
# ==============================================================================
#
# ÜRETİLEN FIGÜRLER (makale sırasıyla):
#   Figure 2: Violin plot — 7 ilacın IC50 dağılımı
#   Figure 3: Heatmap — 7 ilaç × 6 model, test RMSE
#   Figure 4: Bar chart — her ilacın en iyi modeli (RMSE + R², ayrı panel)
#   Figure 5: Scatter — tahmin vs gerçek (her ilacın en iyi modeli)
#   Figure 6: Fold boxplot — nested CV RMSE dağılımı (model kararlılığı)
#   [Figure 7 ve 8 → 09_kegg_go_enrichment.R üretir]
#
# ÜRETİLEN TABLOLAR:
#   Table_1_data_summary.csv
#   Table_3_best_model_per_drug.csv
#   Table_S1_full_test_performance.csv
#   Table_S2_gene_selection.csv
#
# Bağımlı: 00_setup.R + 01–06 çıktıları
# ==============================================================================
setwd(SCRIPTS_PATH)

source("00_setup.R")
cat("\n[07_visualize] Görseller üretiliyor...\n")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(RColorBrewer)
  library(pheatmap)
  library(patchwork)
})

# ==============================================================================
# RENK PALETİ — tüm figürlerde kullanılır
# ==============================================================================
# Model renkleri (Figure 3, 4, 5, 6'da)
MODEL_COLORS <- c(
  "KNN"           = "#A3C4BC",   # Pastel Blue
  "Elastic Net"   = "#F5CAC3",   # Soft Pink
  "Random Forest" = "#D9EAD3",   # Mint Green
  "SVM"           = "#C4A8B8",   # Light Lavender (koyulaştırılmış)
  "XGBoost"       = "#FFE3B3",   # Pale Lemon
  "Neural Net"    = "#D4A5A5"    # Dusty Rose
)

# İlaç renkleri (Figure 2 violin'de)
DRUG_COLORS <- c(
  "RVX-208"   = "#A3C4BC",   # Pastel Blue
  "PFI-1"     = "#F5CAC3",   # Soft Pink
  "JQ1"       = "#D9EAD3",   # Mint Green
  "I-BRD9"    = "#EAD7D1",   # Light Lavender
  "OTX015"    = "#FFE3B3",   # Pale Lemon
  "I-BET-762" = "#C8DEB8",   # Deeper Mint
  "AZD5153"   = "#F0B8B8"    # Deeper Pink
)

# ==============================================================================
# YAYIN TEMASI — tüm ggplot figürlerinde ortak
# ==============================================================================
pub_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor   = element_blank(),
    strip.background   = element_rect(fill = "#F5F0EB", color = "#C8BFB5"),
    strip.text         = element_text(face = "bold", size = 10),
    legend.position    = "right",
    axis.text          = element_text(color = "black"),
    plot.title         = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle      = element_text(size = 10, hjust = 0.5, color = "gray40")
  )


# ==============================================================================
# ÖNCE: ALL_DRUGS_TEST_SUMMARY
# ==============================================================================

cat("\n[Hazırlık] Test sonuçları toplanıyor...\n")

combined_test <- list()
for (drug in DRUG_FOLDER_NAMES) {
  fp <- paste0(BASE_PATH, drug, "/results/test_performances.csv")
  if (file.exists(fp)) {
    df <- read.csv(fp)
    df$Drug <- drug
    combined_test[[drug]] <- df
  }
}
all_test <- bind_rows(combined_test) %>%
  dplyr::select(Drug, Model, Genes, MAE, RMSE, R2)

write.csv(all_test,
          paste0(BASE_PATH, "ALL_DRUGS_TEST_SUMMARY.csv"),
          row.names = FALSE)

all_test$Drug  <- factor(all_test$Drug,  levels = DRUG_FOLDER_NAMES)
all_test$Model <- factor(all_test$Model, levels = MODEL_ORDER)


# ==============================================================================
# FIGURE 2: VİOLİN PLOT
# ==============================================================================
# AMACI: İlaçlar arasındaki IC50 dağılım farkını göster.
#   Pan-BET inhibitörler (JQ1, OTX015, I-BET-762, AZD5153) daha düşük IC50 →
#   daha güçlü aktivite. Seçici olanlar (RVX-208, I-BRD9) daha yüksek IC50.
# VERİ: data_train_50.csv → IC50 sütunu (log10 dönüşümlü)

cat("\n[Fig 2] Violin plot...\n")

violin_data <- list()
for (drug in DRUG_FOLDER_NAMES) {
  fp <- paste0(BASE_PATH, drug, "/data_train_50.csv")
  if (file.exists(fp)) {
    df <- read.csv(fp, row.names = 1, check.names = FALSE)
    violin_data[[drug]] <- data.frame(
      Drug       = drug,
      log10_IC50 = as.numeric(df$IC50)
    )
  }
}
violin_df <- bind_rows(violin_data)
violin_df$Drug <- factor(violin_df$Drug, levels = DRUG_FOLDER_NAMES)

n_labels <- violin_df %>%
  dplyr::group_by(Drug) %>%
  dplyr::summarise(n = n(), .groups = "drop")
y_min <- min(violin_df$log10_IC50, na.rm = TRUE)

fig2_violin <- ggplot(violin_df, aes(x = Drug, y = log10_IC50, fill = Drug)) +
  geom_violin(trim = FALSE, alpha = 0.85, color = "gray40", linewidth = 0.3) +
#  geom_boxplot(width = 0.12, outlier.size = 0.8, fill = "white",
#               color = "gray30", linewidth = 0.4) +  
  geom_text(data = n_labels,
            aes(x = Drug, y = y_min - 0.18, label = paste0("n=", n)),
            inherit.aes = FALSE, size = 3, color = "gray50") +
  scale_fill_manual(values = DRUG_COLORS) +
  labs(x = NULL,
       y = expression(log[10] ~ IC[50] ~ "(µM)")) +
  pub_theme +
  theme(axis.text.x   = element_text(angle = 30, hjust = 1),
        legend.position = "none")

ggsave(paste0(FIGURES_PATH, "Figure2_violin_IC50.png"),
       fig2_violin, width = 9, height = 6, dpi = 300, bg = "white")
ggsave(paste0(FIGURES_PATH, "Figure2_violin_IC50.pdf"),
       fig2_violin, width = 9, height = 6)
cat("  ✓ Figure2_violin_IC50.png/pdf\n")


# ==============================================================================
# FIGURE 3: HEATMAP
# ==============================================================================
# AMACI: 7 ilaç × 6 model RMSE matrisine kuşbakışı bakmak.
#   Yeşil = düşük RMSE = iyi; Kırmızı = yüksek RMSE = kötü.
#   RVX-208 satırı genel olarak en yeşil (en kolay tahmin edilebilir IC50 aralığı).
#   NNet sütunu genel olarak en kırmızı.
# VERİ: all_test → her Drug × Model kombinasyonu için test RMSE

cat("\n[Fig 3] RMSE heatmap...\n")

heatmap_df <- all_test %>%
  dplyr::select(Drug, Model, RMSE) %>%
  pivot_wider(names_from = Model, values_from = RMSE)

heatmap_mat <- as.matrix(heatmap_df[, -1])
rownames(heatmap_mat) <- as.character(heatmap_df$Drug)
col_order   <- order(colMeans(heatmap_mat, na.rm = TRUE))
heatmap_mat <- heatmap_mat[, col_order]
display_mat <- round(heatmap_mat, 3)

# Heatmap renkleri: paletimizle uyumlu pastel → koyu ton
heatmap_colors <- colorRampPalette(c("#A3C4BC", "#FFF8E1", "#D4A5A5"))(100)

png(paste0(FIGURES_PATH, "Figure3_heatmap.png"),
    width = 10, height = 6, units = "in", res = 300, bg = "white")
pheatmap(heatmap_mat,
         cluster_rows    = FALSE,
         cluster_cols    = FALSE,
         display_numbers = display_mat,
         number_color    = "black",
         fontsize_number = 10,
         color           = heatmap_colors,
         border_color    = "white",
         main   = "Test RMSE: 7 BET Inhibitors × 6 ML Models\n(teal = lower RMSE = better; pink = higher RMSE = worse)",
         fontsize    = 11,
         cellwidth   = 70,
         cellheight  = 38)
dev.off()
cat("  ✓ Figure3_heatmap.png\n")


# ==============================================================================
# FIGURE 4: İKİ PANEL BAR CHART — RMSE (üst) + R² (alt)
# ==============================================================================
# AMACI: Her ilaç için kazanan modeli ve performansını iki ayrı metrikle göster.
#   Üst panel: RMSE — düşük = iyi (çubuk kısa = iyi)
#   Alt panel: R²  — yüksek = iyi (çubuk uzun = iyi)
#   Renkler = kazanan model ailesi (hangi algoritma kazandı?)
# VERİ: all_test → her Drug için slice_min(RMSE)

cat("\n[Fig 4] İki panel bar chart (RMSE + R²)...\n")

best_per_drug <- all_test %>%
  dplyr::group_by(Drug) %>%
  dplyr::slice_min(RMSE, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(factor(Drug, levels = DRUG_FOLDER_NAMES))

p_rmse <- ggplot(best_per_drug, aes(x = Drug, y = RMSE, fill = Model)) +
  geom_bar(stat = "identity", width = 0.6,
           color = "gray40", linewidth = 0.3) +
  geom_text(aes(label = round(RMSE, 3)),
            vjust = -0.4, size = 3.2, color = "gray20") +
  scale_fill_manual(values = MODEL_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14))) +
  labs(x = NULL, y = "RMSE (log10 IC50)", fill = "Best model") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

p_r2 <- ggplot(best_per_drug, aes(x = Drug, y = R2, fill = Model)) +
  geom_bar(stat = "identity", width = 0.6,
           color = "gray40", linewidth = 0.3) +
  geom_text(aes(label = round(R2, 3)),
            vjust = -0.4, size = 3.2, color = "gray20") +
  scale_fill_manual(values = MODEL_COLORS) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.14)),
                     limits = c(0, 0.55)) +
  labs(x = NULL, y = expression(R^2), fill = "Best model") +
  pub_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig4_two_panel <- (p_rmse / p_r2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

ggsave(paste0(FIGURES_PATH, "Figure4_bar_best_model.png"),
       fig4_two_panel, width = 10, height = 9, dpi = 300, bg = "white")
ggsave(paste0(FIGURES_PATH, "Figure4_bar_best_model.pdf"),
       fig4_two_panel, width = 10, height = 9)
cat("  ✓ Figure4_bar_best_model.png/pdf\n")


# ==============================================================================
# FIGURE 5: SCATTER — tahmin vs gerçek
# ==============================================================================
# AMACI: Her ilacın kazanan modelinin test setindeki tahmin kalitesini göster.
#   Her panel = bir ilaç. Nokta = bir hücre hattı.
#   Kesikli çizgi = mükemmel tahmin (y=x).
#   Mavi çizgi = gerçek regresyon fit'i.
#   R² ve RMSE panel başlığında.
# VERİ: modeller disk'ten yükleniyor, test setine uygulanıyor.

cat("\n[Fig 5] Scatter — tahmin vs gerçek...\n")

scatter_data <- list()

for (drug in as.character(DRUG_FOLDER_NAMES)) {
  drug_path   <- paste0(BASE_PATH, drug, "/")
  models_path <- paste0(drug_path, "models/")

  best_row        <- best_per_drug %>% dplyr::filter(Drug == drug)
  best_model_name <- as.character(best_row$Model)
  n_genes         <- as.numeric(best_row$Genes)

  test_df <- read.csv(paste0(drug_path, "data_test_", n_genes, ".csv"),
                      row.names = 1, check.names = FALSE)
  test_df$IC50 <- as.numeric(test_df$IC50)

  if (best_model_name == "KNN") {
    load(paste0(models_path, "model_knn.RData"))
    preds <- predict(knn_final, newdata = test_df); rm(knn_final)
  } else if (best_model_name == "Elastic Net") {
    load(paste0(models_path, "model_enet.RData"))
    preds <- predict(enet_final, newdata = test_df); rm(enet_final)
  } else if (best_model_name == "Random Forest") {
    load(paste0(models_path, "model_rf.RData"))
    preds <- predict(rf_final, newdata = test_df); rm(rf_final)
  } else if (best_model_name == "SVM") {
    load(paste0(models_path, "model_svm.RData"))
    preds <- predict(svm_final, newdata = test_df); rm(svm_final)
  } else if (best_model_name == "Neural Net") {
    load(paste0(models_path, "model_nnet.RData"))
    preds <- predict(nnet_final, newdata = test_df); rm(nnet_final)
  } else if (best_model_name == "XGBoost") {
    xgb_final <- xgb.load(paste0(models_path, "model_xgb.json"))
    dt    <- xgb.DMatrix(data = as.matrix(
               test_df[, setdiff(colnames(test_df), "IC50")]))
    preds <- predict(xgb_final, dt); rm(xgb_final)
  }

  r2_val   <- round(cor(test_df$IC50, preds)^2, 3)
  rmse_val <- round(sqrt(mean((test_df$IC50 - preds)^2)), 3)

  scatter_data[[drug]] <- data.frame(
    Drug      = drug,
    Model     = best_model_name,
    Actual    = test_df$IC50,
    Predicted = as.numeric(preds),
    Label     = paste0(drug, "\n", best_model_name,
                       "  R²=", r2_val, "  RMSE=", rmse_val)
  )
  gc(verbose = FALSE)
}

scatter_df       <- bind_rows(scatter_data)
scatter_df$Drug  <- factor(scatter_df$Drug,  levels = DRUG_FOLDER_NAMES)
scatter_df$Model <- factor(scatter_df$Model, levels = MODEL_ORDER)

fig5_scatter <- ggplot(scatter_df,
                       aes(x = Actual, y = Predicted, color = Model)) +
  geom_point(alpha = 0.6, size = 1.6) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "gray40", linewidth = 0.6) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7,
              alpha = 0.15) +
  facet_wrap(~ Label, scales = "free", ncol = 4) +
  scale_color_manual(values = MODEL_COLORS) +
  labs(x = expression("Observed" ~ log[10] ~ IC[50]),
       y = expression("Predicted" ~ log[10] ~ IC[50])) +
  pub_theme +
  theme(legend.position = "none",
        strip.text = element_text(size = 8, face = "bold"))

ggsave(paste0(FIGURES_PATH, "Figure5_scatter_best_models.png"),
       fig5_scatter, width = 14, height = 10, dpi = 300, bg = "white")
ggsave(paste0(FIGURES_PATH, "Figure5_scatter_best_models.pdf"),
       fig5_scatter, width = 14, height = 10)
cat("  ✓ Figure5_scatter_best_models.png/pdf\n")


# ==============================================================================
# FIGURE 6: FOLD BOXPLOT
# ==============================================================================
# AMACI: Nested CV'deki dış fold RMSE'lerini göster → model kararlılığı.
#   Kutu dar = model tutarlı, fold'dan fold'a değişmiyor.
#   Kutu geniş = model unstable.
#   NNet her ilaçta en geniş + en yüksek → tutarsız ve kötü.
# VERİ: fold_details_all_models.csv (03_run_nested_cv üretir)

cat("\n[Fig 6] Fold boxplot...\n")

fold_all <- list()
for (drug in DRUG_FOLDER_NAMES) {
  fold_file     <- paste0(BASE_PATH, drug, "/results/fold_details_all_models.csv")
  decision_file <- paste0(BASE_PATH, drug, "/results/final_gene_selection.csv")
  if (file.exists(fold_file) && file.exists(decision_file)) {
    fold_data <- read.csv(fold_file)
    decision  <- read.csv(decision_file)
    fold_data$Drug <- drug

    fold_filtered <- list()
    for (m in unique(fold_data$Model)) {
      n_chosen <- decision$Final_Choice[decision$Model == m]
      if (length(n_chosen) > 0) {
        sub <- fold_data %>% dplyr::filter(Model == m, Genes == n_chosen)
        fold_filtered[[m]] <- sub
      }
    }
    fold_all[[drug]] <- bind_rows(fold_filtered)
  }
}
fold_df <- bind_rows(fold_all)

model_name_map <- c("KNN" = "KNN", "ENET" = "Elastic Net",
                    "RF"  = "Random Forest", "SVM"  = "SVM",
                    "XGB" = "XGBoost",       "NNET" = "Neural Net")
fold_df$Model <- model_name_map[fold_df$Model]
fold_df$Drug  <- factor(fold_df$Drug,  levels = DRUG_FOLDER_NAMES)
fold_df$Model <- factor(fold_df$Model, levels = MODEL_ORDER)

fig6_fold <- ggplot(fold_df, aes(x = Model, y = RMSE, fill = Model)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.8, linewidth = 0.35,
               outlier.color = "gray50") +
  geom_jitter(width = 0.15, alpha = 0.5, size = 1.2,
              aes(color = Model)) +
  facet_wrap(~ Drug, scales = "free_y", ncol = 3) +
  scale_fill_manual(values  = MODEL_COLORS) +
  scale_color_manual(values = MODEL_COLORS) +
  labs(x = NULL, y = "RMSE (log10 IC50)") +
  pub_theme +
  theme(axis.text.x    = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "none")

ggsave(paste0(FIGURES_PATH, "Figure6_fold_boxplot.png"),
       fig6_fold, width = 13, height = 10, dpi = 300, bg = "white")
ggsave(paste0(FIGURES_PATH, "Figure6_fold_boxplot.pdf"),
       fig6_fold, width = 13, height = 10)
cat("  ✓ Figure6_fold_boxplot.png/pdf\n")


# ==============================================================================
# TABLOLAR
# ==============================================================================

cat("\n[Tablolar] Üretiliyor...\n")

# Tablo 1: Veri özeti
summary_rows <- list()
for (drug in DRUG_FOLDER_NAMES) {
  train_fp <- paste0(BASE_PATH, drug, "/data_train_50.csv")
  test_fp  <- paste0(BASE_PATH, drug, "/data_test_50.csv")
  val_fp   <- paste0(BASE_PATH, drug, "/data_val_50.csv")
  if (file.exists(train_fp)) {
    tr <- nrow(read.csv(train_fp, row.names = 1))
    te <- nrow(read.csv(test_fp,  row.names = 1))
    va <- nrow(read.csv(val_fp,   row.names = 1))
    summary_rows[[drug]] <- data.frame(
      Drug = drug, N_total = tr + te + va,
      Train_n = tr, Test_n = te, Val_n = va
    )
  }
}
table1 <- bind_rows(summary_rows)
write.csv(table1, paste0(FIGURES_PATH, "Table_1_data_summary.csv"),
          row.names = FALSE)
cat("  ✓ Table_1_data_summary.csv\n")
print(table1, row.names = FALSE)

# Tablo 3: En iyi model per ilaç
table3 <- best_per_drug %>%
  dplyr::select(Drug, Best_Model = Model, Genes, RMSE, MAE, R2) %>%
  dplyr::mutate(across(c(RMSE, MAE, R2), ~ round(., 3)))
write.csv(table3, paste0(FIGURES_PATH, "Table_3_best_model_per_drug.csv"),
          row.names = FALSE)
cat("  ✓ Table_3_best_model_per_drug.csv\n")
print(table3, row.names = FALSE)

# Ek S1: Tüm model performansları
table_s1 <- all_test %>%
  dplyr::mutate(across(c(MAE, RMSE, R2), ~ round(., 4)))
write.csv(table_s1, paste0(FIGURES_PATH, "Table_S1_full_test_performance.csv"),
          row.names = FALSE)
cat("  ✓ Table_S1_full_test_performance.csv\n")

# Ek S2: 1-SE gen seçim özeti
gene_sel_all <- list()
for (drug in DRUG_FOLDER_NAMES) {
  fp <- paste0(BASE_PATH, drug, "/results/final_gene_selection.csv")
  if (file.exists(fp)) {
    df <- read.csv(fp); df$Drug <- drug
    gene_sel_all[[drug]] <- df
  }
}
table_s2 <- bind_rows(gene_sel_all) %>%
  dplyr::select(Drug, Model, Best_Genes, Best_RMSE, SE,
                Threshold_1SE, Final_Choice, Note)
write.csv(table_s2, paste0(FIGURES_PATH, "Table_S2_gene_selection.csv"),
          row.names = FALSE)
cat("  ✓ Table_S2_gene_selection.csv\n")

n_50    <- sum(table_s2$Final_Choice == 50)
n_total <- nrow(table_s2)
cat("\n  1-SE kontrolü:", n_50, "/", n_total, "=",
    round(n_50 / n_total * 100, 1), "% → 50 gen seçildi\n")


# ==============================================================================
# ÖZET
# ==============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("✅ 07_visualize.R TAMAMLANDI\n")
cat(strrep("=", 60), "\n")
cat("\nFigürler (", FIGURES_PATH, "):\n", sep = "")
cat("  Figure2_violin_IC50\n")
cat("  Figure3_heatmap\n")
cat("  Figure4_bar_best_model   (2 panel: RMSE + R²)\n")
cat("  Figure5_scatter_best_models\n")
cat("  Figure6_fold_boxplot     (ncol=3, düzgün layout)\n")
cat("\n📌 Figure 7 + 8 = 09_kegg_go_enrichment.R\n")
