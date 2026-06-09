# ==============================================================================
# 99_data_summary.R — Yayın Methods için detaylı veri özeti
# ==============================================================================
# Bu YARDIMCI script (pipeline'ın bir parçası değil) Methods bölümü için
# gereken veri boyutu istatistiklerini üretir:
#
#   1. Gen ekspresyon matrisi: kaç hücre, kaç gen
#   2. Her ilaç için: kaç IC50 değeri var
#   3. Filtreleme öncesi/sonrası: kaç hücre common, kaç IC50<200
#   4. Train/test/val split sayıları
#   5. Eşsiz ilaç-hücre çiftleri
#
# Pipeline'ın 02_preprocess'in çıktısını TAMAMLAR — daha detaylı bilgi verir.
#
# Bağımlı: 00_setup.R + 01_load_data_and_translate.R çıktıları
# Tahmini süre: ~2 dakika
#
# KULLANIM:
#   setwd("/Users/gizemtutkun/Desktop/Gizem_yl/scripts/")
#   source("99_data_summary.R")
# ==============================================================================

source("00_setup.R")
cat("\n[99_data_summary] Detaylı veri özeti çıkarılıyor...\n")

# Çevrilmiş gen verisini yükle
load(paste0(BASE_PATH, "gene_data_translated.RData"))


# ==============================================================================
# 1. GEN EKSPRESYON MATRISI ÖZETİ
# ==============================================================================

cat("\n--- 1. Gene Expression Matrix ---\n")
n_cells_total <- nrow(gene1)
n_genes_total <- ncol(gene1)
cat("  Total cell lines: ", n_cells_total, "\n", sep = "")
cat("  Total genes (after symbol mapping): ", n_genes_total, "\n", sep = "")
cat("  Successfully mapped to symbols: ",
    sum(!grepl("^ID_", colnames(gene1))), "\n", sep = "")
cat("  Unmapped (kept as ID_xxxx): ",
    sum(grepl("^ID_", colnames(gene1))), "\n", sep = "")


# ==============================================================================
# 2. HER İLAÇ İÇİN DETAYLI ÖZET
# ==============================================================================

cat("\n--- 2. Per-drug data summary (filtering before/after) ---\n")

drug_summary <- list()

for (k in seq_along(DRUG_FILES)) {
  drug_file   <- DRUG_FILES[k]
  drug_folder <- DRUG_FOLDER_NAMES[k]

  # Ham ilaç verisi
  drug <- as.data.frame(readxl::read_excel(drug_file, sheet = 1))
  rownames(drug) <- drug[, 1]
  drug <- drug[, -1, drop = FALSE]
  colnames(drug) <- "IC50"

  n_drug_raw <- nrow(drug)

  # Common cells (gene1 ile)
  common <- intersect(rownames(gene1), rownames(drug))
  n_common <- length(common)

  # IC50 dağılımı
  drug_common <- drug[common, , drop = FALSE]
  ic50_vec <- drug_common$IC50
  ic50_vec <- ic50_vec[!is.na(ic50_vec)]

  ic50_min  <- round(min(ic50_vec), 3)
  ic50_max  <- round(max(ic50_vec), 3)
  ic50_med  <- round(median(ic50_vec), 3)
  ic50_mean <- round(mean(ic50_vec), 3)

  # Filtreleme: IC50 < 200, IC50 > 0
  n_after_200 <- sum(ic50_vec < 200, na.rm = TRUE)
  n_positive  <- sum(ic50_vec > 0, na.rm = TRUE)
  n_kept      <- sum(ic50_vec < 200 & ic50_vec > 0, na.rm = TRUE)

  # Yüzde olarak kaybedildi
  pct_kept <- round(100 * n_kept / n_common, 1)

  # Train/Test/Val split (70/20/10) tahmini
  n_train <- round(n_kept * 0.70)
  n_test_val <- n_kept - n_train
  n_test  <- round(n_test_val * (2/3))
  n_val   <- n_test_val - n_test

  drug_summary[[drug_folder]] <- data.frame(
    Drug = drug_folder,
    Drug_File_N = n_drug_raw,
    Common_With_Gene = n_common,
    IC50_min = ic50_min,
    IC50_median = ic50_med,
    IC50_mean = ic50_mean,
    IC50_max = ic50_max,
    Filtered_IC50_lt200 = n_kept,
    Pct_Kept = pct_kept,
    Train_N = n_train,
    Test_N = n_test,
    Val_N = n_val
  )

  cat("  ", drug_folder, "\n", sep = "")
  cat("    Drug file IC50 values: ", n_drug_raw, "\n", sep = "")
  cat("    Common with gene expression: ", n_common, "\n", sep = "")
  cat("    IC50 range: [", ic50_min, ", ", ic50_max, "] µM, median = ",
      ic50_med, "\n", sep = "")
  cat("    After IC50<200 & >0 filter: ", n_kept, " (", pct_kept, "%)\n", sep = "")
  cat("    Train/Test/Val: ", n_train, " / ", n_test, " / ", n_val, "\n\n", sep = "")
}

summary_df <- bind_rows(drug_summary)


# ==============================================================================
# 3. KAYDET — Methods İçin Hazır Tablo
# ==============================================================================

write.csv(summary_df,
          paste0(BASE_PATH, "data_summary_for_methods.csv"),
          row.names = FALSE)

cat("\n--- 3. Methods bölümü için yazı önerisi ---\n\n")

# Methods için yazı şablonu (İngilizce — yayın için)
cat("📝 Suggested Methods text (English):\n\n")
cat("    \"Gene expression data comprised ", n_cells_total,
    " cancer cell lines × ", n_genes_total,
    " genes (Affymetrix probes mapped to gene symbols using\n",
    "    org.Hs.eg.db v3.18). Drug sensitivity data (IC50, µM) for 7 BET\n",
    "    inhibitors (RVX-208, PFI-1, JQ1, I-BRD9, OTX015, I-BET-762, AZD5153)\n",
    "    were obtained from CancerRxGene/GDSC. Following common-cell\n",
    "    intersection, we retained cell lines with IC50 < 200 µM and IC50 > 0\n",
    "    (filter rationale: Indrayanto et al., 2021), then log10-transformed\n",
    "    IC50 values. Each drug-specific dataset was randomly partitioned into\n",
    "    training (70%), test (20%), and validation (10%) sets using stratified\n",
    "    sampling (caret::createDataPartition) to preserve IC50 distribution.\"\n\n",
    sep = "")

# Sayısal özet için tablo
cat("📊 Per-drug numbers (paste into Table 1 or Methods):\n\n")
print(summary_df, row.names = FALSE)

cat("\n", strrep("=", 60), "\n", sep = "")
cat("✅ 99_data_summary.R TAMAMLANDI\n")
cat(strrep("=", 60), "\n", sep = "")
cat("\nÇıktı: ", BASE_PATH, "data_summary_for_methods.csv\n", sep = "")
cat("\n💡 Bu tablo Methods bölümü için kullanılabilir.\n")
cat("💡 Yayında 'Table 1: Data summary' olarak yer alabilir.\n")
