# ==============================================================================
# 02_preprocess_per_drug.R — Her ilaç için veri ön işleme
# ==============================================================================
# Bu script her ilaç için:
#   1. Gen ekspresyon verisi + ilaç IC50 verisini birleştirir
#   2. IC50 < 200 µM filtresi
#   3. log10 dönüşümü
#   4. Train/Test/Val (70/20/10) bölme — caret::createDataPartition (stratified)
#   5. Min-max normalizasyon (sadece train üzerinden)
#   6. Pearson korelasyonu ile top 50/100/150 gen seçimi
#   7. 12 CSV kaydet (her ilaç için)
#
# Bağımlı: 00_setup.R + 01_load_data_and_translate.R
# Tahmini süre: ~5 dakika (7 ilaç toplam)
# ==============================================================================
setwd(SCRIPTS_PATH)

source("00_setup.R")
cat("\n[02_preprocess] Her ilaç için ön işleme başlıyor...\n")

# Adım 01'den çevrilmiş gen verisini yükle
load(paste0(BASE_PATH, "gene_data_translated.RData"))
cat("✓ gene_data_translated.RData yüklendi\n")
cat("  Boyut:", nrow(gene1), "hücre ×", ncol(gene1), "gen (semboller)\n\n")


# ==============================================================================
# TEK İLAÇ İÇİN PREPROCESSİNG FONKSİYONU
# ==============================================================================

preprocess_drug <- function(drug_file, drug_folder) {

  cat(strrep("-", 60), "\n", sep = "")
  cat("İLAÇ: ", drug_folder, "\n", sep = "")
  cat(strrep("-", 60), "\n", sep = "")

  drug_path <- paste0(BASE_PATH, drug_folder, "/")
  dir.create(drug_path, showWarnings = FALSE, recursive = TRUE)

  # ---- 1. İlaç verisi yükle ----
  drug <- as.data.frame(readxl::read_excel(drug_file, sheet = 1))
  rownames(drug) <- drug[, 1]
  drug <- drug[, -1, drop = FALSE]
  colnames(drug) <- "IC50"

  # ---- 2. Gen + ilaç birleştir ----
  common_indices <- intersect(rownames(gene1), rownames(drug))
  gene_filt <- gene1[common_indices, , drop = FALSE]
  drug_filt <- drug[common_indices, , drop = FALSE]
  merged    <- cbind(gene_filt, drug_filt)

  # ---- 3. IC50 < 200 filtresi + log10 ----
  merged    <- subset(merged, IC50 < 200)
  merged    <- merged[merged$IC50 > 0, ]
  merged$IC50 <- log10(merged$IC50)

  # Tamamen NA olan sütunları kaldır
  merged <- merged[, colSums(is.na(merged)) != nrow(merged)]

  cat("  Birleştirilmiş veri: ", nrow(merged), " hücre × ",
      ncol(merged) - 1, " gen\n", sep = "")

  # ---- 4. Train/Test/Val bölme (70/20/10) ----
  set.seed(1234)
  df <- as.data.frame(merged[sample(1:nrow(merged)), ])

  set.seed(GLOBAL_SEED)
  train_idx <- createDataPartition(df$IC50, p = 0.70, list = FALSE)
  train_data <- df[train_idx, ]
  temp_data  <- df[-train_idx, ]

  test_idx  <- createDataPartition(temp_data$IC50, p = 2/3, list = FALSE)
  test_data <- temp_data[test_idx, ]
  val_data  <- temp_data[-test_idx, ]

  # ---- 5. Min-max normalizasyon (TRAIN üzerinden) ----
  train_y <- train_data$IC50
  test_y  <- test_data$IC50
  val_y   <- val_data$IC50

  train_x <- train_data[, setdiff(colnames(train_data), "IC50"), drop = FALSE]
  test_x  <- test_data[,  setdiff(colnames(test_data),  "IC50"), drop = FALSE]
  val_x   <- val_data[,   setdiff(colnames(val_data),   "IC50"), drop = FALSE]

  mins   <- apply(train_x, 2, min, na.rm = TRUE)
  maxs   <- apply(train_x, 2, max, na.rm = TRUE)
  ranges <- maxs - mins
  ranges[ranges == 0] <- 1

  scale_minmax <- function(x) {
    as.data.frame(sweep(sweep(x, 2, mins, FUN = "-"), 2, ranges, FUN = "/"))
  }

  train_scaled <- cbind(IC50 = train_y, scale_minmax(train_x))
  test_scaled  <- cbind(IC50 = test_y,  scale_minmax(test_x))
  val_scaled   <- cbind(IC50 = val_y,   scale_minmax(val_x))

  # ---- 6. Pearson özellik seçimi (TRAIN üzerinden) ----
  train_x_only <- train_scaled[, -1, drop = FALSE]
  train_y_only <- train_scaled$IC50

  correlations <- apply(train_x_only, 2, function(x) {
    cor(x, train_y_only, method = "pearson", use = "complete.obs")
  })

  cor_table <- data.frame(
    Gene = names(correlations),
    Correlation = correlations,
    stringsAsFactors = FALSE
  )
  cor_table <- cor_table[order(abs(cor_table$Correlation), decreasing = TRUE), ]

  # ---- 7. Her gen sayısı için 12 CSV kaydet ----
  for (n in TOP_N_LIST) {
    top_genes <- cor_table$Gene[1:n]

    train_n <- train_scaled[, c("IC50", top_genes), drop = FALSE]
    test_n  <- test_scaled[,  c("IC50", top_genes), drop = FALSE]
    val_n   <- val_scaled[,   top_genes, drop = FALSE]

    val_ic50_n <- data.frame(
      IC50 = val_scaled[rownames(val_n), "IC50"],
      row.names = rownames(val_n)
    )

    write.csv(train_n,
              paste0(drug_path, "data_train_", n, ".csv"), row.names = TRUE)
    write.csv(test_n,
              paste0(drug_path, "data_test_", n, ".csv"), row.names = TRUE)
    write.csv(val_n,
              paste0(drug_path, "data_val_", n, ".csv"), row.names = TRUE)
    write.csv(val_ic50_n,
              paste0(drug_path, "val", n, "_ic50.csv"), row.names = TRUE)
  }

  cat("  ✓ 12 CSV kaydedildi (50/100/150 × train/test/val/val_ic50)\n")

  return(list(
    drug = drug_folder,
    n_cells = nrow(merged),
    n_genes = ncol(merged) - 1,
    train_n = nrow(train_data),
    test_n  = nrow(test_data),
    val_n   = nrow(val_data)
  ))
}


# ==============================================================================
# 7 İLAÇ İÇİN PREPROCESS
# ==============================================================================

cat("\n", strrep("=", 60), "\n", sep = "")
cat("7 ilaç için preprocessing başlıyor...\n")
cat(strrep("=", 60), "\n", sep = "")

preprocess_summary <- list()

for (k in seq_along(DRUG_FILES)) {
  result <- preprocess_drug(DRUG_FILES[k], DRUG_FOLDER_NAMES[k])
  preprocess_summary[[DRUG_FOLDER_NAMES[k]]] <- result
}

# Özet tablo
summary_df <- bind_rows(preprocess_summary)
write.csv(summary_df,
          paste0(BASE_PATH, "preprocessing_summary.csv"),
          row.names = FALSE)

cat("\n", strrep("=", 60), "\n", sep = "")
cat("✅ 02_preprocess_per_drug.R TAMAMLANDI\n")
cat(strrep("=", 60), "\n", sep = "")
print(summary_df, row.names = FALSE)
cat("\n📌 Sonraki adım: 03_run_nested_cv.R (UZUN — 9-13 saat)\n")
