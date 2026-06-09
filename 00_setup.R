# ==============================================================================
# 00_setup.R — Kütüphaneler, yollar, sabitler
# ==============================================================================
# Bu script tüm pipeline için ortak değişkenleri ve kütüphaneleri yükler.
# Diğer tüm scriptler bu dosyayı en başta çalıştırır:
#     source("00_setup.R")
#
# Tahmini süre: ~30 saniye
# ==============================================================================

cat("\n[00_setup] Pipeline kuruluyor...\n")

# ---- Kütüphaneler ----
suppressPackageStartupMessages({
  library(readxl)         # Excel okuma
  library(dplyr)          # Veri işleme
  library(tidyr)          # pivot_wider
  library(stats)          # İstatistik
  library(caret)          # ML wrapper
  library(class)          # KNN
  library(Metrics)        # MAE, RMSE
  library(elasticnet)     # ENet
  library(nnet)           # Neural Net
  library(kernlab)        # SVM
  library(xgboost)        # XGBoost
  library(kknn)           # k-NN regression
  library(randomForest)   # Random Forest
})

# ---- Bioconductor (gen sembol çevirisi) ----
# Sadece gen isim çevirisi için — varsa yükle, yoksa kur
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  cat("📦 org.Hs.eg.db yükleniyor (~3-5 dk, sadece ilk kez)...\n")
  BiocManager::install("org.Hs.eg.db", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

options(error = NULL)

# ---- Ana yollar ----
# Bu klasör Gizem_yl'in tam yolu — tek değişken, her yerde kullanılır.
# Yıllar sonra Gizem_yl klasörü farklı bir yere taşınırsa SADECE buraya bak.
BASE_PATH <- "/Users/gizemtutkun/Desktop/Gizem_yl/"

# Alt klasörler
SCRIPTS_PATH <- paste0(BASE_PATH, "scripts/")
FIGURES_PATH <- paste0(BASE_PATH, "figures/")

# Klasör yoksa oluştur
dir.create(FIGURES_PATH, showWarnings = FALSE, recursive = TRUE)

# Çalışma dizini = Gizem_yl
setwd(BASE_PATH)

# ---- Sabitler ----
DRUG_FILES <- c("RVX-208.xls", "PFI-1.xls", "JQ1.xls", "I-BRD9.xls",
                "OTX015.xls", "I-BET-762.xls", "AZD5153.xls")

DRUG_FOLDER_NAMES <- gsub("\\.xls$", "", DRUG_FILES)

TOP_N_LIST <- c(50, 100, 150)

OUTER_FOLDS <- 5
INNER_FOLDS <- 5
GLOBAL_SEED <- 123

# Yayın için model isim sırası (Fig 5 vb.)
MODEL_ORDER <- c("KNN", "Elastic Net", "Random Forest",
                 "SVM", "XGBoost", "Neural Net")

# Yayın için model renkleri (her figürde aynı)
MODEL_COLORS <- c(
  "KNN"           = "#A3C4BC",   # Pastel Blue
  "Elastic Net"   = "#F5CAC3",   # Soft Pink
  "Random Forest" = "#D9EAD3",   # Mint Green
  "SVM"           = "#C4A8B8",   # Light Lavender (koyulaştırılmış)
  "XGBoost"       = "#FFE3B3",   # Pale Lemon
  "Neural Net"    = "#D4A5A5"    # Dusty Rose
)

# ---- Yardımcı fonksiyon: Entrez ID → Symbol ----
extract_entrez_id <- function(gene_str) {
  cleaned <- gsub("`", "", gene_str)
  cleaned <- gsub("\\.at$", "", cleaned)
  return(cleaned)
}

entrez_to_symbol <- function(entrez_ids) {
  symbols <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys      = as.character(entrez_ids),
    column    = "SYMBOL",
    keytype   = "ENTREZID",
    multiVals = "first"
  )
  # Çevrilemeyenleri orijinal ID ile bırak
  symbols[is.na(symbols)] <- paste0("ID_", entrez_ids[is.na(symbols)])
  return(symbols)
}

cat("✅ Setup tamamlandı.\n")
cat("   Ana yol: ", BASE_PATH, "\n", sep = "")
cat("   İlaç sayısı: ", length(DRUG_FOLDER_NAMES), "\n", sep = "")
cat("   Gen sayıları: ", paste(TOP_N_LIST, collapse = ", "), "\n", sep = "")
