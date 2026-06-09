# ==============================================================================
# 00_RUN_ALL.R — Tüm pipeline'ı tek tıkla çalıştır
# ==============================================================================
# Bu script tüm scriptleri sırayla çağırır.
#
# Tahmini toplam süre: ~10-14 saat (bütünüyle)
#   01: ~2 dakika
#   02: ~5 dakika
#   03: ~9-13 saat (UZUN)
#   04: ~10 saniye
#   05: ~10 dakika
#   06: ~2 dakika
#   07: ~5 dakika
#
# KULLANIM:
#   1. Mac'te uyku modunu kapat: System Settings → Battery → Prevent sleep
#   2. Bu dosyayı çalıştır: Cmd+A → Cmd+Enter
#   3. ☕ uyu / dışarı çık — sabah hepsi hazır
#
# ALTERNATİF: Her scripti TEK TEK çalıştırmak için
#   source("01_load_data_and_translate.R")  # sadece bunu çalıştır
#
# ⚠️  Eğer pipeline ortasında bir hata çıkarsa, kod durur. Hatayı düzeltip
#     o scriptten sonra başlayabilirsin (sıfırdan başlamak gerekmez):
#     source("03_run_nested_cv.R")
#     source("04_statistical_tests.R")
#     ... vb.
# ==============================================================================

cat("\n", strrep("=", 70), "\n", sep = "")
cat("MASTER PIPELINE BAŞLIYOR\n")
cat(strrep("=", 70), "\n", sep = "")
cat("Toplam tahmini süre: 10-14 saat\n")
cat("Mac uyku modunu kapattınız mı? (System Settings → Battery)\n\n")

master_start <- Sys.time()

# Scripts klasöründe çalıştığımızdan emin ol
# (Eğer scripts/ alt klasöründeysen, BASE_PATH/scripts/ olur)
if (basename(getwd()) == "scripts") {
  setwd("..")
}

# Tüm script dosyalarının olduğunu doğrula
script_files <- c(
  "01_load_data_and_translate.R",
  "02_preprocess_per_drug.R",
  "03_run_nested_cv.R",
  "04_statistical_tests.R",
  "05_train_final_models.R",
  "06_test_and_validate.R",
  "07_visualize.R"
)

# Script dosyaları "scripts/" alt klasöründe mi yoksa burada mı?
if (file.exists("scripts/01_load_data_and_translate.R")) {
  scripts_dir <- "scripts/"
} else if (file.exists("01_load_data_and_translate.R")) {
  scripts_dir <- ""
} else {
  stop("Script dosyaları bulunamadı! scripts/ klasörünü kontrol et.")
}

cat("📂 Scripts klasörü: '", scripts_dir, "'\n\n", sep = "")

BASE_PATH <- "/Users/gizemtutkun/Desktop/Gizem_yl/"

# Alt klasörler
SCRIPTS_PATH <- paste0(BASE_PATH, "scripts/")
FIGURES_PATH <- paste0(BASE_PATH, "figures/")

# Sırayla çalıştır
for (script in script_files) {
  full_path <- paste0(scripts_dir, script)
  cat("\n▶️  Çalıştırılıyor: ", script, "\n", sep = "")

  step_start <- Sys.time()
  tryCatch({
    source(full_path, chdir = FALSE)
  }, error = function(e) {
    cat("\n❌ HATA — ", script, ":\n", sep = "")
    cat("   ", conditionMessage(e), "\n")
    cat("\n⚠️  Bu hatayı düzelttikten sonra şu komutla devam edebilirsin:\n")
    cat("    source('", full_path, "')\n", sep = "")
    stop("Pipeline durdu")
  })
  step_end <- Sys.time()
  step_min <- round(as.numeric(difftime(step_end, step_start, units = "mins")), 1)
  cat("✅ ", script, " bitti (", step_min, " dk)\n", sep = "")
}

master_end <- Sys.time()
total_h <- round(as.numeric(difftime(master_end, master_start, units = "hours")), 2)

cat("\n", strrep("=", 70), "\n", sep = "")
cat("🎉 TÜM PİPELİNE BİTTİ!\n")
cat(strrep("=", 70), "\n", sep = "")
cat("Toplam süre: ", total_h, " saat\n", sep = "")
cat("\nÇıktılar:\n")
cat("  ", BASE_PATH, "ALL_DRUGS_TEST_SUMMARY.csv\n", sep = "")
cat("  ", BASE_PATH, "ALL_DRUGS_GENE_SELECTION.csv\n", sep = "")
cat("  ", BASE_PATH, "figures/  → 5 PNG + tablolar\n", sep = "")
cat("  ", BASE_PATH, "{her ilaç}/results/  → her ilaç için detaylı sonuçlar\n", sep = "")
