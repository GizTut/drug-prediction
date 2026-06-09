# ==============================================================================
# 01_load_data_and_translate.R — Veri yükle ve gen sembollerine çevir
# ==============================================================================
# Bu script:
#   1. GeneExp.csv'yi yükler
#   2. Sütun isimlerindeki Entrez ID'leri (örn `2017.at`) gen sembolüne çevirir
#      (örn `2017.at` → `CELSR1`)
#   3. Tüm 7 ilacın IC50 verisini yükler
#   4. Hepsini bir RData dosyasına kaydeder (sonraki scriptler için)
#
# AVANTAJ: Bundan sonra TÜM CSV'lerde, modellerde, görsellerde gen sembolleri
# kullanılır — Entrez ID değil. Hakemler net görür.
#
# Bağımlı: 00_setup.R (BASE_PATH, kütüphaneler)
# Tahmini süre: ~2 dakika
# ==============================================================================
setwd(SCRIPTS_PATH)

source("00_setup.R")
cat("\n[01_load_data] Veri yükleme ve gen çevirisi başlıyor...\n")


# ==============================================================================
# 1. GEN EKSPRESYON VERİSİNİ YÜKLE
# ==============================================================================

cat("\n[1/4] GeneExp.csv okunuyor...\n")
gene0 <- read.csv("GeneExp.csv", row.names = 1, header = TRUE,
                  check.names = FALSE)
gene1 <- t(gene0)

cat("    Boyut: ", nrow(gene1), " hücre × ", ncol(gene1), " gen\n", sep = "")
cat("    Örnek sütun isimleri (ilk 5):\n")
print(head(colnames(gene1), 5))


# ==============================================================================
# 2. ENTREZ ID → GEN SEMBOLÜ ÇEVİRİSİ (TEK SEFERLİK!)
# ==============================================================================

cat("\n[2/4] Sütun isimleri gen sembollerine çevriliyor...\n")

original_ids <- colnames(gene1)
entrez_ids   <- extract_entrez_id(original_ids)

# Eşsiz ID'leri tek seferde çevir (verimli)
unique_ids <- unique(entrez_ids)
cat("    ", length(unique_ids), " eşsiz Entrez ID çevriliyor...\n", sep = "")

id_to_symbol <- entrez_to_symbol(unique_ids)
names(id_to_symbol) <- unique_ids

# Tüm sütunlar için sembolleri ata
gene_symbols <- id_to_symbol[entrez_ids]

# Çevrilmemiş olanları say
unknown_count <- sum(grepl("^ID_", gene_symbols))
known_count   <- length(gene_symbols) - unknown_count

cat("    ✅ Çevrildi   : ", known_count, " gen\n", sep = "")
cat("    ⚠️  Çevrilemedi: ", unknown_count, " gen (eski/çekilmiş ID'ler)\n", sep = "")

# Yinelenen sembolleri ele al (nadir ama olabilir, örn 2 farklı Entrez ID aynı sembol)
if (any(duplicated(gene_symbols))) {
  dup_count <- sum(duplicated(gene_symbols))
  cat("    ⚠️  Yinelenen sembol: ", dup_count, " (suffix ekleniyor)\n", sep = "")

  # Yinelenenlere _2, _3 gibi suffix ekle
  gene_symbols <- make.unique(gene_symbols, sep = "_")
}

# Sütun isimlerini güncelle
colnames(gene1) <- gene_symbols

cat("    Örnek YENİ sütun isimleri (ilk 5):\n")
print(head(colnames(gene1), 5))


# ==============================================================================
# 3. ID-SYMBOL EŞLEME TABLOSU (REFERANS — kalıcı kayıt)
# ==============================================================================

cat("\n[3/4] ID-symbol eşleme tablosu kaydediliyor...\n")

mapping_table <- data.frame(
  Original_ID = original_ids,
  Entrez_ID   = entrez_ids,
  Symbol      = gene_symbols,
  stringsAsFactors = FALSE
)

write.csv(mapping_table,
          paste0(BASE_PATH, "id_to_symbol_mapping.csv"),
          row.names = FALSE)
cat("    💾 id_to_symbol_mapping.csv\n")


# ==============================================================================
# 4. RData OLARAK KAYDET (SONRAKI SCRIPTLER İÇİN)
# ==============================================================================

cat("\n[4/4] gene_data_translated.RData kaydediliyor...\n")

save(gene1, mapping_table,
     file = paste0(BASE_PATH, "gene_data_translated.RData"))

cat("    💾 gene_data_translated.RData\n")
cat("       İçerik: gene1 (sembolik isimlerle), mapping_table\n")


cat("\n", strrep("=", 60), "\n", sep = "")
cat("✅ 01_load_data_and_translate.R TAMAMLANDI\n")
cat(strrep("=", 60), "\n", sep = "")
cat("\nÜretilen dosyalar:\n")
cat("  ", BASE_PATH, "gene_data_translated.RData\n", sep = "")
cat("  ", BASE_PATH, "id_to_symbol_mapping.csv\n", sep = "")
cat("\n📌 Sonraki adım: 02_preprocess_per_drug.R\n")
