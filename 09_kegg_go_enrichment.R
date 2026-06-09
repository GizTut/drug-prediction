# ==============================================================================
# 09_kegg_go_enrichment.R — KEGG ve GO zenginleştirme analizi
# ==============================================================================
setwd(SCRIPTS_PATH)
source("00_setup.R")
cat("\n[09_kegg_go_enrichment] KEGG/GO analizi başlıyor...\n")

ENRICH_LOW  <- "#7BAAA3"
ENRICH_HIGH <- "#E8A09A"

# ==============================================================================
# 1. PAKET YÜKLEME
# ==============================================================================
bioc_pkgs <- c("clusterProfiler", "enrichplot", "DOSE")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("📦 Yükleniyor:", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}
for (pkg in c("igraph", "ggraph")) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
}
suppressPackageStartupMessages({
  library(clusterProfiler); library(enrichplot)
  library(DOSE); library(ggplot2); library(dplyr)
})
cat("✓ Paketler yüklendi\n")

pub_theme <- theme_bw(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle    = element_text(size = 10, hjust = 0.5, color = "gray40"),
    axis.text.y      = element_text(color = "black", size = 10),
    axis.title       = element_text(face = "bold")
  )

# ==============================================================================
# 2. GEN LİSTESİ
# ==============================================================================
cat("\n[2/5] Gen listesi yükleniyor...\n")
top10_fp <- paste0(BASE_PATH, "permutation_top10_per_drug_model.csv")
if (!file.exists(top10_fp)) stop("08 önce çalıştırılmalı: ", top10_fp)
top10            <- read.csv(top10_fp, stringsAsFactors = FALSE)
all_unique_genes <- unique(top10$Gene)
cat("  Eşsiz gen sayısı:", length(all_unique_genes), "\n")

# ==============================================================================
# 3. SYMBOL → ENTREZ ID
# ==============================================================================
cat("\n[3/5] Sembol → Entrez ID...\n")
conv <- AnnotationDbi::select(
  org.Hs.eg.db, keys = all_unique_genes,
  columns = c("SYMBOL", "ENTREZID"), keytype = "SYMBOL"
)
conv       <- conv[!duplicated(conv$SYMBOL), ]
conv       <- conv[!is.na(conv$ENTREZID), ]
all_entrez <- conv$ENTREZID
cat("  ", length(all_entrez), "/", length(all_unique_genes), " gen çevrildi\n", sep = "")

# ==============================================================================
# 4. KEGG
# ==============================================================================
cat("\n[4/5] KEGG zenginleştirme...\n")
kegg_res <- tryCatch({
  enrichKEGG(gene = all_entrez, organism = "hsa",
             pvalueCutoff = 0.05, qvalueCutoff = 0.10)
}, error = function(e) { cat("  ⚠️ KEGG hatası:", conditionMessage(e), "\n"); NULL })

kegg_readable <- NULL
if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
  kegg_readable <- setReadable(kegg_res, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  kegg_df <- as.data.frame(kegg_readable)
  cat("  ✓", nrow(kegg_df), "anlamlı KEGG yolağı\n")
  write.csv(kegg_df, paste0(BASE_PATH, "kegg_enrichment_all_drugs.csv"), row.names = FALSE)
  cat("\n  TOP 12 KEGG:\n")
  print(head(kegg_df[, c("Description","pvalue","qvalue","Count")], 12), row.names = FALSE)
}

# ==============================================================================
# 5. GO BP
# ==============================================================================
cat("\n[5/5] GO Biological Process...\n")
go_bp <- tryCatch({
  enrichGO(gene = all_entrez, OrgDb = org.Hs.eg.db, ont = "BP",
           pvalueCutoff = 0.05, qvalueCutoff = 0.10, readable = TRUE)
}, error = function(e) { cat("  ⚠️ GO BP hatası:", conditionMessage(e), "\n"); NULL })

go_bp_df <- NULL
if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
  go_bp_df <- as.data.frame(go_bp)
  cat("  ✓", nrow(go_bp_df), "anlamlı GO BP terimi\n")
  write.csv(go_bp_df[1:min(225, nrow(go_bp_df)), ],
            paste0(BASE_PATH, "go_BP_enrichment_all_drugs.csv"), row.names = FALSE)
  cat("\n  TOP 10 GO BP:\n")
  print(head(go_bp_df[, c("Description","pvalue","qvalue","Count")], 10), row.names = FALSE)
}

# ==============================================================================
# YARDIMCI: Manuel barplot fonksiyonu
# ==============================================================================
# barplot() ve dotplot() kendi scale_fill'lerini içeriden tanımlıyor.
# Üzerine scale_fill_gradient ekleyince "already present" uyarısı çıkıyor
# ve bizim renkler görmezden geliniyor.
# Çözüm: enrichResult objesini data.frame'e çevirip sıfırdan ggplot çizmek.

make_enrichment_barplot <- function(enrich_obj, top_n = 15,
                                    low_col = ENRICH_LOW,
                                    high_col = ENRICH_HIGH,
                                    title = "", subtitle = "") {
  df <- as.data.frame(enrich_obj)
  df <- df[order(df$pvalue), ]
  df <- head(df, top_n)
  # GeneRatio string'ini sayıya çevir (örn. "8/127" → 0.063)
  df$GeneRatioNum <- sapply(df$GeneRatio, function(x) {
    parts <- strsplit(x, "/")[[1]]
    as.numeric(parts[1]) / as.numeric(parts[2])
  })
  # Uzun isimleri kes (40 karakter)
  df$ShortDesc <- ifelse(nchar(df$Description) > 40,
                         paste0(substr(df$Description, 1, 38), "…"),
                         df$Description)
  df$ShortDesc <- factor(df$ShortDesc,
                         levels = rev(df$ShortDesc[order(df$pvalue)]))

  ggplot(df, aes(x = Count, y = ShortDesc, fill = p.adjust)) +
    geom_bar(stat = "identity", width = 0.75,
             color = "gray40", linewidth = 0.25) +
    scale_fill_gradient(low = low_col, high = high_col,
                        name = "Adjusted\np-value") +
    labs(title = title, subtitle = subtitle,
         x = "Gene count", y = NULL) +
    pub_theme +
    theme(axis.text.y = element_text(size = 10))
}

# ==============================================================================
# FIGURE 7: KEGG BARPLOT — renkler artık doğru
# ==============================================================================
if (!is.null(kegg_readable) && nrow(as.data.frame(kegg_readable)) > 0) {
  fig7 <- make_enrichment_barplot(
    kegg_readable, top_n = 12,
    title    = "KEGG Pathway Enrichment",
    subtitle = paste0(length(all_unique_genes), " recurrent predictor genes")
  )
  ggsave(paste0(FIGURES_PATH, "Figure7_KEGG_barplot.png"),
         fig7, width = 10, height = 7, dpi = 300, bg = "white")
  ggsave(paste0(FIGURES_PATH, "Figure7_KEGG_barplot.pdf"),
         fig7, width = 10, height = 7)
  cat("  ✓ Figure7_KEGG_barplot.png/pdf\n")
}

# ==============================================================================
# FIGURE 8: GO BP BARPLOT — renkler artık doğru
# ==============================================================================
if (!is.null(go_bp) && !is.null(go_bp_df)) {
  fig8 <- make_enrichment_barplot(
    go_bp, top_n = 20,
    title    = "GO Biological Process Enrichment",
    subtitle = paste0(length(all_unique_genes), " recurrent predictor genes")
  )
  ggsave(paste0(FIGURES_PATH, "Figure8_GO_BP_barplot.png"),
         fig8, width = 10, height = 9, dpi = 300, bg = "white")
  ggsave(paste0(FIGURES_PATH, "Figure8_GO_BP_barplot.pdf"),
         fig8, width = 10, height = 9)
  cat("  ✓ Figure8_GO_BP_barplot.png/pdf\n")
}

# ==============================================================================
# FIGURE S1: GO BP NETWORK
# linewidth uyarısı: enrichplot'un iç kodu — suppress ediyoruz
# ==============================================================================
if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 5) {
  fig_s1 <- tryCatch({
    suppressWarnings(
      emapplot(pairwise_termsim(go_bp), showCategory = 30) +
        labs(title = "GO BP Term Similarity Network") +
        pub_theme
    )
  }, error = function(e) {
    cat("  ⚠️ GO network hatası:", conditionMessage(e), "\n"); NULL
  })
  if (!is.null(fig_s1)) {
    ggsave(paste0(FIGURES_PATH, "FigureS1_GO_BP_network.png"),
           fig_s1, width = 12, height = 10, dpi = 300, bg = "white")
    cat("  ✓ FigureS1_GO_BP_network.png\n")
  }
}

# ==============================================================================
# FIGURE S2: KEGG CNETPLOT
# cex.params argümanı eski enrichplot versiyonlarında çalışmıyor.
# Yeni versiyonda cex_label_gene kullanılıyor; her ikisini de deneyen
# tryCatch bloğu yazıyoruz.
# ==============================================================================
if (!is.null(kegg_readable) && nrow(as.data.frame(kegg_readable)) > 3) {
  fig_s2 <- tryCatch({
    # Yeni enrichplot (≥1.18) API
    suppressWarnings(
      cnetplot(kegg_readable, showCategory = 8,
               cex_label_gene = 0.7, cex_label_category = 0.9)
    )
  }, error = function(e1) {
    tryCatch({
      # Eski API fallback
      suppressWarnings(
        cnetplot(kegg_readable, showCategory = 8)
      )
    }, error = function(e2) {
      cat("  ⚠️ KEGG cnetplot hatası:", conditionMessage(e2), "\n"); NULL
    })
  })

  if (!is.null(fig_s2)) {
    fig_s2 <- fig_s2 +
      labs(title = "KEGG Pathway–Gene Network") +
      pub_theme
    ggsave(paste0(FIGURES_PATH, "FigureS2_KEGG_gene_network.png"),
           fig_s2, width = 13, height = 11, dpi = 300, bg = "white")
    cat("  ✓ FigureS2_KEGG_gene_network.png\n")
  }
}

# ==============================================================================
# ÖZET
# ==============================================================================
cat("\n", strrep("=", 60), "\n", sep = "")
cat("✅ 09_kegg_go_enrichment.R TAMAMLANDI\n")
cat(strrep("=", 60), "\n")

cat("\n📋 MAKALEDEKİ SAYILARI KONTROL ET:\n")
cat("  Eşsiz tahmin geni:", length(all_unique_genes), "\n")

if (!is.null(go_bp_df)) {
  check_term <- function(pattern, label) {
    row <- go_bp_df[grepl(pattern, go_bp_df$Description, ignore.case = TRUE), ]
    if (nrow(row) > 0)
      cat(" ", label, "→ p =",
          formatC(row$pvalue[1], format = "e", digits = 2),
          "| q =", formatC(row$qvalue[1], format = "e", digits = 2),
          "| Count =", row$Count[1], "\n")
  }
  check_term("leukocyte prolif",     "leukocyte proliferation    ")
  check_term("lymphocyte prolif",    "lymphocyte proliferation   ")
  check_term("actin filament organ", "actin filament organization")
  check_term("regulation of B cell", "regulation of B cell prolif")
}
