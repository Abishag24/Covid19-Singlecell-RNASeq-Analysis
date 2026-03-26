# ══════════════════════════════════════════════════════════════════════════════
#  COVID-19 Bronchoalveolar Lavage (BALF) Single-Cell RNA-seq Analysis · GSE149689
#  Tissue : Bronchoalveolar lavage fluid
#  Groups  : Severe COVID · Mild/Asymptomatic COVID · Influenza · Healthy
# ══════════════════════════════════════════════════════════════════════════════
# REQUIREMENTS:
#   install.packages(c("Seurat","dplyr","ggplot2","patchwork","Matrix"))
#   BiocManager::install("GEOquery")
# ══════════════════════════════════════════════════════════════════════════════

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
  library(GEOquery)
})

set.seed(42)

# ── 1. Paths ──────────────────────────────────────────────────────────────────

DATA_DIR    <- "C:/Users/viswa/Documents/Downloads/GSE149689"
RESULTS_DIR <- file.path(DATA_DIR, "results")
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)


# ── 2. Load Count Matrix ──────────────────────────────────────────────────────

mtx      <- readMM(file.path(DATA_DIR, "matrix.mtx"))
features <- read.table(file.path(DATA_DIR, "features.tsv"),
                       header = FALSE, stringsAsFactors = FALSE)
barcodes <- read.table(file.path(DATA_DIR, "barcodes.tsv"),
                       header = FALSE, stringsAsFactors = FALSE)

features[, 2] <- make.unique(features[, 2])
rownames(mtx)  <- features[, 2]
colnames(mtx)  <- barcodes[, 1]

cat("📦 Raw matrix:", nrow(mtx), "genes ×", ncol(mtx), "cells\n")


# ── 3. Create Seurat Object ───────────────────────────────────────────────────

seu <- CreateSeuratObject(
  counts       = mtx,
  project      = "COVID19_BALF_GSE149689",
  min.cells    = 3,     # gene must appear in at least 3 cells
  min.features = 200    # cell must express at least 200 genes
)
rm(mtx); gc()
cat("✅ Seurat object created:", ncol(seu), "cells\n")


# ── 4. Add Condition Metadata from GEO ───────────────────────────────────────

gse  <- getGEO("GSE149689", GSEMatrix = TRUE)
meta <- pData(phenoData(gse[[1]]))

meta_clean <- data.frame(
  sample_num = seq_len(nrow(meta)),
  title      = meta$title,
  condition  = case_when(
    grepl("severe",           meta$`subject status:ch1`, ignore.case = TRUE) ~ "Severe_COVID",
    grepl("mild|asymptomatic",meta$`subject status:ch1`, ignore.case = TRUE) ~ "Mild_COVID",
    grepl("influenza",        meta$`subject group:ch1`,  ignore.case = TRUE) ~ "Influenza",
    grepl("healthy",          meta$`subject group:ch1`,  ignore.case = TRUE) ~ "Healthy",
    TRUE ~ "Other"
  ),
  stringsAsFactors = FALSE
)

cat("\nSample → condition mapping:\n")
print(meta_clean[, c("sample_num", "title", "condition")])
cat("\nCondition counts:\n")
print(table(meta_clean$condition))

# Extract sample ID from barcode suffix (e.g., "AAACCTGAG-3" → sample 3)
# This maps each cell to its original GEO sampleseu$sample_num <- as.integer(sub(".*-", "", colnames(seu)))
seu$condition  <- meta_clean$condition[seu$sample_num]
seu$sample_id  <- paste0("S", seu$sample_num)

cat("\nCell counts per condition:\n")
print(table(seu$condition, useNA = "always"))


# ── 5. QC + Filtering ────────────────────────────────────────────────────────

seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Save QC plot 
p_qc <- VlnPlot(seu,
                features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                ncol = 3, pt.size = 0) &
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave(file.path(RESULTS_DIR, "01_QC.png"), p_qc, width = 12, height = 4, dpi = 150)

# Keep cells with 300–6000 genes and <15% mitochondrial reads
seu <- subset(seu,
              subset = nFeature_RNA > 300 &
                nFeature_RNA < 6000 &
                percent.mt   < 15)

cat("✅ After QC:", ncol(seu), "cells retained\n")


# ── 6. Normalize + Highly Variable Genes ─────────────────────────────────────

seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)

top10 <- head(VariableFeatures(seu), 10)
p_hvg <- LabelPoints(
  plot   = VariableFeaturePlot(seu),
  points = top10, repel = TRUE, xnudge = 0, ynudge = 0
)
ggsave(file.path(RESULTS_DIR, "02_HVG.png"), p_hvg, width = 8, height = 5, dpi = 150)
cat("✅ Normalization + HVG done\n")


# ── 7. Scale + PCA ───────────────────────────────────────────────────────────

# Scale only the 2,000 HVGs — faster than all genes, sufficient for PCA
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu,   features = VariableFeatures(seu), npcs = 30, verbose = FALSE)

p_elbow <- ElbowPlot(seu, ndims = 30) +
  geom_vline(xintercept = 15, linetype = "dashed", colour = "firebrick") +
  ggtitle("Elbow plot — dashed line at PC 15")
ggsave(file.path(RESULTS_DIR, "03_Elbow.png"), p_elbow, width = 7, height = 4, dpi = 150)


# ── 8. Clustering + UMAP ─────────────────────────────────────────────────────

seu <- FindNeighbors(seu, dims = 1:15, verbose = FALSE)
seu <- FindClusters( seu, resolution = 0.5, verbose = FALSE)
seu <- RunUMAP(      seu, dims = 1:15, verbose = FALSE)

# Clean cluster factor — prevents "duplicated factor level" errors downstream
seu$seurat_clusters <- factor(as.integer(as.character(Idents(seu))))
Idents(seu) <- seu$seurat_clusters

cat("✅ Clustering done —", length(levels(seu)), "clusters\n")
print(table(Idents(seu)))

p_umap <- DimPlot(seu, reduction = "umap", label = TRUE, pt.size = 0.3) +
  NoLegend() + ggtitle("UMAP — Seurat clusters")
ggsave(file.path(RESULTS_DIR, "04_UMAP_clusters.png"),
       p_umap, width = 7, height = 6, dpi = 150)


# ── 9. Canonical Marker Plots ─────────────────────────────────────────────────

marker_genes <- c(
  "CD3D", "CD4", "CD8A",    # T cells
  "MS4A1", "CD79A",          # B cells
  "LYZ", "FCN1", "IL1B",    # Monocytes
  "NKG7", "GNLY", "TRDC",   # NK / gamma-delta T cells
  "PF4", "PPBP",             # Platelets
  "SIGLEC1", "IFI44L"        # Plasmacytoid / ISG monocytes
)

p_feature <- FeaturePlot(seu,
                         features  = c("CD3D", "MS4A1", "LYZ",
                                       "NKG7", "GNLY", "SIGLEC1"),
                         ncol      = 3,
                         reduction = "umap",
                         pt.size   = 0.2) &
  theme(axis.text = element_blank(), axis.ticks = element_blank())
ggsave(file.path(RESULTS_DIR, "05_feature_plots.png"),
       p_feature, width = 12, height = 8, dpi = 150)

p_dot <- DotPlot(seu, features = marker_genes) +
  RotatedAxis() +
  ggtitle("Canonical marker expression per cluster")
ggsave(file.path(RESULTS_DIR, "06_dotplot.png"),
       p_dot, width = 14, height = 5, dpi = 150)


# ── 10. Find All Markers ──────────────────────────────────────────────────────

cat("⏳ FindAllMarkers running...\n")

markers <- FindAllMarkers(
  seu,
  only.pos            = TRUE,   # positive markers only
  min.pct             = 0.25,   # expressed in at least 25% of cells
  logfc.threshold     = 0.5,    # minimum log2 fold change
  test.use            = "wilcox",
  max.cells.per.ident = 200,    # subsample for speed — no accuracy loss
  verbose             = TRUE
)

top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>%
  ungroup()

write.csv(markers, file.path(RESULTS_DIR, "markers_all.csv"),  row.names = FALSE)
write.csv(top5,    file.path(RESULTS_DIR, "markers_top5.csv"), row.names = FALSE)
cat("✅ Markers saved —", nrow(markers), "DE genes\n")


# ── 11. Cell Type Annotation ──────────────────────────────────────────────────

# Annotation derived from top marker genes (verified against markers_top5.csv):
#
# Cluster 0  → CD8 T cell            CD8A, CD8B, GZMK     cytotoxic T cells
# Cluster 1  → Inflamed Monocyte     IL1B, PTGS2, OSM     COVID-19 inflammatory signature
# Cluster 2  → Naive CD4 T cell      MAL, LEF1, CD28      naive / resting T cells
# Cluster 3  → NK cell               MYOM2, SH2D1B        natural killer cells
# Cluster 4  → Naive B cell          TCL1A, IGHD, FCER2   naive / transitional B cells
# Cluster 5  → NK/gdT cell           KLRC2, TRDC, KLRF1   gamma-delta T / NK cells
# Cluster 6  → Activated Monocyte    CXCL2, CXCL3, EREG  chemokine-high monocytes
# Cluster 7  → ISG Monocyte          SIGLEC1, IFI44L      interferon-stimulated monocytes
# Cluster 8  → CD16 Monocyte         CDKN1C, C1QA         non-classical patrolling monocytes
# Cluster 9  → Platelet              PTCRA, CLEC1B, ESAM  megakaryocyte / platelet lineage
# Cluster 10 → Plasmablast           TNFRSF13B, SPIB, BLK transitional antibody-secreting B
# Cluster 11 → CD14 Monocyte         S100A9, IFITM3       classical monocytes
# Cluster 12 → Early Monocyte        EGR1, FOS, FLT3      immediate-early response monocytes
# Cluster 13 → CD4 Effector T cell   LAG3, TRAC, LAIR2    exhausted / effector CD4 T cells
# Cluster 14 → Cycling T cell        TYMS, CDC45, PKMYT1  proliferating T cells (S/G2M phase)
# Cluster 15 → Neutrophil            HP, RETN, SERPINB2   low-density granulocytes
# Cluster 16 → ISG cell              IFI27, ANKRD22, CCL2 interferon-stimulated signature
# Cluster 17 → Classical DC (cDC2)   FCER1A, CD1C, CD1E   myeloid dendritic cells
# Cluster 18 → Erythrocyte           GYPA, GYPB, SPTA1    red blood cell contamination
# Cluster 19 → Plasma cell           IGHA1, IGHG1, JCHAIN antibody-secreting plasma cells
# Cluster 20 → Memory B cell         CD19, EBF1, FCRL2    class-switched memory B cells
# Cluster 21 → Regulatory T cell     TNFRSF4, IL2RA       CD25+ FOXP3+ Tregs
# Cluster 22 → Pre-B cell            PAX5, POU2AF1, BLNK  B cell precursors

cell_type_map <- c(
  "0"  = "CD8 T cell",
  "1"  = "Inflamed Monocyte",
  "2"  = "Naive CD4 T cell",
  "3"  = "NK cell",
  "4"  = "Naive B cell",
  "5"  = "NK/gdT cell",
  "6"  = "Activated Monocyte",
  "7"  = "ISG Monocyte",
  "8"  = "CD16 Monocyte",
  "9"  = "Platelet",
  "10" = "Plasmablast",
  "11" = "CD14 Monocyte",
  "12" = "Early Monocyte",
  "13" = "CD4 Effector T cell",
  "14" = "Cycling T cell",
  "15" = "Neutrophil",
  "16" = "ISG cell",
  "17" = "Classical DC",
  "18" = "Erythrocyte",
  "19" = "Plasma cell",
  "20" = "Memory B cell",
  "21" = "Regulatory T cell",
  "22" = "Pre-B cell"
)

seu <- RenameIdents(seu, cell_type_map)
seu$cell_type <- as.character(Idents(seu))

p_annot <- DimPlot(seu, reduction = "umap", label = TRUE,
                   label.size = 3.5, pt.size = 0.3, repel = TRUE) +
  NoLegend() + ggtitle("UMAP — annotated cell types")
ggsave(file.path(RESULTS_DIR, "07_UMAP_annotated.png"),
       p_annot, width = 10, height = 8, dpi = 150)
cat("✅ Annotation done\n")
print(sort(table(seu$cell_type)))


# ── 12. Condition UMAPs + Composition ────────────────────────────────────────

p_condition <- DimPlot(seu, reduction = "umap", group.by = "condition",
                       pt.size = 0.2,
                       cols = c("Severe_COVID" = "#E24B4A",
                                "Mild_COVID"   = "#EF9F27",
                                "Influenza"    = "#7F77DD",
                                "Healthy"      = "#1D9E75")) +
  ggtitle("UMAP — coloured by disease group")
ggsave(file.path(RESULTS_DIR, "08_UMAP_condition.png"),
       p_condition, width = 8, height = 6, dpi = 150)

comp_df <- seu@meta.data %>%
  group_by(condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(pct = n / sum(n) * 100)

p_comp <- ggplot(comp_df, aes(x = condition, y = pct, fill = cell_type)) +
  geom_bar(stat = "identity") +
  labs(title = "Cell type composition per disease group",
       x = "Condition", y = "% of cells", fill = "Cell type") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(RESULTS_DIR, "09_composition_bar.png"),
       p_comp, width = 10, height = 6, dpi = 150)

write.csv(comp_df, file.path(RESULTS_DIR, "cell_composition.csv"), row.names = FALSE)


# ── 13. Differential Expression ──────────────────────────────────────────────

cat("\nCondition labels in object:\n")
print(table(seu$condition))

Idents(seu) <- "condition"

# Comparison 1: Severe COVID vs Healthy controls
cat("⏳ DEG: Severe COVID vs Healthy...\n")
deg_sv_hl <- FindMarkers(
  seu,
  ident.1             = "Severe_COVID",
  ident.2             = "Healthy",
  min.pct             = 0.25,
  logfc.threshold     = 0.5,
  max.cells.per.ident = 500,
  verbose             = FALSE
)
deg_sv_hl$gene      <- rownames(deg_sv_hl)
deg_sv_hl$direction <- ifelse(deg_sv_hl$avg_log2FC > 0, "Up in Severe", "Up in Healthy")
write.csv(deg_sv_hl, file.path(RESULTS_DIR, "DEG_SevereCOVID_vs_Healthy.csv"),
          row.names = TRUE)
cat("✅ Severe vs Healthy:", nrow(deg_sv_hl), "DEGs |",
    sum(deg_sv_hl$avg_log2FC > 0), "up in Severe |",
    sum(deg_sv_hl$avg_log2FC < 0), "up in Healthy\n")

# Comparison 2: Severe COVID vs Mild/Asymptomatic COVID
cat("⏳ DEG: Severe vs Mild COVID...\n")
deg_sv_ml <- FindMarkers(
  seu,
  ident.1             = "Severe_COVID",
  ident.2             = "Mild_COVID",
  min.pct             = 0.25,
  logfc.threshold     = 0.5,
  max.cells.per.ident = 500,
  verbose             = FALSE
)
deg_sv_ml$gene      <- rownames(deg_sv_ml)
deg_sv_ml$direction <- ifelse(deg_sv_ml$avg_log2FC > 0, "Up in Severe", "Up in Mild")
write.csv(deg_sv_ml, file.path(RESULTS_DIR, "DEG_SevereCOVID_vs_MildCOVID.csv"),
          row.names = TRUE)
cat("✅ Severe vs Mild:", nrow(deg_sv_ml), "DEGs |",
    sum(deg_sv_ml$avg_log2FC > 0), "up in Severe |",
    sum(deg_sv_ml$avg_log2FC < 0), "up in Mild\n")

# Comparison 3: Severe COVID vs Influenza
# Key comparison — distinguishes COVID-specific immune response from viral illness
cat("⏳ DEG: Severe COVID vs Influenza...\n")
deg_sv_fl <- FindMarkers(
  seu,
  ident.1             = "Severe_COVID",
  ident.2             = "Influenza",
  min.pct             = 0.25,
  logfc.threshold     = 0.5,
  max.cells.per.ident = 500,
  verbose             = FALSE
)
deg_sv_fl$gene      <- rownames(deg_sv_fl)
deg_sv_fl$direction <- ifelse(deg_sv_fl$avg_log2FC > 0,
                              "Up in Severe COVID", "Up in Influenza")
write.csv(deg_sv_fl, file.path(RESULTS_DIR, "DEG_SevereCOVID_vs_Influenza.csv"),
          row.names = TRUE)
cat("✅ Severe COVID vs Influenza:", nrow(deg_sv_fl), "DEGs |",
    sum(deg_sv_fl$avg_log2FC > 0), "up in COVID |",
    sum(deg_sv_fl$avg_log2FC < 0), "up in Influenza\n")


# ── 14. Volcano Plots ─────────────────────────────────────────────────────────

plot_volcano <- function(deg_df, title, col_up, col_dn, filename) {
  p <- ggplot(deg_df, aes(x = avg_log2FC, y = -log10(p_val_adj),
                          colour = direction)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_colour_manual(values = setNames(c(col_up, col_dn),
                                          unique(deg_df$direction))) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.4) +
    geom_hline(yintercept = -log10(0.05),  linetype = "dashed", alpha = 0.4) +
    labs(title  = title,
         x      = "log2 fold change",
         y      = "-log10 adjusted p-value",
         colour = "") +
    theme_classic() +
    theme(legend.position = "bottom")
  ggsave(file.path(RESULTS_DIR, filename), p, width = 8, height = 6, dpi = 150)
}

plot_volcano(deg_sv_hl, "Severe COVID vs Healthy",
             "#E24B4A", "#1D9E75", "10_volcano_Severe_vs_Healthy.png")
plot_volcano(deg_sv_ml, "Severe COVID vs Mild COVID",
             "#E24B4A", "#EF9F27", "11_volcano_Severe_vs_Mild.png")
plot_volcano(deg_sv_fl, "Severe COVID vs Influenza",
             "#E24B4A", "#7F77DD", "12_volcano_Severe_vs_Influenza.png")


# ── 15. Reset Identity + Save ─────────────────────────────────────────────────

Idents(seu) <- "cell_type"
saveRDS(seu, file = file.path(DATA_DIR, "gse149689_final.rds"))

cat("\n🎉 PIPELINE COMPLETE!\n")
cat("📁 Results saved to:", RESULTS_DIR, "\n\n")
