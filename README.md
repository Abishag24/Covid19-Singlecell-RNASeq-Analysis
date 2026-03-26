# scRNA-seq Analysis of COVID-19 Bronchoalveolar Lavage Fluid (BALF)

Single-cell profiling of the lung immune response across 4 disease groups using Seurat (GSE149689).

[![R](https://img.shields.io/badge/R-4.3-276DC3?style=flat-square&logo=r&logoColor=white)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-5.x-2E86AB?style=flat-square)](https://satijalab.org/seurat/)
[![GEO](https://img.shields.io/badge/GEO-GSE149689-4CAF50?style=flat-square)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689)
[![Paper](https://img.shields.io/badge/Paper-Nature%20Medicine%202020-red?style=flat-square)](https://doi.org/10.1038/s41591-020-0901-9)

---

## Background

I wanted to work on a COVID dataset that wasn't just peripheral blood. Most scRNA-seq COVID papers use PBMCs, which tells you what's happening systemically — but not at the actual site of infection. BALF samples the lung directly, so the immune cells you're looking at are the ones actually dealing with the virus.

This dataset has four groups (severe COVID, mild COVID, influenza, healthy), which means you can ask whether what you're seeing is COVID-specific or just a generic viral response. The severe vs. influenza comparison turned out to be the most interesting one.

Dataset: GSE149689 — 20 donors, 85,144 raw cells.

---

## Dataset

Started with 85,144 cells, kept 56,821 after QC (filtered on nGenes 300–6000 and MT% < 15).

| Condition | Donors | Cells (post-QC) |
|---|---|---|
| Severe COVID-19 | 6 | 68,922 |
| Mild / Asymptomatic COVID-19 | 5 | 20,327 |
| Influenza | 5 | 10,171 |
| Healthy Control | 4 | 17,401 |

---

## Pipeline

| Step | Detail |
|---|---|
| Raw matrix | 85,144 × 33,538 genes |
| QC filter | nGenes 300–6000, MT% < 15 |
| Post-QC cells | 56,821 retained |
| NormalizeData | Corrects for sequencing depth |
| FindVariableFeatures | 2,000 HVGs selected |
| ScaleData | Centers and scales expression |
| PCA | 30 PCs computed; top 15 used |
| FindNeighbors + FindClusters | Resolution 0.5 |
| RunUMAP | Visualization only (not used for clustering) |
| Manual annotation | 23 clusters identified |
| FindAllMarkers | Wilcoxon rank-sum per cluster |
| Condition DEGs | Severe vs Healthy / Mild / Influenza |

Resolution 0.5 was chosen after trying 0.3 (too coarse, merged monocyte subtypes) and 0.8 (over-split T cells in a biologically unsupported way). PC cutoff at 15 from the elbow plot.

---

## Methodological Rationale

| Step | Why |
|---|---|
| QC: nGenes 300–6000, MT% < 15 | Removes empty droplets, dying cells, doublets. MT% kept at 15 (not 10) to preserve cells in the inflammatory COVID context |
| Normalization | Corrects sequencing depth so comparisons reflect biology, not library size |
| 2,000 HVGs | Focuses analysis on genes that vary between cell types; removes housekeeping noise |
| PCA → top 15 PCs | Cutoff from elbow plot — beyond 15, PCs capture noise not signal |
| Clustering res 0.5 | Empirically tuned: 0.3 merged monocyte subtypes; 0.8 over-split T cells |
| UMAP | Visualization only — not used for clustering decisions |
| FindAllMarkers | Wilcoxon rank-sum test per cluster to identify marker genes for annotation |

---

## Results

### Cell Types — 23 Clusters

The monocyte compartment is way more heterogeneous than I expected — at least 6 distinct monocyte/DC populations (Clusters 1, 6, 7, 8, 11, 12, 17). This level of detail only shows up because BALF was used instead of blood.

| Cluster | Cell Type | Key Markers |
|---|---|---|
| 0 | CD8 T cell | CD8A, CD8B, GZMK |
| 1 | Inflamed Monocyte | IL1B, PTGS2, OSM |
| 2 | Naive CD4 T cell | MAL, LEF1, CD28 |
| 3 | NK cell | MYOM2, SH2D1B |
| 4 | Naive B cell | TCL1A, IGHD, FCER2 |
| 5 | NK/gdT cell | KLRC2, TRDC, KLRF1 |
| 6 | Activated Monocyte | CXCL2, CXCL3, EREG |
| 7 | ISG Monocyte | SIGLEC1, IFI44L |
| 8 | CD16 Monocyte | CDKN1C, C1QA |
| 9 | Platelet | PTCRA, CLEC1B, ESAM |
| 10 | Plasmablast | TNFRSF13B, SPIB, BLK |
| 11 | CD14 Monocyte | S100A9, IFITM3 |
| 12 | Early Monocyte | EGR1, FOS, FLT3 |
| 13 | CD4 Effector T cell | LAG3, TRAC, LAIR2 |
| 14 | Cycling T cell | TYMS, CDC45, PKMYT1 |
| 15 | Neutrophil | HP, RETN, SERPINB2 |
| 16 | ISG cell | IFI27, ANKRD22, CCL2 |
| 17 | Classical DC (cDC2) | FCER1A, CD1C, CD1E |
| 18 | Erythrocyte* | GYPA, GYPB, SPTA1 |
| 19 | Plasma cell | IGHA1, IGHG1, JCHAIN |
| 20 | Memory B cell | CD19, EBF1, FCRL2 |
| 21 | Regulatory T cell | TNFRSF4, IL2RA |
| 22 | Pre-B cell | PAX5, POU2AF1, BLNK |

> \* Cluster 18 (erythrocytes) is likely contamination or stress erythropoiesis — kept for completeness but excluded from downstream analysis.

### Condition-Level Distribution

Severe COVID patients are heavily skewed toward monocyte clusters. Healthy donors are predominantly T cells. The shift is pretty striking visually — see `08_UMAP_condition.png` and `09_composition_bar.png`.

### Differential Expression

| Comparison | Total DEGs | Notes |
|---|---|---|
| Severe vs Healthy | 755 | Baseline immune landscape shift |
| Severe vs Mild | 561 | More genes up in Mild (313 vs 248) |
| Severe vs Influenza | 843 | More genes up in Influenza (572 vs 271) — most informative comparison |

---

## What Actually Stood Out

**The ISG Monocytes (Cluster 7)** — defined by SIGLEC1 and IFI44L — are expanded specifically in severe COVID, beyond levels seen in influenza. That points to something SARS-CoV-2-specific in the interferon-driven innate response, not just a generic viral pattern.

The severe vs. influenza comparison had the most DEGs (843), but more genes were actually higher in influenza (572) than in severe COVID (271). COVID severity isn't a bigger immune response — it's a more dysregulated one concentrated in specific innate pathways. The severe vs. mild comparison showed the same pattern.

---

## Limitations

- Sample size is small — 20 donors, 4–6 per group. Rare population results should be interpreted carefully.
- No batch correction applied (no Harmony/scVI). Donor-level effects may be influencing some of the clustering.
- Cell type annotation is marker-based only, not validated experimentally (no flow cytometry or CITE-seq to confirm).
- Cross-sectional design — single time point per donor, so no way to track progression.
- Cluster 18 (erythrocytes) is almost certainly a technical artifact.

---

## What I'd Do Next

- Re-run DEGs within individual clusters (e.g., ISG Monocytes only) to pin down which populations are actually driving the condition-level differences.
- Apply Harmony or scVI for batch correction and see if the key findings hold.
- Trajectory analysis with Monocle or RNA velocity to figure out if ISG monocytes are a transient activation state or terminal.
- GSEA / GO enrichment on the DEG lists — complement, coagulation, and interferon pathways are the obvious candidates.
- Integrate with GSE145926 and other published COVID BALF datasets to improve statistical power.

---

## Reproduce

### Install

```r
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork", "Matrix"))
BiocManager::install("GEOquery")
```

### Data

Download `matrix.mtx`, `features.tsv`, and `barcodes.tsv` from [GSE149689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689) and place them in one folder.

### Run

```r
DATA_DIR <- "path/to/GSE149689"  # set this at the top of the script
source("GSE149689_pipeline_fixed.R")
```

`results/` is created automatically. `FindAllMarkers` takes around 15 min; everything else is fast. Tested on R 4.3, Seurat 5.x, Windows 11.

---

## Output Files

```
results/
├── 01_QC_prefilter.png
├── 02_QC_postfilter.png
├── 03_HVG.png
├── 04_Elbow.png
├── 05_feature_plots.png        ← canonical markers on UMAP
├── 06_dotplot.png
├── 07_UMAP_annotated.png       ← main result
├── 08_UMAP_condition.png
├── 09_composition_bar.png
├── 10_volcano_Severe_vs_Healthy.png
├── 11_volcano_Severe_vs_Mild.png
├── 12_volcano_Severe_vs_Influenza.png
├── markers_all.csv             ← 17,807 DE genes
├── markers_top5.csv
├── cell_composition.csv
├── DEG_SevereCOVID_vs_Healthy.csv
├── DEG_SevereCOVID_vs_MildCOVID.csv
└── DEG_SevereCOVID_vs_Influenza.csv

gse149689_final.rds             ← annotated Seurat object
```

---

## Reference

Liao M, Liu Y, Yuan J, et al. Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19. *Nature Medicine.* 2020;26:842–844. https://doi.org/10.1038/s41591-020-0901-9
