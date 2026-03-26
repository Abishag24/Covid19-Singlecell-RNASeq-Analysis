****Single-cell RNA-seq analysis of bronchoalveolar lavage fluid from COVID-19 patients — profiling the lung immune response at the actual site of infection across 4 disease groups using Seurat (GSE149689).**

[![R](https://img.shields.io/badge/R-4.3-276DC3?style=flat-square&logo=r&logoColor=white)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-5.x-2E86AB?style=flat-square)](https://satijalab.org/seurat/)
[![GEO](https://img.shields.io/badge/GEO-GSE149689-4CAF50?style=flat-square)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689)
[![Paper](https://img.shields.io/badge/Paper-Nature%20Medicine%202020-red?style=flat-square)](https://doi.org/10.1038/s41591-020-0901-9)

 Built in R using Seurat. Data is from Liao et al. 2020 (*Nature Medicine*).

---

## Background

I wanted to work on a COVID dataset that wasn't just peripheral blood — most scRNA-seq COVID papers use PBMCs, which tells you what's happening systemically but not at the actual site of infection. BALF samples the lung directly, so the immune cells you're looking at are the ones actually dealing with the virus.

This dataset has four groups (severe COVID, mild COVID, influenza, healthy), which means you can ask whether what you're seeing is COVID-specific or just a generic viral response. The severe vs influenza comparison turned out to be the most interesting one.

Dataset: [GSE149689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689), 20 donors, 85k raw cells.

---

## Dataset

| Condition | Donors | Cells (post-QC) |
|---|---|---|
| Severe COVID-19 | 6 | 8,922 |
| Mild / Asymptomatic COVID-19 | 5 | 20,327 |
| Influenza | 5 | 10,171 |
| Healthy control | 4 | 17,401 |

Started with 85,144 cells, kept 56,821 after QC (filtered on nGenes 300–6000 and MT% < 15).

---

## Pipeline

Pretty standard Seurat workflow:

```
raw matrix (85,144 × 33,538)
  → QC filter → 56,821 cells
  → NormalizeData → FindVariableFeatures (2000 HVGs)
  → ScaleData → PCA (30 PCs, used top 15)
  → FindNeighbors → FindClusters (res 0.5) → RunUMAP
  → manual annotation (23 clusters)
  → FindAllMarkers + condition DEGs
```

I used resolution 0.5 after trying 0.3 (too coarse, merged monocyte subtypes that should be separate) and 0.8 (split a T cell cluster in a way that didn't make biological sense). PC cutoff at 15 from the elbow plot.

---

## Results

### Cell types

Ended up with 23 clusters. Annotated using marker genes from the literature + checking the top markers from `FindAllMarkers` output.

![UMAP annotated](results/07_UMAP_annotated.png)

The monocyte compartment is heterogeneous — there are at least 6 distinct monocyte/DC populations here (Clusters 1, 6, 7, 8, 11, 12, 17). That's one of the things that makes BALF interesting compared to blood.

Full marker table:

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
| 18 | Erythrocyte | GYPA, GYPB, SPTA1 |
| 19 | Plasma cell | IGHA1, IGHG1, JCHAIN |
| 20 | Memory B cell | CD19, EBF1, FCRL2 |
| 21 | Regulatory T cell | TNFRSF4, IL2RA |
| 22 | Pre-B cell | PAX5, POU2AF1, BLNK |

Cluster 18 (erythrocytes) is almost certainly contamination / stress erythropoiesis rather than a real immune population — kept it in for completeness but would exclude it from any downstream analysis.

---

### Where cells come from

![UMAP by condition](results/08_UMAP_condition.png)

Severe COVID patients are heavily skewed toward the monocyte clusters. Healthy donors are mostly T cells. That shift is pretty striking visually.

![Composition](results/09_composition_bar.png)

---

### Marker expression

![Feature plots](results/05_feature_plots.png)

![Dot plot](results/06_dotplot.png)

---

### What actually stood out

**The ISG monocytes (Cluster 7).** SIGLEC1 and IFI44L are the defining markers — this is an interferon-stimulated monocyte population that's expanded specifically in severe COVID. What made it interesting is that it's elevated *beyond* what you see in influenza, which also triggers a strong interferon response. So it's not just "interferon = bad viral infection" — there's something specific about SARS-CoV-2 here. I don't have a clean mechanistic explanation for why, but the pattern is clear in the data.

**Severe vs Influenza was the most informative comparison (843 DEGs).** More genes are higher in influenza (572) than in severe COVID (271). My read on this: COVID severity isn't about having a bigger immune response, it's about having a more dysregulated one concentrated in specific innate pathways.

The Severe vs Mild comparison had a similar pattern — more genes up in mild (313) than severe (248). I wasn't expecting that going in.

![Volcano severe vs healthy](results/10_volcano_Severe_vs_Healthy.png)
*755 DEGs — severe vs healthy*

![Volcano severe vs mild](results/11_volcano_Severe_vs_Mild.png)
*561 DEGs — severe vs mild*

![Volcano severe vs influenza](results/12_volcano_Severe_vs_Influenza.png)
*843 DEGs — severe vs influenza (most interesting one)*

---

## Reproduce

### Install

```r
install.packages(c("Seurat", "dplyr", "ggplot2", "patchwork", "Matrix"))
BiocManager::install("GEOquery")
```

### Data

Download `matrix.mtx`, `features.tsv`, `barcodes.tsv` from [GSE149689](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689) and put them in one folder.

### Run

```r
# Set your data folder at the top of the script
DATA_DIR <- "path/to/GSE149689"

source("GSE149689_pipeline_fixed.R")
# results/ folder gets created automatically
```

FindAllMarkers takes a while (~15 min on my laptop). Everything else is fast.

Tested on R 4.3, Seurat 5.x, Windows 11. Should work on Mac/Linux, haven't tested.

---

## Output files

```
results/
├── 01_QC_prefilter.png
├── 02_QC_postfilter.png
├── 03_HVG.png
├── 04_Elbow.png
├── 05_feature_plots.png          ← canonical markers on UMAP
├── 06_dotplot.png
├── 07_UMAP_annotated.png         ← main result
├── 08_UMAP_condition.png
├── 09_composition_bar.png
├── 10_volcano_Severe_vs_Healthy.png
├── 11_volcano_Severe_vs_Mild.png
├── 12_volcano_Severe_vs_Influenza.png
├── markers_all.csv               ← 17,807 DE genes
├── markers_top5.csv
├── cell_composition.csv
├── DEG_SevereCOVID_vs_Healthy.csv
├── DEG_SevereCOVID_vs_MildCOVID.csv
└── DEG_SevereCOVID_vs_Influenza.csv

gse149689_final.rds               ← annotated Seurat object
```

---

## Reference

Liao M, Liu Y, Yuan J, et al. Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19. *Nature Medicine*. 2020;26:842–844. https://doi.org/10.1038/s41591-020-0901-9
