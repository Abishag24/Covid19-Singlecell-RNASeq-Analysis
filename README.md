**Single-cell RNA-seq analysis of bronchoalveolar lavage fluid (BALF) from COVID-19 patients — profiling the lung immune response at the actual site of infection across 4 disease groups using Seurat (GSE149689).**

[![R](https://img.shields.io/badge/R-4.3-276DC3?style=flat-square&logo=r&logoColor=white)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-5.x-2E86AB?style=flat-square)](https://satijalab.org/seurat/)
[![GEO](https://img.shields.io/badge/GEO-GSE149689-4CAF50?style=flat-square)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689)
[![Paper](https://img.shields.io/badge/Paper-Nature%20Medicine%202020-red?style=flat-square)](https://doi.org/10.1038/s41591-020-0901-9)

 Built in R using Seurat. Data is from Liao et al. 2020 (*Nature Medicine*).

**Project Summary**
Unlike most COVID scRNA-seq studies that use peripheral blood, this analysis targets the lung directly to capture immune dynamics where infection occurs.
Single-cell RNA-seq analysis of bronchoalveolar lavage fluid from COVID-19 patients — profiling the lung immune response at the actual site of infection across 4 disease groups using Seurat (GSE149689).
Show Image
Show Image
Show Image
Show Image
Built in R using Seurat. Data is from Liao et al. 2020 (Nature Medicine).

**Background**
I wanted to work on a COVID dataset that wasn't just peripheral blood — most scRNA-seq COVID papers use PBMCs, which tells you what's happening systemically but not at the actual site of infection. BALF samples the lung directly, so the immune cells you're looking at are the ones actually dealing with the virus.
This dataset has four groups (severe COVID, mild COVID, influenza, healthy), which means you can ask whether what you're seeing is COVID-specific or just a generic viral response. The severe vs influenza comparison turned out to be the most interesting one.
Dataset: GSE149689, 20 donors, 85k raw cells.

**Dataset**
ConditionDonorsCells (post-QC)Severe COVID-1968,922Mild / Asymptomatic COVID-19520,327Influenza510,171Healthy control417,401
Started with 85,144 cells, kept 56,821 after QC (filtered on nGenes 300–6000 and MT% < 15).

**Pipeline**
Pretty standard Seurat workflow:
raw matrix (85,144 × 33,538)
  → QC filter → 56,821 cells
  → NormalizeData → FindVariableFeatures (2000 HVGs)
  → ScaleData → PCA (30 PCs, used top 15)
  → FindNeighbors → FindClusters (res 0.5) → RunUMAP
  → manual annotation (23 clusters)
  → FindAllMarkers + condition DEGs
I used resolution 0.5 after trying 0.3 (too coarse, merged monocyte subtypes that should be separate) and 0.8 (split a T cell cluster in a way that didn't make biological sense). PC cutoff at 15 from the elbow plot.

** Methodological Rationale**
Each step in the pipeline was chosen for a specific reason — not just convention:

**Quality Control (QC)**: Filtering on nGenes (300–6000) and MT% (<15) removes low-quality cells, empty droplets, and dying cells that would introduce noise into clustering. Keeping MT% threshold at 15 (rather than stricter 10) preserved more cells given the inflammatory context of COVID BALF.
Normalization (NormalizeData): Corrects for sequencing depth differences across cells so that gene expression comparisons reflect biology, not technical variation in library size.
Highly Variable Gene Selection (2000 HVGs): Focuses downstream analysis on genes that actually differ between cell types, reducing noise from housekeeping genes.
PCA → top 15 PCs: Dimensionality reduction before clustering. PC cutoff at 15 was chosen from the elbow plot — beyond that, PCs capture noise rather than meaningful variance.
Clustering (resolution 0.5): Resolution was tuned empirically — 0.3 merged distinct monocyte subtypes; 0.8 over-split T cells in a biologically unsupported way. 0.5 gave 23 clusters that aligned with known immune populations.
UMAP: Used for visualization only (not for clustering). Provides intuitive 2D representation of the high-dimensional transcriptional landscape.
Marker Gene Identification (FindAllMarkers): Wilcoxon rank-sum test to identify genes significantly enriched per cluster, enabling biologically-informed cell type annotation.


**Results**
Cell types
Ended up with 23 clusters. Annotated using marker genes from the literature + checking the top markers from FindAllMarkers output.
Show Image
The monocyte compartment is heterogeneous — there are at least 6 distinct monocyte/DC populations here (Clusters 1, 6, 7, 8, 11, 12, 17). That's one of the things that makes BALF interesting compared to blood.
Full marker table:
ClusterCell TypeKey Markers0CD8 T cellCD8A, CD8B, GZMK1Inflamed MonocyteIL1B, PTGS2, OSM2Naive CD4 T cellMAL, LEF1, CD283NK cellMYOM2, SH2D1B4Naive B cellTCL1A, IGHD, FCER25NK/gdT cellKLRC2, TRDC, KLRF16Activated MonocyteCXCL2, CXCL3, EREG7ISG MonocyteSIGLEC1, IFI44L8CD16 MonocyteCDKN1C, C1QA9PlateletPTCRA, CLEC1B, ESAM10PlasmablastTNFRSF13B, SPIB, BLK11CD14 MonocyteS100A9, IFITM312Early MonocyteEGR1, FOS, FLT313CD4 Effector T cellLAG3, TRAC, LAIR214Cycling T cellTYMS, CDC45, PKMYT115NeutrophilHP, RETN, SERPINB216ISG cellIFI27, ANKRD22, CCL217Classical DC (cDC2)FCER1A, CD1C, CD1E18ErythrocyteGYPA, GYPB, SPTA119Plasma cellIGHA1, IGHG1, JCHAIN20Memory B cellCD19, EBF1, FCRL221Regulatory T cellTNFRSF4, IL2RA22Pre-B cellPAX5, POU2AF1, BLNK
Cluster 18 (erythrocytes) is almost certainly contamination / stress erythropoiesis rather than a real immune population — kept it in for completeness but would exclude it from any downstream analysis.

Where cells come from
Show Image
Severe COVID patients are heavily skewed toward the monocyte clusters. Healthy donors are mostly T cells. That shift is pretty striking visually.
Show Image

Marker expression
Show Image
Show Image

What actually stood out
The ISG monocytes (Cluster 7). SIGLEC1 and IFI44L are the defining markers — this is an interferon-stimulated monocyte population that's expanded specifically in severe COVID. What made it interesting is that it's elevated beyond what you see in influenza, which also triggers a strong interferon response. So it's not just "interferon = bad viral infection" — there's something specific about SARS-CoV-2 here. I don't have a clean mechanistic explanation for why, but the pattern is clear in the data.
Severe vs Influenza was the most informative comparison (843 DEGs). More genes are higher in influenza (572) than in severe COVID (271). My read on this: COVID severity isn't about having a bigger immune response, it's about having a more dysregulated one concentrated in specific innate pathways.
The Severe vs Mild comparison had a similar pattern — more genes up in mild (313) than severe (248). I wasn't expecting that going in.
Show Image
755 DEGs — severe vs healthy
Show Image
561 DEGs — severe vs mild
Show Image
843 DEGs — severe vs influenza (most interesting one)

**Key Biological Findings**
Monocyte expansion in severe COVID: Severe patients show a dramatic shift toward innate immune dominance — monocyte clusters (particularly ISG Monocytes, Inflamed Monocytes, and Activated Monocytes) are disproportionately expanded compared to healthy controls and mild cases.
ISG Monocyte signature is COVID-specific: Cluster 7 (SIGLEC1+, IFI44L+) is elevated in severe COVID beyond what's seen even in influenza — suggesting a SARS-CoV-2-specific interferon-driven innate response rather than a generic viral one.
Dysregulation, not just activation: The severe vs influenza comparison (843 DEGs) shows more genes upregulated in influenza (572) than severe COVID (271), indicating severity in COVID is driven by immune dysregulation in specific innate pathways — not a globally stronger response.
T cell compartment is relatively suppressed: Healthy donors are predominantly T cell-rich; severe COVID patients show relative depletion of this adaptive immune compartment at the site of infection, consistent with known impaired adaptive responses in severe disease.
Innate-adaptive imbalance: Taken together, the compositional shifts and DEG patterns suggest a breakdown in the coordination between innate and adaptive immunity in severe COVID-19.


**Biological Significance**
This analysis captures immune dynamics at the actual site of SARS-CoV-2 infection — the lung — rather than systemic blood markers, giving a more direct window into disease pathophysiology.

The ISG monocyte expansion beyond influenza levels points to a virus-specific dysregulation of interferon signaling that may contribute to the cytokine-driven lung damage seen in severe COVID.
The finding that mild COVID has more upregulated genes than severe disease supports a model where effective immune activation (not just strong activation) drives recovery, and where severity reflects qualitative dysregulation.
The 6 distinct monocyte/DC populations in BALF highlight how cellular heterogeneity at the infection site is lost when only peripheral blood is analyzed — a finding relevant to understanding which cell populations are actually targetable therapeutically.
These patterns are consistent with published work on cytokine storm mechanisms and provide a reproducible baseline for comparing future therapeutic interventions at single-cell resolution.


**Limitations**

Sample size: 20 donors total, with 4–6 per group. Compositional shifts and DEG results should be interpreted with caution given limited statistical power for rare populations.
No batch correction applied: Samples were processed together in the pipeline but no explicit batch correction (e.g., Harmony) was used. Donor-level effects may influence clustering.
Cell type annotation is marker-based: Annotations rely on known marker genes from the literature and have not been experimentally validated (e.g., flow cytometry, CITE-seq protein co-expression).
Cross-sectional design: All samples are single time points. There is no longitudinal data to track how these populations shift during disease progression or recovery.
Cluster 18 (erythrocytes): Almost certainly technical contamination or stress erythropoiesis rather than a true immune population — not biologically meaningful for this analysis.


**Future Directions**
Cell-type-specific DEG analysis: Re-run differential expression within individual clusters (e.g., just ISG Monocytes) rather than across all cells, to resolve which cell populations are driving the condition-level differences.
Batch correction: Apply Harmony or scVI to explicitly account for donor-level variation and test whether key findings hold after integration.
Trajectory / pseudotime analysis: Use Monocle or RNA velocity to model monocyte state transitions — particularly whether the ISG monocyte state is a terminal activation state or a transient one.
Pathway enrichment: Run GSEA or GO enrichment on the condition DEG lists to identify which biological processes (complement, coagulation, interferon pathways) are most significantly altered.
Extend to published datasets: Integrate with other published COVID BALF datasets (e.g., GSE145926) to test reproducibility across cohorts and expand statistical power.


**Reproduce**
Install
rinstall.packages(c("Seurat", "dplyr", "ggplot2", "patchwork", "Matrix"))
BiocManager::install("GEOquery")
Data
Download matrix.mtx, features.tsv, barcodes.tsv from GSE149689 and put them in one folder.
Run
r# Set your data folder at the top of the script
DATA_DIR <- "path/to/GSE149689"

source("GSE149689_pipeline_fixed.R")
# results/ folder gets created automatically
FindAllMarkers takes a while (~15 min on my laptop). Everything else is fast.
Tested on R 4.3, Seurat 5.x, Windows 11. Should work on Mac/Linux, haven't tested.

**Output files**
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

**Reference**
Liao M, Liu Y, Yuan J, et al. Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19. Nature Medicine. 2020;26:842–844. https://doi.org/10.1038/s41591-020-0901-9
