**Single-cell RNA-seq analysis of bronchoalveolar lavage fluid (BALF) from COVID-19 patients — profiling the lung immune response at the actual site of infection across 4 disease groups using Seurat (GSE149689).**

[![R](https://img.shields.io/badge/R-4.3-276DC3?style=flat-square&logo=r&logoColor=white)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-5.x-2E86AB?style=flat-square)](https://satijalab.org/seurat/)
[![GEO](https://img.shields.io/badge/GEO-GSE149689-4CAF50?style=flat-square)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149689)
[![Paper](https://img.shields.io/badge/Paper-Nature%20Medicine%202020-red?style=flat-square)](https://doi.org/10.1038/s41591-020-0901-9)

 Built in R using Seurat. Data is from Liao et al. 2020 (*Nature Medicine*).

Show Image
Show Image
Show Image
Show Image

📌 Project Summary
Unlike most COVID scRNA-seq studies that use peripheral blood (PBMCs), this analysis targets the lung directly — sampling immune cells at the actual site of infection. Built in R using Seurat. Data from Liao et al. 2020 (Nature Medicine).

🧬 Background
I wanted to work on a COVID dataset that wasn't just peripheral blood — most scRNA-seq COVID papers use PBMCs, which tells you what's happening systemically but not at the actual site of infection. BALF samples the lung directly, so the immune cells you're looking at are the ones actually dealing with the virus.
This dataset has four groups (severe COVID, mild COVID, influenza, healthy), which means you can ask whether what you're seeing is COVID-specific or just a generic viral response. The severe vs influenza comparison turned out to be the most interesting one.
Dataset: GSE149689 — 20 donors, 85k raw cells.

📊 Dataset
ConditionDonorsCells (post-QC)Severe COVID-1968,922Mild / Asymptomatic COVID-19520,327Influenza510,171Healthy control417,401
Started with 85,144 cells, kept 56,821 after QC (filtered on nGenes 300–6000 and MT% < 15).

🔧 Pipeline
raw matrix (85,144 × 33,538)
  → QC filter → 56,821 cells
  → NormalizeData → FindVariableFeatures (2000 HVGs)
  → ScaleData → PCA (30 PCs, used top 15)
  → FindNeighbors → FindClusters (res 0.5) → RunUMAP
  → manual annotation (23 clusters)
  → FindAllMarkers + condition DEGs
Resolution 0.5 was chosen after trying 0.3 (too coarse, merged monocyte subtypes) and 0.8 (over-split T cells in a biologically unsupported way). PC cutoff at 15 from the elbow plot.

⚙️ Methodological Rationale
Each step was chosen for a specific reason — not just convention:
StepWhyQC (nGenes 300–6000, MT% < 15)Removes empty droplets, dying cells, and doublets. MT% kept at 15 (not stricter 10) to preserve cells in the inflammatory COVID contextNormalizationCorrects for sequencing depth differences so comparisons reflect biology, not library size2000 HVGsFocuses analysis on genes that vary between cell types; removes housekeeping gene noisePCA → top 15 PCsPC cutoff from elbow plot — beyond 15, PCs capture noise not signalClustering res 0.5Empirically tuned: 0.3 merged monocyte subtypes; 0.8 over-split T cellsUMAPVisualization only — not used for clustering decisionsFindAllMarkersWilcoxon rank-sum test per cluster to identify marker genes for annotation

📈 Results
Cell Types — 23 Clusters
Show Image
The monocyte compartment is highly heterogeneous — at least 6 distinct monocyte/DC populations (Clusters 1, 6, 7, 8, 11, 12, 17). This level of detail is only visible because BALF was used instead of blood.
ClusterCell TypeKey Markers0CD8 T cellCD8A, CD8B, GZMK1Inflamed MonocyteIL1B, PTGS2, OSM2Naive CD4 T cellMAL, LEF1, CD283NK cellMYOM2, SH2D1B4Naive B cellTCL1A, IGHD, FCER25NK/gdT cellKLRC2, TRDC, KLRF16Activated MonocyteCXCL2, CXCL3, EREG7ISG MonocyteSIGLEC1, IFI44L8CD16 MonocyteCDKN1C, C1QA9PlateletPTCRA, CLEC1B, ESAM10PlasmablastTNFRSF13B, SPIB, BLK11CD14 MonocyteS100A9, IFITM312Early MonocyteEGR1, FOS, FLT313CD4 Effector T cellLAG3, TRAC, LAIR214Cycling T cellTYMS, CDC45, PKMYT115NeutrophilHP, RETN, SERPINB216ISG cellIFI27, ANKRD22, CCL217Classical DC (cDC2)FCER1A, CD1C, CD1E18Erythrocyte*GYPA, GYPB, SPTA119Plasma cellIGHA1, IGHG1, JCHAIN20Memory B cellCD19, EBF1, FCRL221Regulatory T cellTNFRSF4, IL2RA22Pre-B cellPAX5, POU2AF1, BLNK
*Cluster 18 (erythrocytes) is likely contamination or stress erythropoiesis — kept for completeness but excluded from downstream analysis.

Condition-Level Distribution
Show Image
Severe COVID patients are heavily skewed toward monocyte clusters. Healthy donors are predominantly T cells. The shift is striking visually.
Show Image

Marker Expression
Show Image
Show Image

What Actually Stood Out
The ISG Monocytes (Cluster 7) — defined by SIGLEC1 and IFI44L — are expanded specifically in severe COVID, beyond levels seen in influenza. This suggests something SARS-CoV-2-specific in the interferon-driven innate response, not just a generic viral pattern.
Severe vs Influenza had the most DEGs (843), but more genes were higher in influenza (572) than in severe COVID (271). COVID severity is not a bigger immune response — it's a more dysregulated one concentrated in specific innate pathways.
The Severe vs Mild comparison showed the same pattern: more genes up in mild (313) than severe (248).
Show Image
755 DEGs — severe vs healthy
Show Image
561 DEGs — severe vs mild
Show Image
843 DEGs — severe vs influenza (most informative comparison)

🔬 Key Biological Findings

Monocyte expansion in severe COVID: Dramatic shift toward innate immune dominance — ISG Monocytes, Inflamed Monocytes, and Activated Monocytes are disproportionately expanded vs healthy and mild cases
ISG Monocyte signature is COVID-specific: Cluster 7 (SIGLEC1+, IFI44L+) exceeds influenza levels — pointing to a SARS-CoV-2-specific interferon dysregulation, not a generic viral response
Dysregulation, not volume: Influenza has more upregulated genes than severe COVID — severity is driven by dysregulation in targeted innate pathways, not a globally stronger response
T cell suppression at infection site: Severe COVID shows relative depletion of the adaptive immune compartment in BALF, consistent with impaired adaptive responses in severe disease
Innate-adaptive imbalance: Compositional shifts and DEG patterns together suggest a breakdown in innate-adaptive coordination in severe COVID-19


💡 Biological Significance
This analysis captures immune dynamics at the actual lung infection site — not systemic blood signals — giving a more direct view of COVID pathophysiology.

ISG monocyte expansion beyond influenza levels points to virus-specific interferon signaling dysregulation that may drive cytokine-mediated lung damage
Mild COVID having more upregulated genes than severe disease supports the idea that effective (not just strong) immune activation drives recovery
6 distinct monocyte/DC populations visible only in BALF highlights what's lost when only peripheral blood is analyzed — relevant for identifying therapeutically targetable cell populations
Results align with cytokine storm literature and provide a reproducible baseline for comparing therapeutic interventions at single-cell resolution


⚠️ Limitations

Sample size: 20 donors, 4–6 per group — interpret rare population results with caution
No batch correction: No explicit Harmony/scVI correction applied; donor-level effects may influence clustering
Annotation is marker-based: Cell types assigned from literature markers, not experimentally validated (no flow cytometry or CITE-seq confirmation)
Cross-sectional: Single time point per donor — no longitudinal tracking of disease progression
Cluster 18: Erythrocytes likely represent technical contamination, not a true immune population


🚀 Future Directions

Cell-type-specific DEGs: Re-run differential expression within individual clusters (e.g., ISG Monocytes only) to identify which populations drive condition-level differences
Batch correction: Apply Harmony or scVI to account for donor variation and validate key findings post-integration
Trajectory analysis: Use Monocle or RNA velocity to model monocyte state transitions — is ISG monocyte a terminal or transient activation state?
Pathway enrichment: GSEA / GO enrichment on DEG lists to identify dysregulated biological processes (complement, coagulation, interferon)
Multi-cohort integration: Extend to GSE145926 and other published COVID BALF datasets to test reproducibility and expand statistical power


🔁 Reproduce
Install
rinstall.packages(c("Seurat", "dplyr", "ggplot2", "patchwork", "Matrix"))
BiocManager::install("GEOquery")
Data
Download matrix.mtx, features.tsv, barcodes.tsv from GSE149689 and place in one folder.
Run
rDATA_DIR <- "path/to/GSE149689"   # set this at the top of the script
source("GSE149689_pipeline_fixed.R")
# results/ folder is created automatically
FindAllMarkers takes ~15 min. Everything else is fast. Tested on R 4.3, Seurat 5.x, Windows 11.

📁 Output Files
results/
├── 01_QC_prefilter.png
├── 02_QC_postfilter.png
├── 03_HVG.png
├── 04_Elbow.png
├── 05_feature_plots.png           ← canonical markers on UMAP
├── 06_dotplot.png
├── 07_UMAP_annotated.png          ← main result
├── 08_UMAP_condition.png
├── 09_composition_bar.png
├── 10_volcano_Severe_vs_Healthy.png
├── 11_volcano_Severe_vs_Mild.png
├── 12_volcano_Severe_vs_Influenza.png
├── markers_all.csv                ← 17,807 DE genes
├── markers_top5.csv
├── cell_composition.csv
├── DEG_SevereCOVID_vs_Healthy.csv
├── DEG_SevereCOVID_vs_MildCOVID.csv
└── DEG_SevereCOVID_vs_Influenza.csv

gse149689_final.rds                ← annotated Seurat object

📚 Reference
Liao M, Liu Y, Yuan J, et al. Single-cell landscape of bronchoalveolar immune cells in patients with COVID-19. Nature Medicine. 2020;26:842–844. https://doi.org/10.1038/s41591-020-0901-9
