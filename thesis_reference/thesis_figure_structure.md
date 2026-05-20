# Chromatin Accessibility and Transcription Factor Regulation in Age-Related Diseases
**Student**: Yanchen Zheng (48692452)
**Supervisor**: Dr Ralph Patrick
**Course**: BIOX7026

---

## Abstract
200–300 words

---

## Introduction

### 1. Age-Related Diseases and Aging
- Age-related disease burden across organ systems
- Kidney disease (CKD), Lung disease (IPF), Aorta, T cells as study models

### 2. Epigenetic and Transcription Factor Control
- Chromatin remodelling in aging
- SIPHON model (Patrick et al., 2024)
- AP-1 as key driver of aging chromatin changes

### 3. ATAC-seq as Key Methodology

**Fig. 1. Overview of scATAC-seq pipeline.**
(A) Transposition with Tn5 transposase. (B) Amplification with barcodes. (C) Sequencing. (D) Alignment. (E) Peak calling. (F) Downstream analysis.

### 4. Analysis Approaches
- Quality control and peak calling
- Clustering and cell type annotation
- DAR calling (pseudo-bulk DESeq2)
- Motif enrichment analysis (HOMER)

### 5. Project Rationale and Hypothesis

**Fig. 2. Project overview.**
Schematic of three aims: (1) Data integration of disease scATAC-seq datasets, (2) TF motif enrichment analysis, (3) Mapping disease vs aging chromatin regions.

---

## Materials and Methods

### Data Details

**Table 1. Datasets included in this project.**

| Tissue | Disease | Mouse Model | Contrast | Reference |
|--------|---------|-------------|----------|-----------|
| Kidney | AKI → CKD transition | Ischemia-reperfusion injury (IRI) | Day42 vs Sham | Muto et al., *Sci. Adv.* 2024 |
| Lung | Chronic obstructive pulmonary disease (COPD) | Cigarette smoke + LPS | Case vs Control | Zhang et al., *PLOS One* 2025 |
| Aorta | Aortic aneurysm and dissection (AAD) | HFD + AngII infusion | Challenge vs Control | Zhang et al., *ATVB* 2023 |
| T cells | Immune aging | In-house (Patrick Lab) | Young chronic vs Young control | Unpublished |

### Data Processing
- Download, Cell Ranger ATAC alignment (mm10)
- Peak calling (MACS2)
- Doublet removal
- Quality filtering

### Universal Peak Set
- Consensus peak merging across samples (`merge_consensus_tissues.py`)
- Universal peak set construction across 4 tissues

### DAR Calling
- Pseudo-bulk aggregation per cell type per sample
- DESeq2 (`padj < 0.05`, `|log2FC| > 0.5`)

### Motif Enrichment Analysis
- HOMER `findMotifsGenome.pl` (mm10)
- Background: stable peaks (cell-type-specific) or non-significant (NS) peaks
- Min peaks: 20
- Comparisons: opening vs stable/NS, closing vs stable/NS (and reciprocals)

### Aging DAR Comparison (Aim 3)
- Reference aging dataset: Science paper (Patrick et al., 2024)
- `bedtools intersect` for peak overlap
- Hypergeometric test for enrichment
- Pearson `r` for directional concordance

---

## Results

---

### Aim 1: Data Integration

**Fig. 3. Dataset composition across four disease tissues.**
(A) Number of mice per experimental group (stacked bar chart): Kidney (Sham n=3, Day14 n=3, Day42 n=3); Lung (Control n=2, Case n=4); Aorta (Control n=3, Challenge n=3); T cells (Young control n=3, Juvenile n=2, Young acute n=3, Young chronic n=3, Aged n=4). (B) Number of cells per experimental group per tissue (stacked bar chart): Kidney total ~73,974; Lung ~8,179; Aorta ~9,552; T cells ~11,958. Together these panels establish the experimental design and scale of each dataset.

> *Script*: `scATAC-seq/` QC outputs

**Fig. 4. scATAC-seq quality control metrics across tissues.**
Violin plots showing per-cell distributions across Kidney, Lung, Aorta, and T cells. (A) Fragments per cell (log10 scale). (B) TSS enrichment score. (C) Nucleosome signal. All four tissues pass standard QC thresholds, confirming high-quality chromatin accessibility data suitable for downstream analysis.

> *Script*: `scATAC-seq/` QC outputs

**Fig. 5. Universal peak set composition across four disease tissues.**
Consensus peaks per tissue and universal peak set (bar chart): Kidney 581,481; Lung 245,604; Aorta 230,473; T cells 145,370; merged universal set 667,473. Tissue-specific and universal bars distinguished by colour with dashed separator. The universal peak set captures the full accessible chromatin landscape across all four disease tissues.

> *Script*: `Fig4_universal_peaks.R`

**Fig. 6. Tissue peak overlap in the universal peak set.**
UpSet plot showing intersections of peaks across the four tissues: largest intersection is Kidney-only peaks (292,438); shared peaks across all four tissues = 1,302; pairwise and three-way overlaps range from 7,064 to 73,520. The majority of peaks are tissue-specific, reflecting the divergent chromatin landscapes of distinct disease contexts and justifying a universal peak set approach for cross-tissue comparison.

> *Script*: `Fig4_universal_peaks.R`

**Fig. 7. Cell type clustering and annotation across four tissues.**
(A) UMAP of Kidney dataset: DCT_CNT, DTL_ATL, EC, IC, Injured_PT, LEUK, PC_URO, Pen, PODO_PEC, PST, TAL — including injured and transitional tubular cell states. (B) UMAP of Lung dataset: AT2, B, Ciliated, EC-vasc, Eosinophils, Fib, Mac, Mac-alv, Mac-inter, Mesothelial, Mo-Ly6c+, NK, Pen, SMCs, T. (C) UMAP of Aorta dataset: Endothelial, Fibroblast, Mac, Pericyte, SMC, T-cell. (D) UMAP of T cell dataset: Activated T cell, CD8+ effector, CD8+ T-cell, Cycling T cell, Effector T cell, Innate-like T cell, Memory T cell, Naïve T cell, NK cell, Treg — annotated via gene marker–based subclustering (FindAllMarkers); B cell contaminant cluster (Ebf1+) excluded from downstream analysis.

> *Script*: `Fig5_UMAP_annotation.R`, `Fig5_Tcells_annotation.R`

**Fig. 8. CoveragePlot at marker gene loci validates cell type annotation.**
Coverage plots showing fragment accessibility signal at tissue-specific marker genes. (A) *Epcam* locus (Kidney) — signal enriched in epithelial cell types (PT, TAL, DCT). (B) *Sftpc* locus (Lung) — signal enriched in AT2 cells. (C) *Acta2* locus (Aorta) — signal enriched in SMC. (D) *Cd3e* locus (T cells) — signal enriched in T cell subtypes. The concordance between chromatin accessibility and known marker gene expression confirms the accuracy of the chromatin-based cell type annotation.

> *Script*: `Fig_CoveragePlot.R`

---

### Aim 2: Transcription Factor Motif Enrichment

**Fig. 9. Differentially accessible region (DAR) quantification across all four disease tissues.**
Multi-panel diverging bar chart overview of opening (red) and closing (blue) DARs identified by pseudo-bulk DESeq2 analysis (padj < 0.05).
(A) Kidney (Day42 vs Sham) — PT and TAL show the largest number of DARs.
(B) Lung (Case vs Control) — AT2, Fib, and Mac-alv show the strongest response.
(C) Aorta (Challenge vs Control) — SMC is the predominantly affected cell type.
(D) T cells (Young_chronic vs Young_control) — chromatin remodelling distributed across T cell subtypes.

> *Script*: `DAR/Fig9_DAR_counts.R`

**Fig. 10. Four-tissue integrated HOMER motif enrichment heatmap.**
Side-by-side ComplexHeatmap panels. (A) Left — NS background: HOMER enrichment comparing opening/closing DARs against non-stable (NS) peaks as background. (B) Right — Stable background: same comparison using cell-type-specific stable peaks. Rows = TF motifs in fixed order; columns = cell types stratified by tissue (Lung, Aorta, Kidney, T cells — tissue identity shown as top colour bar) and by 4 comparison types (Opening vs background, Closing vs background, Background vs Opening, Background vs Closing). Colour = signed −log10(P-value): red = enriched in opening; blue = enriched in closing. AP-1 row group (top) and CTCF row group annotated on the right. Concordance between NS and Stable panels demonstrates that AP-1/bZIP enrichment in disease-associated opening DARs is robust and not background-dependent.

> *Script*: `CMpaper_heatmap/4tissues_integrate_heatmap.R`, `CMpaper_heatmap/assemble_Fig7.R`

**Fig. 11. Cell-type identity TF enrichment establishes chromatin identity baseline.**
Six HOMER motif enrichment barplots (2×3 layout). Each panel shows −log10(P-value) for the top cell-type-specific TF motifs (blue bars; red dashed line = significance threshold). Top row: (A) T cell — TCF3/LEF1/TCF7 lymphoid TFs; (B) Kidney PT — HNF1b, PPARα, HNF4a tubular epithelial TFs; (C) Kidney TAL — HNF family and KLF factors. Bottom row: (D) Lung AT2 — Smad2, HIF-1α alveolar TFs; (E) Lung Fib — KLF/ETS fibroblast TFs; (F) Aorta SMC — WT1, AP-2α, KLF4 vascular smooth muscle TFs. These cell-type identity TFs serve as the reference baseline for interpreting closing DAR motif enrichment in disease (Fig. 13–16).

> *Script*: `CMpaper_heatmap/assemble_Fig6.R`
> *Input*: `HOMER_Plots/DAR_closing_vs_opening/Cell_type_identity.pdf`

**Fig. 12. Shared TF motifs enriched in opening DARs across disease tissues.**
Cross-tissue dotplot of TF motifs present in ≥3 cell types/tissues (opening DARs only). AP-1 family (JunB, Fra2, Fos, Atf3, BATF), C/EBP, and ETS motifs consistently enriched in opening DARs across Kidney, Lung, Aorta, and T cells — indicating a common inflammatory chromatin opening programme across disease contexts. Bubble size = % target sequences; colour = −log10(P-value).

> *Script*: `HOMER/HOMER_Plot/Fig10_11_16_dotplots.R`
> *Output*: `TF_motif_plots/Fig10_Opening_shared_dotplot.png`

**Fig. 13. Kidney TF motif enrichment in disease DARs.**
(A) Opening DARs — amalgamated across all kidney cell types (Day42 vs Sham): AP-1 family (AP-1, JunB, Fos, Fra1, Fra2, Atf3, BATF) uniformly enriched.
(B) Closing DARs — per-cell-type top 25 motifs (7 cell types): PT and TAL show the strongest HNF1b, HNF4a, PPARα, and COUP-TFII loss, indicating progressive loss of tubular epithelial identity. DCT shows a similar but earlier response.
(C) PT cell type focus — Opening (AP-1 gain) vs Closing (HNF family/KLF loss) side-by-side, directly illustrating chromatin remodelling from an epithelial to an inflammatory state in proximal tubule cells.

> *Script*: `HOMER/HOMER_Plot/homer_kidney_plot.R`

**Fig. 14. Lung TF motif enrichment in disease DARs.**
(A) Opening DARs — amalgamated across all lung cell types (Case vs Control): AP-1/bZIP motifs broadly enriched; ETS factors additionally enriched in endothelial and macrophage populations.
(B) AT2 opening and closing DARs: Opening enriches for stress-response TFs (Stat3, Smad2/HIF-1α); closing shows loss of alveolar epithelial identity TFs.
(C) Fibroblast (Fib) opening DARs: KLF/ETS motifs enriched, supporting a pro-fibrotic TF network in COPD fibroblasts.
(D) Macrophage (Mac-alv) opening DARs: AP-1/C/EBP enrichment, consistent with inflammatory macrophage activation.

> *Script*: `HOMER/HOMER_Plot/homer_lung_plot.R`

**Fig. 15. Aorta TF motif enrichment in disease DARs.**
(A) SMC opening DARs (Challenge vs Control): WT1, AP-2α (TFAP2A), KLF4, and AP-1/bZIP motifs enriched — indicating vascular SMC activation and phenotypic switching.
(B) SMC closing DARs: Loss of canonical SMC contractile identity TFs; COUP-TFII/EAR2 nuclear receptor motifs enriched in closing DARs.
(C) Macrophage opening DARs: AP-1/C/EBP enrichment consistent with atherosclerotic plaque-associated inflammation.

> *Script*: `HOMER/HOMER_Plot/homer_aorta_plot.R`

**Fig. 16. T cell TF motif enrichment in disease DARs.**
(A) Opening DARs (Young_chronic vs Young_control): AP-1 family (JunB, Fos, Atf3) and ETS family (ETS1, ETV1, ERG) enriched in opening DARs, consistent with chronic T cell activation and altered effector programming.
(B) Closing DARs: Lymphoid identity TFs (TCF3/E-protein, LEF1, IRF4) enriched in closing DARs, indicating loss of naïve/memory T cell chromatin identity during chronic stimulation — paralleling the loss of epithelial identity seen in Kidney and Lung.

> *Script*: `HOMER/HOMER_Plot/homer_tcells_plot.R`

**Fig. 17. Tissue-specific TF motif enrichment in opening DARs.**
Dotplot of TF motifs unique to a single tissue (opening DARs). ETS/CTCF motifs specific to Aorta SMC (WT1, ERG); homeobox motifs specific to T cells (Hoxc13, HOXB13); STAT/TEAD motifs in Kidney IC and Lung. These tissue-specific opening signatures reveal disease-context-specific TF programmes beyond the shared AP-1 signature.

> *Script*: `HOMER/HOMER_Plot/Fig10_11_16_dotplots.R`
> *Output*: `TF_motif_plots/Fig16_Opening_specific_dotplot.png`


---

### Aim 3: Mapping Disease vs Aging Chromatin Regions

**Fig. 18. Pseudobulk DESeq2 aging DARs from Science paper data.**
(A) Number of aging-associated opening DARs (Aged vs Young, `padj < 0.05`) per cell type in Kidney. (B) Same for Lung.

> *Script*: `sciencepaper_DAR_comparison/02_pesudobulk_DESeq2.sh`

**Fig. 19. Enrichment of disease DARs in aging DARs — Kidney & Lung.**
(A) Enrichment fold (observed/expected) bar chart for opening (red) and closing (blue) disease DARs in aging regions. Significance stars from hypergeometric test. Cell types: Kidney DCT, PT, TAL; Lung AT2, Fib, Mac-alv.
(B) Pearson `r` correlation between disease and aging `log2FC` among aging-significant peaks (`padj < 0.05`). `n` = number of overlapping peaks per cell type.

> *Script*: `sciencepaper_DAR_comparison/03_overlap_and_plots.sh` → `Aim3_DAR_overlap_summary.png`

**Fig. 20. Scatter plots of disease vs aging log2FC.**
8-panel (3 Kidney top row + 3 Lung bottom row) scatter plots. X-axis: disease `log2FC` (`padj < 0.05`, `|log2FC| > 0.5`). Y-axis: aging `log2FC` (all peaks). Red = aging significant (`padj < 0.05`), orange = trending (`padj < 0.2`), grey = ns. Regression line through aging-significant points.

> *Script*: `sciencepaper_DAR_comparison/03_overlap_and_plots.sh` → `DAR_scatter_combined.png`

**Fig. 21. Classification of disease DARs relative to aging.**
(A) Stacked bar chart showing absolute numbers of disease-specific, shared, and aging-specific DARs per cell type (opening and closing).
(B) Proportion of disease DARs shared with aging DARs (%).
(C) Number of shared peaks per cell type.

> *Script*: `sciencepaper_DAR_comparison/04_region_classification.sh`

---

## Discussion

### Universal Peak Set and Data Integration
- Quality of the integrated dataset
- Challenges of cross-tissue peak merging

### TF Network Remodelling in Age-Related Diseases
- AP-1 enrichment across all four tissues (consistent with SIPHON model)
- Shared opening signature vs tissue-specific closing signature
- Loss of cell type identity TFs in closing DARs (HNF in kidney, Smad2/NKX in lung, KLF4 in aorta, TCF7 in T cells)
- Stable vs non-stable DAR classification confirms robustness of results

### Disease vs Aging Chromatin Overlap (Aim 3)
- Kidney shows stronger concordance with aging (`r = 0.45–0.62` in PT, TAL)
- Lung shows weaker but detectable overlap (`r = 0.23–0.31` in AT2, Mac-alv)
- Shared regions as potential aging-disease regulatory hubs
- Cell-type specificity of the overlap

### Limitations
- Mouse model may not fully recapitulate human aging
- Pseudo-bulk DESeq2 requires sufficient samples per group
- Lung Fib not tissue-specific — interpretation requires caution

---

## Supplementary Figures

**Supp. Fig. 1.** Full DAR tables for all cell types and tissues (opening and closing counts per contrast).
**Supp. Fig. 2.** Per-cell-type HOMER bubble plots for all remaining Kidney cell types (Endothelial, IC, Macrophages, PC).
**Supp. Fig. 3.** Per-cell-type HOMER bubble plots for remaining Lung cell types (B, NK, Ciliated, EC-vasc, Mo-Ly6c+, SMCs, T).
**Supp. Fig. 4.** UpSet plots of peak overlap across cell types within each tissue.
**Supp. Fig. 5.** Full scatter plots including EC-vasc and T cells (Lung) for Aim 3 disease vs aging comparison.
**Supp. Fig. 6.** -log10(p-value) bar charts for disease DAR enrichment in aging regions (all cell types).
**Supp. Fig. 7.** Comprehensive cross-tissue TF motif enrichment dotplot (AllTissues). Full dotplot showing top enriched motifs for all cell types across tissues for opening DARs. Supports Fig. 10 and Fig. 12–17 with a complete motif landscape view.

> *Script*: `HOMER/HOMER_Plot/Fig10_11_16_dotplots.R` → `Fig11_AllTissues_Opening_dotplot.png`

---

## Figure Checklist

| Figure | Description | Script | Status |
|--------|-------------|--------|--------|
| Fig. 3 | Dataset composition: mice + cells per group | `scATAC-seq/` QC outputs | ☐ |
| Fig. 4 | QC violin plots: fragments, TSS enrichment, nucleosome signal | `scATAC-seq/` QC outputs | ☐ |
| Fig. 5 | Universal peak set bar chart (counts per tissue) | `Fig4_universal_peaks.R` | ☐ |
| Fig. 6 | UpSet plot — tissue peak overlap in universal set | `Fig4_universal_peaks.R` | ☐ |
| Fig. 7 | UMAP clustering (4 tissues, detailed cell types) | `Fig5_UMAP_annotation.R`, `Fig5_Tcells_annotation.R` | ☐ |
| Fig. 8 | CoveragePlot at marker genes — validates annotation (4 tissues) | `Fig_CoveragePlot.R` | ✅ |
| Fig. 9 | DAR counts — all 4 tissues (diverging bar chart) | `DAR/Fig9_DAR_counts.R` | ✅ |
| Fig. 10 | 4-tissue integrated HOMER heatmap (NS + Stable) | `CMpaper_heatmap/4tissues_integrate_heatmap.R` | ☐ |
| Fig. 11 | Cell identity TF barplots (6 panels, 2×3) — baseline reference | `CMpaper_heatmap/assemble_Fig6.R` | ☐ |
| Fig. 12 | Shared TF motifs — Opening dotplot (≥3 units) | `HOMER/HOMER_Plot/Fig10_11_16_dotplots.R` | ✅ |
| Fig. 13 | Kidney HOMER — amalgamated + PT/TAL closing detail | `HOMER/HOMER_Plot/homer_kidney_plot.R` | ☐ |
| Fig. 14 | Lung HOMER — amalgamated + AT2/Fib/Mac detail | `HOMER/HOMER_Plot/homer_lung_plot.R` | ☐ |
| Fig. 15 | Aorta HOMER — SMC + Macrophage | `HOMER/HOMER_Plot/homer_aorta_plot.R` | ☐ |
| Fig. 16 | T cell HOMER — Opening + Closing | `HOMER/HOMER_Plot/homer_tcells_plot.R` | ☐ |
| Fig. 17 | Tissue-specific TF motifs — Opening dotplot | `HOMER/HOMER_Plot/Fig10_11_16_dotplots.R` | ✅ |
| Fig. 18 | Aging DAR counts | `sciencepaper_DAR_comparison/02` | ☐ |
| Fig. 19 | Aim3 summary (enrichment fold + Pearson r) | `sciencepaper_DAR_comparison/03` | ☐ |
| Fig. 20 | Disease vs aging scatter plots | `sciencepaper_DAR_comparison/03` | ☐ |
| Fig. 21 | Region classification stacked bars | `sciencepaper_DAR_comparison/04` | ☐ |
| Supp. Fig. 7 | AllTissues comprehensive dotplot | `HOMER/HOMER_Plot/Fig10_11_16_dotplots.R` | ✅ |
