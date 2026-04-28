# 🧬 scATAC-seq Analysis Pipeline — Human PBMCs (10k cells)

**Dataset:** 10x Genomics PBMC 10k ATACv2 (Chromium Controller)  
**Genome:** hg38 (GRCh38 — Ensembl v98)  
**Tools:** Signac v1.17 · Seurat v5 · R 4.6  
**Author:** Zain Ul Abiddin

> 📺 **Acknowledgement:** A portion of the code structure and workflow approach in
> this pipeline was adapted from tutorials by the
> **[Bioinformagician YouTube Channel](https://www.youtube.com/@Bioinformagician)**.
> Their clear and practical teaching style helped shape how this pipeline was
> developed and documented. Highly recommended for anyone learning bioinformatics in R.

---

## 📋 Table of Contents

1. [What is scATAC-seq? — A Beginner's Guide](#1-what-is-scatac-seq--a-beginners-guide)
2. [What Does This Pipeline Produce?](#2-what-does-this-pipeline-produce)
3. [Requirements — What to Install](#3-requirements--what-to-install)
4. [Data Downloads](#4-data-downloads)
5. [Part 1 — Standard Workflow (Steps 1–6)](#5-part-1--standard-workflow)
   - [Step 1 — Load Data](#step-1--load-data)
   - [Step 2 — Filter Non-Standard Chromosomes](#step-2--filter-non-standard-chromosomes)
   - [Step 3 — Gene Annotations](#step-3--gene-annotations-hg38)
   - [Step 4 — Quality Control](#step-4--quality-control)
   - [Step 5 — TF-IDF Normalisation & LSI](#step-5--tf-idf-normalisation--lsi)
   - [Step 6 — UMAP & Clustering](#step-6--umap--leiden-clustering)
6. [Part 2 — Downstream Workflow (Steps 7–11)](#6-part-2--downstream-workflow)
   - [Step 7 — Gene Activity Scoring](#step-7--gene-activity-scoring)
   - [Step 8 — Marker Gene Visualisation](#step-8--marker-gene-visualisation)
   - [Step 9 — Label Transfer from scRNA-seq](#step-9--label-transfer-from-scrna-seq)
   - [Step 10 — Differential Accessibility](#step-10--differential-accessibility-analysis)
   - [Step 11 — Coverage Plots](#step-11--coverage-genome-browser-plots)
7. [Statistical Methods Explained](#7-statistical-methods-explained)
8. [Tunable Parameters](#8-tunable-parameters)
9. [Common Errors & Fixes](#9-common-errors--fixes)
10. [Project Structure](#10-project-structure)
11. [References & Further Reading](#11-references--further-reading)

---

## 1. What is scATAC-seq? — A Beginner's Guide

### The central question
Standard scRNA-seq asks *"what genes is this cell expressing right now?"*  
scATAC-seq asks a deeper question: *"which parts of this cell's genome are physically
accessible — open and available to be switched on?"*

### The biology in plain English

Your DNA is about 2 metres long, packed into a nucleus smaller than a full stop.
It is wound around protein spools called **nucleosomes** to fit:

```
Closed chromatin (gene OFF):  DNA tightly wound → inaccessible
                               ██████████████████████████████
                               ↑ nucleosome  ↑ nucleosome

Open chromatin (gene ON):     DNA unwound → transcription factors can bind
                               ···[TF]·····[TF]···················
                               ↑ open/accessible region
```

Regions that are **open** are available for transcription factors to bind and
activate nearby genes. **scATAC-seq maps these open regions at single-cell resolution.**

### How the experiment works

The **Tn5 transposase** enzyme is added to cell nuclei. It acts like molecular
scissors that preferentially cut open (accessible) chromatin and tag both cut
ends with sequencing adapters:

```
Step 1 — Cells are permeabilised (nuclei isolated)
Step 2 — Tn5 transposase added:
         Open chromatin:   ···Tn5···Tn5···Tn5···  ← many cuts → many fragments
         Closed chromatin: ████████████████████   ← no cuts → few fragments
Step 3 — Cut fragments are sequenced
Step 4 — Fragment positions are mapped to the genome
         More fragments at a position = more accessible at that position
```

### scATAC-seq vs scRNA-seq — what's different?

| Property | scRNA-seq | scATAC-seq |
|----------|-----------|------------|
| What is measured | Gene expression (mRNA) | Chromatin accessibility (DNA) |
| Data | Integer counts per gene | Sparse 0/1 per genomic peak |
| Sparsity | ~90% zeros | ~95–99% zeros |
| Matrix size | cells × ~20,000 genes | cells × ~100,000–200,000 peaks |
| Normalisation | Log-normalise | **TF-IDF** (from text mining) |
| Dim. reduction | PCA | **LSI** (Latent Semantic Indexing) |
| Critical issue | None | **LSI component 1 = depth artefact — must skip** |
| Main R toolkit | Seurat | **Signac** (extends Seurat) |

---

## 2. What Does This Pipeline Produce?

Running this script from start to finish produces:

| Output | Description |
|--------|-------------|
| QC plots | Density scatter, fragment histogram, violin plots |
| UMAP coloured by clusters | 15 chromatin-based clusters of PBMCs |
| Gene activity feature plots | 6 canonical PBMC marker genes on UMAP |
| Label transfer UMAP | Cell types predicted from matched scRNA-seq |
| Differential accessibility | Peaks more open in CD4 Naive vs CD14+ Monocytes |
| Coverage browser plots | Genome-browser tracks for each cell type |
| `pbmc_atac_final.rds` | Fully annotated Seurat object |
| `pbmc_atac_clusters.csv` | Cell barcodes with cluster and cell type labels |

**Cell types identified in this dataset:**

| Cell type | Biological role |
|-----------|-----------------|
| CD4 Naive | Helper T cells — resting, never encountered antigen |
| CD4 Memory | Helper T cells — previously activated, long-lived |
| CD8 Naive | Cytotoxic T cells — resting |
| CD8 effector | Cytotoxic T cells — activated, killing infected cells |
| Double negative T cell | Rare T cell subset (CD4⁻ CD8⁻) |
| NK dim / NK bright | Natural killer cells — innate immune effectors |
| CD14+ Monocytes | Classical monocytes — phagocytosis, inflammation |
| CD16+ Monocytes | Non-classical monocytes — patrol blood vessels |
| Dendritic cell | Antigen presenting cells |
| pDC | Plasmacytoid dendritic cells — anti-viral interferon |
| B cell progenitor | Immature B cells |
| pre-B cell | B cell development stage |
| Platelet | Thrombocytes — blood clotting |

---

## 3. Requirements — What to Install

**R 4.2 or higher** required.

```r
# CRAN packages
install.packages(c("remotes", "Matrix", "irlba", "ggplot2", "patchwork", "httr"))

# Signac (latest from GitHub)
remotes::install_github("stuart-lab/signac", ref = "develop")

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("AnnotationHub", "GenomicRanges", "GenomeInfoDb", "biovizBase"))
```

> ⚠️ **`biovizBase` is required** — without it, `GetGRangesFromEnsDb()` fails.
> This is not mentioned in the standard Signac documentation.

---

## 4. Data Downloads

Run these commands in your terminal inside your working directory:

```bash
# Peak-barcode matrix (~162 MB)
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5 \
     -O filtered_peak_bc_matrix.h5

# Per-cell QC metadata (~31 MB)
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv \
     -O singlecell.csv

# Fragment file (~2.5 GB) + index
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz \
     -O fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz.tbi \
     -O fragments.tsv.gz.tbi

# scRNA-seq reference for label transfer (~164 MB)
wget https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds
```

> ⚠️ **All files download into your working directory** (not a `data/` subfolder).
> The script reads them from the working directory. Only saved outputs go into `data/`.

---

## 5. Part 1 — Standard Workflow

---

### Step 1 — Load Data

```r
counts      <- Read10X_h5(filename = "filtered_peak_bc_matrix.h5")
metadata    <- read.csv("singlecell.csv", header=TRUE, row.names=1)
chrom_assay <- CreateChromatinAssay(counts=counts, sep=c(":","-"),
                                     fragments="fragments.tsv.gz",
                                     min.cells=10, min.features=200)
pbmc <- CreateSeuratObject(counts=chrom_assay, assay="peaks", meta.data=metadata)
```

**What is a fragment file?**  
The largest file (~2.5 GB) records every DNA cut position made by Tn5:
```
chr1    10073    10498    ACGTACGT-1    1
chrom   start    end      cell barcode  duplicates
```
More fragments at a genomic location = more accessible chromatin there.

**What is the peak matrix?**  
CellRanger ATAC counts how many fragments from each cell overlap each called
peak. The result is a sparse **cells × peaks** integer matrix (like RNA's
cells × genes matrix).

**Why `assay = "peaks"` not `"ATAC"`?**  
Current Signac convention names the assay `"peaks"`. This means all metadata
columns are named `nCount_peaks` and `nFeature_peaks` — not `nCount_ATAC`.
Using the wrong name in downstream functions will cause errors.

---

### Step 2 — Filter Non-Standard Chromosomes

```r
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc       <- pbmc[as.vector(peaks.keep), ]
```

CellRanger includes peaks from scaffold sequences like `KI270713.1` (chromosome
patches). These are not real chromosomes, add noise, and can cause warnings.
This step keeps only `chr1–chr22`, `chrX`, `chrY`, `chrMT`.

---

### Step 3 — Gene Annotations (hg38)

```r
ah        <- AnnotationHub()
ensdb_v98 <- ah[["AH75011"]]          # Ensembl v98 for hg38
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)
seqlevelsStyle(annotations) <- "UCSC" # 1,2,X → chr1,chr2,chrX
genome(annotations) <- "hg38"
Annotation(pbmc) <- annotations
```

**Why does this step matter so much?**  
Without gene annotations, Signac cannot:
- Calculate TSS enrichment scores for QC
- Score gene activity from peaks
- Show gene names on coverage plots

**The chromosome naming pitfall:**  
Ensembl databases use bare numbers (`1`, `2`, `X`) while 10x data uses
UCSC style (`chr1`, `chr2`, `chrX`). `seqlevelsStyle(annotations) <- "UCSC"`
converts automatically. Getting this wrong causes all annotation lookups to
return empty results, silently corrupting downstream analysis.

---

### Step 4 — Quality Control

#### The five QC metrics

| Metric | Column name | Good range | What it measures |
|--------|-------------|-----------|-----------------|
| Total fragments | `nCount_peaks` | 9,771–52,418 (5th–95th %ile) | Total DNA cuts per cell |
| TSS enrichment | `TSS.enrichment` | 4.59–6.55 (5th–95th %ile) | Signal at gene start sites |
| Nucleosome signal | `nucleosome_signal` | < 4 | Chromatin structural quality |
| Reads in peaks | `pct_reads_in_peaks` | > FRiP cutoff | Signal quality vs background |
| Blacklist ratio | `blacklist_ratio` | ≈ 0 | Artefact region contamination |

#### 📊 Plot 1 — Density Scatter: nCount_peaks vs TSS enrichment

![Density scatter plot showing nCount_peaks vs TSS enrichment with quantile lines](Rplot00_density_scatter.png)

**How to read this plot:**  
Each dot is one cell. The colour shows density (yellow = many cells, purple = few).
The red lines mark the 5th, 10th, 90th, and 95th percentiles.

- The **dense yellow cloud** (top-right) = high-quality cells with many fragments AND strong TSS signal
- **Isolated black dots** (bottom-left, low counts) = empty droplets or dead cells
- **Outliers** (extreme right) = potential doublets

The quantiles shown in the title (nCount: 5%:9,771 ... 95%:52,418 | TSS: 5%:4.59 ... 95%:6.55) are used to set the dynamic filtering thresholds.

---

#### 📊 Plot 2 — Fragment Size Histogram

![Fragment histogram showing NS<4 and NS>4 groups](Rplot01_fragment_histogram.png)

**What is nucleosome banding and why does it matter?**  

DNA is wound around nucleosomes every ~200 bp. Healthy ATAC-seq cuts
preferentially in the gaps between nucleosomes:

```
Nucleosome positions on DNA:
─────[NUCL]────────[NUCL]────────[NUCL]──────
     │200bp│       │200bp│       │200bp│
Tn5 cuts in the gaps (accessible linker DNA)

Resulting fragment sizes:
<200bp  = sub-nucleosomal (cuts within one open region)
~200bp  = mono-nucleosomal
~400bp  = di-nucleosomal
~600bp  = tri-nucleosomal
```

**Left panel (NS < 4 — good cells):** Clear sub-nucleosomal enrichment (tall peak at <200bp) with decreasing secondary peaks at 200bp and 400bp — the expected healthy banding pattern.

**Right panel (NS > 4 — bad cells):** The banding pattern is dominated by the mono-nucleosomal peak (~200bp) with very few sub-nucleosomal fragments. These cells have poor Tn5 insertion efficiency.

---

#### 📊 Plot 3 — QC Violin Plots

![Violin plots of all 5 QC metrics](Rplot02_violin_qc.png)

**How to read violin plots:**  
Each violin shows the distribution of cells for one metric. Wider = more cells
at that value. The centre line = median.

Key observations from this data:
- **nCount_peaks:** Wide range (0–100,000+) with a concentration around 20,000–50,000 for good cells
- **TSS.enrichment:** Clean bell shape centred around 5–6 (healthy signal)
- **blacklist_ratio:** Nearly all cells at 0 — very clean dataset
- **nucleosome_signal:** Most cells < 2, with a small tail of poor-quality cells above 4
- **pct_reads_in_peaks:** Concentration at 60–80% indicating high signal quality

**Dynamic thresholds used (computed from the data's percentiles):**
```r
ncount_min = max(200,  quantile(nCount_peaks,       0.02))
ncount_max = min(1e6,  quantile(nCount_peaks,       0.98))
frip_min   = max(10,   quantile(pct_reads_in_peaks, 0.05))
ns_max     = min(10,   quantile(nucleosome_signal,  0.95))
bl_max     = min(0.10, quantile(blacklist_ratio,    0.95))
```
Unlike fixed thresholds, these adapt to your specific dataset and always
retain a sensible proportion of cells.

---

### Step 5 — TF-IDF Normalisation & LSI

```r
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)
DepthCor(pbmc)
```

#### Why not use log-normalisation like scRNA-seq?

ATAC-seq peaks are **near-binary**: a peak is either open (1) or closed (0) in
a given cell. Log-normalisation is designed for continuous count data and works
poorly with binary data.

**TF-IDF** (Term Frequency — Inverse Document Frequency) is borrowed from
text mining where documents are bags of words (also binary: word either
appears or not):

```
TF  = (cuts in this peak for this cell) / (total cuts for this cell)
      → normalises for sequencing depth

IDF = log(1 + total cells / cells with this peak open)
      → upweights RARE, cell-type-specific peaks
      → downweights COMMON peaks open in all cells

TF-IDF = TF × IDF  →  the final normalised value
```

#### What is LSI?

**LSI (Latent Semantic Indexing)** is the ATAC-seq equivalent of PCA.
Applied to the TF-IDF matrix, it compresses ~165,000 peaks into ~30
latent dimensions capturing the main patterns of chromatin accessibility.

#### 📊 Plot 4 — Depth Correlation Plot (Critical!)

![Depth correlation plot showing LSI1 correlates with sequencing depth](Rplot03_depth_cor.png)

**This is the most important diagnostic plot in scATAC-seq.**

The plot shows the correlation between each LSI component and total
sequencing depth (total fragments per cell).

- **Component 1:** Correlation ≈ +0.88 → **This is a technical artefact, not biology**
- **Components 2–10:** Correlation ≈ 0 → **These capture biology**

**If you include Component 1 in UMAP and clustering, cells will cluster by
how deeply they were sequenced — not by cell type.** This is why all
downstream analysis uses `dims = 2:30`.

---

### Step 6 — UMAP & Leiden Clustering

```r
pbmc <- RunUMAP(pbmc, reduction="lsi", dims=2:30)
pbmc <- FindNeighbors(pbmc, reduction="lsi", dims=2:30)
pbmc <- FindClusters(pbmc, algorithm=3, resolution=0.5)
```

#### What is UMAP?

UMAP (Uniform Manifold Approximation and Projection) projects the 29 LSI
dimensions (2–30) into 2D for visualisation. Similar cells end up close
together.

#### What is Leiden clustering?

Leiden clustering (`algorithm=3`) builds a graph where each cell is connected
to its nearest neighbours, then finds densely connected communities (clusters).

#### 📊 Plot 5 — UMAP with Leiden Clusters

![UMAP showing 15 Leiden clusters of PBMCs](Rplot04_umap_clusters.png)

**15 clusters (0–14) are identified.** At this stage the clusters have no
biological labels — just numbers. Notice:
- Large upper group (clusters 1, 2, 3) = likely T cells (most abundant in blood)
- Isolated cluster (6, right) = likely B cells
- Dense lower cluster (0) = likely monocytes
- Small isolated clusters (11, 14) = rare populations

Clusters are assigned biological identities in Step 9 using label transfer.

---

## 6. Part 2 — Downstream Workflow

---

### Step 7 — Gene Activity Scoring

```r
gene.activities <- GeneActivity(pbmc)
pbmc[["RNA"]] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(pbmc, assay="RNA", normalization.method="LogNormalize",
                      scale.factor = median(pbmc$nCount_RNA))
```

**What is gene activity scoring?**  
ATAC-seq peaks tell us about chromatin accessibility — not gene expression.
But we can **estimate** how active a gene is by summing all fragments
overlapping its gene body and 2 kb upstream promoter:

```
Gene activity score =  Σ (fragment counts in peaks overlapping gene body + promoter)
                       ────────────────────────────────────────────────────────────
                       Example: 15 peaks overlap the CD3D gene body
                                Those peaks collectively have 847 fragments
                                → Gene activity score for CD3D = 847
```

> ⚠️ **Important caveat:** Open chromatin ≠ active transcription. Gene activity
> is an **approximation** used for exploratory visualisation. Use label transfer
> (Step 9) for reliable cell type assignment.

---

### Step 8 — Marker Gene Visualisation

```r
DefaultAssay(pbmc) <- "RNA"
FeaturePlot(pbmc, features=c("MS4A1","CD3D","LEF1","NKG7","TREM1","LYZ"),
            max.cutoff="q95", ncol=3)
```

#### 📊 Plot 6 — Gene Activity Feature Plots

![Feature plots showing gene activity for 6 canonical PBMC markers](Rplot05_marker_genes.png)

**How to read these plots:**  
Each plot shows one gene. Blue = high gene activity (open chromatin), grey = low.
Clusters where a gene is blue = likely that cell type.

**What each marker tells us:**

| Gene | What it marks | What we see |
|------|--------------|-------------|
| `MS4A1` (CD20) | B cells | Isolated right cluster (cluster 6) — confirms B cells |
| `CD3D` | All T cells | Large upper clusters — confirms T cells |
| `LEF1` | Naive T cells | Concentrated in top-left (clusters 1, 3) |
| `NKG7` | NK cells | Small isolated left cluster (cluster 5) |
| `TREM1` | Monocytes | Large lower cluster (cluster 0) |
| `LYZ` | Monocytes | Same lower cluster |

`max.cutoff = "q95"` caps the colour scale at the 95th percentile so that
rare outlier cells don't make all other cells appear grey.

---

### Step 9 — Label Transfer from scRNA-seq

```r
pbmc_rna <- readRDS("pbmc_10k_v3.rds")
transfer.anchors <- FindTransferAnchors(reference=pbmc_rna, query=pbmc,
                                         reduction="cca")
predicted.labels <- TransferData(anchorset=transfer.anchors,
                                  refdata=pbmc_rna$celltype,
                                  weight.reduction=pbmc[["lsi"]], dims=2:30)
pbmc <- AddMetaData(pbmc, metadata=predicted.labels)
pbmc <- subset(pbmc, subset=prediction.score.max > 0.5)
```

**What is label transfer?**  
We have a **matched scRNA-seq dataset** from the same tissue (PBMCs from the
same donors) where cell types are already known.

Label transfer finds "anchor" pairs — one ATAC cell and one RNA cell that look
similar in a shared embedding space — and borrows the RNA cell's label:

```
scRNA-seq cell with known label "CD4 Naive"
         ↕ (similar in CCA embedding space)
scATAC-seq cell with unknown identity
         → Predicted label: "CD4 Naive" (score: 0.87)
```

**Why CCA (Canonical Correlation Analysis)?**  
CCA finds dimensions that maximise correlation **between** the two datasets.
This is appropriate for cross-modality integration (ATAC → RNA). CCA is
different from `rpca` which is used for same-modality batch correction.

**The confidence filter** `prediction.score.max > 0.5` removes cells where
the model could not confidently assign a cell type.

#### 📊 Plot 7 — scATAC-seq vs scRNA-seq Before Transfer

![Side-by-side UMAP of ATAC clusters and RNA cell types before transfer](Rplot06_atac_vs_rna.png)

The left panel shows the ATAC-seq UMAP (unlabelled clusters). The right panel
shows the RNA reference with known cell type labels. Note that the overall
topology is similar — both show a large T cell cloud, isolated B cell cluster,
and a monocyte cluster — confirming the integration will work well.

#### 📊 Plot 8 — Label Transfer Results

![Side-by-side UMAP showing predicted ATAC labels vs true RNA labels](Rplot07_label_transfer.png)

After transfer, 14 distinct cell types are predicted in the ATAC-seq data.
Compare left (ATAC predicted) with right (RNA ground truth) — the cluster
shapes and relative sizes closely match, confirming the transfer was accurate.

Key observations:
- CD4 Naive and CD4 Memory correctly separate in the upper cluster region
- CD8 effector and CD8 Naive cluster together as expected (similar chromatin)
- Monocyte subtypes (CD14+ and CD16+) correctly identify the lower large cluster
- Rare populations (pDC, Platelet) are identified in small isolated clusters

---

### Step 10 — Differential Accessibility Analysis

```r
DefaultAssay(pbmc) <- "peaks"
da_peaks <- FindMarkers(pbmc, ident.1="CD4 Naive", ident.2="CD14+ Monocytes",
                         test.use="LR", latent.vars="nCount_peaks")
```

**What are differentially accessible (DA) peaks?**  
Genomic regions that are significantly **more open** in one cell type compared
to another. These correspond to regulatory elements (promoters, enhancers)
that are active specifically in that cell type.

**Why logistic regression (`test.use="LR"`) not Wilcoxon?**  
- Wilcoxon test is rank-based and works well for continuous RNA counts
- ATAC data is near-binary (0/1) — logistic regression handles this better
- `latent.vars="nCount_peaks"` removes depth-driven false positives: a peak
  might appear "more open" in one group simply because that group was more
  deeply sequenced

**The results table columns:**

| Column | Meaning |
|--------|---------|
| `avg_log2FC > 0` | Peak more open in CD4 Naive |
| `avg_log2FC < 0` | Peak more open in CD14+ Monocytes |
| `p_val_adj` | Adjusted p-value (Bonferroni) — use this for significance |
| `pct.1` | % of CD4 Naive cells with this peak open |
| `pct.2` | % of CD14+ Monocytes cells with this peak open |

#### 📊 Plot 9 — Top DA Peak Violin + UMAP

![Violin plot and feature plot for the top differentially accessible peak](Rplot08_da_peak.png)

The top DA peak is `chr12-119988511-119989430` (overlapping the **BICDL1** gene).

**Left (violin):** The peak is highly accessible in CD4 Naive T cells (teal, wide
violin) but closed in CD14+ Monocytes (pink, flat at 0). This is a stark
contrast.

**Right (UMAP):** Blue = high accessibility. The upper T cell clusters light up
strongly, while the lower monocyte region is grey — confirming this is a
T-cell-specific regulatory element.

---

### Step 11 — Coverage (Genome Browser) Plots

```r
CoveragePlot(pbmc, region=rownames(da_peaks)[1],
             extend.upstream=40000, extend.downstream=20000)
CoveragePlot(pbmc, region="CD8A")
```

Coverage plots are the single-cell equivalent of a **genome browser track**.
They show fragment density along a genomic region, one track per cell type.

**Reading a coverage plot:**
```
Track structure (top to bottom):
  ┌─────────────────────────────────┐
  │ Cell type 1: ─────▁▂████▂▁───── │  ← fragment density
  │ Cell type 2: ─────────────────── │  ← closed = flat
  │ Cell type 3: ───▁▂███▂▁───────── │  ← accessible
  ├─────────────────────────────────┤
  │ Genes:      →→→[GENE]→→→→→→→→→  │  ← gene structure
  ├─────────────────────────────────┤
  │ Peaks:      ████    ████         │  ← called peaks
  └─────────────────────────────────┘
```

#### 📊 Plot 10 — Coverage at Top DA Peak (BICDL1 region)

![Coverage plot for top DA peak on chr12 in the BICDL1 region](Rplot09_coverage_da_peak.png)

**Region:** chr12:119,960,000–120,010,000 (near the BICDL1 gene)

Key observations:
- **CD4 Naive, CD4 Memory, CD8 cells, NK cells** all show a clear peak at the same position — a T/NK lineage-specific regulatory element
- **CD14+ Monocytes, CD16+ Monocytes, Dendritic cells** show flat signal — this region is closed in myeloid cells
- **B cell progenitor** also has a small peak — may share some lymphoid regulatory programme

The grey bars at the bottom (Peaks track) show the called peaks — the sharp
signal aligns with the called peak, confirming CellRanger correctly identified
this accessible region.

---

#### 📊 Plot 11 — Coverage at CD8A Gene

![Coverage plot for the CD8A gene on chr2](Rplot10_coverage_cd8a.png)

**Gene:** CD8A (chr2:86,785,000–86,808,000)

CD8A encodes the alpha chain of the CD8 co-receptor on cytotoxic T cells.

Key observations:
- **CD8 effector and CD8 Naive** show the strongest peaks — as expected, they express CD8A
- **Double negative T cells and NK dim** also show signal (these can express CD8-like receptors)
- **CD4 Naive, CD14+ Monocytes, B cells** are flat — confirming cell-type specificity
- Multiple peaks visible across the gene body — some at the promoter (left), some in introns
- The **Genes track** shows the CD8A gene structure (exons as filled boxes, introns as arrows)

---

## 7. Statistical Methods Explained

### TF-IDF (Term Frequency — Inverse Document Frequency)

Normalises the peak-cell matrix accounting for both cell depth and peak rarity:

```
TF(peak p, cell c)  =  count(p,c) / Σ_peaks count(p,c)
IDF(peak p)          =  log(1 + N_cells / N_cells_with_peak)
TF-IDF(p,c)         =  TF × IDF
```

### LSI (Latent Semantic Indexing)

Singular Value Decomposition (SVD) applied to the TF-IDF matrix:
```
TF-IDF matrix (cells × peaks) = U × Σ × V^T
```
The columns of `U` are the LSI components (like PCs in PCA).
**Component 1 always captures depth** — always skip it.

### Leiden Clustering

Optimises a **modularity** function to find densely-connected communities
in the k-nearest-neighbour graph. The `resolution` parameter controls
community size: higher = more/smaller clusters.

### Logistic Regression (DA testing)

Models peak accessibility as a binary outcome:
```
logit(P(peak open)) = β₀ + β₁(cell_type) + β₂(log(nCount_peaks))
```
`β₁` = the effect of cell type, controlling for depth via `β₂`.
The test statistic follows a likelihood ratio test (LR test).

### CCA (Canonical Correlation Analysis) for Label Transfer

Finds linear combinations of features in ATAC and RNA that are maximally
correlated across the two datasets. Cells are embedded in this shared space
and labelled by k-nearest-neighbour voting from the RNA reference.

---

## 8. Tunable Parameters

| Parameter | Default | Location | Effect |
|-----------|---------|----------|--------|
| `min.cells` | `10` | `CreateChromatinAssay` | Min cells a peak must appear in |
| `min.features` | `200` | `CreateChromatinAssay` | Min peaks a cell must have |
| QC filter percentiles | 2nd/98th/5th/95th | Step 4 | Adjust stringency of cell filtering |
| `min.cutoff` | `"q0"` | `FindTopFeatures` | `"q75"` for top 25% peaks only (faster) |
| LSI dims | `2:30` | `RunUMAP`, `FindNeighbors` | Always start at 2; check `DepthCor()` |
| `resolution` | `0.5` | `FindClusters` | Higher = more/smaller clusters |
| Prediction threshold | `0.5` | `subset()` | Raise for stricter label transfer |
| DA comparison | CD4 Naive vs CD14+ Mono | Step 10 | Any two cell types |
| Coverage window | 40,000 / 20,000 bp | `CoveragePlot` | Show more/less genomic context |

---

## 9. Common Errors & Fixes

| Error message | Cause | Fix |
|---|---|---|
| `Please install biovizBase` | Package not installed | `BiocManager::install("biovizBase")` |
| `chrchrX` in seqlevels | Double chr prefix from EnsDb | Use `seqlevelsStyle(annotations) <- "UCSC"` |
| `No cells found` after QC | Fixed thresholds too strict | Use percentile-based dynamic thresholds |
| `No matching chromosomes` in GeneActivity | Annotation chr mismatch | Check `seqlevels(Annotation(pbmc))` starts with `chr`, not `chrchr` |
| `'fast' argument deprecated` | Old Signac API | Remove `fast=FALSE` from `TSSEnrichment()` |
| `'pcaproject' not valid` | Seurat v5 removed it | Use `reduction="cca"` for cross-modality transfer |
| `'ATAC' not found` | Assay named `"peaks"` | Use `assay="peaks"` and `nCount_peaks` throughout |
| `cannot open connection data/` | `data/` directory missing | Add `if (!dir.exists("data")) dir.create("data")` |

---

## 10. Project Structure

```
your-project/
│
├── scATACSeq_complete_pipeline_fixed.R  ← Main analysis script
├── README.md                            ← This file
│
├── [working directory — data files]
│   ├── filtered_peak_bc_matrix.h5       ← Downloaded
│   ├── singlecell.csv                   ← Downloaded
│   ├── fragments.tsv.gz                 ← Downloaded
│   ├── fragments.tsv.gz.tbi             ← Downloaded
│   └── pbmc_10k_v3.rds                  ← Downloaded (RNA reference)
│
├── [plots/]
│   ├── Rplot00_density_scatter.png
│   ├── Rplot01_fragment_histogram.png
│   ├── Rplot02_violin_qc.png
│   ├── Rplot03_depth_cor.png
│   ├── Rplot04_umap_clusters.png
│   ├── Rplot05_marker_genes.png
│   ├── Rplot06_atac_vs_rna.png
│   ├── Rplot07_label_transfer.png
│   ├── Rplot08_da_peak.png
│   ├── Rplot09_coverage_da_peak.png
│   └── Rplot10_coverage_cd8a.png
│
└── data/                                ← Created automatically by script
    ├── pbmc_atac_part1.rds              ← Saved after Part 1
    ├── pbmc_atac_final.rds              ← Saved after Part 2
    └── pbmc_atac_clusters.csv          ← Cell type assignments
```

> ⚠️ **Never commit large files to GitHub.** Add to `.gitignore`:
> ```
> *.h5
> *.tsv.gz
> *.tbi
> *.rds
> data/
> ```

---

## 11. References & Further Reading

### Tools used
- [Signac — Stuart Lab](https://stuartlab.org/signac/) — the main scATAC-seq analysis toolkit
- [Seurat v5 — Satija Lab](https://satijalab.org/seurat/) — single-cell framework
- [AnnotationHub — Bioconductor](https://bioconductor.org/packages/AnnotationHub/)

### Key papers
- [Stuart et al. 2021 — Signac (Nature Methods)](https://www.nature.com/articles/s41592-021-01282-5)
- [Cusanovich et al. 2018 — TF-IDF for scATAC-seq (Science)](https://www.science.org/doi/10.1126/science.aab1601)
- [Stuart et al. 2019 — Label transfer across modalities (Cell)](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8)
- [Heumos et al. 2023 — Best practices for single-cell analysis (Nature Reviews Genetics)](https://www.nature.com/articles/s41576-023-00586-w)

### Learning resources
- 📺 **[Bioinformagician YouTube](https://www.youtube.com/@Bioinformagician)** — Excellent tutorials on scRNA-seq and scATAC-seq in R (portions of this pipeline adapted from their content)
- [Official Signac PBMC vignette](https://stuartlab.org/signac/articles/pbmc_vignette)
- [10x Genomics ATAC-seq documentation](https://support.10xgenomics.com/single-cell-atac)
- [ENCODE blacklist regions](https://github.com/Boyle-Lab/Blacklist)

### Data source
- [10x Genomics PBMC 10k ATACv2 dataset](https://www.10xgenomics.com/datasets/10-k-human-pbm-cs-atac-v-1-1-chromium-x-1-standard-2-0-0)

---

*Pipeline completed April 2026 — Signac v1.17, Seurat v5, R 4.6.0 (Windows x64)*
