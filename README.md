# 🧬 scATAC-seq Analysis Pipeline & Interactive EDA App
## Human PBMCs — 10x Genomics ATACv2 (10,000 cells)
### A Beginner's Complete Guide to the Methods, Results & Code



---

## 📋 Table of Contents

1. [Repository Contents & File Guide](#1-repository-contents--file-guide)
2. [What is scATAC-seq? A Plain English Explanation](#2-what-is-scatac-seq-a-plain-english-explanation)
3. [How to Read This Readme](#3-how-to-read-this-readme)
4. [Part 1 — Standard Workflow: Methods & Results](#4-part-1--standard-workflow-methods--results)
   - [Step 1: Loading the Data](#step-1-loading-the-data)
   - [Step 2: Filtering Non-Standard Chromosomes](#step-2-filtering-non-standard-chromosomes)
   - [Step 3: Gene Annotations](#step-3-gene-annotations)
   - [Step 4: Quality Control — The 5 QC Metrics](#step-4-quality-control--the-5-qc-metrics)
   - [Step 5: TF-IDF Normalisation](#step-5-tf-idf-normalisation)
   - [Step 6: LSI Dimensionality Reduction](#step-6-lsi-latent-semantic-indexing)
   - [Step 7: UMAP + Leiden Clustering](#step-7-umap--leiden-clustering)
5. [Part 2 — Downstream Workflow: Methods & Results](#5-part-2--downstream-workflow-methods--results)
   - [Step 8: Gene Activity Scoring](#step-8-gene-activity-scoring)
   - [Step 9: Marker Gene Visualisation](#step-9-marker-gene-visualisation)
   - [Step 10: Label Transfer from scRNA-seq](#step-10-label-transfer-from-scrna-seq)
   - [Step 11: Differential Accessibility Analysis](#step-11-differential-accessibility-analysis)
   - [Step 12: Coverage Plots](#step-12-coverage-plots)
6. [Interactive Shiny App — User Guide](#6-interactive-shiny-app--user-guide)
7. [How to Adapt This Pipeline for Your Own Data](#7-how-to-adapt-this-pipeline-for-your-own-data)
8. [Requirements & Data Downloads](#8-requirements--data-downloads)
9. [Common Errors & Fixes](#9-common-errors--fixes)
10. [Further Reading & References](#10-further-reading--references)

---

## 1. Repository Contents & File Guide

```
your-project/
│
├── ATACSeq_complete_pipeline.R   ← Run this FIRST — full analysis
├── Shiney_App_For_EDA.R          ← Run this SECOND — interactive app
├── README.md                     ← This file
│
├── Rplot.png      ← QC Plot 1:  Density scatter (nCount vs TSS)
├── Rplot01.png    ← QC Plot 2:  Fragment size histogram
├── Rplot02.png    ← QC Plot 3:  5-metric violin plots
├── Rplot03.png    ← QC Plot 4:  LSI depth correlation
├── Rplot04.png    ← UMAP 1:     Leiden clusters (before labelling)
├── Rplot05.png    ← UMAP 2:     6 marker genes activity
├── Rplot06.png    ← Transfer 1: ATAC vs RNA side-by-side (before)
├── Rplot07.png    ← Transfer 2: Predicted cell types (after)
├── Rplot08.png    ← DA Plot:    Top DA peak violin + UMAP
├── Rplot09.png    ← Coverage 1: chr12 — BICDL1 gene region
├── Rplot10.png    ← Coverage 2: chr2  — CD8A gene
│
└── data/
    ├── pbmc_atac_final.rds         ← Full annotated object (needed by app)
    ├── pbmc_atac_part1.rds         ← Part 1 checkpoint
    ├── da_peaks_results.rds        ← DA results table
    └── pbmc_atac_clusters.csv      ← Cell barcode → cell type table
```

> ⚠️ **Never commit .rds, .h5, .tsv.gz, or .tbi files to GitHub** — they are
> hundreds of megabytes. Add `data/` and `*.h5` to your `.gitignore`.

---

## 2. What is scATAC-seq? A Plain English Explanation

### The Big Picture

Every cell in your body has the same DNA but behaves differently — liver
cells make liver proteins, T cells fight infection. How? The answer lies in
which parts of the DNA are **physically accessible** to the machinery that
reads genes.

Think of DNA like a book in a library:
- **Open chromatin** = the book is on the desk, open, being read
- **Closed chromatin** = the book is locked in a cabinet

**scATAC-seq** (single-cell Assay for Transposase-Accessible Chromatin using
sequencing) takes a snapshot of which pages are "open" in each individual cell.

### How the Experiment Works

```
Individual cells
       ↓
  Add Tn5 transposase enzyme (like molecular scissors)
       ↓
  Tn5 cuts wherever chromatin is OPEN (accessible)
  Tn5 CANNOT cut where chromatin is CLOSED (wrapped in nucleosomes)
       ↓
  Sequence the cut fragments
       ↓
  More fragments at a location = more accessible chromatin there
```

### The Data We Get

The result is a **cells × peaks matrix** — similar to scRNA-seq's cells × genes
matrix, but instead of gene expression counts, we have **chromatin accessibility**
at 165,000+ genomic regions (peaks) for each of 10,000 cells:

```
            Peak_1          Peak_2          Peak_3    ...
Cell_1        0               1               0
Cell_2        0               0               2
Cell_3        3               1               0
...
```

- `0` = closed chromatin (Tn5 did not cut here in this cell)
- `1` or `2` = open chromatin (Tn5 cut here — region is accessible)

---

## 3. How to Read This README

Each section follows the same structure:

```
WHAT IT DOES    — plain English explanation, no jargon
WHY WE DO IT    — the biological / statistical reason
THE MATHS       — the formula (shown simply, with explanation)
THE CODE        — the R function used (one line, no detail needed)
THE RESULT      — what the output plot shows
WHAT TO CHANGE  — which numbers to adjust for your own data
```

---

## 4. Part 1 — Standard Workflow: Methods & Results

---

### Step 1: Loading the Data

**WHAT IT DOES**

Reads three input files into R:
1. **Peak-barcode matrix** (`.h5`) — the cells × peaks count table
2. **Metadata CSV** — per-cell QC statistics pre-computed by CellRanger ATAC
3. **Fragment file** (`.tsv.gz`) — the raw list of every DNA cut site

**WHY WE DO IT**

The `.h5` matrix alone tells you *how many* fragments fell on each peak per
cell, but to compute QC metrics like TSS enrichment (Step 4), we need to
look up exact fragment positions on the genome — that requires the fragment file.

**THE CODE**

```r
counts      <- Read10X_h5("filtered_peak_bc_matrix.h5")
metadata    <- read.csv("singlecell.csv")
chrom_assay <- CreateChromatinAssay(counts = counts, fragments = "fragments.tsv.gz")
pbmc        <- CreateSeuratObject(counts = chrom_assay, assay = "peaks")
```

**WHAT TO CHANGE**

```r
# 👇 Update these three filenames to your own data files
counts   <- Read10X_h5("YOUR_matrix.h5")
metadata <- read.csv("YOUR_singlecell.csv")
CreateChromatinAssay(counts = counts, fragments = "YOUR_fragments.tsv.gz",
                     min.cells    = 10,   # lower if dataset is small
                     min.features = 200)  # lower if cells are sparse
```

---

### Step 2: Filtering Non-Standard Chromosomes

**WHAT IT DOES**

Removes any peaks not on the 22 autosomes (chr1–chr22) plus chrX, chrY.

**WHY WE DO IT**

CellRanger ATAC calls peaks on all sequences in the reference genome,
including unplaced scaffold sequences (e.g. `KI270713.1`) that represent
incomplete or patch regions of the human genome. These are not real
chromosomes, add noise to the analysis, and cause warnings in downstream
tools.

**THE CODE**

```r
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc       <- pbmc[as.vector(peaks.keep), ]
```

**RESULT:** Reduces ~165,434 peaks → ~165,376 peaks (removes ~58 scaffold peaks).

---

### Step 3: Gene Annotations

**WHAT IT DOES**

Downloads the human gene annotation database (hg38 / GRCh38, Ensembl v98)
and attaches it to the Seurat object. This maps every peak to the nearest gene.

**WHY WE DO IT**

Peaks are just coordinates like `chr1:100000-100500`. Without annotations,
we cannot know which gene they belong to. Gene annotations are required for:
- TSS enrichment scoring (Step 4)
- Gene activity scoring (Step 8)
- Coverage plots with gene track (Step 12)

**CRITICAL DETAIL — chromosome name matching**

The 10x data uses UCSC-style names: `chr1`, `chr2`, `chrX`...
The Ensembl database uses bare numbers: `1`, `2`, `X`...
Without converting, none of the peaks would match any genes.

```r
# Convert 1, 2, X → chr1, chr2, chrX
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations)         <- "hg38"
```

**WHAT TO CHANGE**

```r
# For hg19 data (older 10x datasets):
BiocManager::install("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# For mouse data (mm10):
ah <- AnnotationHub()
ensdb_mouse <- ah[["AH57770"]]   # Mus musculus GRCm38
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_mouse)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"
```

---

### Step 4: Quality Control — The 5 QC Metrics

This is the most important step before any analysis. We remove low-quality
cells using five metrics. Here is every metric explained from scratch.

---

#### QC Metric 1 — nCount_peaks (Total Fragment Count)

**WHAT IT IS**

The total number of DNA fragments detected in one cell.

**WHY IT MATTERS**

- Too few fragments (< 3,000–9,000) → likely an **empty droplet**
  (the sequencer captured a bead with no actual cell)
- Too many fragments (> 30,000–100,000) → likely a **doublet**
  (two cells accidentally captured together and counted as one)

**THE THRESHOLD (ATACv2 dataset)**

```r
nCount_peaks > 9000 & nCount_peaks < 100000
```

---

#### QC Metric 2 — TSS Enrichment Score

**WHAT IT IS**

A signal-to-noise ratio measuring how much the ATAC-seq signal is concentrated
at **Transcription Start Sites (TSS)** — the points on the genome where genes
begin to be read.

**THE BIOLOGY**

Chromatin near active gene promoters is almost always open — the cell needs
to be able to start transcribing. A good ATAC-seq experiment should therefore
show much higher fragment density right at TSS positions compared to the
surrounding background.

**THE FORMULA**

```
TSS Enrichment = (signal at TSS region) / (background flanking signal)

Where:
  signal at TSS    = average fragment count in ±50 bp window around TSS
  background       = average fragment count in flanking regions (±900–1000 bp)

Good cell:   TSS enrichment > 4    (signal is at least 4× background)
Bad cell:    TSS enrichment < 4    (signal barely above background noise)
```

**THE VISUAL**

```
Good cell (high TSS enrichment):
                          ████████████████
            ──────────────                ──────────────
                         -100bp  TSS  +100bp

Dead cell (low TSS enrichment):
            ─────────────────────────────────────────────
                         -100bp  TSS  +100bp   (flat — no signal)
```

**THE THRESHOLD**

```r
TSS.enrichment > 4      # ATACv2; use > 2 or > 3 for older protocols
```

---

#### QC Metric 3 — Nucleosome Signal

**WHAT IT IS**

The ratio of mono-nucleosomal fragments (~200 bp) to sub-nucleosomal
fragments (< 147 bp) in a single cell.

**THE BIOLOGY**

DNA in the nucleus is wound around protein spools called **nucleosomes**,
each wrapping ~147 base pairs of DNA. In a healthy ATAC-seq experiment, Tn5
cuts in three main places:
1. **Open chromatin** (< 147 bp fragments) — lots of these = good
2. **Between nucleosomes** (~200 bp mono-nucleosomal fragments)
3. **Between two nucleosomes** (~400 bp di-nucleosomal fragments)

The fragment size histogram of a healthy cell shows a clear banding pattern:

```
Fragment length histogram (healthy cell):
Count  ██████
       ███████
       ████████████
       ████████████████
       █████████████████  ←← sub-nucleosomal peak (~100 bp)
       ██████████
       ██████
       ████    ████████   ←← mono-nucleosomal peak (~200 bp)
       ██        ██
       █           █      ←← di-nucleosomal peak (~400 bp)
       ─────────────────────────────
       0   100  200  300  400  500  600 bp
```

A **high nucleosome signal** means the cell has too few short fragments
relative to nucleosomal ones — suggesting the Tn5 enzyme did not penetrate
the nucleus properly, or the cell is damaged.

**THE FORMULA**

```
Nucleosome Signal = (fragments between 147–294 bp) / (fragments < 147 bp)

Good cell:    < 4    (mostly sub-nucleosomal = lots of open chromatin)
Bad cell:     > 4    (too many nucleosomal fragments)
```

**THE THRESHOLD**

```r
nucleosome_signal < 4
```

---

#### QC Metric 4 — pct_reads_in_peaks (FRiP Score)

**WHAT IT IS**

The percentage of a cell's total fragments that overlap a called ATAC-seq
peak (a region identified by CellRanger as reproducibly accessible).

**WHY IT MATTERS**

In a high-quality cell, most of the Tn5 cuts happen in **real open chromatin
regions** — these become the called peaks. If many fragments fall *outside*
peaks, they represent background noise (random cuts in closed chromatin).

**THE FORMULA**

```
pct_reads_in_peaks = (fragments overlapping peaks / total fragments) × 100%

Good cell:    > 40%   (ATACv2; use > 15% for older v1 protocols)
Bad cell:     < 15%   (most fragments are noise)
```

**THE THRESHOLD**

```r
pct_reads_in_peaks > 40    # ATACv2; change to > 15 for ATACv1
```

---

#### QC Metric 5 — Blacklist Ratio

**WHAT IT IS**

The fraction of a cell's fragments that fall in **ENCODE blacklist regions** —
parts of the genome known to produce false, artefactual ATAC-seq signal.

**WHY THESE REGIONS EXIST**

Certain genomic regions (centromeres, telomeres, satellite repeats) are
physically unstable, have high copy numbers, or have sequence features
that cause sequencing artefacts. Reads mapping there are almost always
noise, not real chromatin accessibility.

**THE FORMULA**

```
Blacklist Ratio = (fragments in blacklist regions) / (fragments in peaks)

Good cell:    < 0.01  (< 1% of reads in artefact zones)
Bad cell:     > 0.05  (too many artefact reads)
```

**THE THRESHOLD**

```r
blacklist_ratio < 0.01    # strict; relax to < 0.05 for noisier datasets
```

---

#### Applying All 5 Filters Together

```r
pbmc <- subset(pbmc,
  subset = nCount_peaks         > 9000   &
           nCount_peaks         < 100000 &
           pct_reads_in_peaks   > 40     &
           blacklist_ratio      < 0.01   &
           nucleosome_signal    < 4      &
           TSS.enrichment       > 4
)
```

**WHAT TO CHANGE FOR YOUR DATA**

> 💡 Never use fixed thresholds blindly. Always look at the violin plots
> and density scatter first, then set thresholds at where the distributions
> show natural breaks.

| Your situation | What to change |
|---|---|
| Fewer than expected cells remain | Relax thresholds (lower min, raise max) |
| Old ATACv1 protocol | `pct_reads_in_peaks > 15`, `TSS.enrichment > 2` |
| Heart / muscle tissue | `nucleosome_signal < 6` (naturally higher) |
| Mouse data | Change `"^MT-"` to `"^mt-"` for mitochondrial genes |
| Very deep sequencing (> 100k reads/cell) | Raise `nCount_peaks < 200000` |

---

#### QC Plot 1 — Density Scatter

![QC Density Scatter](Rplot.png)

**How to read this plot:**
- Every dot = one cell
- x-axis = total fragments (log scale)
- y-axis = TSS enrichment score
- Red lines = 5th, 10th, 90th, 95th percentile thresholds
- The **dense yellow-orange cloud** in the upper right = your high-quality cells
- Dots far from the main cloud are poor-quality cells to remove

In this dataset, the main cloud sits between ~10,000–45,000 fragments and
TSS enrichment of 4.6–6.6. This is excellent ATACv2 quality.

---

#### QC Plot 2 — Fragment Size Histogram

![Fragment Size Histogram](Rplot01.png)

**How to read this plot:**
- Left panel: cells with nucleosome signal < 4 (good cells)
- Right panel: cells with nucleosome signal > 4 (poor cells)
- The left panel shows clear peaks at ~100 bp and ~200 bp — the classic
  nucleosomal banding pattern confirming high-quality chromatin accessibility

---

#### QC Plot 3 — Violin Plots (5 Metrics)

![QC Violins](Rplot02.png)

**How to read this plot:**
- Each violin shows the distribution of that metric across all cells
- Wider violin = more cells at that value
- Use these shapes to choose your threshold cutoffs:
  - `nCount_peaks`: cut off the thin tails at low and high ends
  - `blacklist_ratio`: should be very narrow and centred near 0
  - `pct_reads_in_peaks`: should peak above 40% for good ATACv2 data

---

### Step 5: TF-IDF Normalisation

**WHAT IT DOES**

Converts the raw fragment counts in the cells × peaks matrix into normalised
values that are comparable across cells with different sequencing depths.

**WHY WE DO NOT USE LOG-NORMALISATION HERE**

In scRNA-seq, we use `log(x + 1)` normalisation. ATAC-seq data is different:
- Values are almost always 0, 1, or 2 (very few 3+)
- The data is much more sparse (95–99% zeros vs ~90% in RNA-seq)
- Log-normalisation was designed for continuous count data — it does not
  work well on near-binary accessibility data

**THE TF-IDF FORMULA**

TF-IDF was originally invented for **text search engines** (searching documents
for words). In our case: cells = documents, peaks = words.

```
Step 1 — Term Frequency (TF): normalise each cell for sequencing depth
   TF(peak p, cell c) = fragments at peak p in cell c
                        ─────────────────────────────────────
                        total fragments in cell c

Step 2 — Inverse Document Frequency (IDF): upweight rare peaks
   IDF(peak p) = log( 1 + total_cells / cells_where_peak_p_is_open )

Step 3 — Multiply:
   TF-IDF(peak p, cell c) = TF × IDF
```

**In plain English:**
- A peak open in **many cells** (e.g. a housekeeping gene promoter) gets a
  **low IDF** → down-weighted (not informative for cell type differences)
- A peak open in only **a few cells** (e.g. a T cell-specific enhancer) gets
  a **high IDF** → up-weighted (tells you something specific about that cell)

**THE CODE**

```r
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")  # keep all peaks
```

**WHAT TO CHANGE**

```r
# Use only the top 25% most variable peaks (faster for large datasets):
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q75")
```

---

### Step 6: LSI (Latent Semantic Indexing)

**WHAT IT DOES**

Compresses the 165,000-peak TF-IDF matrix into ~30 meaningful summary
dimensions (components) that capture the main patterns of variation between cells.

**WHY WE CANNOT USE PCA**

Principal Component Analysis (PCA) is the standard approach in scRNA-seq,
but PCA assumes data is approximately normally distributed. ATAC-seq data
is near-binary and extremely sparse — PCA performs poorly on it.

**LSI is the ATAC-seq equivalent of PCA.** It applies **Singular Value
Decomposition (SVD)** to the TF-IDF matrix, which is mathematically better
suited for sparse binary data.

**THE CRITICAL DEPTH CORRELATION CHECK**

```r
DepthCor(pbmc)   # always run this before anything else
```

**LSI component 1 almost always captures sequencing depth, not biology.**

The reason: cells with more total fragments (deeper sequencing) systematically
show more accessible peaks everywhere — not because their chromatin is actually
more open, but because they were sequenced more. LSI component 1 picks this
up as the biggest pattern of variation.

```
Expected output from DepthCor():
  Component 1: correlation ≈ +0.88–0.98   ← depth artefact → SKIP
  Component 2: correlation ≈  0.0         ← real biology starts here
  Component 3: correlation ≈  0.05
  ...
```

**We always use `dims = 2:30` in all downstream steps, never `dims = 1:30`.**

---

#### QC Plot 4 — Depth Correlation

![Depth Correlation](Rplot03.png)

**How to read this plot:**
- x-axis = LSI component number
- y-axis = Pearson correlation with total fragment count per cell
- Component 1 (≈ +0.88) captures depth → **excluded from all downstream steps**
- Components 2–10 hover near zero → they capture **real biological variation**

This plot confirmed that our pipeline correctly excludes component 1 by
using `dims = 2:30`.

---

### Step 7: UMAP + Leiden Clustering

#### UMAP (Uniform Manifold Approximation and Projection)

**WHAT IT DOES**

Takes the 30 LSI components and projects them into a 2D scatter plot that
can be displayed on screen.

**WHY WE USE UMAP (NOT PCA OR t-SNE)**

- PCA is linear — it cannot show the curved structure of biological cell types
- t-SNE is slow and loses global structure (how far apart clusters are)
- UMAP preserves both local structure (cells that are similar are close)
  *and* some global structure (clusters that are biologically related are
  nearer to each other)

**THE KEY PRINCIPLE**

```
Cells close together on UMAP  →  similar chromatin accessibility  ✅
Distance BETWEEN clusters      →  not always biologically meaningful ⚠️
```

Always look at both the UMAP *and* the gene markers to interpret clusters.

**THE CODE**

```r
pbmc <- RunUMAP(pbmc, reduction = "lsi", dims = 2:30)
```

#### Leiden Clustering

**WHAT IT DOES**

Groups cells into clusters by building a **neighbourhood graph** (connecting
each cell to its most similar neighbours in LSI space) and then finding
dense communities in that graph.

**THE ALGORITHM**

```
Step 1: FindNeighbors()
        → For each cell, find its k nearest neighbours in LSI space
        → Build a graph where cells are nodes and similar cells are connected

Step 2: FindClusters() with algorithm = 3 (Leiden)
        → Find groups of cells that are more densely connected to each other
          than to the rest of the graph
        → The 'resolution' parameter controls how many groups are found
```

**LEIDEN vs LOUVAIN**

We use Leiden (`algorithm = 3`) rather than the default Louvain algorithm.
Leiden is an improvement that guarantees well-connected clusters (no
"disconnected" clusters where cells in the same group are not actually near
each other in the graph).

**THE RESOLUTION PARAMETER**

```
Low resolution (0.1–0.3)  → fewer, broader clusters  (coarse cell types)
Medium (0.5)              → good default starting point
High resolution (0.8–1.5) → many fine-grained clusters (rare subtypes)
```

**WHAT TO CHANGE**

```r
pbmc <- FindNeighbors(pbmc, reduction = "lsi", dims = 2:30,
                      k.param = 15)  # 👇 lower = tighter clusters

pbmc <- FindClusters(pbmc, algorithm = 3,
                     resolution = 0.5)  # 👇 change this to tune cluster count
```

---

#### UMAP Plot 1 — Leiden Clusters

![UMAP Leiden](Rplot04.png)

**How to read this plot:**
- Every dot = one cell
- Dots of the same colour = cells in the same Leiden cluster (0–14)
- The algorithm found **15 clusters** at resolution = 0.5
- Large central mass = T cells; isolated islands = monocytes, B cells, NK cells

---

## 5. Part 2 — Downstream Workflow: Methods & Results

---

### Step 8: Gene Activity Scoring

**WHAT IT DOES**

Converts the chromatin accessibility data into an *estimate* of gene
expression by summing fragment counts in each gene's regulatory region.

**THE FORMULA**

```
Gene Activity Score for gene G in cell C =
    Σ (all fragments in cell C that overlap gene G's body + 2kb upstream)
```

The 2 kb upstream region captures the **promoter** — the DNA sequence that
directly controls whether a gene is turned on or off. An open promoter usually
means the gene is being expressed.

**IMPORTANT CAVEAT**

```
Open chromatin ≠ active transcription
```

Gene activity is an *approximation*. A gene's promoter may be open without
the gene being actively expressed at that moment. Treat gene activity scores
as a rough guide for cluster interpretation — not as precise expression values.

**NORMALISATION**

After computing gene activity, we normalise using log-normalisation with a
scale factor equal to the median total gene activity per cell:

```
Normalised activity = log( (raw activity / total activity in cell) × scale_factor + 1 )
```

Using the **median** (rather than the standard 10,000) as the scale factor
is better here because gene activity values are naturally smaller in magnitude
than RNA expression counts.

**THE CODE**

```r
gene.activities <- GeneActivity(pbmc)
pbmc[["RNA"]]   <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(pbmc, assay = "RNA",
                      scale.factor = median(pbmc$nCount_RNA))
```

---

### Step 9: Marker Gene Visualisation

**WHAT IT DOES**

Projects the gene activity score for known **canonical marker genes** onto
the UMAP to help identify what biological cell type each cluster represents.

**HOW TO INTERPRET**

```
Dark blue region on UMAP  →  cells with HIGH activity at this gene locus
                               (chromatin is open near this gene)
Light/white region         →  cells with LOW activity
                               (chromatin is closed near this gene)
```

**KNOWN PBMC MARKER GENES**

| Gene | What it marks | Why it is a marker |
|------|--------------|-------------------|
| `MS4A1` (CD20) | B cells | Encodes a B cell surface protein; only expressed in B lineage |
| `CD3D` | All T cells | Part of the T cell receptor complex |
| `LEF1` | Naive T cells | Transcription factor active in resting T cells |
| `NKG7` | NK cells | Natural killer cell granule protein |
| `TREM1` | Monocytes | Inflammatory monocyte marker |
| `LYZ` | Monocytes | Lysozyme — expressed in myeloid cells |

**THE `max.cutoff = "q95"` PARAMETER**

A few cells may have very high activity (outliers). Without a cap, these
outliers pull the colour scale so high that most cells look white.
`max.cutoff = "q95"` caps the colour scale at the 95th percentile value,
so the colour variation is meaningful for the majority of cells.

---

#### UMAP Plot 2 — Marker Gene Activity

![Marker Genes](Rplot05.png)

**How to read this plot:**
- Six panels, one per marker gene
- Dark blue = high chromatin accessibility near that gene
- `MS4A1` (top left) is restricted to a small separate cluster → B cells
- `CD3D` and `LEF1` are co-active in the large upper clusters → T cells
- `LYZ` and `TREM1` mark the isolated left-side cluster → monocytes

---

### Step 10: Label Transfer from scRNA-seq

**WHAT IT DOES**

Borrows **known cell type labels** from a matched scRNA-seq experiment and
assigns them to ATAC-seq cells based on molecular similarity.

**WHY THIS IS NEEDED**

ATAC-seq data alone is hard to annotate. We can identify clusters, but
confirming that "cluster 3 is CD4 Naive T cells" requires knowing which
genes are expressed — and ATAC-seq only tells us which regions are
accessible, not which genes are actively transcribed.

**THE ALGORITHM — CCA (Canonical Correlation Analysis)**

```
Step 1: Take the scATAC-seq object and the scRNA-seq reference object

Step 2: FindTransferAnchors(reduction = "cca")
        → CCA finds a shared mathematical space where both datasets
          (ATAC and RNA) can be compared despite measuring different things
        → It identifies "anchor" pairs: one ATAC cell + one RNA cell that
          are highly similar in this shared space

Step 3: TransferData()
        → For each ATAC cell, look at its RNA anchors
        → Assign the most common cell type label among those anchors
        → Also compute a confidence score (prediction.score.max)

Step 4: Filter cells with low confidence (prediction.score.max < 0.5)
        → Only keep cells where the model is reasonably certain
```

**WHY CCA AND NOT RPCA?**

There are two main methods in Seurat for finding anchors:
- `rpca` (Reciprocal PCA) — fast, works well for **same modality** (RNA → RNA)
- `cca` (Canonical Correlation Analysis) — slower, works better for
  **cross-modality** (ATAC → RNA), because it finds correlations between
  two completely different types of measurements

**THE CODE**

```r
transfer.anchors <- FindTransferAnchors(reference = pbmc_rna,
                                         query = pbmc, reduction = "cca")
predicted.labels <- TransferData(anchorset = transfer.anchors,
                                  refdata = pbmc_rna$celltype,
                                  weight.reduction = pbmc[["lsi"]],
                                  dims = 2:30)
pbmc <- AddMetaData(pbmc, predicted.labels)
pbmc <- subset(pbmc, prediction.score.max > 0.5)  # high confidence only
```

---

#### Transfer Plot 1 — Before Transfer (Side-by-Side)

![Label Transfer](Rplot06.png)

**How to read this plot:**
- Left: ATAC-seq UMAP with raw Leiden cluster numbers
- Right: RNA-seq UMAP with known cell type labels
- The similar topology of the two UMAPs confirms that the two datasets share
  the same biological structure, validating that label transfer will work

---

#### Transfer Plot 2 — After Transfer (Predicted Cell Types)

![Predicted Labels](Rplot07.png)

**How to read this plot:**
- The ATAC-seq UMAP now has meaningful biological labels on every cluster
- **Cell types identified in this PBMC dataset:**

| Cell Type | Full Name | Key Biology |
|---|---|---|
| CD4 Naive | Naive CD4+ helper T cell | Resting, never activated |
| CD4 Memory | Memory CD4+ T cell | Previously activated; faster response |
| CD8 Naive | Naive CD8+ cytotoxic T cell | Resting killer T cell |
| CD8 effector | Effector CD8+ T cell | Actively killing infected cells |
| CD14+ Monocytes | Classical monocytes | Inflammatory, phagocytic |
| CD16+ Monocytes | Patrolling monocytes | Non-classical; tissue surveillance |
| NK bright | Bright NK cells | Cytotoxic natural killer cells |
| NK dim | Dim NK cells | Regulatory natural killer cells |
| B cell progenitor | Immature B cells | Early B lineage development |
| pre-B cell | Pre-B cells | B cell precursors |
| Dendritic cell | Conventional dendritic cells | Antigen presentation |
| pDC | Plasmacytoid dendritic cells | Interferon production |
| Double negative T cell | DN T cells | Unconventional T lineage |
| Platelet | Thrombocytes | Blood clotting |

---

### Step 11: Differential Accessibility Analysis

**WHAT IT DOES**

Finds **which peaks are significantly more open** in one cell type compared
to another — identifying the genomic regulatory elements that define each
cell type's identity.

**THE STATISTICAL TEST — Logistic Regression (LR)**

We use logistic regression rather than a simple t-test because:

1. ATAC-seq values are near-binary (0/1/2), not normally distributed
2. We need to control for a **confounding variable** — total fragment count
3. Logistic regression allows us to control for confounders via `latent.vars`

**THE MODEL**

```
For each peak P:

  Is peak P open in this cell? (0 or 1)  ~  Cell type + Total fragment count

  ↑ outcome variable                        ↑ predictor    ↑ confounder

If "Cell type" is a significant predictor after controlling for fragment count,
then peak P is differentially accessible.
```

**WHY WE CONTROL FOR TOTAL FRAGMENT COUNT**

Imagine Cell Type A has 20,000 total fragments and Cell Type B has 5,000.
Without controlling for depth, every peak will look more accessible in
Type A — not because it's truly more open, but because we sequenced it more.
`latent.vars = "nCount_peaks"` tells the model to account for this.

**FOLD CHANGE INTERPRETATION**

```
avg_log2FC > 0  →  peak more open in ident.1 (CD4 Naive)
avg_log2FC < 0  →  peak more open in ident.2 (CD14+ Monocytes)

|log2FC| = 1   means the peak is 2× more accessible in one group
|log2FC| = 2   means 4× more accessible
```

**THE CODE**

```r
da_peaks <- FindMarkers(pbmc,
  ident.1     = "CD4 Naive",
  ident.2     = "CD14+ Monocytes",
  test.use    = "LR",
  latent.vars = "nCount_peaks"
)
```

**WHAT TO CHANGE**

```r
# Compare any two cell types:
da_peaks <- FindMarkers(pbmc,
  ident.1 = "YOUR_CELL_TYPE_1",   # 👇 any cell type from the UMAP labels
  ident.2 = "YOUR_CELL_TYPE_2",
  test.use = "LR",
  latent.vars = "nCount_peaks"
)

# Find markers for ONE cell type vs ALL others:
da_peaks <- FindMarkers(pbmc, ident.1 = "CD4 Naive",
                         test.use = "LR", latent.vars = "nCount_peaks")
```

---

#### DA Plot — Violin + UMAP for Top DA Peak

![DA Peaks](Rplot08.png)

**How to read this plot:**
- **Left violin:** accessibility distribution of the top DA peak
  `chr12-119988511-119989430` in CD4 Naive (teal) vs CD14+ Monocytes (red)
  → The peak is clearly more accessible in CD4 Naive T cells
- **Right UMAP:** the same peak's accessibility projected across all cells
  → Blue (open) regions correspond exactly to the T cell clusters

---

### Step 12: Coverage Plots

**WHAT IT DOES**

Creates a genome-browser-style visualisation showing the density of DNA
fragments along a chromosomal region, one track per cell type.

**HOW TO READ COVERAGE PLOTS**

```
Track height at any position =
  number of fragment centres in that ~100 bp window
  ÷ total fragments in that cell type (normalised)

Tall peak  →  many fragments cutting here  →  accessible chromatin
Flat line  →  few fragments               →  closed chromatin
```

**THE ANATOMY OF A COVERAGE PLOT**

```
[Cell type 1 track] ─────▁▂████▂▁──────────────
[Cell type 2 track] ──────────────────────────── (closed)
[Cell type 3 track] ─────▁▂████▂▁──────────────
                                    ...
[Gene track]         ←──────── GENE_NAME ────────►   (gene model)
[Peaks track]                ▬ ▬▬  ▬         (called peaks)
                    chromosome position (bp)
```

**THE CODE**

```r
CoveragePlot(pbmc, region = "CD8A",    # gene name
             extend.upstream   = 3000,  # bp upstream of gene
             extend.downstream = 3000)  # bp downstream of gene
```

**WHAT TO CHANGE**

```r
# Show a wider region (more genomic context):
CoveragePlot(pbmc, region = "CD8A",
             extend.upstream = 20000, extend.downstream = 10000)

# Use coordinates instead of gene name:
CoveragePlot(pbmc, region = "chr2-86780000-86815000")

# Change which cell types are shown (set before calling):
Idents(pbmc) <- pbmc$predicted.id
# Then subset to only the cell types you want to see
sub <- subset(pbmc, idents = c("CD8 Naive", "CD8 effector", "CD4 Naive"))
CoveragePlot(sub, region = "CD8A")
```

---

#### Coverage Plot 1 — chr12 Region (BICDL1 gene)

![Coverage chr12 BICDL1](Rplot09.png)

**How to read this plot:**
- 14 cell type tracks shown
- The sharp accessibility peak at ~chr12:119,983,000 is specifically open
  in **T cell subtypes** (CD4 Naive, CD4 Memory, CD8 Naive, CD8 effector,
  Double negative T cell) but closed in monocytes and B cells
- Grey bars at the bottom = called ATAC-seq peaks
- Gene model below shows the BICDL1 gene structure (exons = thick boxes,
  introns = arrows, direction = transcription direction)

---

#### Coverage Plot 2 — CD8A Gene (chr2)

![Coverage CD8A](Rplot10.png)

**How to read this plot:**
- CD8A chromatin is highly accessible specifically in **CD8 Naive** and
  **CD8 effector** T cells (tall, sharp peaks)
- NK cells show moderate signal (NK cells weakly express CD8A)
- CD4 T cells, Monocytes, and B cells show minimal accessibility
- This perfectly matches known CD8A biology and validates the whole pipeline

---

## 6. Interactive Shiny App — User Guide

The Shiny app (`Shiney_App_For_EDA.R`) lets you explore any gene and
cell type combination without writing any code.

### How to Launch

```r
# Prerequisite: data/pbmc_atac_final.rds must exist (created by pipeline)
shiny::runApp("Shiney_App_For_EDA.R")
```

### Three-Step Workflow

```
Step 1 — Pick a Chromosome
         Dropdown shows chr1–chr22, chrX, chrY
         Selecting it automatically filters the gene list to only genes
         on that chromosome

Step 2 — Pick a Gene
         Type to search (e.g. type "CD8" to see CD8A, CD8B)
         Only genes present in the gene activity matrix are shown

Step 3 — Tick Cell Types (>= 2)
         Selected types: shown in colour on UMAP, included in violin/coverage
         Unselected types: shown in grey on UMAP, excluded from violin/coverage
         ALL cluster labels always visible on UMAP regardless of selection
```

### Five Output Panels

| Tab | What You See | How to Interpret |
|-----|-------------|-----------------|
| **Gene Activity** | Violin + box plot per cell type | Wider violin = more cells at that accessibility level; Wilcoxon brackets show statistical significance |
| **UMAP** | Full UMAP with cluster labels; selected types coloured by gene activity | Yellow label background = selected cluster; blue intensity = accessibility |
| **Coverage** | Genome browser tracks for selected cell types | Tall peaks = open chromatin at that gene locus |
| **Gene Info** | Table: chromosome, coordinates, gene size | Verify which genomic region you are analysing |
| **Summary Stats** | N cells, mean/median/SD per cell type | Coloured bar shows relative accessibility across cell types |

### Adjusting the Coverage Window

In the Coverage tab, two numeric inputs control how much genomic context
is shown around the gene:

```
Extend upstream (bp):   default 3,000  → raise to 20,000 to see enhancers
Extend downstream (bp): default 3,000  → raise to see downstream elements
```

---

## 7. How to Adapt This Pipeline for Your Own Data

### Different Organism

```r
# Mouse (Mus musculus, mm10):
ensdb_mouse <- ah[["AH57770"]]   # search with query(ah, "EnsDb.Mmusculus")
annotations <- GetGRangesFromEnsDb(ensdb = ensdb_mouse)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

# Mitochondrial genes — lowercase in mouse:
# No change needed; the pipeline doesn't filter mt genes for ATAC
```

### Different Genome Version

```r
# hg19 (GRCh37) — older 10x ATAC datasets:
BiocManager::install("EnsDb.Hsapiens.v75")
library(EnsDb.Hsapiens.v75)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "hg19"
```

### Different Number of Clusters

```r
# Too few clusters (major cell types merged):
pbmc <- FindClusters(pbmc, algorithm = 3, resolution = 0.8)  # increase

# Too many clusters (over-splitting):
pbmc <- FindClusters(pbmc, algorithm = 3, resolution = 0.3)  # decrease

# Try multiple resolutions at once and compare:
pbmc <- FindClusters(pbmc, algorithm = 3, resolution = c(0.2, 0.5, 0.8, 1.2))
```

### Different DA Comparison

```r
# Any two cell types from your predicted.id labels:
da_peaks <- FindMarkers(pbmc,
  ident.1     = "NK bright",     # 👈 change these two lines
  ident.2     = "NK dim",
  test.use    = "LR",
  latent.vars = "nCount_peaks",
  min.pct     = 0.05   # lower this if your cell types are small
)
```

### Changing the Shiny App Working Directory

```r
# At the top of Shiney_App_For_EDA.R, update this line:
setwd("C:/YOUR/PATH/TO/PROJECT")   # 👈 change to where your data/ folder is
```

---

## 8. Requirements & Data Downloads

### Install Packages (run once)

```r
if (!requireNamespace("remotes",     quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

remotes::install_github("stuart-lab/signac", ref = "develop")

install.packages(c("shiny", "shinydashboard", "shinycssloaders",
                   "ggplot2", "ggrepel", "dplyr", "ggsignif", "DT"))

BiocManager::install(c("Seurat", "Signac", "AnnotationHub",
                        "GenomicRanges", "GenomeInfoDb", "biovizBase"))
```

### Download Data Files (run in terminal)

```bash
# Create a project folder and download all required files
mkdir my_atac_project && cd my_atac_project

# Peak matrix
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5 \
     -O filtered_peak_bc_matrix.h5

# Cell metadata
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv \
     -O singlecell.csv

# Fragment file + index (BOTH required)
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz \
     -O fragments.tsv.gz
wget https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz.tbi \
     -O fragments.tsv.gz.tbi

# scRNA-seq reference for label transfer
wget https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds
```

### Run Order

```r
# 1. Set your working directory and run the pipeline
setwd("C:/YOUR/PATH/TO/my_atac_project")
source("ATACSeq_complete_pipeline.R")    # ~30–90 min

# 2. After pipeline finishes, launch the app
shiny::runApp("Shiney_App_For_EDA.R")
```

---

## 9. Common Errors & Fixes

| Error Message | Cause | Fix |
|---|---|---|
| `Cannot find data/pbmc_atac_final.rds` | Pipeline not completed | Run `ATACSeq_complete_pipeline.R` first |
| `No matching chromosomes found in fragment file` | Double chr prefix (chrchrX) in annotations | The pipeline auto-repairs this; if it persists, run `seqlevelsStyle(annotations) <- "UCSC"` |
| `All cells have same TSS.enrichment value` | Chromosome name mismatch between annotation and fragment file | Direct consequence of the chr prefix issue above |
| `No cells found after QC subset` | Thresholds too strict for your dataset | Use the dynamic threshold section of the corrected pipeline |
| `Error: 'pcaproject' not a valid reduction` | Old Seurat v4 code in Seurat v5 | Change to `reduction = "cca"` for ATAC→RNA label transfer |
| `Coverage plot: attempt to set colnames` | Signac region format wrong | Use `"chr1-1000-2000"` not `"chr1:1000-2000"` |
| `Error: biovizBase not found` | Missing dependency for GetGRangesFromEnsDb | `BiocManager::install("biovizBase"); library(biovizBase)` |
| `App: gene dropdown empty` | No genes on selected chromosome in lookup | Try chr1 or chr2 which have the most genes |
| `Coverage tab blank` | Fragment file moved from original path | Place `fragments.tsv.gz` and `.tbi` in the same directory as your `.R` scripts |

---

## 10. Further Reading & References

### Papers You Should Read

- **Stuart et al. 2021 — Signac** (the tool we use):
  [Nature Methods](https://www.nature.com/articles/s41592-021-01282-5)

- **Cusanovich et al. 2018 — TF-IDF + LSI for scATAC-seq** (the core method):
  [Science](https://www.science.org/doi/10.1126/science.aab1601)

- **Stuart et al. 2019 — Label Transfer / Integration** (CCA method):
  [Cell](https://www.cell.com/cell/fulltext/S0092-8674(19)30559-8)

- **Heumos et al. 2023 — Best Practices for scATAC-seq** (comprehensive review):
  [Nature Reviews Genetics](https://www.nature.com/articles/s41576-023-00586-w)

### Documentation & Tutorials

- [Signac official PBMC vignette](https://stuartlab.org/signac/articles/pbmc_vignette) — the tutorial this pipeline is based on
- [Signac function reference](https://stuartlab.org/signac/reference/) — every function explained
- [ENCODE ATAC-seq standards](https://www.encodeproject.org/atac-seq/) — QC thresholds from the consortium
- [10x Genomics PBMC ATACv2 dataset](https://www.10xgenomics.com/datasets) — the data used here

### Video Resources

- [Bioinformagician YouTube](https://www.youtube.com/@Bioinformagician) — code tutorial for this pipeline
- [StatQuest — UMAP Explained](https://www.youtube.com/watch?v=eN0wFzBA4Sc) — plain English explanation of UMAP
- [StatQuest — PCA / SVD](https://www.youtube.com/watch?v=FgakZw6K1QQ) — background for understanding LSI

---

## Acknowledgements

- **[Bioinformagician (YouTube)](https://www.youtube.com/@Bioinformagician)** —
  ATAC seq explanation and basic concept.
- **[Signac / Stuart Lab](https://stuartlab.org/signac/)** — core toolkit and methods
- **[10x Genomics](https://www.10xgenomics.com)** — PBMC ATACv2 dataset
- **[Ensembl / EMBL-EBI](https://www.ensembl.org)** — hg38 gene annotations (v98)
- **[ENCODE Project](https://www.encodeproject.org)** — QC standards and blacklist regions
- **[Bioconductor](https://bioconductor.org)** — AnnotationHub, GenomicRanges, biovizBase

---

*Pipeline: Signac v1.17+ · Seurat v5 · R 4.6 · Genome: hg38 (GRCh38) · Ensembl v98*
*Dataset: 10x Genomics PBMC 10k ATACv2 Chromium Controller*
