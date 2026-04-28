# ============================================================
# scATAC-seq Complete Pipeline — Human PBMCs (10k cells)
# Dataset  : 10x Genomics PBMC 10k ATACv2 (Chromium Controller)
# Genome   : hg38 (GRCh38 — Ensembl v98 via AnnotationHub)
# Tools    : Signac v1.17+ + Seurat v5 + R 4.6
#


# ── Install packages (run once, then comment out) ───────────
# if (!requireNamespace("remotes",     quietly = TRUE)) install.packages("remotes")
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# remotes::install_github("stuart-lab/signac", ref = "develop")
# install.packages(c("Matrix", "irlba", "ggplot2", "patchwork", "httr"))
# BiocManager::install(c("AnnotationHub", "GenomicRanges", "GenomeInfoDb",
#                         "biovizBase"))   # FIX 1: biovizBase required by
#                                          # GetGRangesFromEnsDb() — was missing


# ── Load libraries ───────────────────────────────────────────
library(Signac)
library(Seurat)
library(GenomicRanges)
library(GenomeInfoDb)
library(AnnotationHub)
library(biovizBase)      
library(patchwork)
library(ggplot2)
# ============================================================
# PART 1 — STANDARD WORKFLOW
# ============================================================


# ── Set working directory ────────────────────────────────────
setwd("C:/ZAIN/CODES/R/Projects")

# Create output directory for saved objects
if (!dir.exists("data")) dir.create("data")


# ── 1. Load Data ─────────────────────────────────────────────
# Files should already be downloaded. If not, run the download block below.
# ── Optional: download files (skip if already done) ──────────
# library(httr)
# my_ua <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
# GET("https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5",
#     write_disk("filtered_peak_bc_matrix.h5", overwrite=TRUE), user_agent(my_ua), timeout(600))
# GET("https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv",
#     write_disk("singlecell.csv", overwrite=TRUE), user_agent(my_ua))
# GET("https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz",
#     write_disk("fragments.tsv.gz", overwrite=TRUE), user_agent(my_ua), timeout(30000))
# GET("https://cf.10xgenomics.com/samples/cell-atac/2.1.0/10k_pbmc_ATACv2_nextgem_Chromium_Controller/10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz.tbi",
#     write_disk("fragments.tsv.gz.tbi", overwrite=TRUE), user_agent(my_ua))
# GET("https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds",
#     write_disk("pbmc_10k_v3.rds", overwrite=TRUE), user_agent(my_ua), timeout(300))

counts <- Read10X_h5(filename = "filtered_peak_bc_matrix.h5")

metadata <- read.csv(
  file      = "singlecell.csv",
  header    = TRUE,
  row.names = 1
)
head(metadata)

chrom_assay <- CreateChromatinAssay(
  counts       = counts,
  sep          = c(":", "-"),
  fragments    = "fragments.tsv.gz",
  min.cells    = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts    = chrom_assay,
  assay     = "peaks",
  meta.data = metadata
)

pbmc[["peaks"]]
pbmc


# ── 2. Filter Non-Standard Chromosomes ───────────────────────
peaks.keep <- seqnames(granges(pbmc)) %in% standardChromosomes(granges(pbmc))
pbmc       <- pbmc[as.vector(peaks.keep), ]
cat("Features after chromosome filter:", nrow(pbmc), "\n")


# ── 3. Add Gene Annotations (hg38 — Ensembl v98) ─────────────
ah <- AnnotationHub()
query(ah, "EnsDb.Hsapiens.v98")
ensdb_v98 <- ah[["AH75011"]]

annotations <- GetGRangesFromEnsDb(ensdb = ensdb_v98)

# EnsDb already stores some names differently in newer Bioconductor versions.
# Solution: use seqlevelsStyle() which handles this automatically and safely.
# If seqlevelsStyle() fails (some systems), the manual fallback below is used.
tryCatch({
  seqlevelsStyle(annotations) <- "UCSC"
  cat("seqlevelsStyle() succeeded — chromosome names set to UCSC style\n")
}, error = function(e) {
  # Fallback: strip any existing "chr" prefix first, then add it cleanly
  current_levels <- seqlevels(annotations)
  # Remove any "chr" prefix that may already be there (avoids chrchr double)
  stripped <- sub("^chr", "", current_levels)
  seqlevels(annotations) <<- paste0("chr", stripped)
  cat("Manual chr prefix applied\n")
})
genome(annotations) <- "hg38"

# Verify no double prefix — should show chr1, chr2 ... chrX, chrY, chrMT
cat("First 5 seqlevels:", head(seqlevels(annotations), 5), "\n")
cat("Last 5 seqlevels:",  tail(seqlevels(annotations), 5), "\n")

# Attach to the Seurat object
Annotation(pbmc) <- annotations
Annotation(pbmc)


# ── 4. Quality Control ───────────────────────────────────────

# FIX 3: NucleosomeSignal() and TSSEnrichment() are deprecated in Signac
# v1.17 — both now warn "Use 'ATACqc' instead."
# However they still work and produce correct output. We keep them for
# compatibility and suppress the deprecation warning so it doesn't alarm you.
# The ATACqc() replacement API is shown in comments below for future reference.

suppressWarnings(
  pbmc <- NucleosomeSignal(object = pbmc)
)
# Future ATACqc equivalent:
# pbmc <- ATACqc(pbmc, metric = "nucleosome_signal")

suppressWarnings(
  pbmc <- TSSEnrichment(object = pbmc)
)
# Future ATACqc equivalent:
# pbmc <- ATACqc(pbmc, metric = "TSS")

# Fraction of reads in peaks
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

# FIX 4: FractionCountsInRegion() is deprecated in Signac v1.17.
# The error log shows it still works but with a warning.
# Fall back to the simple ratio from CellRanger metadata — more reliable
# and does not require downloading additional AnnotationHub resources.
# The CellRanger metadata already has blacklist_region_fragments pre-computed.
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments /
                        pbmc$peak_region_fragments

# Guard against NaN from division by zero (cells with 0 peak fragments)
pbmc$blacklist_ratio[is.nan(pbmc$blacklist_ratio)] <- 0

cat("QC metrics computed. Column names:\n")
print(colnames(pbmc@meta.data))

# ....Diagnostic: check if TSS enrichment has variation ......
# FIX 5: The error log shows VlnPlot warns "All cells have the same value
# of TSS.enrichment." This means TSSEnrichment() ran but produced a flat
# value — usually because the chromosome names in annotations did not match
# the fragment file at run time. If you see this, check that Annotation(pbmc)
# shows chrX/chrY style names (not chrchrX). FIX 2 above addresses this.
tss_range <- range(pbmc$TSS.enrichment, na.rm = TRUE)
cat("TSS.enrichment range: min =", tss_range[1], "max =", tss_range[2], "\n")

if (tss_range[1] == tss_range[2]) {
  cat("WARNING: TSS enrichment has no variation — annotation chromosome names\n")
  cat("         likely did not match fragment file. Check seqlevels(Annotation(pbmc)).\n")
  cat("         Proceeding without TSS filter — use other QC metrics only.\n")
}


# ....Visualise QC ............................................
DensityScatter(pbmc, x = "nCount_peaks", y = "TSS.enrichment",
               log_x = TRUE, quantiles = TRUE)

pbmc$nucleosome_group <- ifelse(pbmc$nucleosome_signal > 4, "NS > 4", "NS < 4")
FragmentHistogram(object = pbmc, group.by = "nucleosome_group")

VlnPlot(
  object   = pbmc,
  features = c("nCount_peaks", "TSS.enrichment", "blacklist_ratio",
               "nucleosome_signal", "pct_reads_in_peaks"),
  pt.size  = 0.1,
  ncol     = 5
)


# ....Apply QC filter .........................................
# FIX 6: "No cells found" error occurred because strict ATACv2 thresholds
# were applied without first checking actual data distributions. The fix is
# to compute thresholds dynamically from the data's own percentiles, which
# ensures the filter always retains a sensible proportion of cells.
#
# Print the actual quantiles first so you can see your data
cat("\n=== QC Metric Distributions ===\n")
cat("nCount_peaks quantiles (1%, 5%, 95%, 99%):\n")
print(quantile(pbmc$nCount_peaks, c(0.01, 0.05, 0.95, 0.99), na.rm = TRUE))
cat("pct_reads_in_peaks quantiles:\n")
print(quantile(pbmc$pct_reads_in_peaks, c(0.01, 0.05, 0.95), na.rm = TRUE))
cat("nucleosome_signal quantiles:\n")
print(quantile(pbmc$nucleosome_signal, c(0.95, 0.99), na.rm = TRUE))
cat("TSS.enrichment quantiles:\n")
print(quantile(pbmc$TSS.enrichment, c(0.01, 0.05, 0.1), na.rm = TRUE))
cat("blacklist_ratio quantiles:\n")
print(quantile(pbmc$blacklist_ratio, c(0.95, 0.99), na.rm = TRUE))

# Dynamic thresholds: use percentile-based cutoffs that adapt to this dataset
# These replace the fixed ATACv2 vignette thresholds that caused "No cells found"
ncount_min    <- max(200,   quantile(pbmc$nCount_peaks,         0.02,  na.rm = TRUE))
ncount_max    <- min(1e6,   quantile(pbmc$nCount_peaks,         0.98,  na.rm = TRUE))
frip_min      <- max(10,    quantile(pbmc$pct_reads_in_peaks,   0.05,  na.rm = TRUE))
ns_max        <- min(10,    quantile(pbmc$nucleosome_signal,    0.95,  na.rm = TRUE))
bl_max        <- min(0.10,  quantile(pbmc$blacklist_ratio,      0.95,  na.rm = TRUE))

# TSS enrichment threshold: only apply if the metric has real variation
tss_has_variation <- (tss_range[2] - tss_range[1]) > 0.1
tss_min <- if (tss_has_variation) max(1, quantile(pbmc$TSS.enrichment, 0.05, na.rm=TRUE)) else 0

cat("\nApplying QC thresholds:\n")
cat(sprintf("  nCount_peaks:        %d – %d\n",   round(ncount_min), round(ncount_max)))
cat(sprintf("  pct_reads_in_peaks:  > %.1f%%\n",  frip_min))
cat(sprintf("  nucleosome_signal:   < %.2f\n",    ns_max))
cat(sprintf("  blacklist_ratio:     < %.4f\n",    bl_max))
cat(sprintf("  TSS.enrichment:      > %.2f  (applied: %s)\n", tss_min, tss_has_variation))

before <- ncol(pbmc)

if (tss_has_variation) {
  pbmc <- subset(
    x      = pbmc,
    subset = nCount_peaks       > ncount_min  &
             nCount_peaks       < ncount_max  &
             pct_reads_in_peaks > frip_min    &
             blacklist_ratio    < bl_max      &
             nucleosome_signal  < ns_max      &
             TSS.enrichment     > tss_min
  )
} else {
  # TSS enrichment flat — skip that filter
  pbmc <- subset(
    x      = pbmc,
    subset = nCount_peaks       > ncount_min  &
             nCount_peaks       < ncount_max  &
             pct_reads_in_peaks > frip_min    &
             blacklist_ratio    < bl_max      &
             nucleosome_signal  < ns_max
  )
}

after <- ncol(pbmc)
cat(sprintf("Cells: %d → %d  (%.1f%% retained)\n",
            before, after, after / before * 100))


# ── 5. Normalisation & Dimensionality Reduction (LSI) ────────
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = "q0")
pbmc <- RunSVD(pbmc)
DepthCor(pbmc)    # LSI1 should be ~1.0 correlation — always skip it below


# ── 6. UMAP & Leiden Clustering ──────────────────────────────
pbmc <- RunUMAP(
  object    = pbmc,
  reduction = "lsi",
  dims      = 2:30
)

pbmc <- FindNeighbors(
  object    = pbmc,
  reduction = "lsi",
  dims      = 2:30
)

pbmc <- FindClusters(
  object     = pbmc,
  algorithm  = 3,
  resolution = 0.5
)

DimPlot(object = pbmc, label = TRUE) + NoLegend()

# FIX 9: create data/ directory before saving
if (!dir.exists("data")) dir.create("data")
saveRDS(pbmc, file = "data/pbmc_atac_part1.rds")
cat("Part 1 saved → data/pbmc_atac_part1.rds\n")


# ============================================================
# PART 2 — DOWNSTREAM WORKFLOW
# ============================================================

# Reload Part 1 if starting here fresh:
pbmc <- readRDS("data/pbmc_atac_part1.rds")

pbmc
DimPlot(object = pbmc, label = TRUE) + NoLegend()


# ── 7. Gene Activity Scoring ─────────────────────────────────
# FIX 7: GeneActivity() failed with "No matching chromosomes found in fragment
# file." This happens when the annotations attached to the object have double
# chr prefix (chrchrX) from FIX 2 not yet applied, OR when the annotation
# object stored inside pbmc still has wrong chromosome names.
#
# We verify and repair the chromosome names inside the annotation before
# calling GeneActivity().

anno_check <- Annotation(pbmc)
current_seqlevels <- seqlevels(anno_check)
has_double_chr    <- any(grepl("^chrchr", current_seqlevels))

if (has_double_chr) {
  cat("WARNING: Double chr prefix detected in annotations — fixing now...\n")
  fixed_levels <- sub("^chrchr", "chr", current_seqlevels)
  seqlevels(anno_check) <- fixed_levels
  Annotation(pbmc)      <- anno_check
  cat("Fixed. First 5 seqlevels:", head(seqlevels(Annotation(pbmc)), 5), "\n")
} else {
  cat("Annotation chromosome names OK:", head(current_seqlevels, 5), "\n")
}

gene.activities <- GeneActivity(pbmc)
cat("Gene activity matrix dimensions:", dim(gene.activities), "\n")
gene.activities[1:5, 1:5]

pbmc[["RNA"]] <- CreateAssayObject(counts = gene.activities)

scale_factor <- median(pbmc$nCount_RNA)
cat("Gene activity normalisation scale factor:", round(scale_factor), "\n")

pbmc <- NormalizeData(
  object               = pbmc,
  assay                = "RNA",
  normalization.method = "LogNormalize",
  scale.factor         = scale_factor
)

cat("Assays in object:\n")
print(pbmc@assays)


# ── 8. Marker Gene Visualisation ─────────────────────────────
DefaultAssay(pbmc) <- "RNA"

# Check which markers are actually present in the gene activity matrix
desired_markers <- c("MS4A1", "CD3D", "LEF1", "NKG7", "TREM1", "LYZ")
available_markers <- desired_markers[desired_markers %in% rownames(pbmc)]
cat("Marker genes found:", paste(available_markers, collapse = ", "), "\n")

if (length(available_markers) > 0) {
  FeaturePlot(
    object     = pbmc,
    features   = available_markers,
    pt.size    = 0.1,
    max.cutoff = "q95",
    ncol       = min(3, length(available_markers))
  )
} else {
  cat("No marker genes found. Check that gene activity scoring ran correctly.\n")
}


# ── 9. Label Transfer from scRNA-seq ─────────────────────────
# FIX 11: the RNA reference was downloaded to the working directory, not data/
# The error log shows it as "pbmc_10k_v3.rds" in the working dir listing.
rna_path <- if (file.exists("pbmc_10k_v3.rds")) {
  "pbmc_10k_v3.rds"
} else if (file.exists("data/pbmc_10k_v3.rds")) {
  "data/pbmc_10k_v3.rds"
} else {
  stop("Cannot find pbmc_10k_v3.rds. Download it with:\n",
       "  download.file('https://signac-objects.s3.amazonaws.com/pbmc_10k_v3.rds', 'pbmc_10k_v3.rds')")
}
cat("Loading RNA reference from:", rna_path, "\n")

pbmc_rna <- readRDS(rna_path)
pbmc_rna <- UpdateSeuratObject(pbmc_rna)

head(pbmc_rna@meta.data)

p1 <- DimPlot(pbmc,     reduction = "umap") +
      NoLegend() + ggtitle("scATAC-seq")
p2 <- DimPlot(pbmc_rna, reduction = "umap", group.by = "celltype",
              label = TRUE, repel = TRUE) +
      NoLegend() + ggtitle("scRNA-seq")
p1 | p2

# CCA is correct for cross-modality (ATAC → RNA) label transfer
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query     = pbmc,
  reduction = "cca"
)

predicted.labels <- TransferData(
  anchorset        = transfer.anchors,
  refdata          = pbmc_rna$celltype,
  weight.reduction = pbmc[["lsi"]],
  dims             = 2:30
)
head(predicted.labels)

pbmc <- AddMetaData(object = pbmc, metadata = predicted.labels)

plot1 <- DimPlot(pbmc,     reduction = "umap", group.by = "predicted.id",
                 label = TRUE, repel = TRUE) +
         NoLegend() + ggtitle("scATAC-seq — Predicted")
plot2 <- DimPlot(pbmc_rna, reduction = "umap", group.by = "celltype",
                 label = TRUE, repel = TRUE) +
         NoLegend() + ggtitle("scRNA-seq — True")
plot1 | plot2

# Keep only high-confidence predictions (score > 0.5)
pbmc <- subset(pbmc, subset = prediction.score.max > 0.5)
cat("Cells after confidence filter:", ncol(pbmc), "\n")


# ── 10. Differential Accessibility Analysis ──────────────────
DefaultAssay(pbmc) <- "peaks"
Idents(pbmc) <- factor(pbmc$predicted.id,
                       levels = sort(unique(pbmc$predicted.id)))

# Show available identities so you can confirm the group names
cat("Available cell type identities:\n")
print(table(Idents(pbmc)))

# Pick the two most common cell types if CD4 Naive / CD14+ Mono not present
all_ids    <- as.character(unique(Idents(pbmc)))
ident1_use <- if ("CD4 Naive"        %in% all_ids) "CD4 Naive"        else all_ids[1]
ident2_use <- if ("CD14+ Monocytes"  %in% all_ids) "CD14+ Monocytes"  else all_ids[2]

cat("Comparing:", ident1_use, "vs", ident2_use, "\n")

da_peaks <- FindMarkers(
  object      = pbmc,
  ident.1     = ident1_use,
  ident.2     = ident2_use,
  test.use    = "LR",
  latent.vars = "nCount_peaks"
)
head(da_peaks)
saveRDS(da_peaks, file = "da_peaks_results.rds", compress = "gzip")
fc <- FoldChange(pbmc, ident.1 = ident1_use, ident.2 = ident2_use)
fc <- fc[order(fc$avg_log2FC, decreasing = TRUE), ]
head(fc)

da_plot1 <- VlnPlot(pbmc, features = rownames(da_peaks)[1], pt.size = 0.1,
                    idents = c(ident1_use, ident2_use))
da_plot2 <- FeaturePlot(pbmc, features = rownames(da_peaks)[1], pt.size = 0.1)
da_plot1 | da_plot2


# ── 11. Coverage Plots ───────────────────────────────────────
# FIX 8: CoveragePlot(region = rownames(da_peaks)[1]) can fail with
# "attempt to set colnames on object with less than two dimensions" when
# the region string is not properly formatted or when the gene annotation
# chromosome names don't match the fragment file.
# Safe approach: use a known valid gene name string instead.

# Try top DA peak first; fall back to a known gene if it fails
tryCatch({
  CoveragePlot(
    object            = pbmc,
    region            = rownames(da_peaks)[1],
    extend.upstream   = 40000,
    extend.downstream = 20000
  )
}, error = function(e) {
  cat("CoveragePlot with peak coordinates failed:", conditionMessage(e), "\n")
  cat("Trying with gene name instead...\n")
})

# Coverage plot using a gene name (more robust than raw coordinates)
# CD8A is a canonical CD8+ T cell marker that should be present
tryCatch({
  CoveragePlot(pbmc, region = "CD8A")
}, error = function(e) {
  cat("CoveragePlot with CD8A failed:", conditionMessage(e), "\n")
  cat("This usually means annotation chr names still mismatch fragment file.\n")
  cat("Check: head(seqlevels(Annotation(pbmc)))\n")
  cat("It should show chr1, chr2 etc — not chrchr1.\n")
})


# ── Save Final Object ────────────────────────────────────────
if (!dir.exists("data")) dir.create("data")

saveRDS(pbmc, file = "data/pbmc_atac_final.rds")
cat("Final object saved → data/pbmc_atac_final.rds\n")

# Export cluster labels as CSV
# Build safely — only include columns that exist
cluster_cols <- list(
  barcode        = colnames(pbmc),
  leiden_cluster = as.character(pbmc$seurat_clusters)
)
if ("predicted.id" %in% colnames(pbmc@meta.data)) {
  cluster_cols$predicted_id     <- as.character(pbmc$predicted.id)
  cluster_cols$prediction_score <- pbmc$prediction.score.max
}
cluster_df <- as.data.frame(cluster_cols, row.names = colnames(pbmc))
write.csv(cluster_df, "pbmc_atac_clusters.csv", row.names = TRUE)
cat("Cluster table saved → pbmc_atac_clusters.csv\n")
