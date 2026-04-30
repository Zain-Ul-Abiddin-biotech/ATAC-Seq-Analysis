# ============================================================
# Shiny EDA App — scATAC‑seq Chromatin Accessibility Explorer
# Dataset  : 10x Genomics PBMC 10k ATACv2
# Requires : data/pbmc_atac_final.rds (from ATACSeq_complete_pipeline.R)
# ============================================================
# Workflow  : Select Chromosome → Gene → Cell Types (≥2) → Plot
# Outputs   : Violin + box, UMAP feature, Coverage plot,
#             Gene info table, DA summary table
# ============================================================

# ── Libraries ────────────────────────────────────────────────
library(shiny)
library(shinydashboard)
library(Signac)
library(Seurat)
library(GenomicRanges)
library(ggplot2)
library(dplyr)
library(ggsignif)
library(patchwork)
library(DT)
library(shinycssloaders)   

# ── 1. Load the saved Seurat object ──────────────────────────
setwd("C:/ZAIN/CODES/R")
message("Loading pbmc_atac_final.rds …")
pbmc <- readRDS("data/pbmc_atac_final.rds")
message("Object loaded: ", ncol(pbmc), " cells")

# ── 1b. Attach fragment file if not already stored ───────────
# CHANGE THIS PATH TO MATCH YOUR REAL FILE NAME
fragment_path <- "C:/ZAIN/CODES/R/pbmc_atac_fragments.tsv.gz"

if (is.null(Fragments(pbmc)) || length(Fragments(pbmc)) == 0) {
  if (file.exists(fragment_path)) {
    pbmc <- SetFragments(pbmc, fragments = fragment_path)
    message("Fragment file attached from: ", fragment_path)
  } else {
    warning("Fragment file NOT FOUND at: ", fragment_path,
            "\nCoverage plot will not work. Update 'fragment_path' to the correct file.")
  }
}

# ── 2. Repair chromosome names if double‑prefix exists ───────
anno <- Annotation(pbmc[["peaks"]])
if (!is.null(anno)) {
  lvls <- seqlevels(anno)
  if (any(grepl("^chrchr", lvls))) {
    seqlevels(anno) <- sub("^chrchr", "chr", lvls)
    Annotation(pbmc[["peaks"]]) <- anno
    message("Double chr prefix repaired in stored annotations.")
  }
}

# ── 3. Determine cell type column ─────────────────────────────
cell_type_col <- if ("predicted.id" %in% colnames(pbmc@meta.data)) {
  "predicted.id"
} else {
  "seurat_clusters"
}
all_cell_types <- sort(unique(as.character(pbmc[[cell_type_col]][, 1])))
message("Cell types found: ", paste(all_cell_types, collapse = ", "))

# ── 4. Ensure RNA (gene activity) assay exists ───────────────
if (!"RNA" %in% Assays(pbmc)) {
  message("RNA assay missing — computing gene activity …")
  DefaultAssay(pbmc) <- "peaks"
  ga <- GeneActivity(pbmc)
  pbmc[["RNA"]] <- CreateAssayObject(counts = ga)
  pbmc <- NormalizeData(pbmc, assay = "RNA",
                        normalization.method = "LogNormalize",
                        scale.factor = median(pbmc$nCount_peaks))
}

# ── 5. Build chromosome → gene lookup from annotation ─────────
gene_annot <- Annotation(pbmc[["peaks"]])

if (is.null(gene_annot) || length(gene_annot) == 0)
  stop("No gene annotation found in the peaks assay. ",
       "Re-run ATACSeq_complete_pipeline.R.")

mcol_names <- colnames(mcols(gene_annot))
gene_col <- if ("gene_name" %in% mcol_names) "gene_name" else
  if ("symbol"    %in% mcol_names) "symbol"    else
    stop("Cannot find gene_name / symbol column in annotation mcols.")

gene_df <- data.frame(
  gene   = as.character(mcols(gene_annot)[[gene_col]]),
  chr    = as.character(seqnames(gene_annot)),
  start  = start(gene_annot),
  end    = end(gene_annot),
  strand = as.character(strand(gene_annot)),
  stringsAsFactors = FALSE
) |>
  dplyr::filter(!is.na(gene), gene != "") |>
  dplyr::filter(gene %in% rownames(pbmc[["RNA"]])) |>
  dplyr::filter(grepl("^chr[0-9XYM]", chr)) |>
  dplyr::distinct(gene, chr, .keep_all = TRUE) |>
  dplyr::arrange(chr, start)

chromosomes <- sort(unique(gene_df$chr))

message("Chromosomes in lookup: ", paste(chromosomes, collapse = ", "))
message("Total genes in lookup: ", nrow(gene_df))

if (nrow(gene_df) == 0)
  stop("gene_df is empty — no annotated genes found in the RNA assay.")

# ── 5b. Pre‑compute UMAP cluster centroids for labels (ALL clusters) ──
umap_mat <- as.data.frame(Embeddings(pbmc, "umap"))
colnames(umap_mat) <- c("UMAP_1", "UMAP_2")
umap_mat$celltype  <- pbmc[[cell_type_col]][, 1]

cluster_centroids <- umap_mat |>
  dplyr::group_by(celltype) |>
  dplyr::summarise(
    cx = median(UMAP_1),
    cy = median(UMAP_2),
    .groups = "drop"
  )

message("Startup complete.")

# ── 6. UI ─────────────────────────────────────────────────────
ui <- dashboardPage(
  skin = "blue",
  
  dashboardHeader(title = "scATAC‑seq Explorer"),
  
  dashboardSidebar(
    width = 290,
    tags$div(
      style = "padding: 12px 15px 0 15px;",
      
      # Step 1 — Chromosome
      tags$h5(tags$b("Step 1 — Chromosome"),
              style = "color:#ecf0f1; margin-bottom:4px;"),
      selectInput("chrom", NULL,
                  choices  = chromosomes,
                  selected = chromosomes[1],
                  width    = "100%"),
      
      tags$hr(style = "border-color:#4a6278; margin:10px 0;"),
      
      # Step 2 — Gene
      tags$h5(tags$b("Step 2 — Gene on selected chromosome"),
              style = "color:#ecf0f1; margin-bottom:4px;"),
      selectizeInput("gene", NULL,
                     choices  = NULL,
                     options  = list(placeholder = "Choose chromosome first …"),
                     width    = "100%"),
      
      tags$hr(style = "border-color:#4a6278; margin:10px 0;"),
      
      # Step 3 — Cell types
      tags$h5(tags$b("Step 3 — Cell types (select ≥ 2)"),
              style = "color:#ecf0f1; margin-bottom:4px;"),
      tags$small("Only cells belonging to ticked types are shown.",
                 style = "color:#bdc3c7;"),
      tags$br(), tags$br(),
      checkboxGroupInput("celltypes", NULL,
                         choices  = all_cell_types,
                         selected = all_cell_types[seq_len(min(2, length(all_cell_types)))]),
      tags$br(),
      
      # Buttons
      actionButton("plotBtn",   "▶  Update Plots",
                   class = "btn-primary btn-block",
                   style = "font-weight:bold;"),
      tags$br(), tags$br(),
      actionButton("selectAll",  "Select All Types",  class = "btn-default btn-block btn-sm"),
      actionButton("clearAll",   "Clear All Types",   class = "btn-default btn-block btn-sm"),
      tags$br()
    )
  ),
  
  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper, .right-side { background-color: #f4f6f9; }
      .nav-tabs-custom .nav-tabs li.active a { color: #2c3e50; font-weight: bold; }
      .info-box { min-height: 75px; }
      .info-box-number { font-size: 18px; }
    "))),
    
    # Info boxes at the top
    fluidRow(
      infoBoxOutput("infoGene",      width = 4),
      infoBoxOutput("infoChrom",     width = 4),
      infoBoxOutput("infoCellCount", width = 4)
    ),
    
    # Validation message
    fluidRow(
      column(12,
             uiOutput("validationMsg")
      )
    ),
    
    # Main tab panels
    tabBox(
      width = 12,
      title = NULL,
      
      # Tab 1 — Gene Activity Violin
      tabPanel(
        title = tagList(icon("chart-bar"), " Gene Activity"),
        fluidRow(
          column(12,
                 plotOutput("accessPlot", height = "480px") |>
                   withSpinner(color = "#2980b9")
          )
        ),
        fluidRow(
          column(12,
                 tags$small(tags$b("How to read: "),
                            "Each violin shows the distribution of log‑normalised gene ",
                            "activity scores (estimated from chromatin accessibility) for ",
                            "that cell type. Wider = more cells at that expression level. ",
                            "* / ** / *** = Wilcoxon p < 0.05 / 0.01 / 0.001.",
                            style = "color:#7f8c8d; padding: 6px 0;")
          )
        )
      ),
      
      # Tab 2 — UMAP Feature Plot
      tabPanel(
        title = tagList(icon("project-diagram"), " UMAP"),
        fluidRow(
          column(12,
                 plotOutput("umapPlot", height = "480px") |>
                   withSpinner(color = "#2980b9")
          )
        ),
        fluidRow(
          column(12,
                 tags$small(tags$b("How to read: "),
                            "Each dot is one cell. Only cells from the selected cell types ",
                            "are highlighted in colour; all others are shown in grey. ",
                            "Colour intensity = gene activity level. White labels show ",
                            "the cell‑type name at its median UMAP centre.",
                            style = "color:#7f8c8d; padding: 6px 0;")
          )
        )
      ),
      
      # Tab 3 — Coverage Plot
      tabPanel(
        title = tagList(icon("wave-square"), " Coverage"),
        fluidRow(
          column(12,
                 plotOutput("coveragePlot", height = "600px") |>
                   withSpinner(color = "#2980b9")
          )
        ),
        fluidRow(
          column(12,
                 tags$small(tags$b("How to read: "),
                            "Each horizontal track shows fragment density (normalised) ",
                            "along the chromosome for one cell type. ",
                            "Tall peaks = accessible chromatin. Grey bars at the bottom = called peaks.",
                            style = "color:#7f8c8d; padding: 6px 0;")
          )
        )
      ),
      
      # Tab 4 — Gene Info Table
      tabPanel(
        title = tagList(icon("table"), " Gene Info"),
        fluidRow(
          column(12, br(),
                 DTOutput("geneTable")
          )
        )
      ),
      
      # Tab 5 — Summary Stats Table
      tabPanel(
        title = tagList(icon("list-alt"), " Summary Stats"),
        fluidRow(
          column(12, br(),
                 DTOutput("summaryTable")
          )
        )
      )
    )
  )
)

# ── 7. Server ─────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ── 7a. Update gene choices when chromosome changes ──────────
  observeEvent(input$chrom, {
    genes_on_chr <- gene_df$gene[gene_df$chr == input$chrom]
    genes_on_chr <- sort(unique(genes_on_chr))
    
    if (length(genes_on_chr) == 0) {
      updateSelectizeInput(session, "gene",
                           choices  = c("No genes on this chromosome" = ""),
                           server   = TRUE)
    } else {
      updateSelectizeInput(session, "gene",
                           choices  = genes_on_chr,
                           selected = genes_on_chr[1],
                           server   = TRUE)
    }
  }, ignoreInit = FALSE)
  
  # ── 7b. Select / Clear all cell types ────────────────────────
  observeEvent(input$selectAll, {
    updateCheckboxGroupInput(session, "celltypes",
                             selected = all_cell_types)
  })
  observeEvent(input$clearAll, {
    updateCheckboxGroupInput(session, "celltypes", selected = character(0))
  })
  
  # ── 7c. Validation helper ─────────────────────────────────────
  validation_ok <- reactive({
    gene_ok <- !is.null(input$gene) && nchar(input$gene) > 0 &&
      input$gene != "No genes on this chromosome" &&
      input$gene %in% rownames(pbmc[["RNA"]])
    cell_ok  <- length(input$celltypes) >= 2
    list(gene = gene_ok, cell = cell_ok, ok = gene_ok && cell_ok)
  })
  
  output$validationMsg <- renderUI({
    v <- validation_ok()
    msgs <- character(0)
    if (!v$gene) msgs <- c(msgs, "⚠ Please select a valid gene.")
    if (!v$cell) msgs <- c(msgs, "⚠ Please select at least 2 cell types.")
    if (length(msgs)) {
      tags$div(class = "alert alert-warning",
               style = "margin: 6px 15px;",
               HTML(paste(msgs, collapse = "<br/>")))
    }
  })
  
  # ── 7d. Subset pbmc to selected cell types ────────────────────
  pbmc_sub <- eventReactive(input$plotBtn, {
    req(validation_ok()$ok)
    keep <- pbmc[[cell_type_col]][, 1] %in% input$celltypes
    pbmc[, keep]
  })
  
  # Gene activity data frame for the chosen gene + cell types
  activity_df <- eventReactive(input$plotBtn, {
    req(validation_ok()$ok)
    sub_obj <- pbmc_sub()
    DefaultAssay(sub_obj) <- "RNA"
    df <- FetchData(sub_obj, vars = input$gene)
    colnames(df) <- "activity"
    df$celltype <- factor(sub_obj[[cell_type_col]][, 1],
                          levels = input$celltypes)
    df
  })
  
  # ── 7e. Info boxes ────────────────────────────────────────────
  output$infoGene <- renderInfoBox({
    gene_sel <- if (!is.null(input$gene) && nchar(input$gene) > 0) input$gene else "—"
    infoBox("Selected Gene", gene_sel, icon = icon("dna"),
            color = "blue", fill = TRUE)
  })
  
  output$infoChrom <- renderInfoBox({
    chr_sel <- if (!is.null(input$chrom)) input$chrom else "—"
    info <- gene_df[gene_df$gene == input$gene & gene_df$chr == chr_sel, ]
    loc <- if (nrow(info) > 0)
      paste0(chr_sel, ":", format(info$start[1], big.mark = ","),
             "–", format(info$end[1], big.mark = ","))
    else chr_sel
    infoBox("Genomic Location", loc, icon = icon("map-marker-alt"),
            color = "green", fill = TRUE)
  })
  
  output$infoCellCount <- renderInfoBox({
    n <- if (!is.null(input$celltypes) && length(input$celltypes) >= 1) {
      sum(pbmc[[cell_type_col]][, 1] %in% input$celltypes)
    } else 0
    infoBox("Cells Selected", format(n, big.mark = ","),
            icon = icon("circle-notch"), color = "yellow", fill = TRUE)
  })
  
  # ── 7f. Tab 1 — Gene Activity Violin ─────────────────────────
  output$accessPlot <- renderPlot({
    input$plotBtn
    req(validation_ok()$ok)
    df <- activity_df()
    
    p <- ggplot(df, aes(x = celltype, y = activity, fill = celltype)) +
      geom_violin(scale = "width", trim = FALSE, alpha = 0.75, linewidth = 0.4) +
      geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white",
                   alpha = 0.85, linewidth = 0.5) +
      scale_fill_brewer(palette = "Set2") +
      labs(
        title    = paste0("Chromatin Accessibility — ", input$gene,
                          "  (", input$chrom, ")"),
        subtitle = paste0("Gene activity proxy from fragment counts over gene body + 2kb promoter"),
        x        = "Cell Type",
        y        = "Log‑normalised Gene Activity Score",
        fill     = "Cell Type"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        axis.text.x      = element_text(angle = 35, hjust = 1, size = 11),
        plot.title       = element_text(face = "bold", size = 14),
        plot.subtitle    = element_text(colour = "grey50", size = 10),
        legend.position  = "none",
        panel.grid.minor = element_blank()
      )
    
    if (length(input$celltypes) == 2) {
      ct_levels <- levels(df$celltype)
      n1 <- sum(df$celltype == ct_levels[1])
      n2 <- sum(df$celltype == ct_levels[2])
      if (n1 > 1 && n2 > 1) {
        p <- p + ggsignif::geom_signif(
          comparisons      = list(ct_levels),
          test             = "wilcox.test",
          map_signif_level = function(pv)
            if (pv < 0.001) "***" else if (pv < 0.01) "**" else if (pv < 0.05) "*" else "ns",
          textsize  = 4.5,
          tip_length = 0.01,
          vjust      = -0.2
        )
      }
    } else {
      ct_levels <- levels(df$celltype)
      pairs <- combn(ct_levels, 2, simplify = FALSE)
      valid_pairs <- Filter(function(pr)
        sum(df$celltype == pr[1]) > 1 && sum(df$celltype == pr[2]) > 1,
        pairs)
      if (length(valid_pairs) > 0 && length(valid_pairs) <= 6) {
        p <- p + ggsignif::geom_signif(
          comparisons      = valid_pairs,
          test             = "wilcox.test",
          step_increase    = 0.08,
          map_signif_level = function(pv)
            if (pv < 0.001) "***" else if (pv < 0.01) "**" else if (pv < 0.05) "*" else "ns",
          textsize  = 3.5,
          tip_length = 0.01
        )
      }
    }
    p
  })
  
  # ── 7g. Tab 2 — UMAP Feature Plot (with cluster labels) ──────
  output$umapPlot <- renderPlot({
    input$plotBtn
    req(validation_ok()$ok)
    
    DefaultAssay(pbmc) <- "RNA"
    
    expr_all <- FetchData(pbmc, vars = input$gene)
    highlight_mask <- pbmc[[cell_type_col]][, 1] %in% input$celltypes
    expr_all[!highlight_mask, 1] <- NA
    
    pbmc_tmp <- pbmc
    pbmc_tmp[["plot_expr"]] <- expr_all[, 1]
    
    umap_coords <- as.data.frame(Embeddings(pbmc_tmp, "umap"))
    colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
    umap_coords$expr     <- pbmc_tmp$plot_expr
    umap_coords$selected <- highlight_mask
    umap_coords$celltype <- pbmc_tmp[[cell_type_col]][, 1]
    
    expr_vals   <- umap_coords$expr[!is.na(umap_coords$expr)]
    cap_val     <- if (length(expr_vals) > 0) quantile(expr_vals, 0.95) else 1
    
    ggplot() +
      # Layer 1: all unselected cells in grey
      geom_point(
        data  = umap_coords[!umap_coords$selected, ],
        aes(x = UMAP_1, y = UMAP_2),
        colour = "grey85", size = 0.3, alpha = 0.5
      ) +
      # Layer 2: selected cells coloured by gene activity
      geom_point(
        data  = umap_coords[umap_coords$selected, ],
        aes(x = UMAP_1, y = UMAP_2, colour = pmin(expr, cap_val)),
        size = 0.6, alpha = 0.85
      ) +
      scale_colour_gradient(
        low     = "#d0e8f7",
        high    = "#08306b",
        na.value = "grey85",
        name    = "Gene\nActivity",
        limits  = c(0, cap_val)
      ) +
      # Layer 3: cell‑type labels at median centre (ALL clusters)
      geom_label(
        data    = cluster_centroids,
        aes(x = cx, y = cy, label = celltype),
        fill    = "white",
        alpha   = 0.75,
        size    = 3.5,
        fontface = "bold",
        color   = "black",
        label.size = 0.2,
        show.legend = FALSE
      ) +
      labs(
        title    = paste0("UMAP — ", input$gene, "  Gene Activity"),
        subtitle = paste0("Highlighted: ", paste(input$celltypes, collapse = ", ")),
        x        = "UMAP 1",
        y        = "UMAP 2"
      ) +
      theme_minimal(base_size = 13) +
      theme(
        plot.title      = element_text(face = "bold", size = 14),
        plot.subtitle   = element_text(colour = "grey50", size = 10),
        panel.grid      = element_blank(),
        legend.position = "right"
      ) +
      guides(colour = guide_colorbar(barheight = 8, barwidth = 1.2))
  })
  
  # ── 7h. Tab 3 — Coverage Plot (FIXED) ─────────────────────────
  output$coveragePlot <- renderPlot({
    input$plotBtn
    req(validation_ok()$ok)
    req(input$gene %in% gene_df$gene)   # gene must be in lookup
    
    info <- gene_df[gene_df$gene == input$gene, ]
    if (nrow(info) == 0) {
      return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                                 label = "Gene not found in annotation.") + theme_void())
    }
    
    region_str <- paste0(
      info$chr[1], "-",
      max(1, info$start[1] - 3000), "-",
      info$end[1] + 3000
    )
    
    sub_pbmc <- pbmc_sub()
    DefaultAssay(sub_pbmc) <- "peaks"
    Idents(sub_pbmc) <- sub_pbmc[[cell_type_col]][, 1]
    
    tryCatch({
      CoveragePlot(
        object            = sub_pbmc,
        region            = region_str,
        extend.upstream   = 0,   # already included
        extend.downstream = 0
      )
    }, error = function(e) {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, size = 5, colour = "grey40",
                 label = paste0("Coverage plot unavailable.\n",
                                conditionMessage(e))) +
        theme_void()
    })
  })
  
  # ── 7i. Tab 4 — Gene Info Table ──────────────────────────────
  output$geneTable <- renderDT({
    req(input$gene, nchar(input$gene) > 0)
    info <- gene_df[gene_df$gene == input$gene, ]
    if (nrow(info) == 0) {
      info <- data.frame(Message = "Gene not found in annotation.")
    } else {
      info <- info |>
        dplyr::transmute(
          Gene       = gene,
          Chromosome = chr,
          Start      = format(start, big.mark = ","),
          End        = format(end,   big.mark = ","),
          Strand     = strand,
          Length_bp  = format(end - start + 1, big.mark = ",")
        ) |>
        dplyr::distinct()
    }
    datatable(info,
              rownames  = FALSE,
              options   = list(dom = "t", pageLength = 5),
              class     = "compact stripe hover")
  })
  
  # ── 7j. Tab 5 — Summary Stats Table ──────────────────────────
  output$summaryTable <- renderDT({
    input$plotBtn
    req(validation_ok()$ok)
    df <- activity_df()
    
    summary_tbl <- df |>
      dplyr::group_by(celltype) |>
      dplyr::summarise(
        N_cells      = dplyr::n(),
        Mean_activity = round(mean(activity, na.rm = TRUE), 4),
        Median       = round(median(activity, na.rm = TRUE), 4),
        SD           = round(sd(activity, na.rm = TRUE), 4),
        Min          = round(min(activity, na.rm = TRUE), 4),
        Max          = round(max(activity, na.rm = TRUE), 4),
        Pct_nonzero  = round(mean(activity > 0, na.rm = TRUE) * 100, 1),
        .groups      = "drop"
      ) |>
      dplyr::rename(
        `Cell Type`          = celltype,
        `N Cells`            = N_cells,
        `Mean Activity`      = Mean_activity,
        `Median`             = Median,
        `SD`                 = SD,
        `Min`                = Min,
        `Max`                = Max,
        `% Cells Active (>0)` = Pct_nonzero
      )
    
    datatable(summary_tbl,
              rownames  = FALSE,
              options   = list(dom = "t", pageLength = 20),
              class     = "compact stripe hover") |>
      DT::formatStyle(
        "Mean Activity",
        background = DT::styleColorBar(range(summary_tbl$`Mean Activity`), "#2980b9"),
        backgroundSize   = "98% 88%",
        backgroundRepeat = "no-repeat",
        backgroundPosition = "center"
      )
  })
}

# ── 8. Launch ─────────────────────────────────────────────────
shinyApp(ui, server)

