# --- WebGestaltR KEGG ORA on your DE outputs ---
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(Seurat)
  library(WebGestaltR)
  library(ggplot2)
})

seurat_obj_path <- paste0(
  "output/annotation/marker_annotation/data/",
  "seurat_filtered_annoted.rds"
)
de_dir <- "output/differential_expression"
out_dir_enrich <- file.path(de_dir, "WebGestalt_KEGG_ORA")
dir.create(out_dir_enrich, showWarnings = FALSE, recursive = TRUE)

# Load Seurat to define a data-driven background (all expressed genes)
obj <- LoadSeuratRds(seurat_obj_path)
DefaultAssay(obj) <- "RNA"
bg_genes <- rownames(GetAssayData(obj, assay = "RNA", slot = "data"))
bg_genes <- unique(na.omit(bg_genes))

# Helper: run ORA for one gene set
run_kegg_ora <- function(genes, project, outdir) {
  genes <- unique(na.omit(genes))
  if (length(genes) < 10) {
    message(sprintf(
      "[%s] Skipped: too few genes (%d).",
      project, length(genes)
    ))
    return(invisible(NULL))
  }
  res <- WebGestaltR(
    organism = "hsapiens",
    enrichMethod = "ORA",
    enrichDatabase = "pathway_KEGG",
    interestGene = genes,
    interestGeneType = "genesymbol",
    referenceGene = bg_genes,
    referenceGeneType = "genesymbol",
    projectName = project,
    fdrMethod = "BH",
    sigMethod = "fdr",
    fdrThr = 0.05,
    minNum = 5,
    isOutput = FALSE
  )
  # Save if we got something back
  if (is.null(res) || !is.data.frame(res) || nrow(res) == 0) {
    message(sprintf("[%s] No significant terms (empty result).", project))
    return(invisible(NULL))
  }

  # Normalize names → lower-case, then sort safely
  res <- tibble::as_tibble(res)
  names(res) <- tolower(names(res))

  fdr_col <- if ("fdr" %in% names(res)) "fdr" else NULL
  p_col <- if ("pvalue" %in% names(res)) "pvalue" else NULL
  er_col <- dplyr::first(
    intersect(
      c("enrichmentratio", "richfactor", "ratio", "es"),
      names(res)
    )
  )

  if (!is.null(fdr_col) && !is.null(er_col)) {
    res <- dplyr::arrange(res, .data[[fdr_col]], dplyr::desc(.data[[er_col]]))
  } else if (!is.null(fdr_col) && !is.null(p_col)) {
    res <- dplyr::arrange(res, .data[[fdr_col]], .data[[p_col]])
  } else if (!is.null(fdr_col)) {
    res <- dplyr::arrange(res, .data[[fdr_col]])
  }

  res$project <- project
  readr::write_csv(res, file.path(outdir, paste0(project, "_KEGG_ORA.csv")))
  saveRDS(res, file.path(outdir, paste0(project, "_KEGG_ORA.rds")))
  invisible(res)
}

# Find all DE CSVs produced by your script
de_files <- list.files(de_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)

for (f in de_files) {
  de <- readr::read_csv(f, show_col_types = FALSE)
  # pick the correct LFC column your script emits
  lfc_col <- if ("avg_log2FC" %in% names(de)) "avg_log2FC" else "avg_logFC"
  if (!("p_val_adj" %in% names(de))) next

  # define gene sets. Use an FDR filter; optional LFC magnitude if you want
  de_up <- de %>%
    filter(p_val_adj < 0.05, !!sym(lfc_col) > 1.5) %>%
    pull(gene) %>%
    unique()
  de_down <- de %>%
    filter(p_val_adj < 0.05, !!sym(lfc_col) < -1.5) %>%
    pull(gene) %>%
    unique()

  # project labels from filename
  base <- basename(f) %>% str_remove("\\.csv$")
  run_kegg_ora(
    de_up,
    project = paste0(base, "_KEGG_UP"),
    outdir = out_dir_enrich
  )
  run_kegg_ora(
    de_down,
    project = paste0(base, "_KEGG_DOWN"),
    outdir = out_dir_enrich
  )
}
message("WebGestaltR KEGG ORA done. Results under: ", out_dir_enrich)


seu <- obj
# Subset cells of the desired type
gene <- "PAEP"

Idents(seu) <- "celltype_general"
cells_use <- WhichCells(seu, ident = "Stroma HT") # if Idents(seu) is celltype
# or explicitly: cells_use <- rownames(seu@meta.data)[seu$celltype == "Stromal"]

# Get counts from the RNA assay
expr_vals <- GetAssayData(
  seu,
  assay = "RNA", slot = "counts"
)[gene, c(cells_use)]

# Sum across those cells
sum_expr <- sum(expr_vals)
sum_expr

n_detected <- sum(expr_vals > 0)

cells_use_c <- WhichCells(seu, ident = "Stroma No-HT")
expr_vals_c <- GetAssayData(
  seu,
  assay = "RNA", slot = "counts"
)[gene, c(cells_use_c)]
sum_expr_c <- sum(expr_vals_c)
sum_expr_c

n_detected_c <- sum(expr_vals_c > 0)


length(cells_use_c)
n_detected_c
sum_expr_c


n_detected
sum_expr

sum_expr / sum_expr_c
pct1 <- n_detected / length(cells_use)
pct2 <- n_detected_c / length(cells_use_c)


dp <- DotPlot(
  seu,
  features = gene,
  group.by = "celltype_general",
  scale = FALSE
)


extra_text <- paste0(
  "pct Stroma HT = ", round(pct1, digits = 4),
  "\npct Stroma No-HT = ", round(pct2, digits = 4),
  "\nexpr Stroma HT = ", round(sum_expr),
  "\nexpr Stroma C = ", round(sum_expr_c)
)

# Add annotation at bottom
dp + labs(caption = extra_text)

ggsave("output/differential_expression/dotplot_paep.png",
  bg = "white",
  width = 6
)
