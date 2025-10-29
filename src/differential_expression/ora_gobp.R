# --- WebGestaltR KEGG ORA on your DE outputs ---
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(Seurat)
  library(WebGestaltR)
})

seurat_obj_path <- paste0(
  "output/annotation/marker_annotation/data/",
  "seurat_filtered_annoted.rds"
)
de_dir <- "output/differential_expression"
out_dir_enrich <- file.path(de_dir, "WebGestalt_GOBP_ORA")
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
    message(sprintf("[%s] Skipped: too few genes (%d).", project, length(genes)))
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

  if (is.null(res) || !is.data.frame(res) || nrow(res) == 0) {
    message(sprintf("[%s] No significant terms (empty result).", project))
    return(invisible(NULL))
  }

  res <- tibble::as_tibble(res)
  names(res) <- tolower(names(res))

  fdr_col <- if ("fdr" %in% names(res)) "fdr"
  p_col <- if ("pvalue" %in% names(res)) "pvalue"
  er_col <- dplyr::first(intersect(
    c("enrichmentratio", "richfactor", "ratio", "es"),
    names(res)
  ))

  # merged sorting logic into one call
  if (!is.null(fdr_col)) {
    sort_vars <- c(fdr_col, er_col, p_col)
    sort_vars <- sort_vars[!is.na(sort_vars)]
    res <- dplyr::arrange(res, !!!rlang::syms(sort_vars))
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
    project = paste0(base, "_GOBP_UP"), outdir = out_dir_enrich
  )
  run_kegg_ora(
    de_down,
    project = paste0(base, "_GOBP_DOWN"), outdir = out_dir_enrich
  )
}
message("WebGestaltR GOBP ORA done. Results under: ", out_dir_enrich)
