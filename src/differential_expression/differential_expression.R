#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(MAST)
  library(dplyr)
  library(tidyverse)
})

# ---- User parameters ----
seurat_obj_path <- paste0(
  "output/annotation/marker_annotation/data/",
  "seurat_filtered_annoted.rds"
)
cluster_contrasts <- list(
  c("Stroma HT", "Stroma No-HT"),
  c("Epithelium HT", "Epithelium No-HT"),
  c("Cycling epithelium HT", "Cycling epithelium No-HT")
)

output_dir <- "output/differential_expression"
assay_use <- "RNA" # use log-normalized RNA for MAST
min_pct <- 0.1 # detectability filter
logfc_thr <- 0.00 # effect-size filter (ln FC)
latent_vars <- c("chip", "nCount_RNA", "mitoRatio", "nFeature_RNA")
max_cells <- NULL # e.g., 5000 to cap per ident
set.seed(42)

# ---- Load & prep ----
obj <- LoadSeuratRds(seurat_obj_path)
stopifnot(is(obj, "Seurat"))

DefaultAssay(obj) <- assay_use
if (!"data" %in% Layers(obj)) {
  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    verbose = FALSE
  )
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

log_path <- file.path(output_dir, "runlog.txt")
cat(format(Sys.time()), "Starting DE run\n", file = log_path)



obj$celltype_general <- obj$celltype

obj$celltype_general[
  which(obj$celltype == "Decidual stroma")
] <- "Stroma HT"
obj$celltype_general[
  which(obj$celltype == "Non-decidual stroma")
] <- "Stroma HT"
obj$celltype_general[
  which(obj$celltype == "Stromal")
] <- "Stroma No-HT"

obj$celltype_general[
  which(obj$celltype == "HT epithelium I")
] <- "Epithelium HT"
obj$celltype_general[
  which(obj$celltype == "HT epithelium II")
] <- "Epithelium HT"
obj$celltype_general[
  which(obj$celltype == "Luminal-like epithelium")
] <- "Epithelium No-HT"

obj$celltype_general[
  which(obj$celltype == "Myofibroblast" & obj$group == "HT")
] <- "Myofibroblast HT"
obj$celltype_general[
  which(obj$celltype == "Myofibroblast" & obj$group == "CNT")
] <- "Myofibroblast No-HT"

obj$celltype_general[
  which(obj$celltype == "Cycling epithelium" & obj$group == "HT")
] <- "Cycling epithelium HT"
obj$celltype_general[
  which(obj$celltype == "Cycling epithelium" & obj$group == "CNT")
] <- "Cycling epithelium No-HT"

DimPlot(
  obj,
  group.by = "celltype_general",
  reduction = "umap.chip.harmony",
  label = TRUE
)
ggsave(
  filename = "output/differential_expression/plots/umap_general_celltypes.png",
  dpi = 400,
  width = 7,
  height = 7,
  bg = "white"
)


Idents(obj) <- "celltype_general"

group_comb <- TRUE

# ---- Iterate contrasts ----
for (contrast in cluster_contrasts) {
  c1 <- contrast[1]
  c2 <- contrast[2]
  out_csv <- file.path(output_dir, sprintf("DE_%s_vs_%s.csv", c1, c2))

  cells_1 <- WhichCells(obj, idents = c1)
  cells_2 <- WhichCells(obj, idents = c2)
  if (length(cells_1) < 20 || length(cells_2) < 20) {
    cat(
      sprintf(
        "%s Skipping %s vs %s (n=%d,%d)\n",
        format(Sys.time()), c1, c2, length(cells_1), length(cells_2)
      ),
      file = log_path, append = TRUE
    )
    next
  }

  obj_sub <- subset(obj, cells = c(cells_1, cells_2))
  Idents(obj_sub) <- factor(Idents(obj_sub), levels = c(c1, c2))

  # Optional balancing
  if (!is.null(max_cells)) {
    obj_sub <- subset(
      obj_sub,
      downsample = pmin(max_cells, min(table(Idents(obj_sub))))
    )
  }

  BiocParallel::register(MulticoreParam(workers = 8))

  de <- FindMarkers(
    object = obj_sub,
    ident.1 = c1,
    ident.2 = c2,
    test.use = "MAST",
    latent.vars = intersect(latent_vars, colnames(obj_sub@meta.data)),
    min.pct = min_pct,
    logfc.threshold = logfc_thr,
    only.pos = FALSE
  ) %>%
    tibble::rownames_to_column("gene") %>%
    arrange(p_val_adj, desc(avg_log2FC))

  # annotate sample sizes
  de$N_cells_c1 <- length(WhichCells(obj_sub, idents = c1))
  de$N_cells_c2 <- length(WhichCells(obj_sub, idents = c2))

  readr::write_csv(de, out_csv)
  cat(
    sprintf(
      "%s %s vs %s done. n=%d/%d → %s\n",
      format(Sys.time()), c1, c2,
      de$N_cells_c1[1], de$N_cells_c2[1], out_csv
    ),
    file = log_path, append = TRUE
  )
}


# You can later read all CSVs,
# combine p-values, and recompute an across-contrast BH.
sessionInfo()
