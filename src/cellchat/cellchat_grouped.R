suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(CellChat)
  library(patchwork)
})

# -------------------- CLI options --------------------
option_list <- list(
  make_option("--seurat_rds",
    type = "character",
    default = paste0(
      "output/annotation/marker_annotation/data/",
      "seurat_filtered_annoted.rds"
    ),
    help = "Path to Seurat .rds input [default: %default]"
  ),
  make_option("--outdir",
    type = "character", default = "output/ccc_pooled",
    help = "Directory to write outputs [default: %default]"
  ),
  make_option("--species",
    type = "character", default = "human",
    help = "Species: human|mouse [default: %default]"
  ),
  make_option("--celltype_col",
    type = "character", default = "celltype_common",
    help = "Metadata column with cell identities [default: %default]"
  ),
  make_option("--min_cells",
    type = "integer", default = 10,
    help = "Min cells per celltype to keep edges [default: %default]"
  ),
  make_option("--ncores",
    type = "integer", default = 8,
    help = "Cores for computeCommunProb (if supported) [default: %default]"
  )
)
opt <- parse_args(OptionParser(option_list = option_list))

dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(opt$outdir, "log.txt")
sink(logf, split = TRUE)

message("# ==== Args ====")
print(opt)

# -------------------- Parallel config --------------------
n_workers <- opt$ncores
if (requireNamespace("future", quietly = TRUE)) {
  options(future.globals.maxSize = 8 * 1024^3) # 8 GiB
  if (.Platform$OS.type == "windows") {
    future::plan(future::multisession, workers = n_workers)
  } else {
    future::plan(future::multicore, workers = n_workers)
  }
  on.exit(
    {
      try(future::plan(future::sequential), silent = TRUE)
    },
    add = TRUE
  )
}

# -------------------- Load data --------------------
stopifnot(file.exists(opt$seurat_rds))
seu <- readRDS(opt$seurat_rds) # Using base readRDS for portability

# Optional: harmonize some fine-grained labels into common ones
if (!("celltype_common" %in% colnames(seu@meta.data))) {
  seu$celltype_common <- seu@meta.data[[opt$celltype_col]]
  seu$celltype_common[seu$celltype %in% c(
    "Stromal",
    "Non-decidual stroma",
    "Decidual stroma"
  )] <- "Stroma"
  seu$celltype_common[seu$celltype %in% c(
    "Luminal-like epithelium",
    "HT epithelium I",
    "HT epithelium II"
  )] <- "Epithelium"
}

meta <- seu@meta.data
stopifnot(opt$celltype_col %in% colnames(meta))

# -------------------- Species DBs --------------------
if (tolower(opt$species) == "human") {
  cellchatDB <- CellChatDB.human # nolint
  PPI <- PPI.human # nolint
} else if (tolower(opt$species) == "mouse") {
  cellchatDB <- CellChatDB.mouse # nolint
  PPI <- PPI.mouse # nolint
} else {
  stop("species must be 'human' or 'mouse'")
}

# -------------------- CellChat pooled analysis --------------------
Idents(seu) <- seu@meta.data[[opt$celltype_col]]
seu$samples <- seu$sample_id %||% NULL # optional if present

message("# Running CellChat for ALL samples pooled")
cellchat <- createCellChat(object = seu, group.by = opt$celltype_col)
cellchat@DB <- cellchatDB

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- smoothData(cellchat, adj = PPI)

set.seed(1)
cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = opt$min_cells)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat)

# -------------------- Save objects & tables --------------------
saveRDS(cellchat, file.path(opt$outdir, "cellchat_all.rds"))
cellchat <- readRDS(file.path(opt$outdir, "cellchat_all.rds"))

df_comm <- subsetCommunication(cellchat) # all edges (ligand–receptor level)
write.csv(
  df_comm,
  file.path(opt$outdir, "communication_table_all.csv"),
  row.names = FALSE
)

df_path <- subsetCommunication(
  cellchat,
  slot.name = "netP"
) # pathway-level edges if desired
if (!is.null(df_path) && nrow(df_path) > 0) {
  write.csv(
    df_path,
    file.path(opt$outdir, "communication_table_pathway_all.csv"),
    row.names = FALSE
  )
}

# -------------------- Key plots (optional but handy) --------------------
try(
  {
    png(
      file.path(
        opt$outdir,
        "bubble_pathways_all.png"
      ),
      width = 1600,
      height = 1200,
      res = 200
    )
    print(netAnalysis_signalingRole_network(cellchat, slot.name = "netP"))
    dev.off()
  },
  silent = TRUE
)

try(
  {
    png(
      file.path(
        opt$outdir,
        "circle_overall_all.png"
      ),
      width = 3200,
      height = 4000,
      res = 400
    )
    print(netVisual_circle(
      cellchat@net$count,
      vertex.weight = as.numeric(table(cellchat@idents)),
      weight.scale = TRUE,
      label.edge = FALSE,
      edge.width.max = 20,
      title.name = "Number of interactions - ALL"
    ))
    dev.off()

    svg(file.path(opt$outdir, "circle_overall_all.svg"), width = 8, height = 10)
    print(netVisual_circle(
      cellchat@net$count,
      vertex.weight = as.numeric(table(cellchat@idents)),
      weight.scale = TRUE,
      label.edge = FALSE,
      edge.width.max = 20,
      title.name = "Number of interactions - ALL"
    ))
    dev.off()
  },
  silent = TRUE
)

try(
  {
    topn <- 10
    pathways_show <- head(cellchat@netP$pathways, topn)
    png(
      file.path(
        opt$outdir,
        sprintf(
          "heatmap_top%02d_pathways_all.png",
          topn
        )
      ),
      width = 1600, height = 1200, res = 200
    )
    print(
      netAnalysis_signalingRole_heatmap(
        cellchat,
        pattern = "all",
        signaling = pathways_show
      )
    )
    dev.off()
  },
  silent = TRUE
)

# -------------------- Reproducibility --------------------
writeLines(
  c(capture.output(sessionInfo())),
  con = file.path(opt$outdir, "sessionInfo.txt")
)

message("# Done. Outputs in: ", opt$outdir)
sink()
