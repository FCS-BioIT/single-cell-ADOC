suppressPackageStartupMessages({
  library(optparse)
  library(Seurat)
  library(CellChat)
  library(patchwork)
})


# ---- Parallel config (do this ONCE near the top) ----
n_workers <- opt$ncores
if (requireNamespace("future", quietly = TRUE)) {
  options(future.globals.maxSize = 8 * 1024^3) # 8 GiB; adjust to your RAM
  if (.Platform$OS.type == "windows") {
    future::plan(future::multisession, workers = n_workers)
  } else {
    # multicore is a bit lighter on *nix; falls back if not supported
    future::plan(future::multicore, workers = n_workers)
  }
  on.exit(
    {
      try(future::plan(future::sequential), silent = TRUE)
    },
    add = TRUE
  )
}


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
    type = "character", default = "output/ccc",
    help = "Directory to write outputs [default: %default]"
  ),
  make_option("--species",
    type = "character", default = "human",
    help = "Species: human|mouse [default: %default]"
  ),
  make_option("--group_col",
    type = "character", default = "group",
    help = "Metadata column with conditions [default: %default]"
  ),
  make_option("--celltype_col",
    type = "character", default = "celltype",
    help = "Metadata column with cell identities [default: %default]"
  ),
  make_option("--min_cells",
    type = "integer", default = 10,
    help = "Min cells per group/celltype to keep edges [default: %default]"
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

# ---- Load data ----
stopifnot(file.exists(opt$seurat_rds))
seu <- LoadSeuratRds(file = opt$seurat_rds)

seu$group <- factor(seu$group, levels = c("CNT", "HT"))

unique(seu$celltype)

# Celltypes in common
seu$celltype_common <- seu$celltype

seu$celltype_common[which(
  seu$celltype == "Stromal" |
    seu$celltype == "Non-decidual stroma" |
    seu$celltype == "Decidual stroma"
)] <- "Stroma"
seu$celltype_common[which(
  seu$celltype == "Luminal-like epithelium" |
    seu$celltype == "HT epithelium I" |
    seu$celltype == "HT epithelium II"
)] <- "Epithelium"

unique(seu$celltype_common)
opt$celltype_col <- "celltype_common"
# Checks
meta <- seu@meta.data
for (col in c(opt$group_col, opt$celltype_col)) {
  if (!col %in% colnames(meta)) {
    stop(sprintf("Metadata column '%s' not found.", col))
  }
}
if (length(unique(meta[[opt$group_col]])) < 2) {
  stop(
    sprintf(
      "'%s' must have at least two levels (e.g., Control/Hormone).",
      opt$group_col
    )
  )
}
meta$group <- factor(meta$group, levels = c("CNT", "HT"))
# ---- Species DBs ----
if (tolower(opt$species) == "human") {
  cellchatdb <- CellChatDB.human
  PPI <- PPI.human # nolint
} else if (tolower(opt$species) == "mouse") {
  cellchatdb <- CellChatDB.mouse
  PPI <- PPI.mouse # nolint
} else {
  stop("species must be 'human' or 'mouse'")
}

# ---- Split by condition ----
Idents(seu) <- factor(meta[[opt$celltype_col]])
seu_list <- SplitObject(seu, split.by = opt$group_col)
cond_names <- names(seu_list)
message("# Conditions detected: ", paste(cond_names, collapse = ", "))

# ---- Core pipeline for one condition ----
analyze_one <- function(obj, cond, cellchatdb, PPI, min_cells, ncores) {
  Idents(obj) <- obj@meta.data[[opt$celltype_col]]
  obj$samples <- obj$sample_id
  cellchat <- createCellChat(object = obj, group.by = opt$celltype_col)
  cellchat@DB <- cellchatdb

  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- smoothData(cellchat, adj = PPI)


  set.seed(1)
  cellchat <- computeCommunProb(cellchat, raw.use = FALSE)
  cellchat <- filterCommunication(cellchat, min.cells = min_cells)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  cellchat <- netAnalysis_computeCentrality(cellchat)

  # Save per-condition outputs
  saveRDS(cellchat, file.path(opt$outdir, sprintf("cellchat_%s.rds", cond)))
  # Communication table
  df <- subsetCommunication(cellchat)
  write.csv(
    df,
    file.path(opt$outdir, sprintf("communication_table_%s.csv", cond)),
    row.names = FALSE
  )

  # Key plots
  try(
    {
      png(
        file.path(opt$outdir, sprintf("bubble_pathways_%s.png", cond)),
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
        file.path(opt$outdir, sprintf("circle_overall_%s.png", cond)),
        width = 1600,
        height = 1600,
        res = 200
      )
      print(netVisual_circle(cellchat@net$count,
        vertex.weight = as.numeric(table(cellchat@idents)),
        weight.scale = TRUE,
        label.edge = FALSE,
        title.name = paste0("Number of interactions - ", cond)
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
          sprintf("heatmap_top%02d_pathways_%s.png", topn, cond)
        ),
        width = 1600,
        height = 1200,
        res = 200
      )
      print(
        netAnalysis_signalingRole_heatmap(
          cellchat,
          pattern = "all", signaling = pathways_show
        )
      )
      dev.off()
    },
    silent = TRUE
  )

  cellchat
}

# ---- Run per condition ----
cellchat_list <- list()
for (cond in cond_names) {
  message("# Running CellChat for: ", cond)
  cellchat_list[[cond]] <- analyze_one(
    seu_list[[cond]],
    cond,
    cellchatdb,
    PPI,
    opt$min_cells, opt$ncores
  )
}

# ---- Comparison ----
message("# Merging for comparison")
cellchat_merged <- mergeCellChat(
  cellchat_list,
  add.names = (names(cellchat_list))
)
saveRDS(
  cellchat_merged,
  file.path(opt$outdir, "cellchat_merged_comparison.rds")
)

cellchat_list <- list(
  CNT = readRDS("output/ccc/cellchat_CNT.rds"),
  HT = readRDS("output/ccc/cellchat_HT.rds")
)
# Circle plot of celltype interaction by condition
for (cond in c("HT", "CNT")) {
  svg(
    file.path(opt$outdir, sprintf("circle_overall_%s.svg", cond)),
    width = 8, height = 10
  )
  print(netVisual_circle(
    cellchat_list[[cond]]@net$count,
    vertex.weight = as.numeric(table(cellchat_list[[cond]]@idents)),
    weight.scale = TRUE,
    label.edge = FALSE,
    title.name = paste0("Number of interactions - ", cond)
  ))
  dev.off()
}


# Global interaction comparison
try(
  {
    gg1 <- compareInteractions(
      cellchat_merged,
      show.legend = FALSE,
      group = 1:2
    )
    gg2 <- compareInteractions(
      cellchat_merged,
      show.legend = FALSE,
      group = 1:2,
      measure = "weight"
    )
    g <- gg1 + gg2
    png(
      file.path(opt$outdir, "compare_interaction_count_weight.png"),
      width = 2000,
      height = 900,
      res = 200
    )
    print(g)
    dev.off()
  },
  silent = TRUE
)

# Differential edges & pathways
try(
  {
    png(
      file.path(opt$outdir, "diff_interactions_matrix.png"),
      width = 1600,
      height = 1400,
      res = 200
    )
    print(
      netVisual_diffInteraction(cellchat_merged,
        weight.scale = TRUE,
        measure = "count"
      )
    )
    dev.off()
  },
  silent = TRUE
)

try(
  {
    png(
      file.path(opt$outdir, "diff_signaling_pathways.png"),
      width = 1600,
      height = 1400,
      res = 200
    )
    print(
      netVisual_diffSignaling(cellchat_merged,
        weight.scale = TRUE
      )
    )
    dev.off()
  },
  silent = TRUE
)

# Rank pathways by change
try(
  {
    png(
      file.path(opt$outdir, "rank_pathways_comparison.png"),
      width = 2000,
      height = 5400,
      res = 400
    )
    print(rankNet(cellchat_merged, mode = "comparison", stacked = TRUE))
    dev.off()
  },
  silent = TRUE
)

# Export tables: per-condition + delta
message("# Exporting communication tables")
for (cond in names(cellchat_list)) {
  df <- subsetCommunication(cellchat_list[[cond]])
  write.csv(
    df,
    file.path(opt$outdir, sprintf("communication_table_%s.csv", cond)),
    row.names = FALSE
  )
}

# Pathway-level comparison table (if available)
# Note: CellChat stores comparison in merged object; extract deltas if present
delta_list <- list()
for (path in cellchat_merged@netP$pathways) {
  # per pathway matrices by condition
  mats <- cellchat_merged@netP$interaction[[path]]
  if (is.list(mats) && length(mats) >= 2) {
    # Assume first two correspond to the two groups in merge order
    mat_delta <- mats[[2]] - mats[[1]]
    df_delta <- as.data.frame(as.table(mat_delta))
    colnames(df_delta) <- c("source", "target", "delta_weight")
    df_delta$pathway <- path
    delta_list[[path]] <- df_delta
  }
}
if (length(delta_list) > 0) {
  delta_df <- do.call(rbind, delta_list)
  write.csv(
    delta_df, file.path(opt$outdir, "delta_interactions_by_pathway.csv"),
    row.names = FALSE
  )
}

# Centrality changes per pathway (summary)
try(
  {
    png(
      file.path(opt$outdir, "centrality_diff_contribution.png"),
      width = 1800,
      height = 1200,
      res = 200
    )
    print(
      netAnalysis_contribution(
        cellchat_merged,
        signaling = NULL,
        measure = "overall"
      )
    )
    dev.off()
  },
  silent = TRUE
)

# Session info for reproducibility
writeLines(
  c(capture.output(sessionInfo())),
  con = file.path(opt$outdir, "sessionInfo.txt")
)

message("# Done. Outputs in: ", opt$outdir)



# ---- Bundle outputs for downstream plotting scripts ----
# Order the two groups explicitly (edit to your desired order):
desired_order <- intersect(c("CNT", "HT"), names(cellchat_list))
if (length(desired_order) != 2) desired_order <- names(cellchat_list)[1:2]

cc_list <- cellchat_list[desired_order]
multiple_cc <- cellchat_merged

# Build a flat 'net' data.frame per pathway with dataset labels
mk_net <- function(obj_list) {
  out <- lapply(names(obj_list), function(grp) {
    df <- subsetCommunication(
      obj_list[[grp]],
      slot.name = "net"
    ) # edges at pathway level
    if (!is.null(df) && nrow(df) > 0) {
      df$datasets <- grp
    }
    df
  })
  out <- do.call(rbind, out)
  # standardize a few expected column names used later
  if (
    !"pathway_name" %in% colnames(out) && "pathway_name_2" %in% colnames(out)) {
    out$pathway_name <- out$pathway_name_2
  }
  out
}

net_df <- tryCatch(mk_net(cc_list), error = function(e) NULL)

ccc_bundle <- list(
  CC_list = cc_list, # named list of two CellChat objects
  multiple_CC = multiple_cc, # merged comparison object
  net = net_df # optional flat table of edges
)

saveRDS(ccc_bundle, file.path(opt$outdir, "ccc_bundle_HT_vs_CNT.rds"))
message("# Saved bundle: ", file.path(opt$outdir, "ccc_bundle_HT_vs_CNT.rds"))


sink()
