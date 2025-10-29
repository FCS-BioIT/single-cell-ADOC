library(Seurat)
library(dittoSeq)
library(clustree)
library(ggplot2)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(pheatmap)
library(dplyr)
library(tibble)
library(stringr)

# Loading the dataset to query

out_dir <- "output/reference_projection/"
out_data <- paste0(out_dir, "data/")
out_plots <- paste0(out_dir, "plots/")
dir.create(out_data)
dir.create(out_plots)
seu <- LoadSeuratRds(file = "output/annotation/marker_annotation/data/seurat_filtered_annoted.rds") # nolint

# Labelling query (organ on chip sc rnaseq)

seu$celltype_general <- seu$celltype

seu$celltype_general[
  which(
    seu$celltype == "Decidual stroma" | seu$celltype == "Non-decidual stroma"
  )
] <- "HT stroma"


# Load preprocessed reference:
reference <- LoadSeuratRds(
  file = paste0(
    "references/sc-rnaseq/",
    "endometriumAtlas_cells_with_counts_integrated.rds"
  )
) # nolint

# First approach
# find anchors
anchors <- FindTransferAnchors(
  reference = reference,
  query = seu,
  reduction = "cca",
  dims = 1:30
)

saveRDS(anchors, file = paste0(out_data, "anchors_cca.Rds"))
anchors <- readRDS(paste0(out_data, "anchors_cca.Rds"))


# transfer labels
predictions <- TransferData(
  anchorset = anchors,
  refdata = reference$celltype,
  weight.reduction = "cca",
  dims = 1:30
)

seu <- AddMetaData(object = seu, metadata = predictions)

# Graphic representations
DimPlot(
  seu,
  group.by = "predicted.id",
  reduction = "umap.chip.harmony",
  label = TRUE
)
ggsave(
  filename = paste0(out_plots, "umap_heca_celltypes.png"),
  width = 10
)

pred_scores_lineages_colnames <- c(
  "prediction.score.Mesenchymal",
  "prediction.score.Epithelial",
  "prediction.score.Immune",
  "prediction.score.Endothelial"
)
pred_scores <- c(
  # glandular
  "prediction.score.SOX9_basalis",

  # glandular proliferative (No HT)
  "prediction.score.SOX9_functionalis_I",
  "prediction.score.SOX9_functionalis_II",
  # glandular secretory (HT)
  "prediction.score.preGlandular",
  "prediction.score.Glandular",
  "prediction.score.Glandular_secretory",
  "prediction.score.Glandular_secretory_FGF7",

  # luminal proliferative (No HT)
  "prediction.score.SOX9_luminal",
  # luminal secretory (HT)
  "prediction.score.preLuminal",
  "prediction.score.Luminal",

  # ciliated
  "prediction.score.preCiliated",
  "prediction.score.Ciliated",

  # MUC5B KRT5
  "prediction.score.MUC5B",
  "prediction.score.KRT5",

  # Stroma prolif (No HT)
  "prediction.score.eStromal_MMPs",
  "prediction.score.eStromal",
  "prediction.score.eStromal_cycling",
  # Stroma gland (HT)
  "prediction.score.dStromal_early",
  "prediction.score.dStromal_mid",
  "prediction.score.dStromal_late"
)


celltypes <- c(
  "Cycling epithelium",
  "Luminal-like epithelium",
  "HT epithelium I",
  "HT epithelium II",
  "Stromal",
  "Non-decidual stroma",
  "Decidual stroma",
  "Myofibroblast"
)

celltype_general <- c(
  "Cycling epithelium",
  "Luminal-like epithelium",
  "HT epithelium I",
  "HT epithelium II",
  "Stromal",
  "HT stroma",
  "Myofibroblast"
)


seu$celltype_general <- factor(seu$celltype_general, levels = celltype_general)


pred_cell <- seu@meta.data[, grep(
  colnames(seu@meta.data),
  pattern = "prediction.score"
)]

openxlsx::write.xlsx(
  pred_cell,
  file = paste0(out_data, "scores_cca.xlsx"), rowNames = TRUE
)



# Transfer also the celltypes from the main lineages in Vento-Tormo dataset
# transfer labels
predictions_lineage <- TransferData(
  anchorset = anchors,
  refdata = reference$lineage,
  weight.reduction = "cca",
  dims = 1:30
)

colnames(predictions_lineage) <- paste0(
  colnames(predictions_lineage),
  "_lineage"
)

seu <- AddMetaData(object = seu, metadata = predictions_lineage)

pred_scores_lineages <- colnames(predictions_lineage)[-c(1, 6)]

library(ComplexHeatmap)
mat_lineages <- as.matrix(seu@meta.data[, c(pred_scores_lineages)])
mat <- as.matrix(seu@meta.data[, c(pred_scores_lineages, pred_scores)])
celltypes <- seu$celltype_general

avg_mat_lineages <- aggregate(
  mat_lineages,
  by = list(Celltype = celltypes), FUN = mean
)
avg_mat <- aggregate(mat, by = list(CellType = celltypes), FUN = mean)

rownames(avg_mat) <- avg_mat$CellType
avg_mat <- avg_mat[, -1]

rownames(avg_mat_lineages) <- avg_mat_lineages$Celltype
avg_mat_lineages <- avg_mat_lineages[, -1]

openxlsx::write.xlsx(
  avg_mat_lineages,
  paste0(out_data, "main_lineages_heca_scores.xlsx"),
  rowNames = TRUE
)

openxlsx::write.xlsx(
  avg_mat,
  paste0(out_data, "all_celltypes_heca_scores.xlsx"),
  rowNames = TRUE
)


ht_mat <- t(avg_mat[, colSums(avg_mat) > 0.1])
ht_mat <- t(avg_mat[, colSums(avg_mat) >= 0.0])

# ---- row groups ----
n <- nrow(ht_mat)
lineage_idx <- which(grepl("_lineage$", rownames(ht_mat))) # first 4 expected
last4_idx <- (n - 3):n

row_group <- rep("Epithelial cell states", n)
row_group[last4_idx] <- "Stromal cell states"
if (length(lineage_idx) > 0) row_group[lineage_idx] <- "Main lineages"

ann_row <- data.frame(
  Category = factor(row_group,
    levels = c(
      "Main lineages",
      "Epithelial cell states",
      "Stromal cell states"
    )
  ),
  row.names = rownames(ht_mat)
)

# gaps after last lineage row, and before the last 4 rows
gaps_row <- sort(unique(c(
  ifelse(length(lineage_idx) > 0, max(lineage_idx), 4),
  n - 4
)))

# ---- column groups (4 epi, 4 stromal) ----
stopifnot(ncol(ht_mat) == 7) # as per your spec
ann_col <- data.frame(
  Compartment = factor(
    c(rep("Epithelial cells", each = 4), rep("Stromal cells", each = 3)),
    levels = c("Epithelial cells", "Stromal cells")
  ),
  row.names = colnames(ht_mat)
)

ann_colors <- list(
  Category = c(
    "Main lineages" = "#023047",
    "Epithelial cell states" = "#fb8500",
    "Stromal cell states" = "#219ebc"
  ),
  Compartment = c(
    "Epithelial cells" = "#fb5607",
    "Stromal cells" = "#8338ec"
  )
)

svglite::svglite(
  file = paste0(out_plots, "heatmap_cca_heca_ct_general_nofilt.svg"),
  height = 8, width = 16,
  fix_text_size = FALSE,
  system_fonts = list(sans = "Arial") # map generic 'sans' -> an installed font
)

pheatmap(
  mat = t(ht_mat),
  cluster_rows = FALSE, cluster_cols = FALSE,
  color = c("white", "#f94144"),
  display_numbers = TRUE, number_color = "black",
  cellwidth = 25, cellheight = 25,
  gaps_col = gaps_row,
  gaps_row = 4,
  annotation_col = ann_row,
  annotation_row = ann_col,
  annotation_colors = ann_colors,
  border_color = "grey",
  angle_col = "45"
)

dev.off()


# Separating lineages from specific celltypes



ht_mat <- t(avg_mat[, colSums(avg_mat) > 0.1])
ht_mat <- t(avg_mat[, colSums(avg_mat) >= 0.0])
ht_mat_specific_ct <- ht_mat[str_detect(
  rownames(ht_mat),
  pattern = "_lineage",
  negate = TRUE
), ]


# ---- row groups ----
n <- nrow(ht_mat_specific_ct)
last4_idx <- (n - 3):n

row_group <- rep("Epithelial cell states", n)
row_group[last4_idx] <- "Stromal cell states"


ann_row <- data.frame(
  Category = factor(row_group,
    levels = c("Epithelial cell states", "Stromal cell states")
  ),
  row.names = rownames(ht_mat_specific_ct)
)

# gaps after last lineage row, and before the last 4 rows
gaps_row <- n - 4

# ---- column groups (4 epi, 4 stromal) ----
stopifnot(ncol(ht_mat_specific_ct) == 7) # as per your spec
ann_col <- data.frame(
  Compartment = factor(
    c(
      rep("Epithelial cells", each = 4),
      rep("Stromal cells", each = 3)
    ),
    levels = c("Epithelial cells", "Stromal cells")
  ),
  row.names = colnames(ht_mat_specific_ct)
)

ann_colors <- list(
  Category = c(
    "Epithelial cell states" = "#fb8500",
    "Stromal cell states" = "#219ebc"
  ),
  Compartment = c(
    "Epithelial cells" = "#fb5607",
    "Stromal cells" = "#8338ec"
  )
)

svglite::svglite(
  file = paste0(out_plots, "heatmap_cca_heca_ct_no_lineage_nofilt.svg"),
  height = 12, width = 16,
  fix_text_size = FALSE,
  system_fonts = list(sans = "Arial") # map generic 'sans' -> an installed font
)

pheatmap(
  mat = t(ht_mat_specific_ct),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = c("white", "#f94144"),
  display_numbers = TRUE,
  number_color = "black",
  cellwidth = 25,
  cellheight = 25,
  gaps_row = 4,
  annotation_row = ann_col,
  annotation_colors = ann_colors,
  border_color = "grey",
  angle_col = "45",
  annotation_names_row = FALSE
)

dev.off()

# Only with main lineages

ht_mat_lineages <- t(avg_mat_lineages)

svglite::svglite(
  file = paste0(out_plots, "heatmap_cca_heca_ct_lineage.svg"),
  height = 8,
  width = 8,
  fix_text_size = FALSE,
  system_fonts = list(sans = "Arial") # map generic 'sans' -> an installed font
)

pheatmap(
  mat = t(ht_mat_lineages),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = c("white", "#f94144"),
  display_numbers = TRUE,
  number_color = "black",
  cellwidth = 25,
  cellheight = 25,
  gaps_row = 4,
  annotation_row = ann_col,
  annotation_colors = ann_colors,
  border_color = "grey",
  angle_col = "45",
  annotation_names_row = FALSE
)

dev.off()
