library(Seurat)
library(dittoSeq)
library(clustree)
library(ggplot2)
library(SeuratObject)
library(openxlsx)


# Loading the dataset to query

out_dir <- "output/annotation/marker_annotation/"
out_data <- paste0(out_dir, "data/")
out_plots <- paste0(out_dir, "plots/")

seu <- LoadSeuratRds(file = "output/sample_integration/data/chip/seurat_filtered.rds") # nolint


seu$celltype <- as.character(seu$RNA_snn_res.0.5)


seu$celltype[which(seu$RNA_snn_res.0.5 == "0")] <- "Luminal-like epithelium"
seu$celltype[which(seu$RNA_snn_res.0.5 == "1")] <- "HT epithelium I"
seu$celltype[which(seu$RNA_snn_res.0.5 == "2")] <- "Stromal"
seu$celltype[which(seu$RNA_snn_res.0.5 == "3")] <- "Cycling epithelium"
seu$celltype[which(seu$RNA_snn_res.0.5 == "4")] <- "Luminal-like epithelium"
seu$celltype[which(seu$RNA_snn_res.0.5 == "5")] <- "HT epithelium II"
seu$celltype[which(seu$RNA_snn_res.0.5 == "6")] <- "Decidual stroma"
seu$celltype[which(seu$RNA_snn_res.0.6 == "6")] <- "Decidual stroma"
seu$celltype[which(seu$RNA_snn_res.0.6 == "7")] <- "Non-decidual stroma"
seu$celltype[which(seu$RNA_snn_res.0.5 == "7")] <- "Luminal-like epithelium"
seu$celltype[which(seu$RNA_snn_res.0.5 == "8")] <- "Myofibroblast"
seu$celltype[which(seu$RNA_snn_res.0.5 == "10")] <- "Non-decidual stroma"
seu <- subset(seu, subset = !RNA_snn_res.0.5 == "9")

Idents(seu) <- "celltype"

# finding markers for each cell type
future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = seu) <- "RNA"
sc_celltypes <- FindAllMarkers(object = seu, assay = "RNA", min.pct = 0.3)
future::plan("sequential")

write.csv(
  sc_celltypes,
  file = paste0("output/sample_integration/data/chip/markers_annotation.csv")
)

write.xlsx(
  sc_celltypes,
  "output/sample_integration/data/chip/markers_annotation.xlsx"
)


SaveSeuratRds(
  seu,
  file = paste0(out_data, "seurat_filtered_annoted.rds")
)


# Loading the dataset to query

out_dir <- "output/annotation/marker_annotation_v2/"
out_data <- paste0(out_dir, "data/")
out_plots <- paste0(out_dir, "plots/")

seu_dir <- paste0(
  "output/annotation/marker_annotation/data/",
  "seurat_filtered_annoted.rds"
)

seu <- LoadSeuratRds(file = seu_dir)


seu$celltype_v2 <- as.character(seu$celltype)

seu$celltype_v2[which(seu$celltype == "Stromal")] <- "Stroma"
seu$celltype_v2[which(seu$celltype == "Decidual stroma")] <- "HT stroma"
seu$celltype_v2[which(seu$celltype == "Non-decidual stroma")] <- "HT stroma"





Idents(seu) <- "celltype_v2"

# finding markers for each cell type
future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(object = seu) <- "RNA"
sc_celltypes <- FindAllMarkers(
  object = seu,
  assay = "RNA",
  min.pct = 0.3
)
future::plan("sequential")

write.csv(
  sc_celltypes,
  file = paste0(out_data, "markers_annotation.csv")
)

write.xlsx(sc_celltypes, paste0(out_data, "markers_annotation.xlsx"))

cols <- c(
  "Cycling epithelium" = "#F8766D",
  "HT epithelium I" = "#7CAE00",
  "HT epithelium II" = "#00BE67",
  "Luminal-like epithelium" = "#00BFC4",
  "Myofibroblast" = "#00A9FF",
  "HT stroma" = "#C77CFF",
  "Stroma" = "#FF61CC"
)

DimPlot(
  seu,
  group.by = "celltype_v2",
  reduction = "umap.chip.harmony",
  cols = cols
) + theme(legend.position = "none")


ggsave(
  paste0(out_plots, "umap_celltypes.png"),
  width = 7,
  height = 7,
  dpi = 400
)
ggsave(
  paste0(out_plots, "umap_celltypes.svg"),
  width = 7,
  height = 7,
  device = "svg",
  dpi = 400,
  fix_text_size = FALSE
)


SaveSeuratRds(
  seu,
  file = paste0(out_data, "seurat_filtered_annoted.rds")
)

seu <- LoadSeuratRds(
  file = paste0(out_data, "seurat_filtered_annoted.rds")
)


DimPlot(
  seu,
  group.by = "celltype_v2",
  reduction = "umap.chip.harmony"
) +
  DimPlot(
    seu,
    group.by = "celltype_v2",
    reduction = "umap.harmony"
  )

DimPlot(
  seu,
  group.by = "celltype_v2",
  reduction = "umap.chip.harmony",
  split.by = "group"
) +
  DimPlot(
    seu,
    group.by = "celltype_v2",
    reduction = "umap.harmony",
    split.by = "group"
  )

DimPlot(
  seu,
  reduction = "umap.chip.harmony",
  cols = c("CNT" = "#0f99b2", "HT" = "#ce0665"),
  group.by = "group"
) +
  theme(legend.position = "none")
ggsave(
  filename = paste0(out_plots, "umap_group.svg"),
  device = "svg",
  dpi = 400,
  height = 7,
  width = 7
)


# Subsets with a more general annotation
unique(seu$celltype_general)
seu_ss <- subset(
  seu,
  subset = celltype_general == "Epithelium HT" | celltype_general == "Epithelium No-HT" # nolint
)

seu_ss$celltype_general <- factor(
  seu_ss$celltype_general,
  levels = c("Epithelium No-HT", "Epithelium HT")
)
Idents(seu_ss) <- "celltype_general"

VlnPlot(seu_ss,
  features = c("PAEP", "SPP1"), pt.size = 0,
  cols = c(
    "Epithelium No-HT" = "#0f99b2",
    "Epithelium HT" = "#ce0665"
  )
)
ggsave(
  filename = paste0(
    out_plots, "vln_paep_spp1_epi_groups.svg"
  ),
  device = "svg",
  dpi = 400,
  height = 7,
  width = 8
)


## BOXPLOT OF GENES

# Extract expression data for PAEP and SPP1
df <- FetchData(
  seu_ss,
  vars = c("PAEP", "SPP1", "celltype_general")
)

# Convert to long format for ggplot
df_long <- df %>%
  pivot_longer(
    cols = c("PAEP", "SPP1"),
    names_to = "gene",
    values_to = "expression"
  )

# Make boxplot
p <- ggplot(
  df_long,
  aes(
    x = celltype_general,
    y = expression,
    fill = celltype_general
  )
) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~gene, scales = "free_y") +
  scale_fill_manual(values = c(
    "Epithelium No-HT" = "#0f99b2",
    "Epithelium HT" = "#ce0665"
  )) +
  theme_minimal(base_size = 14) +
  labs(x = "Cell type", y = "Expression")

# Save
ggsave(
  paste0(out_plots, "boxplot_paep_spp1_epi_groups.svg"),
  p,
  device = "svg",
  dpi = 400,
  height = 7,
  width = 8,
  bg = "white"
)

dittoBarPlot(
  seu,
  group.by = "sample_id",
  var = "celltype_v2",
  scale = "count",
  color.panel = cols,
  x.reorder = c(1, 3, 5, 2, 4, 6)
)
ggsave(
  filename = paste0(out_plots, "barplot_sample_celltype.svg"),
  device = "svg",
  dpi = 400,
  bg = "white",
  width = 5,
  height = 4,
  fix_text_size = FALSE
)
ggsave(
  filename = paste0(out_plots, "barplot_sample_celltype.png"),
  device = "png",
  dpi = 400,
  bg = "white",
  width = 5,
  height = 4
)
