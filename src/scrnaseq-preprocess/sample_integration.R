# Title: single cell exploratory analysis & integration
# Author: Diego Amoros & Jaime Llera
# Date: 15-10-2024
# Inputs:
#   - CellBender H5 results files, h5
#   - metadata, csv
# Outputs:
#   - merged dataset, .rds
#   - QC plots, jpg


# libraries and utils and constants
suppressPackageStartupMessages({
  library(patchwork)
  library("Seurat")
  library("scCustomize")
  library("knitr")
  library("DT")
  library("stringr")
  library("RColorBrewer")
  library("harmony")
  library("ggpubr")
  library("grid")
  library("ggplot2")
  library("ggrepel")
  library("clustree")
  library("dittoSeq")
  library("scuttle")
  library("dplyr")
})

source("src/scrnaseq-preprocess/utils_sample_integration.R")

# set.seed
set.seed(47)

## Packages and versions
paste0("Cell Ranger 9.1.0")
paste0("Seurat ", packageVersion("Seurat"))
paste0("harmony ", packageVersion("harmony"))


## Folders
dir_in <- "output/qc/seurat_qc/"
dir_out <- "output/sample_integration"
dir_plots <- paste0(dir_out, "/plots/")
dir_tabl <- paste0(dir_out, "/data/")

dir.create(
  path = dir_out,
  showWarnings = FALSE
)
dir.create(
  path = dir_plots,
  showWarnings = FALSE
)
dir.create(
  path = dir_tabl,
  showWarnings = FALSE
)



## Merge the data
### RDSSeurat files to load
files <- list.files(
  path = dir_in,
  full.names = TRUE,
  pattern = "\\.rds$"
)
files_list <- as.list(files)

### Define Sample Names
names(files_list) <- gsub(
  pattern = "-hq-cells.rds",
  replacement = "",
  basename(files)
)

### Load RDSSeurat files into list of Seurat Objects
sc_multiple_list <- lapply(files_list, LoadSeuratRds)

### Combine to unique seurat object
sc_multiple <- merge(
  x = sc_multiple_list[[1]],
  y = sc_multiple_list[2:length(sc_multiple_list)],
  add.cell.ids = names(sc_multiple_list)
)
rm(sc_multiple_list)
gc(verbose = FALSE, full = TRUE)

## Add metadata
sc_multiple@meta.data$Sample <- sc_multiple@meta.data$sample_id

## Extract group_tissue Groups
sc_multiple@meta.data$group <- substr(
  sc_multiple@meta.data$sample_id,
  start = 8,
  stop = 10
)

# Save combined object
# ========================
SaveSeuratRds(
  sc_multiple,
  paste0(
    dir_out,
    "/seurat_merged.rds"
  )
)

gc(verbose = FALSE, full = TRUE)

## MAD and min.diff stats
outlier_mito <- scuttle::isOutlier(
  sc_multiple$mitoRatio,
  type = "higher",
  nmads = 3,
  batch = sc_multiple$sample_id
)

sc_multiple <- AddMetaData(
  sc_multiple,
  metadata = outlier_mito,
  col.name = "isoutlier_mito"
)

## PLOT
df <- sc_multiple@meta.data %>%
  mutate(
    outlier_color = ifelse(isoutlier_mito, "red", "black")
  )

p1_outlayers <- ggplot(
  df,
  aes(x = sample_id, y = mitoRatio)
) +
  geom_boxplot(
    fill = "lightblue",
    color = "black",
    alpha = 0.7
  ) + # Boxplot
  geom_jitter(
    aes(color = outlier_color),
    width = 0.2,
    size = 0.5,
    alpha = 0.6
  ) +
  scale_color_identity() +
  labs(
    title = "mitoRatio by protocol",
    y = "mitoRatio",
    x = "Protocol"
  ) +
  theme_minimal()

p1_outlayers

ggsave(
  plot = p1_outlayers,
  filename = paste0(dir_plots, "/boxplot_outlayers.jpg"),
  width = 8,
  height = 6,
  device = "png",
  dpi = 100,
  bg = "white"
)


p1_outlayers5 <- ggplot(df, aes(x = sample_id, y = mitoRatio)) +
  geom_boxplot(
    fill = "lightblue",
    color = "black",
    alpha = 0.7
  ) + # Boxplot
  geom_jitter(
    aes(
      color = outlier_color
    ),
    width = 0.2,
    size = 0.1,
    alpha = 0.6
  ) +
  scale_color_identity() +
  labs(
    title = "mitoRatio",
    y = "mitoRatio",
    x = "Protocol"
  ) +
  theme_minimal() +
  ylim(c(0, 25))

ggsave(
  plot = p1_outlayers5,
  filename = paste0(
    dir_plots, "/boxplot_outlayers_nocolor.png"
  ),
  width = 8,
  height = 6,
  device = "png",
  dpi = 100,
  bg = "white"
)




# Filters already applied at 500 features, 15% mito and 1000 counts


# Normalize
# ========================
sc_multiple <- NormalizeData(
  sc_multiple,
  verbose = TRUE,
  assay = "RNA",
  normalization.method = "LogNormalize",
  scale.factor = 10000
)


# Find Variable features
# ========================
# To create lists of mitochondrial, ribosomal and cell cycle and Hemoglobin
#  genes to remove them previous to High Variable Fetures (HVF) calculation

mito_genes <- row.names(
  x = sc_multiple
)[grep(
  pattern = "^MT-",
  x = row.names(sc_multiple)
)]
ribo_genes <- grep(
  pattern = "^RP[SL][[:digit:]]|^RP[[:digit:]]|^RPSA",
  rownames(sc_multiple),
  value = TRUE
)
cc_genes <- c(
  "UBR7", "RFC2", "RAD51", "MCM2", "TIPIN", "MCM6", "UNG", "POLD3", "WDR76",
  "CLSPN", "CDC45", "CDC6", "MSH2", "MCM5", "POLA1", "MCM4", "RAD51AP1",
  "GMNN", "RPA2", "CASP8AP2", "HELLS", "E2F8", "GINS2", "PCNA", "NASP",
  "BRIP1", "DSCC1", "DTL", "CDCA7", "CENPU", "ATAD2", "CHAF1B", "USP1",
  "SLBP", "RRM1", "FEN1", "RRM2", "EXO1", "CCNE2", "BLM", "PRIM1", "UHRF1",
  "NCAPD2", "ANLN", "TACC3", "HMMR", "GTSE1", "NDC80", "AURKA", "TPX2",
  "BIRC5", "G2E3", "CBX5", "RANGAP1", "CTCF", "CDCA3", "TTK", "SMC4",
  "ECT2", "CENPA", "CDC20", "NEK2", "CENPF", "HJURP", "CKS2", "DLGAP5",
  "PIMREG", "TOP2A", "PSRC1", "CDCA8", "CKAP2", "NUSAP1", "KIF23", "KIF11",
  "KIF20B", "CENPE", "GAS2L3", "KIF2C", "NUF2", "ANP32E", "LBR", "MKI67",
  "CCNB2", "CDC25C", "HMGB2", "CKAP2L", "BUB1", "CDK1", "CKS1B", "UBE2C",
  "CKAP5", "AURKB", "CDCA2", "TUBB4B", "JPT1"
)
hbb_genes <- c("HBA1", "HBA2", "HBB")

genes_to_include_hvf <- setdiff(
  rownames(sc_multiple),
  c(mito_genes, ribo_genes, cc_genes, hbb_genes)
)


# To find variable features (4000) of the whole seurat object

sc_multiple <- FindVariableFeatures(
  sc_multiple,
  selection.method = "vst",
  nfeatures = 4000,
  verbose = TRUE
)


## Remove mito genes, ribo genes, cell cycle genes and HBB genes
## from variable features
variable_features_filt <- intersect(
  VariableFeatures(sc_multiple),
  genes_to_include_hvf
)
length(variable_features_filt)

## Calculate cycle score
# Temporary merged object for dataset cell cycle calculation
sc_multiple_join <- JoinLayers(sc_multiple)

# Scoring the cell cycle with sets from seurat
s_genes <- Seurat::cc.genes$s.genes
g2m_genes <- Seurat::cc.genes$g2m.genes

sc_multiple_join <- CellCycleScoring(
  sc_multiple_join,
  s.features = s_genes,
  g2m.features = g2m_genes
)

cycleScore <- data.frame(
  Phase = sc_multiple_join$Phase,
  G2M.Score = sc_multiple_join$G2M.Score,
  S.score = sc_multiple_join$S.Score
)

cycleScore <- cycleScore[colnames(sc_multiple), ]

sc_multiple <- AddMetaData(
  sc_multiple,
  metadata = cycleScore
)

# remove temporary joined object
rm(sc_multiple_join)
gc()


# Scale Data and dim reduction
# ========================
# set.seed
set.seed(47)

sc_multiple <- ScaleData(
  object = sc_multiple,
  features = variable_features_filt,
  verbose = TRUE
)

sc_multiple <- RunPCA(sc_multiple,
  features = variable_features_filt,
  reduction.name = "pca.unintegrated",
  verbose = TRUE
)

sc_multiple <- RunUMAP(sc_multiple,
  dims = 1:30,
  reduction = "pca.unintegrated",
  reduction.name = "umap.unintegrated",
  verbose = TRUE
)

## Input Atlas UMAP by celltype
umap_unintegrated <- DimPlot(
  sc_multiple,
  reduction = "umap.unintegrated",
  group.by = "sample_id", raster = FALSE
) +
  theme_classic(base_size = 10) +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("Unintegrated dataset")

umap_unintegrated

ggsave(
  plot = umap_unintegrated,
  filename = paste0(dir_plots, "/plot_umap_unintegrated.png"),
  width = 10,
  height = 6,
  device = "png",
  dpi = 400,
  bg = "white"
)


umap_unintegrated_splitted <- DimPlot(
  sc_multiple,
  reduction = "umap.unintegrated",
  group.by = "sample_id",
  split.by = "group",
  raster = FALSE
) +
  theme_classic(base_size = 8) +
  labs(x = "UMAP1", y = "UMAP2")

umap_unintegrated_splitted

ggsave(
  plot = umap_unintegrated_splitted,
  filename = paste0(dir_plots, "/plot_umap_unintegrated_sp.png"),
  width = 14,
  height = 5,
  device = "png",
  dpi = 400,
  bg = "white"
)


# Integration
# #####################################
# set.seed
set.seed(47)
gc()

# Integration by sample
sc_multiple <- IntegrateLayers(
  object = sc_multiple,
  method = HarmonyIntegration,
  orig.reduction = "pca.unintegrated",
  new.reduction = "pca.harmony",
  features = variable_features_filt,
  verbose = TRUE,
  theta = 2
)

sc_multiple <- RunUMAP(sc_multiple,
  dims = 1:30,
  reduction = "pca.harmony",
  reduction.name = "umap.harmony",
  verbose = TRUE
)

## Input Atlas UMAP by celltype
umap_harmony <- DimPlot(
  sc_multiple,
  reduction = "umap.harmony",
  group.by = "sample_id", raster = FALSE
) +
  theme_classic(base_size = 10) +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("Harmony by sample integrated dataset")

umap_harmony

ggsave(
  plot = umap_harmony,
  filename = paste0(dir_plots, "/plot_umap_harmony.png"),
  width = 6,
  height = 6,
  device = "png",
  dpi = 400,
  bg = "white"
)


umap_harmony_splitted <- DimPlot(
  sc_multiple,
  reduction = "umap.harmony",
  group.by = "sample_id",
  split.by = "group",
  raster = FALSE
) +
  theme_classic(base_size = 8) +
  labs(x = "UMAP1", y = "UMAP2")

umap_harmony_splitted

ggsave(
  plot = umap_harmony_splitted,
  filename = paste0(dir_plots, "/plot_umap_harmony_sp.png"),
  width = 14,
  height = 7,
  device = "png",
  dpi = 400,
  bg = "white"
)



# Plotting metadata and genes
feat_plot <- plot_features(sc_multiple, "umap.harmony")

ggsave(
  plot = feat_plot,
  filename = paste0(dir_plots, "/plot_markers.harmony_beforeFilt.jpg"),
  width = 14,
  height = 24,
  device = "jpg",
  dpi = 150,
  bg = "white"
)

# Save integrated object
# ========================
SaveSeuratRds(sc_multiple, paste0(dir_out, "/integration_sample_id.rds"))
gc(verbose = FALSE, full = TRUE)


# First annotation
############################################
# FindClusters 0.2 (resolution parameter)
sc_multiple$SampleIdents <- Idents(sc_multiple)
sc_multiple <- FindNeighbors(
  sc_multiple,
  reduction = "pca.harmony",
  dims = 1:30,
  verbose = TRUE
)

res <- 0.2
sc_multiple <- FindClusters(
  sc_multiple,
  resolution = res,
  verbose = TRUE
)

p6 <- DimPlot(
  sc_multiple,
  reduction = "umap.harmony",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) + NoLegend()

p6

## Gene markers of each cluster
###############################

sc_multiple_join <- JoinLayers(sc_multiple)

future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
DefaultAssay(sc_multiple_join) <- "RNA"
all_markers <- FindAllMarkers(
  sc_multiple_join,
  min.pct = 0.25
)
future::plan("sequential")

all_markers <- all_markers[, c(7, 1:6)]
write.table(all_markers,
  file = paste0(dir_tabl, "/markers_res02_preliminary.csv"),
  row.names = FALSE, quote = FALSE, sep = ","
)

rm(sc_multiple_join)
gc(full = TRUE, reset = TRUE)


# Clustree
########################################
sc_multiple <- JoinLayers(sc_multiple)

n_clus <- c()
n_res <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
for (res in n_res) {
  sc_multiple <- FindClusters(
    sc_multiple,
    resolution = res,
    verbose = FALSE,
    random.seed = 42
  )
  nclusters <- eval(parse(text = paste0(
    "max(as.numeric(unique(sc_multiple$RNA_snn_res.", res, ")))"
  )))
  n_clus <- c(n_clus, nclusters)
  print(paste0(
    "Find Clusters with ",
    res,
    " resolution: clusters Detected ", nclusters
  ))
}


clustreeplot <- clustree(sc_multiple, prefix = "RNA_snn_res.")
print(clustreeplot)

ggsave(
  plot = clustreeplot,
  filename = paste0(dir_plots, "/clustreeplot_harmony.png"),
  width = 8,
  height = 10,
  device = "png",
  dpi = 400,
  bg = "white"
)


plot_list <- list()
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
for (res in resolutions) {
  plot_list[[paste0("res_", res)]] <- DimPlot(
    sc_multiple,
    reduction = "umap.harmony",
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    group.by = paste0(paste0("RNA_snn_res.", res))
  ) +
    ggtitle(paste0("Res:", res)) + theme_classic(base_size = 8) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
}

ggsave(
  plot =
    plot_list[["res_0.1"]] + plot_list[["res_0.2"]] + plot_list[["res_0.3"]] + plot_list[["res_0.4"]], # nolint
  filename = paste0(dir_plots, "/clustree_A.png"),
  width = 12,
  height = 10,
  device = "png",
  dpi = 400,
  bg = "white"
)

ggsave(
  plot =
    plot_list[["res_0.5"]] + plot_list[["res_0.6"]] + plot_list[["res_0.7"]] + plot_list[["res_0.8"]], # nolint
  filename = paste0(dir_plots, "/clustree_B.png"),
  width = 12,
  height = 10,
  device = "png",
  dpi = 400,
  bg = "white"
)



gc(verbose = FALSE, full = TRUE)
# Second annotation
############################################
# FindClusters 0.3 (resolution parameter)

res <- 0.3
sc_multiple <- FindClusters(
  sc_multiple,
  resolution = res,
  verbose = TRUE
)

umap_res03 <- DimPlot(
  sc_multiple,
  reduction = "umap.harmony",
  label = TRUE,
  repel = TRUE,
  raster = TRUE,
  cols = colorRampPalette(brewer.pal(8, "Set2"))(14)
)

ggsave(
  plot = umap_res03,
  filename = paste0(dir_plots, "/umap_cluster03.png"),
  width = 6,
  height = 5,
  device = "png",
  dpi = 150,
  bg = "white"
)

umap_res03_split <- DimPlot(
  sc_multiple,
  reduction = "umap.harmony",
  split.by = "group",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) + labs(x = "UMAP1", y = "UMAP2")

ggsave(
  plot = umap_res03_split,
  filename = paste0(dir_plots, "/umap_res03_split.png"),
  width = 10,
  height = 5,
  device = "png",
  dpi = 400,
  bg = "white"
)



HLplots <- list()
# for (clus in 0:length(unique(sc_multiple$seurat_clusters)) - 1) {
for (i in seq_along(unique(sc_multiple$seurat_clusters))) {
  CellstoHighlight <- WhichCells(sc_multiple, idents = i - 1)
  HLplots[[i]] <- DimPlot(
    sc_multiple,
    reduction = "umap.harmony",
    label = TRUE,
    repel = TRUE,
    raster = TRUE,
    group.by = "seurat_clusters",
    cells.highlight = CellstoHighlight,
    cols.highlight = "red",
    cols = "grey90"
  ) +
    theme_classic(base_size = 8) +
    labs(x = "UMAP1", y = "UMAP2") +
    ggtitle(
      paste0(
        "Cluster ", i - 1, ", ", length(CellstoHighlight), " cells"
      )
    ) +
    NoLegend()
}

umap_harmony_clusters_highlight <- ggarrange(plotlist = HLplots)


ggsave(
  plot = umap_harmony_clusters_highlight,
  filename = paste0(dir_plots, "/umap_cluster03_highlight.png"),
  width = 13,
  height = 10,
  device = "png",
  dpi = 400,
  bg = "white"
)




## Get all cluster markers

DefaultAssay(sc_multiple) <- "RNA"

future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
all_markers_03 <- FindAllMarkers(sc_multiple, min.pct = 0.25)

future::plan("sequential")


write.table(all_markers_03[, c(7, 1:6)],
  file = paste0(dir_tabl, "/markers_res03_preliminary.csv"),
  row.names = FALSE,
  quote = FALSE,
  sep = ","
)
all_markers_03 <- read.table(
  file = paste0(dir_tabl, "/markers_res03_preliminary.csv"),
  sep = ",",
  col.names = 1,
  row.names = NULL
)
# Save integrated object
# ========================
SaveSeuratRds(
  sc_multiple,
  paste0(dir_out, "/integration_sample_id.rds")
)

sc_multiple <- LoadSeuratRds(paste0(dir_out, "/integration_sample_id.rds"))

gc(verbose = FALSE, full = TRUE)



## SOME PLOTS FOR REPORTS
# Extract metadata
metadata <- sc_multiple@meta.data
cell_counts <- as.data.frame(table(metadata$Sample, metadata$group))
colnames(cell_counts) <- c("Sample", "Group", "CellCount")


t_barplot <- ggplot(cell_counts, aes(x = Sample, y = CellCount, fill = Group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("coral1", "royalblue1")) +
  theme_minimal() +
  labs(
    title = "Number of Cells by Sample and Group",
    x = "Sample",
    y = "Number of Cells",
    fill = "Group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

t_barplot

ggsave(
  plot = t_barplot,
  filename = paste0(dir_plots, "/sample_group_barplot.png"),
  width = 8,
  height = 4,
  device = "png",
  dpi = 100,
  bg = "white"
)





# Count the number of cells per sample
cell_counts <- table(metadata$sample_id)
cell_counts_df <- as.data.frame(cell_counts)
colnames(cell_counts_df) <- c("Sample", "CellCount")

# Calculate the average number of cells
average_cells <- mean(cell_counts_df$CellCount)

# Create the barplot and add the average line
s_barplot <- ggplot(cell_counts_df, aes(x = Sample, y = CellCount)) +
  geom_bar(stat = "identity", fill = "royalblue1") +
  geom_hline(yintercept = average_cells, linetype = "dashed", color = "red") +
  theme_minimal() +
  labs(
    title = "Number of Cells per Sample",
    x = "Sample",
    y = "Number of Cells",
    caption = paste("Average number of cells:", round(average_cells))
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot
ggsave(
  plot = s_barplot,
  filename = paste0(dir_plots, "/sample_nocond_barplot.png"),
  width = 14,
  height = 5,
  device = "png",
  dpi = 100,
  bg = "white"
)


## Markers plot

df <- read.table(
  paste0(dir_tabl, "/markers_res03_preliminary.csv"),
  sep = ",", header = TRUE
)

df <- df[df$p_val_adj < 0.05, ]

# Add a new column to classify genes as positive or negative based on avg_log2FC
df <- df %>%
  mutate(gene_class = ifelse(avg_log2FC >= 0, "positive", "negative"))

# Count the number of positive and negative genes per cluster
df_summary <- df %>%
  group_by(cluster, gene_class) %>%
  summarise(count = n(), .groups = "drop")

# Create the bar plot
nmarkers_barplot <- ggplot(
  df_summary,
  aes(x = factor(cluster), y = count, fill = gene_class)
) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Cluster", y = "Number of DE Genes", fill = "Gene Class") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("negative" = "blue", "positive" = "red"))

nmarkers_barplot

ggsave(
  plot = nmarkers_barplot,
  filename = paste0(dir_plots, "/nMarkers_barplot.png"),
  width = 6,
  height = 3,
  device = "png",
  dpi = 400,
  bg = "white"
)



## Cells by cluster

# Extract the cluster assignments
cluster_counts <- table(sc_multiple$seurat_clusters)

# Convert to data frame for plotting
cluster_df <- as.data.frame(cluster_counts)
colnames(cluster_df) <- c("Cluster", "CellCount")

# Create the barplot
ncells_barplot <- ggplot(cluster_df, aes(x = Cluster, y = CellCount)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Cluster",
    y = "Number of Cells",
    title = "Number of Cells in Each Cluster"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_viridis_c() +
  scale_y_continuous(breaks = seq(0, max(cluster_df$CellCount), by = 10000))



ggsave(
  plot = ncells_barplot,
  filename = paste0(dir_plots, "/nCells_barplot.png"),
  width = 6,
  height = 3,
  device = "png",
  dpi = 400,
  bg = "white"
)


# Some configurations of plots:
layout <- "AABB
           #CC#"


VlnPlot(
  sc_multiple,
  features = c("nFeature_RNA", "mitoRatio"),
  split.by = "RNA_snn_res.0.3",
  pt.size = 0
) +
  DimPlot(
    sc_multiple,
    group.by = "RNA_snn_res.0.3",
    reduction = "umap.harmony",
    label = TRUE
  ) +
  patchwork::plot_layout(design = layout)

ggsave(
  filename = paste0(dir_plots, "vln_qc_cluster_03.png"),
  device = "png",
  dpi = 400,
  bg = "white"
)

top_5 <- all_markers_03 %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)
DotPlot(
  sc_multiple,
  features = unique(top_5$gene),
  group.by = "RNA_snn_res.0.3"
) + coord_flip()
ggsave(
  filename = paste0(dir_plots, "dotplot_cluster_03.png"),
  device = "png",
  dpi = 400,
  height = 8,
  width = 6,
  bg = "white"
)

DotPlot(
  sc_multiple,
  features = unique(top_5$gene),
  group.by = "RNA_snn_res.0.3"
) + coord_flip() +
  DimPlot(
    sc_multiple,
    group.by = "RNA_snn_res.0.3",
    reduction = "umap.harmony",
    label = TRUE
  ) +
  patchwork::plot_layout(design = "ABB\nABB")

ggsave(
  filename = paste0(
    dir_plots,
    "markers_cluster_dotplot_umap_res03.png"
  ),
  device = "png",
  dpi = 400,
  width = 12,
  height = 6.2,
  scale = 1.3,
  bg = "white"
)


# Eliminating clusters that provide little to no information
# (low quality cells remaining)

sc_filt <- sc_multiple[
  , -which(
    sc_multiple$RNA_snn_res.0.3 == "2" |
      sc_multiple$RNA_snn_res.0.3 == "5" |
      sc_multiple$RNA_snn_res.0.3 == "6"
  )
]


# Integration without removing HT
# #####################################
# set.seed
set.seed(47)
gc()
sc_filt_join <- JoinLayers(sc_filt)

# Split dataset
sc_filt_join@meta.data[["chip"]] <-
  stringr::str_extract(sc_filt_join$sample_id, pattern = "[0-9]{2}")

sc_filt[["RNA"]] <- split(
  sc_filt_join[["RNA"]],
  sc_filt_join$chip
)

sc_filt <- NormalizeData(sc_filt)


sc_filt <- FindVariableFeatures(
  sc_filt,
  selection.method = "vst",
  nfeatures = 4000,
  verbose = TRUE
)


## Remove mito genes, ribo genes, cell cycle genes and HBB genes
## from variable features
variable_features_filt_chip <- intersect(
  VariableFeatures(sc_filt),
  genes_to_include_hvf
)
length(variable_features_filt_chip)

set.seed(47)

sc_filt <- ScaleData(
  object = sc_filt,
  features = variable_features_filt_chip,
  verbose = TRUE
)

sc_filt <- RunPCA(sc_filt,
  features = variable_features_filt_chip,
  reduction.name = "pca.chip.unintegrated",
  verbose = TRUE
)

# Integration by sample
sc_filt <- IntegrateLayers(
  object = sc_filt,
  method = HarmonyIntegration,
  orig.reduction = "pca.chip.unintegrated",
  new.reduction = "pca.chip.harmony",
  features = variable_features_filt_chip,
  verbose = TRUE,
  theta = 2
)

sc_filt <- RunUMAP(sc_filt,
  dims = 1:30,
  reduction = "pca.chip.harmony",
  reduction.name = "umap.chip.harmony",
  verbose = TRUE
)

## Input Atlas UMAP by celltype
umap_chip_harmony <- DimPlot(
  sc_filt,
  reduction = "umap.chip.harmony",
  group.by = "sample_id", raster = FALSE
) +
  theme_classic(base_size = 10) +
  labs(x = "UMAP1", y = "UMAP2") +
  ggtitle("Harmony by chip integrated dataset")

umap_chip_harmony

ggsave(
  plot = umap_chip_harmony,
  filename = paste0(dir_plots, "/plot_umap_chip_harmony.png"),
  width = 6,
  height = 6,
  device = "png",
  dpi = 400,
  bg = "white"
)


umap_harmony_splitted <- DimPlot(
  sc_filt,
  reduction = "umap.chip.harmony",
  group.by = "sample_id",
  split.by = "group",
  raster = FALSE
) +
  theme_classic(base_size = 8) +
  labs(x = "UMAP1", y = "UMAP2")

umap_harmony_splitted

ggsave(
  plot = umap_harmony_splitted,
  filename = paste0(dir_plots, "/plot_umap_chip_harmony_sp.png"),
  width = 14,
  height = 7,
  device = "png",
  dpi = 400,
  bg = "white"
)

# Add chip metadata
sc_filt$chip <- str_extract(string = sc_filt$sample_id, pattern = "[0-9]{2}")

dittoBarPlot(sc_filt, var = "chip", group.by = "group")
ggsave(
  filename = paste0(
    dir_plots, "/chip/barplot_chip_group.png"
  ),
  device = "png",
  width = 3,
  height = 5,
  dpi = 400
)

# Clustree
########################################
sc_filt <- JoinLayers(sc_filt)

sc_filt <- FindNeighbors(
  sc_filt,
  reduction = "pca.chip.harmony",
  features = variable_features_filt_chip,
  dims = 1:30,
  verbose = TRUE
)

n_clus <- c()
n_res <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0)
for (res in n_res) {
  sc_filt <- FindClusters(
    sc_filt,
    resolution = res,
    verbose = FALSE,
    random.seed = 42
  )
  nclusters <- eval(parse(text = paste0(
    "max(as.numeric(unique(sc_filt$RNA_snn_res.", res, ")))"
  )))
  n_clus <- c(n_clus, nclusters)
  print(paste0(
    "Find Clusters with ",
    res,
    " resolution: clusters Detected ", nclusters
  ))
}


clustreeplot <- clustree(sc_filt, prefix = "RNA_snn_res.")
print(clustreeplot)

ggsave(
  plot = clustreeplot,
  filename = paste0(dir_plots, "/chip/clustreeplot_chip_harmony.png"),
  width = 8,
  height = 10,
  device = "png",
  dpi = 400,
  bg = "white"
)


plot_list <- list()
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
for (res in resolutions) {
  plot_list[[paste0("res_", res)]] <- DimPlot(
    sc_filt,
    reduction = "umap.chip.harmony",
    label = TRUE,
    repel = TRUE,
    raster = FALSE,
    group.by = paste0(paste0("RNA_snn_res.", res))
  ) +
    ggtitle(paste0("Res:", res)) + theme_classic(base_size = 8) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 2))
}

ggsave(
  plot =
    plot_list[["res_0.1"]] + plot_list[["res_0.2"]] + plot_list[["res_0.3"]] + plot_list[["res_0.4"]], # nolint
  filename = paste0(dir_plots, "/chip/clustree_A.png"),
  width = 12,
  height = 10,
  device = "png",
  dpi = 400,
  bg = "white"
)

ggsave(
  plot =
    plot_list[["res_0.5"]] + plot_list[["res_0.6"]] + plot_list[["res_0.7"]] + plot_list[["res_0.8"]], # nolint
  filename = paste0(dir_plots, "/chip/clustree_B.png"),
  width = 12,
  height = 10,
  device = "png",
  dpi = 400,
  bg = "white"
)


DefaultAssay(sc_filt) <- "RNA"

sc_filt <- SetIdent(sc_filt, value = "RNA_snn_res.0.5")

Idents(sc_filt)

future::plan("multisession", workers = 8)
options(future.globals.maxSize = 4000 * 1024^2)
allmarkers_05 <- FindAllMarkers(sc_filt, min.pct = 0.25)
future::plan("sequential")


write.table(allmarkers_05[, c(7, 1:6)],
  file = paste0(dir_tabl, "/chip/markers_res05.csv"),
  row.names = FALSE,
  quote = FALSE,
  sep = ","
)

openxlsx::write.xlsx(
  allmarkers_05,
  file = paste0(dir_tabl, "/chip/markers_res05.xlsx")
)

dittoBarPlot(
  sc_filt,
  group.by = "RNA_snn_res.0.5",
  var = "sample_id",
  split.by = "group"
)
ggsave(filename = paste0(dir_plots, "chip/barplot_group_sample.png"), dpi = 400)

SaveSeuratRds(sc_filt, file = paste0(dir_tabl, "chip/seurat_filtered.rds"))
sc_filt <- LoadSeuratRds(file = paste0(dir_tabl, "chip/seurat_filtered.rds"))

top_5_05 <- allmarkers_05 %>%
  group_by(cluster) %>%
  filter(p_val_adj < 0.05) %>%
  slice_max(order_by = avg_log2FC, n = 5)

DotPlot(
  sc_filt,
  features = top_5_05$gene,
  cols = "Spectral"
) +
  coord_flip() +
  DimPlot(
    sc_filt,
    group.by = "RNA_snn_res.0.5",
    reduction = "umap.chip.harmony",
    pt.size = 1,
    label = TRUE
  ) +
  plot_layout(
    design = "ABBB"
  )

ggsave(
  filename = paste0(dir_plots, "/chip/markers_umap_chip.png"),
  width = 12,
  height = 10
)

DimPlot(
  sc_filt,
  group.by = "RNA_snn_res.0.5",
  split.by = "group",
  reduction = "umap.chip.harmony"
)
ggsave(
  filename = paste0(dir_plots, "/chip/umap_chip_split.png"),
  width = 12,
  height = 10
)

FeaturePlot(
  sc_filt,
  features = c(
    "EPCAM",
    "KRT8",
    "DCN",
    "COL1A1",
    "LUM",
    "PRL",
    "ESR1",
    "PRL",
    "FOXO1"
  ),
  min.cutoff = "q5",
  max.cutoff = "q95",
  reduction = "umap.chip.harmony"
)
