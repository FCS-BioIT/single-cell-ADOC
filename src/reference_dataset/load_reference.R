library(Seurat)
library(dittoSeq)
library(harmony)

reference <- LoadSeuratRds(
  "references/sc-rnaseq/endometriumAtlasV2_cells_with_counts.rds"
)


## Keep controls only
reference <- reference[, which(reference$Group != "Endo_Superficial")]
reference <- reference[, which(reference$Hormonal.treatment == "nan")]

reference

# Processing the healthy normalized data
reference <- reference %>%
  NormalizeData(
    verbose = TRUE,
    assay = "RNA",
    normalization.method = "LogNormalize",
    scale.factor = 10000
  ) %>%
  FindVariableFeatures(
    selection.method = "vst",
    nfeatures = 2000,
    verbose = TRUE
  )


reference <- ScaleData(
  reference,
  features = VariableFeatures(object = reference),
  verbose = TRUE
)

reference <- RunPCA(
  reference,
  features = VariableFeatures(reference),
  verbose = TRUE
)


reference <- RunHarmony(
  reference,
  group.by.vars = c("dataset", "sample"),
  theta = c(2, 2),
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony.dataset.sample"
)


reference <- RunUMAP(
  reference,
  dims = 1:30,
  reduction = "harmony.dataset.sample",
  reduction.name = "umap.harmony.dataset.sample"
)

DimPlot(
  reference,
  reduction = "umap.harmony.dataset.sample",
  group.by = c("dataset"), raster = FALSE
)

reference <- FindNeighbors(reference, reduction = "harmony", dims = 1:30)
reference <- FindClusters(
  reference,
  resolution = 2,
  cluster.name = "harmony_clusters"
)


SaveSeuratRds(
  reference,
  file = paste0(
    "references/sc-rnaseq/",
    "endometriumAtlas_cells_with_counts_integrated.rds"
  )
)
reference <- LoadSeuratRds(
  file = paste0(
    "references/sc-rnaseq/",
    "endometriumAtlas_cells_with_counts_integrated.rds"
  )
)

# Loading the dataset to query
sc <- LoadSeuratRds(
  file = "output/annotation/marker_annotation/data/seurat_filtered_annoted.rds"
) # nolint
