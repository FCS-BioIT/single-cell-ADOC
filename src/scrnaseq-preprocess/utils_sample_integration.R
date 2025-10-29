library("Seurat")
library("ggpubr")
library("patchwork")
library("Matrix")
library("ggplot2")
library("ggpubr")
## Utils
`%ni%` <- Negate(`%in%`)



#' #############################################################################
#' plot_features
#' This function is intended to be part of data integration
#'
#' @param sc_multiple Seurat object
#' @param dimred Dim reduction
#'
#' @return Plots arranged in single plot
#'
#' @examples
#' \dontrun{
#' plot_QCsample_RNA(sc_multiple, dimred)
#' }
#'
#' @export

plot_features <- function(sc_multiple,
                          dimred) {
  DefaultAssay(sc_multiple) <- "RNA"

  umap1 <- FeaturePlot(
    object = sc_multiple,
    features = c("nCount_RNA", "nFeature_RNA"),
    min.cutoff = "q1",
    max.cutoff = "q99",
    raster = FALSE,
    reduction = dimred,
    ncol = 2
  ) & labs(x = "UMAP1", y = "UMAP2")
  umap2 <- FeaturePlot(
    object = sc_multiple,
    features = c("mitoRatio", "scDblFinder_RNA.score"),
    min.cutoff = "q1",
    max.cutoff = "q99",
    raster = FALSE,
    reduction = dimred,
    ncol = 2
  ) & labs(x = "UMAP1", y = "UMAP2")
  umap3 <- FeaturePlot(
    object = sc_multiple,
    features = c("S.score", "G2M.Score"),
    min.cutoff = "q1",
    max.cutoff = "q99",
    raster = FALSE,
    reduction = dimred,
    ncol = 2
  ) & labs(
    x = "UMAP1",
    y = "UMAP2"
  )

  feat1 <- FeaturePlot(
    sc_multiple,
    reduction = dimred,
    features = c("EPCAM", "PAX8", "WFDC2"),
    min.cutoff = "q5",
    max.cutoff = "q95",
    ncol = 3,
    raster = FALSE
  ) & labs(
    x = "UMAP1",
    y = "UMAP2"
  )

  feat2 <- FeaturePlot(
    sc_multiple,
    reduction = dimred,
    features = c("DCN", "LUM", "MECOM"),
    min.cutoff = "q5",
    max.cutoff = "q95",
    ncol = 3,
    raster = FALSE
  ) & labs(
    x = "UMAP1",
    y = "UMAP2"
  )

  plot_int <- ggarrange(umap1, umap2, umap3, feat1, feat2, nrow = 5, ncol = 1)
  plot_int
}


#' #############################################################################
#' harmony_integration
#' To integrate data
#'
#' @param sc_obj Seurat object
#'
#' @return Integrated Seurat Object
#' @files Plot files of Integrated Seurat Object
#' @files RDS file of Integrated Seurat Object
#'
#'
#' @export

harmony_integration <- function(
    sc_obj,
    integration_vars,
    scale.vars,
    genes_to_include_hvf,
    output_dir) {
  # Set seed and free some memory
  gc(full = TRUE, reset = TRUE)
  set.seed(47)

  ## Create output folder
  dir.create(
    path = output_dir,
    showWarnings = FALSE
  )

  dir.create(
    path = paste0(output_dir, "/Plots"),
    showWarnings = FALSE
  )

  # Split dataset
  sc_obj@meta.data[["integrationCovariates"]] <-
    Reduce(function(x, y) paste(x, y, sep = "_"), integration_vars)

  sc_obj[["RNA"]] <- split(
    sc_obj[["RNA"]],
    sc_obj$integrationCovariates
  )

  # FindVariableFeatures, Scale and PCA
  sc_obj <- FindVariableFeatures(
    sc_obj,
    selection.method = "vst",
    nfeatures = 5000,
    verbose = TRUE
  )

  VariableFeatures(sc_obj) <- intersect(
    VariableFeatures(sc_obj), genes_to_include_hvf
  )[1:4000]

  sc_obj <- ScaleData(
    object = sc_obj,
    verbose = TRUE,
    features = VariableFeatures(sc_obj),
    vars.to.regress = scale.vars
  )

  sc_obj <- RunPCA(
    sc_obj,
    features = VariableFeatures(sc_obj),
    verbose = TRUE
  )

  # Integrate
  gc(full = TRUE, reset = TRUE)
  sc_obj <- IntegrateLayers(
    object = sc_obj,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = paste0(
      "pca.harmony.", paste(names(integration_vars), collapse = ".")
    ),
    verbose = TRUE
  )

  sc_obj <- RunUMAP(
    sc_obj,
    dims = 1:30,
    reduction = paste0(
      "pca.harmony.", paste(names(integration_vars), collapse = ".")
    ),
    reduction.name = paste0(
      "umap.harmony.", paste(names(integration_vars), collapse = ".")
    ),
    verbose = TRUE
  )


  ## PLOTS A
  umap1_integrated <- DimPlot(
    sc_obj,
    reduction = paste0(
      "umap.harmony.", paste(names(integration_vars), collapse = ".")
    ),
    group.by = "Sample",
    raster = FALSE
  ) +
    theme_classic(base_size = 10) +
    labs(x = "UMAP1", y = "UMAP2")

  umap2_integrated <- DimPlot(
    sc_obj,
    reduction = paste0(
      "umap.harmony.", paste(names(integration_vars), collapse = ".")
    ),
    group.by = "tissue",
    raster = FALSE
  ) +
    theme_classic(base_size = 15) +
    labs(x = "UMAP1", y = "UMAP2")

  umap_integrated <- wrap_plots(
    umap1_integrated,
    umap2_integrated,
    ncol = 2
  )

  ggsave(
    plot = umap_integrated,
    filename = paste0(
      output_dir, "/Plots/",
      "plot_umap.integrated.",
      paste(names(integration_vars), collapse = "."), ".jpg"
    ),
    width = 18,
    height = 6,
    device = "jpg",
    dpi = 100,
    bg = "white"
  )


  ## PLOTS B
  umap1_integrated_sp1 <- DimPlot(
    sc_obj,
    reduction = paste0(
      "umap.harmony.", paste(names(integration_vars), collapse = ".")
    ),
    group.by = "protocol",
    split.by = "tissue",
    raster = FALSE
  ) +
    theme_classic(base_size = 15) +
    labs(x = "UMAP1", y = "UMAP2")

  ggsave(
    plot = umap1_integrated_sp1,
    filename = paste0(
      output_dir, "/Plots/",
      "plot_umap.integrated.",
      paste(names(integration_vars), collapse = "."), ".sp1.jpg"
    ),
    width = 12,
    height = 4,
    device = "jpg",
    dpi = 150,
    bg = "white"
  )


  ## PLOTS C
  umap1_integrated_sp2 <- DimPlot(
    sc_obj,
    reduction = paste0(
      "umap.harmony.", paste(names(integration_vars), collapse = ".")
    ),
    group.by = "tissue",
    split.by = "protocol",
    raster = FALSE
  ) +
    theme_classic(base_size = 15) +
    labs(x = "UMAP1", y = "UMAP2")

  ggsave(
    plot = umap1_integrated_sp2,
    filename = paste0(
      output_dir, "/Plots/",
      "plot_umap.integrated.",
      paste(names(integration_vars), collapse = "."), ".sp2.jpg"
    ),
    width = 18,
    height = 4,
    device = "jpg",
    dpi = 150,
    bg = "white"
  )


  ## PLOTS D
  umap1_integrated_sp3 <- DimPlot(
    sc_obj,
    reduction = paste0(
      "umap.harmony.", paste(names(integration_vars), collapse = ".")
    ),
    group.by = "protocol",
    split.by = "group",
    raster = FALSE
  ) +
    theme_classic(base_size = 15) +
    labs(x = "UMAP1", y = "UMAP2")

  ggsave(
    plot = umap1_integrated_sp3,
    filename = paste0(
      output_dir, "/Plots/",
      "plot_umap.integrated.",
      paste(names(integration_vars), collapse = "."), ".sp3.jpg"
    ),
    width = 15,
    height = 5,
    device = "jpg",
    dpi = 100,
    bg = "white"
  )

  ## PLOTS E
  feat_plot1 <- plot_features(
    sc_obj,
    paste0(
      "umap.harmony.", paste(names(integration_vars), collapse = ".")
    )
  )

  ggsave(
    plot = feat_plot1,
    filename = paste0(
      output_dir, "/Plots/",
      "plot_markers.integrated.",
      paste(names(integration_vars), collapse = "."), ".jpg"
    ),
    width = 12,
    height = 24,
    device = "jpg",
    dpi = 150,
    bg = "white"
  )

  sc_obj <- JoinLayers(sc_obj)

  SaveSeuratRds(
    object = sc_obj,
    file = paste0(
      output_dir,
      "/", "integration.",
      paste(names(integration_vars), collapse = "."), ".rds"
    )
  )

  gc(full = TRUE, reset = TRUE)

  sc_obj
}
