library("Seurat")
library("scDblFinder")
library("ggpubr")
library("patchwork")
library("Matrix")
library("ggplot2")
library("RColorBrewer")
library("ggpubr")
## Utils
`%ni%` <- Negate(`%in%`)


#' ###########################################################################
#' plot_qcsample_rna
#' QC Plots
#' This function is intended to be part of data preprocess and filtering
#'
#' @param sc_obj Seurat object
#' @param nfeature_rna_threshold
#' @param mitoratio_rna_threshold
#' @param ncounts_rna_threshold
#'
#' @return Filtered Seurat Object
#'
#' @examples
#' \dontrun{
#' plot_qcsample_rna(
#'   seurat_obj,
#'   nfeature_rna_threshold,
#'   mitoratio_rna_threshold,
#'   ncounts_rna_threshold
#' )
#' }
#'
#' @export

plot_qcsample_rna <- function(sc_obj,
                              nfeature_rna_threshold,
                              mitoratio_rna_threshold,
                              ncounts_rna_threshold) {
  if (length(Cells(sc_obj)) < 50) {
    runpca_npcs <- length(Cells(sc_obj)) - 1
  } else {
    runpca_npcs <- 50
  }

  DefaultAssay(sc_obj) <- "RNA"
  vln <- VlnPlot(
    object = sc_obj,
    features = c("nCount_RNA", "nFeature_RNA", "mitoRatio"),
    pt.size = 0,
    ncol = 3
  )
  sc_obj <- NormalizeData(
    object = sc_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  sc_obj <- FindVariableFeatures(
    object = sc_obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  sc_obj <- ScaleData(
    object = sc_obj,
    features = Seurat::VariableFeatures(object = sc_obj),
    verbose = FALSE
  )
  sc_obj <- RunPCA(
    object = sc_obj,
    npcs = runpca_npcs,
    features = Seurat::VariableFeatures(object = sc_obj),
    verbose = FALSE
  )
  sc_obj <- RunUMAP(sc_obj, dims = 1:10, verbose = FALSE)
  umap1 <- FeaturePlot(
    object = sc_obj,
    features = c("nCount_RNA"),
    min.cutoff = "q5",
    max.cutoff = "q90",
    raster = FALSE
  )
  umap2 <- FeaturePlot(
    object = sc_obj,
    features = c("nFeature_RNA"),
    min.cutoff = "q5",
    max.cutoff = "q90",
    raster = FALSE
  )
  umap3 <- FeaturePlot(
    object = sc_obj[, which(!is.na(sc_obj$mitoRatio))],
    features = c("mitoRatio"),
    min.cutoff = "q5",
    max.cutoff = "q90",
    raster = FALSE
  )

  hq_cells <- (sc_obj$nFeature_RNA >= nfeature_rna_threshold) &
    (!is.na(sc_obj$mitoRatio)) &
    (sc_obj$mitoRatio < mitoratio_rna_threshold) &
    (sc_obj$nCount_RNA >= ncounts_rna_threshold
    )

  sc_obj <- AddMetaData(
    object = sc_obj,
    metadata = hq_cells,
    col.name = "hq_cells"
  )
  umap_fail_1 <- DimPlot(object = sc_obj, group.by = "hq_cells")
  umap5 <- FeaturePlot(
    object = sc_obj,
    features = c("LUM"),
    min.cutoff = "q5",
    max.cutoff = "q90",
    raster = FALSE
  )
  umap6 <- FeaturePlot(
    object = sc_obj,
    features = c("EPCAM"),
    min.cutoff = "q5",
    max.cutoff = "q90",
    raster = FALSE
  )
  umap7 <- FeaturePlot(
    object = sc_obj,
    features = c("ACTG2"),
    min.cutoff = "q5",
    max.cutoff = "q90",
    raster = FALSE
  )
  umap8 <- FeaturePlot(
    object = sc_obj,
    features = c("VWF"),
    min.cutoff = "q5",
    max.cutoff = "q90",
    raster = FALSE
  )
  umap9 <- FeaturePlot(
    object = sc_obj,
    features = c("PTPRC"),
    min.cutoff = "q5",
    max.cutoff = "q90",
    raster = FALSE
  )
  umap10 <- FeaturePlot(
    object = sc_obj,
    features = c("RGS5"),
    min.cutoff = "q5",
    max.cutoff = "q90",
    raster = FALSE
  )

  p1 <- VlnPlot(
    object = sc_obj,
    features = "EPCAM"
  ) + NoLegend() + ggtitle("EPCAM") & theme(axis.title.x = element_blank())
  p2 <- VlnPlot(
    object = sc_obj, features = "PAX8"
  ) + NoLegend() +
    ggtitle("PAX8\n(Epithelial)") & theme(axis.title.x = element_blank())
  p3 <- VlnPlot(
    object = sc_obj,
    features = "WFDC2"
  ) + NoLegend() +
    ggtitle("WFDC2") & theme(axis.title.x = element_blank())
  vln2 <- wrap_plots(p1, p2, p3, ncol = 3)

  p1 <- VlnPlot(
    object = sc_obj,
    features = "ACTG2"
  ) + NoLegend() +
    ggtitle("ACTG2") & theme(axis.title.x = element_blank())
  p2 <- VlnPlot(
    object = sc_obj,
    features = "MYH11"
  ) +
    NoLegend() +
    ggtitle(
      "MYH11\n(Smooth Muscle Cells)"
    ) & theme(axis.title.x = element_blank())
  p3 <- VlnPlot(
    object = sc_obj,
    features = "DES"
  ) +
    NoLegend() +
    ggtitle("DES") & theme(axis.title.x = element_blank())
  vln3 <- wrap_plots(p1, p2, p3, ncol = 3)

  p1 <- VlnPlot(
    object = sc_obj,
    features = "TPPP3"
  ) +
    NoLegend() +
    ggtitle("TPPP3\n(Ciliated)") & theme(axis.title.x = element_blank())
  p2 <- VlnPlot(
    object = sc_obj,
    features = "PTGS1"
  ) +
    NoLegend() +
    ggtitle("PTGS1\n(Luminal)") & theme(axis.title.x = element_blank())
  p3 <- VlnPlot(
    object = sc_obj,
    features = "PAEP"
  ) +
    NoLegend() +
    ggtitle("PAEP\n(Glandular)") & theme(axis.title.x = element_blank())
  vln4 <- wrap_plots(p1, p2, p3, ncol = 3)

  p1 <- VlnPlot(
    object = sc_obj,
    features = "LUM"
  ) +
    NoLegend() +
    ggtitle("LUM\n(Stromal)") & theme(axis.title.x = element_blank())
  p2 <- VlnPlot(
    object = sc_obj,
    features = "VWF"
  ) +
    NoLegend() +
    ggtitle("VWF\n(Endothelial)") & theme(axis.title.x = element_blank())
  p3 <- VlnPlot(
    object = sc_obj,
    features = "RGS5"
  ) +
    NoLegend() +
    ggtitle(
      "RGS5\n(Perivascular)"
    ) & theme(axis.title.x = element_blank())
  vln5 <- wrap_plots(p1, p2, p3, ncol = 3)

  ggarrange <- ggarrange(
    umap1,
    umap2,
    umap3,
    umap_fail_1,
    umap5,
    umap6,
    umap7,
    umap8,
    umap9,
    umap10,
    nrow = 5,
    ncol = 2
  )
  ggarrange2 <- ggarrange(vln2, vln3, vln4, vln5, nrow = 4, ncol = 1)
  plot_qc <- vln / ggarrange / ggarrange2 + plot_layout(heights = c(1, 10, 10))
  plot_qc
}


#' ###########################################################################
#' remove_bad_cells_rna
#' To remove bad cells using RNA data
#' This function is intended to be part of data preprocess and filtering
#'
#' @param sc_obj Seurat object
#' @param nfeature_rna_threshold
#' @param mitoratio_rna_threshold
#' @param ncounts_rna_threshold
#'
#' @return Filtered Seurat Object
#'
#' @examples
#' \dontrun{
#' remove_bad_cells_rna(
#'   seurat_obj,
#'   nfeature_rna_threshold,
#'   mitoratio_rna_threshold,
#'   ncounts_rna_threshold
#' )
#' }
#'
#' @export

remove_bad_cells_rna <- function(
    sc_obj,
    nfeature_rna_threshold,
    mitoratio_rna_threshold,
    ncounts_rna_threshold) {
  # Removing cells with NA mitoratio (they are almost empty cells)
  sc_obj <- sc_obj[, which(is.finite(sc_obj$mitoRatio))]

  # Removing cells using thresholds
  sc_obj <- subset(
    x = sc_obj,
    subset = nFeature_RNA >= nfeature_rna_threshold &
      mitoRatio < mitoratio_rna_threshold &
      nCount_RNA >= ncounts_rna_threshold
  )
  sc_obj
}


#' ###########################################################################
#' remove_multiplets_rna
#' To remove multiplets using RNA data
#' This function is intended to be part of data preprocess and filtering
#'
#' @param sc_obj Seurat object
#'
#' @return Filtered Seurat Object
#'
#' @examples
#' \dontrun{
#' remove_multiplets_rna(seurat_obj)
#' }
#'
#' @export

remove_multiplets_rna <- function(sc_obj) {
  # Removing doublets using scDblFinder
  sce <- scDblFinder(as.SingleCellExperiment(sc_obj), verbose = FALSE)
  sc_obj$scDblFinder_RNA.class <- sce$scDblFinder.class
  sc_obj$scDblFinder_RNA.score <- sce$scDblFinder.score
  sc_obj <- sc_obj[, which(sc_obj$scDblFinder_RNA.class != "doublet")]
  sc_obj$scDblFinder_RNA.class <- NULL
  sc_obj
}




#' ###########################################################################
#' plot_qcsample_rna
#' QC Plots
#' This function is intended to be part of data preprocess and filtering
#'
#' @param sc_obj Seurat object
#' @param reference nfeature_rna_threshold

#'
#' @return Filtered Seurat Object
#'
#' @examples
#' \dontrun{
#' transfer_labels(
#'   seurat_obj,
#'   reference
#' )
#' }
#'
#' @export

transfer_labels <- function(sc_obj,
                            reference) {
  if (length(Cells(sc_obj)) < 50) {
    runpca_npcs <- length(Cells(sc_obj)) - 1
  } else {
    runpca_npcs <- 50
  }

  DefaultAssay(sc_obj) <- "RNA"
  sc_obj <- NormalizeData(
    object = sc_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  sc_obj <- FindVariableFeatures(
    object = sc_obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  sc_obj <- ScaleData(
    object = sc_obj,
    features = Seurat::VariableFeatures(object = sc_obj),
    verbose = FALSE
  )
  sc_obj <- RunPCA(
    object = sc_obj,
    npcs = runpca_npcs,
    features = Seurat::VariableFeatures(object = sc_obj),
    verbose = FALSE
  )
  sc_obj <- RunUMAP(sc_obj, dims = 1:10, verbose = FALSE)

  atlasdata_seurat_anchors <- FindTransferAnchors(
    reference = reference,
    query = sc_obj,
    dims = 1:30,
    normalization.method = "LogNormalize",
    reference.reduction = "pca"
  )
  predictions_lineage <- TransferData(
    anchorset = atlasdata_seurat_anchors,
    refdata = reference$lineage,
    dims = 1:30
  )
  predictions_celltype <- TransferData(
    anchorset = atlasdata_seurat_anchors,
    refdata = reference$celltype,
    dims = 1:30
  )

  predictions_lineage$predicted.id[
    predictions_lineage$prediction.score.max < 0.75
  ] <- "weak_assignment"
  predictions_celltype$predicted.id[
    predictions_celltype$prediction.score.max < 0.75
  ] <- "weak_assignment"

  colnames(predictions_lineage) <- paste(
    colnames(predictions_lineage),
    "lineage",
    sep = "."
  )
  colnames(predictions_celltype) <- paste(
    colnames(predictions_celltype),
    "celltype",
    sep = "."
  )

  sc_obj <- AddMetaData(sc_obj, metadata = predictions_lineage)
  sc_obj <- AddMetaData(sc_obj, metadata = predictions_celltype)

  wapc_lineage <- round(
    table(predictions_lineage$predicted.id.lineage)["weak_assignment"] /
      sum(table(predictions_lineage$predicted.id.lineage)) * 100, 2
  )
  wapc_celltype <- round(
    table(predictions_celltype$predicted.id.celltype)["weak_assignment"] /
      sum(table(predictions_celltype$predicted.id.celltype)) * 100, 2
  )

  mypalette <- colorRampPalette(
    brewer.pal(8, "Set2")
  )(ncol(predictions_lineage))
  mypalette[which(
    unique(
      sc_obj$predicted.id.lineage
    )[order(unique(sc_obj$predicted.id.lineage))] == "weak_assignment"
  )] <- "darkgray"
  umap1 <- DimPlot(
    object = sc_obj,
    group.by = "predicted.id.lineage",
    label = TRUE,
    cols = mypalette
  ) +
    ggtitle(
      paste0(
        "Predicted Lineages (weak assignment = ",
        wapc_lineage, "%)"
      )
    )

  mypalette <- colorRampPalette(
    brewer.pal(8, "Set2")
  )(ncol(predictions_celltype))
  mypalette[which(unique(
    sc_obj$predicted.id.celltype
  )[order(
    unique(sc_obj$predicted.id.celltype)
  )] == "weak_assignment")] <- "darkgray"

  umap2 <- DimPlot(
    object = sc_obj,
    group.by = "predicted.id.celltype",
    label = TRUE,
    cols = mypalette
  ) +
    ggtitle(paste0(
      "Predicted Cell Types (weak assignment = ",
      wapc_celltype, "%)"
    ))

  cells_lineages <- as.list(rep(0, 5))
  names(
    cells_lineages
  ) <- c(
    "Epithelial",
    "Mesenchymal",
    "Endothelial",
    "Immune",
    "weak_assignment"
  )

  lineages <- as.list(table(predictions_lineage$predicted.id.lineage))

  if (any(names(lineages) == "Epithelial")) {
    cells_lineages$Epithelial <- lineages$Epithelial
  } else {
    cells_lineages$Epithelial <- 0
  }

  if (any(names(lineages) == "Mesenchymal")) {
    cells_lineages$Mesenchymal <- lineages$Mesenchymal
  } else {
    cells_lineages$Mesenchymal <- 0
  }

  if (any(names(lineages) == "Endothelial")) {
    cells_lineages$Endothelial <- lineages$Endothelial
  } else {
    cells_lineages$Endothelial <- 0
  }

  if (any(names(lineages) == "Immune")) {
    cells_lineages$Immune <- lineages$Immune
  } else {
    cells_lineages$Immune <- 0
  }

  if (any(names(lineages) == "weak_assignment")) {
    cells_lineages$weak_assignment <- lineages$weak_assignment
  } else {
    cells_lineages$weak_assignment <- 0
  }



  plot_qc <- ggarrange(umap1, umap2, nrow = 2, ncol = 1)
  results <- list(plot_qc, cells_lineages)

  results
}
