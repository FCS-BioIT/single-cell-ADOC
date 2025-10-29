# Title: single cell/single nuclei RNA Quality Contol and Preprocessing
# Author: Diego Amoros modified by Jaime Llera Oyola
# Date: 06-08-2024
# Inputs:
#   - CellBender H5 results files, h5
#   - metadata, csv
# Outputs:
#   - filtering results, .rds
#   - QC plots, jpg
#   - QC stats, csv
# Dependencies:
#   - libraries: seurat v5, scCustomize
#   - custom: preprocess.R


# libraries and utils and constants
library("Seurat")
library("scCustomize")
library("knitr")
library("DT")
library("stringr")
library("RColorBrewer")
library("janitor")
library("optparse")




option_list <- list(
  make_option("--metadata", type = "character", help = "Path to metadata CSV"),
  make_option("--input_dir", type = "character", help = "Root input directory with CellBender/DecontX results"),
  make_option("--output_dir", type = "character", help = "Root output directory for QC results"),
  make_option("--atlas", type = "character", help = "Path to reference atlas .rds file"),
  make_option("--utils", type = "character", default = "src/scrnaseq-preprocess/utils_qc.R", help = "Path to utils_qc.R")
)

opt <- parse_args(OptionParser(option_list = option_list))

source(opt$utils)

set.seed(47)




## QC Filter Parameters
nfeature_rna_threshold <- 500
mitoratio_rna_threshold <- 15
ncounts_rna_threshold <- 1000




# Quality control

# QC Analysis
## Packages and versions
paste0("Seurat ", packageVersion("Seurat"))
paste0("scCustomize ", packageVersion("scCustomize"))
paste0("scDblFinder ", packageVersion("scDblFinder"))



## Load Roser-Vento Lab Cell Atlas
atlasdata_cell <- readRDS(
  opt$atlas
)
## Keep controls only
atlasdata_cell <- atlasdata_cell[
  , which(atlasdata_cell$Group != "Endo_Superficial")
]
atlasdata_cell <- atlasdata_cell[
  , which(atlasdata_cell$Hormonal.treatment == "nan")
]

atlasdata_cell <- FindVariableFeatures(
  atlasdata_cell,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)
atlasdata_cell <- ScaleData(
  atlasdata_cell,
  features = VariableFeatures(object = atlasdata_cell),
  verbose = FALSE
)
atlasdata_cell <- RunPCA(
  atlasdata_cell,
  features = VariableFeatures(object = atlasdata_cell),
  verbose = TRUE
)
gc()




# Each CellBender h5 file must be inside a folder named with
# the name of the sample.
# Folder name of sample folders and bioinformatics_id colname
# of metadata need to be the same
metadata <- read.csv(file = opt$metadata)

root_input_folders <- c(
  opt$input_dir
)
output_rds_folder <- file.path(opt$output_dir, "seurat_qc")
root_output_qc_folder <- opt$output_dir

# List samples to process
list_dir_input_folders <- list.files(
  path = root_input_folders,
  pattern = "_decontx\\.h5$",
  full.names = TRUE
)
sample_names <- basename(tools::file_path_sans_ext(list_dir_input_folders))
sample_numbers <- as.numeric(substring(sample_names, 3, 4))

list_dir_input_folders <- list_dir_input_folders[order(sample_numbers)]



## Create matrices for store cell type composition
df_lineages_bf <- data.frame(
  row.names = sample_names,
  Epithelial = rep(NA, length(sample_names)),
  Mesenchymal = rep(NA, length(sample_names)),
  Endothelial = rep(NA, length(sample_names)),
  Immune = rep(NA, length(sample_names)),
  weak_assignment = rep(NA, length(sample_names))
)

df_lineages_af <- data.frame(
  row.names = sample_names,
  Epithelial = rep(NA, length(sample_names)),
  Mesenchymal = rep(NA, length(sample_names)),
  Endothelial = rep(NA, length(sample_names)),
  Immune = rep(NA, length(sample_names)),
  weak_assignment = rep(NA, length(sample_names))
)

# Inicialize QC list
if (exists("qc_stats_df")) {
  rm("qc_stats_df")
}

qc_stats <- list()

# Loop over samples
for (id in seq_along(list_dir_input_folders)) {
  sample_file <- list_dir_input_folders[id]
  sample <- basename(sample_file)
  sample <- sub("_decontx.*$", "", sample)

  qc_stats$sample <- sample
  if (!any(metadata$sample_id == sample, na.rm = TRUE)) {
    stop(
      paste0(
        "Sample ",
        sample,
        " does not exist in metadata column sample_id"
      )
    )
  }

  # Create directory for QC plots
  dir_qc <- paste0(root_output_qc_folder, "/", sample)
  dir.create(path = dir_qc, showWarnings = FALSE, recursive = TRUE)
  dir_filtrds <- paste0(output_rds_folder)
  dir.create(path = dir_filtrds, showWarnings = FALSE, recursive = TRUE)

  # Create seurat object with scRNA-seq data
  # read decontx converted matrices
  counts <- Read10X_h5(filename = sample_file)


  sc_obj <- CreateSeuratObject(counts = counts, project = sample)

  # Add a metadata info to each cell in the sample
  for (metadata_column in seq_along(ncol(metadata))) {
    sc_obj <- AddMetaData(
      object = sc_obj,
      metadata = metadata[
        which(metadata$sample_id == sample),
        metadata_column
      ],
      col.name = colnames(metadata)[metadata_column]
    )
  }

  # Calculate the proportion of transcripts mapping to mitochondrial genes
  sc_obj[["mitoRatio"]] <- PercentageFeatureSet(sc_obj, pattern = "^MT-")


  # Prefilter statistics
  ########################################
  qc_stats$ncell_pre <- ncol(sc_obj)
  qc_stats$gene_mean_pre <- mean(sc_obj$nFeature_RNA)
  qc_stats$gene_median_pre <- median(sc_obj$nFeature_RNA)
  qc_stats$gene_sd_pre <- sd(sc_obj$nFeature_RNA)
  qc_stats$mitoratio_mean_pre <- mean(
    sc_obj$mitoRatio[which(!is.na(sc_obj$mitoRatio))]
  )
  qc_stats$mitoratio_median_pre <- median(
    sc_obj$mitoRatio[which(!is.na(sc_obj$mitoRatio))]
  )
  qc_stats$mitoratio_sd_pre <- sd(
    sc_obj$mitoRatio[which(!is.na(sc_obj$mitoRatio))]
  )
  qc_stats$counts_mean_pre <- mean(sc_obj$nCount_RNA)
  qc_stats$counts_median_pre <- median(sc_obj$nCount_RNA)
  qc_stats$counts_sd_pre <- sd(sc_obj$nCount_RNA)

  # Draw PRE plots
  plot_pre <- plot_qcsample_rna(
    sc_obj,
    nfeature_rna_threshold,
    mitoratio_rna_threshold,
    ncounts_rna_threshold
  )
  ggsave(
    plot = plot_pre,
    filename = paste0(dir_qc, "/vln_umap_pre.png"),
    width = 10,
    height = 48,
    device = "png",
    dpi = 100,
    bg = "white"
  )

  ## Transfer reference labels to query
  tranf_labels <- transfer_labels(sc_obj, atlasdata_cell)

  ggsave(
    plot = tranf_labels[[1]],
    filename = paste0(dir_qc, "/vln_umap_pre_labels.png"),
    width = 10,
    height = 15,
    device = "png",
    dpi = 100,
    bg = "white"
  )

  df_lineages_bf[
    sample,
    names(tranf_labels[[2]])
  ] <- unlist(tranf_labels[[2]])


  # Filtering
  ########################################
  # Filter Bad cells
  sc_obj <- remove_bad_cells_rna(
    sc_obj, nfeature_rna_threshold,
    mitoratio_rna_threshold,
    ncounts_rna_threshold
  )
  qc_stats$n_cells_before_doublets <- ncol(sc_obj)

  # Filter Multiplets
  sc_obj <- remove_multiplets_rna(sc_obj)

  # Postfilter statistics
  ########################################
  qc_stats$n_doublets <- qc_stats$n_cells_before_doublets - ncol(sc_obj)
  qc_stats$doublet_ratio <- qc_stats$n_doublets / ncol(sc_obj)
  qc_stats$n_cells_hq <- ncol(sc_obj)
  qc_stats$filtered_cells_ratio <- 1 - (qc_stats$n_cells_hq / qc_stats$ncell_pre) # nolint
  qc_stats$gene_mean_filt <- mean(sc_obj$nFeature_RNA)
  qc_stats$gene_median_filt <- median(sc_obj$nFeature_RNA)
  qc_stats$gene_sd_filt <- sd(sc_obj$nFeature_RNA)
  qc_stats$mitoratio_mean_filt <- mean(sc_obj$mitoRatio)
  qc_stats$mitoratio_median_filt <- median(sc_obj$mitoRatio)
  qc_stats$mitoratio_sd_filt <- sd(sc_obj$mitoRatio)
  qc_stats$counts_mean_filt <- mean(sc_obj$nCount_RNA)
  qc_stats$counts_median_filt <- median(sc_obj$nCount_RNA)
  qc_stats$counts_sd_filt <- sd(sc_obj$nCount_RNA)

  # Add decontx fraction to qc_stats
  decontx_csv_path <- sub("_decontx\\.h5$", "_decontx_metrics.csv", sample_file)

  if (!file.exists(decontx_csv_path)) {
    warning(paste("Metrics file not found for", sample, ":", decontx_csv_path))
    qc_stats$avg_decontx_fraction_counts_removed <- NA
  } else {
    decontx_csv <- read.csv(decontx_csv_path, header = TRUE)
    qc_stats$avg_decontx_fraction_counts_removed <- mean(decontx_csv$decontx_contamination)
  }
  qc_stats$avg_decontx_fraction_counts_removed <-
    mean(decontx_csv[, "decontx_contamination"])



  # Draw POST plots
  plot_post <- plot_qcsample_rna(
    sc_obj,
    nfeature_rna_threshold,
    mitoratio_rna_threshold,
    ncounts_rna_threshold
  )
  ggsave(
    plot = plot_post,
    filename = paste0(dir_qc, "/vln_umap_post.png"),
    width = 10,
    height = 48,
    device = "png",
    dpi = 72,
    bg = "white"
  )
  print("Transfering labels to filtered object:")
  ## Transfer reference labels to query
  tranf_labels <- transfer_labels(sc_obj, atlasdata_cell)

  print("Label transfer complete!")

  ggsave(
    plot = tranf_labels[[1]],
    filename = paste0(dir_qc, "/vln_umap_post_labels.png"),
    width = 10,
    height = 15,
    device = "png",
    dpi = 100,
    bg = "white"
  )

  print("Filling dataframe with celltypes from transfer")
  df_lineages_af[sample, names(tranf_labels[[2]])] <-
    unlist(tranf_labels[[2]])

  # Add to qc_stats dataframe
  print("Adding qc stats to dataframe")
  if (!exists("qc_stats_df")) {
    qc_stats_df <- as.data.frame(qc_stats)
  } else {
    qc_stats_df <- rbind(qc_stats_df, qc_stats)
  }

  # Save filtered seurat obj
  ########################################
  print("Saving seurat object")
  SaveSeuratRds(sc_obj,
    file = paste0(output_rds_folder, "/", sample, "-hq-cells.rds")
  )
  rm(sc_obj)
  invisible(gc(full = TRUE, verbose = FALSE))
}


# Save data
########################################
qc_stats_df_round <- data.frame(lapply(qc_stats_df, function(x) {
  if (is.numeric(x)) {
    round(x, 4)
  } else {
    x
  }
}))
write.csv(
  x = qc_stats_df_round,
  file = paste0(root_output_qc_folder, "/", "QC_scRNA_stats.csv"),
  row.names = FALSE
)

write.csv(
  x = cbind(rownames(df_lineages_bf), df_lineages_bf),
  file = paste0(root_output_qc_folder, "/", "lineages_pre.csv"),
  row.names = FALSE
)

write.csv(
  x = cbind(rownames(df_lineages_af), df_lineages_af),
  file = paste0(root_output_qc_folder, "/", "lineages_post.csv"),
  row.names = FALSE
)



##### SessionInfo and Enviroment
# system(paste("conda export >", "./envs/main.yaml"))
# writeLines(
#   capture.output(sessionInfo()),
#   "./envs/R_session_info-sc-sn.txt"
# )
