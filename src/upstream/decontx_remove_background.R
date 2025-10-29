#!/usr/bin/env R
suppressPackageStartupMessages({
  library(optparse)
  library(Matrix)
  library(DelayedArray)
  library(celda) # decontX
  library(DropletUtils)
  library(rhdf5)
})

# ---- CLI ----
option_list <- list(
  make_option("--input", type = "character"),
  make_option("--input_raw", type = "character"),
  make_option("--output_h5", type = "character"),
  make_option("--output_csv", type = "character")
)
opt <- parse_args(OptionParser(option_list = option_list))
stopifnot(file.exists(opt$input))


# ---- Read 10x H5, run DecontX ----
sce <- read10xCounts(opt$input, col.names = TRUE) # filtered cellranger HDF5
sce_raw <- read10xCounts(opt$input_raw, col.names = TRUE) # raw cellranger cells
# Remove zero counts
print("Cells before removing cells with 0 counts")
print(dim(counts(sce)))

sce <- sce[, colSums(counts(sce)) > 0]

print("Cells after removing cells with 0 counts")
print(dim(counts(sce)))



# Run decontx
sce <- decontX(sce, background = sce_raw)


# Mathing barcodes
x <- assays(sce)$decontXcounts

# Make sure barcodes are present and aligned
if (is.null(colnames(x))) colnames(x) <- colnames(sce)
stopifnot(identical(colnames(x), colnames(sce)))

# Coerce to sparse if needed (after size is manageable from A/B)
if (!inherits(x, "dgCMatrix")) {
  x <- as(x, "dgCMatrix") # (works once ncol is reasonable)
}

dim(x) # expect genes x cells
sum(x@x == 0) # should be < length(x@x)
stopifnot(identical(colnames(x), colnames(sce)))


# Row metadata fallbacks
rd <- as.data.frame(rowData(sce))
gene_id <- if ("ID" %in% names(rd)) rd$ID else rownames(sce)
gene_name <- if ("Symbol" %in% names(rd)) rd$Symbol else gene_id
gene_type <- if ("Type" %in% names(rd)) {
  rd$Type
} else {
  rep(
    "Gene Expression",
    length(gene_id)
  )
}
genome <- if ("Genome" %in% names(rd)) rd$Genome else "unknown"


# Fallback if gene_id is still NULL/NA (should be rare)
if (is.null(gene_id)) gene_id <- paste0("gene_", seq_len(nrow(x)))

stopifnot(length(gene_id) == nrow(x))

print(opt$output_h5)

out_path <- opt$output_h5

# ---- Write 10x **HDF5** (feature-barcode v3 schema under /matrix) ----
write10xCounts(
  path        = opt$output_h5,
  x           = x, # dgCMatrix (genes x cells)
  barcodes    = colnames(x),
  gene.id     = gene_id,
  gene.symbol = gene_name,
  gene.type   = gene_type,
  genome      = genome,
  version     = "3", # 10x feature-barcode schema used by current Cell Ranger
  type        = "HDF5",
  overwrite   = TRUE,
  chemistry   = "Custom"
)


# ---- Metrics CSV ----
metrics <- data.frame(
  barcode = colnames(sce),
  nUMI_raw = Matrix::colSums(assays(sce)$counts),
  nUMI_decontx = Matrix::colSums(assays(sce)$decontXcounts),
  decontx_contamination = sce$decontX_contamination,
  stringsAsFactors = FALSE
)
write.csv(metrics, opt$output_csv, row.names = FALSE)
