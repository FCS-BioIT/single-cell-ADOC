#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Matrix)
  library(DropletUtils)
  library(optparse)
})

option_list <- list(
  make_option("--matrix", type = "character", help = "Matrix .mtx file"),
  make_option("--features", type = "character", help = "Features .tsv file"),
  make_option("--barcodes", type = "character", help = "Barcodes .tsv file"),
  make_option("--output", type = "character", help = "Output HDF5 file (.h5)")
)
opt <- parse_args(OptionParser(option_list = option_list))

# ---- Validate inputs ----
if (!file.exists(opt$matrix)) stop("Matrix file not found: ", opt$matrix)
if (!file.exists(opt$features)) stop("Features file not found: ", opt$features)
if (!file.exists(opt$barcodes)) stop("Barcodes file not found: ", opt$barcodes)

# ---- Read input ----
cat("\n[INFO] Reading input files...\n")
mat <- readMM(opt$matrix)
genes <- read.delim(opt$features, header = FALSE, stringsAsFactors = FALSE)[, 1]
barcodes <- read.delim(opt$barcodes, header = FALSE, stringsAsFactors = FALSE)[, 1]

if (nrow(mat) != length(genes)) {
  stop("Mismatch between gene count and matrix rows: ", nrow(mat), " vs ", length(genes))
}
if (ncol(mat) != length(barcodes)) {
  stop("Mismatch between barcode count and matrix columns: ", ncol(mat), " vs ", length(barcodes))
}

rownames(mat) <- genes
colnames(mat) <- barcodes

# ---- Prepare output ----
out_dir <- dirname(opt$output)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("[INFO] Writing 10x HDF5 to:", opt$output, "\n")

# ---- Direct write (no temp dir) ----
tryCatch(
  {
    DropletUtils::write10xCounts(
      path = opt$output,
      x = mat,
      barcodes = colnames(mat),
      gene.id = rownames(mat),
      gene.symbol = rownames(mat),
      version = "3",
      type = "HDF5",
      overwrite = TRUE
    )
  },
  error = function(e) {
    cat("[ERROR] DropletUtils::write10xCounts failed:\n", e$message, "\n")
    q(status = 1)
  }
)

# ---- Check that file exists ----
if (!file.exists(opt$output)) {
  cat("[ERROR] Expected output not found:", opt$output, "\n")
  q(status = 1)
}

cat("[SUCCESS] HDF5 successfully written:", opt$output, "\n")
