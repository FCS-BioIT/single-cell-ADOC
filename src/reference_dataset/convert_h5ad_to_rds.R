
# Usage:
# Access the dataset published in https://cellgeni.cog.sanger.ac.uk/vento/reproductivecellatlas/endometriumAtlasV2_cells_with_counts.h5ad 
# The first argument will be an h5ad matrix that will be converted to rds format

args = commandArgs(trailingOnly=TRUE)

# Load libraries
library(Seurat)
library(anndata)
library(Matrix)



annFile<-args[1]
output<-sub(".h5ad$",".rds",annFile)


# Read in the dataset
anndataObj <- read_h5ad(annFile,backed = "r")

invisible(gc())


if (length(anndataObj$raw) != 0) {
 if (class(anndataObj$raw) == "dgRMatrix") {
  counts = as(t(anndataObj$raw), "CsparseMatrix")
 } else {
  counts = t(anndataObj$raw)
 }

} else {

 if (class(anndataObj$X) == "dgRMatrix") {
  counts = as(t(anndataObj$X), "CsparseMatrix")
 } else {
  counts = t(anndataObj$X)
 }
}

seuratObj <- CreateSeuratObject(counts = counts, meta.data = anndataObj$obs,assay = "RNA")

invisible(gc())
rm(anndataObj)
rm(counts)
invisible(gc())
saveRDS(seuratObj, file=output)

