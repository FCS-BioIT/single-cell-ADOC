library(dittoSeq)
library(Seurat)

seu <- LoadSeuratRds(
  file = "output/annotation/marker_annotation/data/seurat_filtered_annoted.rds"
)


VlnPlot(
  seu,
  group.by = "sample_id",
  features = c("nFeature_RNA"),
  pt.size = 0,
  log = TRUE
)
ggsave(
  filename = "output/qc/plots/vln_qc_nfeature.png",
  width = 5,
  height = 4,
  dpi = 400
)
ggsave(
  filename = "output/qc/plots/vln_qc_nfeature.svg",
  fix_text_size = FALSE,
  width = 5,
  height = 4,
  dpi = 400
)


VlnPlot(
  seu,
  group.by = "sample_id",
  features = c("nCount_RNA"),
  pt.size = 0,
  log = TRUE
)
ggsave(
  filename = "output/qc/plots/vln_qc_ncount.png",
  width = 5,
  height = 4,
  dpi = 400
)
ggsave(
  filename = "output/qc/plots/vln_qc_ncount.svg",
  fix_text_size = FALSE,
  width = 5,
  height = 4,
  dpi = 400
)


VlnPlot(
  seu,
  group.by = "sample_id",
  features = c("mitoRatio"),
  pt.size = 0
)
ggsave(
  filename = "output/qc/plots/vln_qc_mitoratio.png",
  width = 5,
  height = 4,
  dpi = 400
)
ggsave(
  filename = "output/qc/plots/vln_qc_mitoratio.svg",
  fix_text_size = FALSE,
  width = 5,
  height = 4,
  dpi = 400
)
