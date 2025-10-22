library(ggplot2)
library(tidyverse)
plots_dir <- "output/differential_expression/WebGestalt_KEGG_ORA/plots/"

oras <- list.files(
  "output/differential_expression/WebGestalt_KEGG_ORA",
  full.names = TRUE,
  pattern = ".rds"
)

for (ora_file in oras) {
  ora <- readRDS(ora_file)
  contrast_name <- basename(ora_file)
  contrast_name <- gsub(contrast_name, pattern = ".rds", replacement = "")

  ora_filt <- ora %>%
    slice_max(enrichmentratio, n = 10) %>%
    arrange(desc(enrichmentratio))

  ggplot(
    ora_filt,
    aes(
      y = reorder(description, enrichmentratio),
      x = enrichmentratio,
      color = fdr,
      size = size
    )
  ) +
    geom_point() +
    scale_color_distiller(palette = "Spectral", direction = 1)

  ggsave(
    filename = paste0(plots_dir, contrast_name, ".png"),
    device = "png",
    dpi = 400,
    height = 4
  )
}
