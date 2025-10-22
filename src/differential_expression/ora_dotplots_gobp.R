library(ggplot2)
library(tidyverse)
plots_dir <- paste0(
  "output/differential_expression/",
  "WebGestalt_GOBP_ORA/plots_selection/"
)

oras <- list.files(
  "output/differential_expression/WebGestalt_GOBP_ORA",
  full.names = TRUE,
  pattern = ".rds"
)

for (ora_file in oras) {
  ora <- readRDS(ora_file)
  ora$geneset_description <- paste0(ora$geneset, " - ", ora$description)
  contrast_name <- basename(ora_file)
  contrast_name <- gsub(contrast_name, pattern = ".rds", replacement = "")

  ora_filt <- ora %>%
    slice_max(enrichmentratio, n = 10) %>%
    arrange(desc(enrichmentratio))


  ggplot(
    ora_filt,
    aes(
      y = reorder(geneset_description, enrichmentratio),
      x = enrichmentratio,
      color = fdr
    )
  ) +
    # segment inherits color mapping
    geom_segment(
      aes(
        x = 0,
        xend = enrichmentratio,
        y = geneset_description,
        yend = geneset_description
      ),
      linewidth = 1
    ) +
    geom_point(size = 4) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "Enrichment ratio", y = NULL, color = "FDR") +
    theme_minimal(base_size = 14)

  ggsave(
    filename = paste0(plots_dir, contrast_name, ".png"),
    device = "png",
    dpi = 400,
    height = 4,
    width = 8,
    bg = "white"
  )
  ggsave(
    filename = paste0(plots_dir, contrast_name, ".svg"),
    device = "svg",
    dpi = 400,
    height = 4,
    width = 8,
    bg = "white",
    fix_text_size = FALSE
  )
}
