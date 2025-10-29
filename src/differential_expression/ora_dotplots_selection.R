library(ggplot2)
library(tidyverse)

plots_dir <- paste0(
  "output/differential_expression/WebGestalt_GOBP_ORA/",
  "plots_selection/"
)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

# <-— edit this vector
pathways <- c(
  "GO:0007165",
  "positive regulation of transcription by RNA polymerase II",
  "GO:0006915 - apoptotic process"
)
sel <- tolower(pathways)

oras <- list.files(
  "output/differential_expression/WebGestalt_GOBP_ORA",
  full.names = TRUE,
  pattern = "\\.rds$"
)

for (ora_file in oras) {
  ora <- readRDS(ora_file)
  ora <- ora %>%
    mutate(
      geneset_description = paste0(geneset, " - ", description),
      key_id = tolower(geneset),
      key_desc = tolower(description),
      key_both = tolower(geneset_description)
    )

  # keep only selected pathways (match by ID, description, or combined label)
  ora_sel <- ora %>%
    filter(key_id %in% sel | key_desc %in% sel | key_both %in% sel)

  if (nrow(ora_sel) == 0) {
    message("No selected pathways found in: ", basename(ora_file))
    next
  }

  # order y-axis by input order (not by enrichment), using first match per input
  # build a level vector in the order provided by `pathways`
  matched_levels <- map_chr(sel, \(s) {
    hit <- ora_sel$geneset_description[
      ora_sel$key_id == s | ora_sel$key_desc == s | ora_sel$key_both == s
    ]
    if (length(hit)) hit[1] else NA_character_
  }) %>%
    discard(is.na) %>%
    unique()

  ora_sel <- ora_sel %>%
    mutate(
      geneset_description = factor(geneset_description,
        levels = matched_levels
      )
    ) %>%
    arrange(geneset_description)

  contrast_name <- sub("\\.rds$", "", basename(ora_file))

  p <- ggplot(
    ora_sel,
    aes(
      y = geneset_description,
      x = enrichmentratio,
      color = fdr
    )
  ) +
    geom_segment(
      aes(x = 0, xend = enrichmentratio, yend = geneset_description),
      linewidth = 1
    ) +
    geom_point(size = 4) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = "Enrichment ratio", y = NULL, color = "FDR") +
    theme_minimal(base_size = 14)

  ggsave(paste0(plots_dir, contrast_name, ".png"), p,
    device = "png", dpi = 400, height = 4, width = 8, bg = "white"
  )
  ggsave(paste0(plots_dir, contrast_name, ".svg"), p,
    device = "svg", dpi = 400, height = 4, width = 8, bg = "white",
    fix_text_size = FALSE
  )
}
