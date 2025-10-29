suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(readr)
})

de_dir <- "output/differential_expression"
plot_dir <- file.path(de_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

de_files <- list.files(de_dir, pattern = "^DE_.*\\.csv$", full.names = TRUE)

for (f in de_files) {
  de <- readr::read_csv(f, show_col_types = FALSE)
  lfc_col <- if ("avg_log2FC" %in% names(de)) "avg_log2FC" else "avg_logFC"
  stopifnot(all(c("gene", "p_val_adj", lfc_col) %in% names(de)))

  df <- de |>
    mutate(
      log2FC = .data[[lfc_col]],
      negLog10FDR = -log10(p_val_adj + 1e-300),
      direction = case_when(
        p_val_adj < 0.05 & log2FC > 1.5 ~ "Up in ident.1",
        p_val_adj < 0.05 & log2FC < -1.5 ~ "Down in ident.1",
        TRUE ~ "NS"
      )
    )

  lab <- bind_rows(
    df |>
      filter(direction == "Up in ident.1") |>
      arrange(desc(log2FC)) |>
      slice_head(n = 10),
    df |>
      filter(direction == "Down in ident.1") |>
      arrange(desc(log2FC)) |>
      slice_tail(n = 10)
  )

  if (grepl("DE_Epithelium", f)) {
    genes <- c(
      "PAEP", "SCGB2A1", "CEACAM5", "CEACAM6", "PI3", "TCN1", "ZG16B",
      "RNASE1", "STC1", "STC2", "THSD7A", "SLCO2B1", "IL24", "MUC5B",
      "SPP1", "UPK3B", "TM4SF4", "CLIC5", "CDH6", "VGLL1", "SOSTDC1",
      "PLCE1", "BNC2", "RIMS1", "NYAP2"
    )
    lab <- df %>% filter(gene %in% genes)
  } else if (grepl("DE_Stroma", f)) {
    genes <- c(
      "PRL", "PTHLH", "PTGS2", "IL24", "CSF3",
      "GDF15", "NRG1", "CXCR4", "EDNRB",
      "HLA-DMB", "MT1G", "MT1H", "AGPAT5",
      "FABP4", "KRTAP1-5", "KRT34", "TAGLN",
      "IFI44L", "OXTR", "IFI44", "IFIT1",
      "CNN1", "ADAMTS17", "WNT2", "MKI67"
    )
    lab <- df %>% filter(gene %in% genes)
  }

  print(lab)

  p <- ggplot(df, aes(x = log2FC, y = negLog10FDR, color = direction)) +
    geom_point(size = 0.9, alpha = 0.6) +
    geom_vline(xintercept = c(-1.5, 1.5), linetype = "dotted") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
    scale_color_manual(
      values = c(
        "Up in ident.1" = "#ce0665",
        "Down in ident.1" = "#0f99b2",
        "NS" = "grey75"
      )
    ) +
    labs(
      title = basename(f),
      x = "avg_log2FC (ident.1 / ident.2)",
      y = expression(-log[10]("FDR"))
    ) +
    theme_classic(base_size = 11) +
    geom_text_repel(
      data = lab,
      aes(label = gene),
      size = 2.6,
      max.overlaps = 100,
      show.legend = FALSE
    )

  outfile_png <- file.path(
    plot_dir, paste0(tools::file_path_sans_ext(basename(f)), "_volcano.png")
  )
  outfile_svg <- file.path(
    plot_dir, paste0(tools::file_path_sans_ext(basename(f)), "_volcano.svg")
  )
  ggsave(
    outfile_png, p,
    width = 7.2, height = 5.2, dpi = 300, bg = "white"
  )
  ggsave(
    outfile_svg, p,
    width = 7.2, height = 5.2, dpi = 300, bg = "white", fix_text_size = FALSE
  )
}
