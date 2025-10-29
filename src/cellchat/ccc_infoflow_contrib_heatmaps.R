suppressPackageStartupMessages({
  library(optparse)
  library(CellChat)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(forcats)
  library(readr)
  library(ComplexHeatmap)
})

# ------------------------- CLI -------------------------
option_list <- list(
  make_option("--ccc_rds",
    type = "character", default = "output/ccc/cellchat_merged_comparison.rds",
    help = paste(
      "RDS with list: cc_list (two CellChat objects),",
      " multiple_CC (merged comparison), net"
    )
  ),
  make_option("--outdir", type = "character", default = "output/ccc"),
  make_option("--group1", type = "character", default = "CNT"),
  make_option("--group2", type = "character", default = "HT"),
  make_option("--color1", type = "character", default = "#5F9ED2"), # for group1
  make_option("--color2", type = "character", default = "#D92F02"), # for group2
  make_option(
    "--groupname1",
    type = "character", default = "CNT"
  ), # display name for group1
  make_option(
    "--groupname2",
    type = "character", default = "HT"
  ), # display name for group2
  make_option("--alpha", type = "double", default = 0.05),
  make_option("--tol", type = "double", default = 0.05),
  make_option("--img_w", type = "double", default = 12),
  make_option("--img_h", type = "double", default = 10),
  make_option("--dpi", type = "integer", default = 400),
  make_option(
    "--format",
    type = "character",
    default = "svg"
  ) # png or svg or pdf
)
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

# ------------------------- CLI (add these two options) -----------------------
option_list <- append(option_list, list(
  make_option("--paths",
    type = "character",
    default = "PROGESTERONE",
    help = "Comma-separated signaling pathways to plot as circle plots"
  ),
  make_option("--circle_w", type = "double", default = 10),
  make_option("--circle_h", type = "double", default = 8)
))
opt <- parse_args(OptionParser(option_list = option_list))
dir.create(
  file.path(opt$outdir, "circle_plots"),
  showWarnings = FALSE,
  recursive = TRUE
)


# ------------------------- I/O -------------------------
ccc_merged <- readRDS(opt$ccc_rds)
ccc_ht <- readRDS("output/ccc/cellchat_HT.rds")
ccc_cnt <- readRDS("output/ccc/cellchat_CNT.rds")
cc_list <- list("HT" = ccc_ht, "CNT" = ccc_cnt)

ccc <- list("cc_list" = cc_list, "multiple_CC" = ccc_merged)


stopifnot("cc_list" %in% names(ccc), "multiple_CC" %in% names(ccc))

# Ensure cc_list contains the requested groups and order them as group1, group2
cc_names <- names(ccc$cc_list)
if (is.null(cc_names)) {
  stop(
    "ccc$cc_list must have names corresponding to groups in your comparison."
  )
}
needed <- c(opt$group1, opt$group2)
if (!all(needed %in% cc_names)) {
  stop(sprintf(
    "Requested groups not found. Have: {%s}; need: {%s}",
    paste(cc_names, collapse = ", "), paste(needed, collapse = ", ")
  ))
}
ccc$cc_list <- ccc$cc_list[needed]
names(ccc$cc_list) <- needed

# ------------------------- Helpers -------------------------
height_from_n <- function(n) {
  h <- (16 * n) / 99
  if (h < 4) h <- 4
  h
}
display_map <- function(x) {
  x <- as.character(x)
  x[x == opt$group1] <- opt$groupname1
  x[x == opt$group2] <- opt$groupname2
  factor(x,
    levels = c(opt$groupname1, opt$groupname2)
  )
}

# ------------------------- Info Flow -------------------------
data_rel <- rankNet(ccc$multiple_CC,
  mode = "comparison",
  stacked = TRUE,
  do.stat = TRUE,
  return.data = TRUE,
  cutoff.pvalue = opt$alpha,
  thresh = opt$alpha,
  tol = opt$tol
)
data_abs <- rankNet(ccc$multiple_CC,
  mode = "comparison",
  stacked = FALSE,
  do.stat = TRUE,
  return.data = TRUE,
  cutoff.pvalue = opt$alpha,
  thresh = opt$alpha
)

df_rel <- data_rel$signaling.contribution %>%
  filter(pvalues < opt$alpha) %>%
  mutate(group = factor(as.character(group),
    levels = c(opt$group1, opt$group2)
  ))

df_abs <- data_abs$signaling.contribution %>%
  filter(name %in% unique(df_rel$name)) %>%
  mutate(group = factor(as.character(group),
    levels = c(opt$group1, opt$group2)
  ))

# Defining the plots of contributions in each condition
cols_use <- c(opt$color1, opt$color2)
# color labels by dominant group (using
# contribution.relative.1 = grp2/grp1 convention in CellChat)

colors.text <- ifelse(
  (df_rel$contribution.relative.1 < 1 - opt$tol) &
    (df_rel$pvalues < opt$alpha),
  cols_use[1],
  ifelse(
    (df_rel$contribution.relative.1 > 1 + opt$tol) &
      (df_rel$pvalues < opt$alpha),
    cols_use[2],
    "black"
  )
)

p_rel <- ggplot(df_rel, aes(x = name, y = contribution, fill = group)) +
  geom_bar(stat = "identity", width = 0.75, position = "fill") +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(x = NULL, y = "Relative information flow") +
  scale_fill_manual(values = cols_use, name = "") +
  coord_flip() +
  CellChat_theme_opts() +
  theme_classic() +
  theme(axis.text.y = element_text(colour = colors.text))

p_abs <- ggplot(df_abs, aes(x = name, y = contribution.scaled, fill = group)) +
  geom_bar(stat = "identity", width = 0.75, position = position_dodge(0.8)) +
  labs(x = NULL, y = "Information flow") +
  scale_fill_manual(values = cols_use, name = "") +
  coord_flip() +
  CellChat_theme_opts() +
  theme_classic() +
  theme(axis.text.y = element_text(colour = colors.text))

p_combo <- ggarrange(
  p_rel,
  p_abs,
  ncol = 2,
  common.legend = TRUE,
  legend = "top"
)

ggsave(file.path(opt$outdir, paste0("infoflow_rel_abs.", opt$format)),
  p_combo,
  width = opt$img_w, height = opt$img_h, dpi = opt$dpi, bg = "white"
)

write_csv(df_rel, file.path(opt$outdir, "infoflow_relative_significant.csv"))
write_csv(df_abs, file.path(opt$outdir, "infoflow_absolute_matched.csv"))

# ------------------------- Contribution plots -------------------------
sig_paths <- unique(df_rel$name)
for (pathway in sig_paths) {
  c1 <- tryCatch(
    netAnalysis_contribution(ccc$cc_list[[opt$group1]],
      signaling = pathway, return.data = TRUE
    ),
    error = function(e) NULL
  )
  c2 <- tryCatch(
    netAnalysis_contribution(ccc$cc_list[[opt$group2]],
      signaling = pathway, return.data = TRUE
    ),
    error = function(e) NULL
  )
  if (is.null(c1) && is.null(c2)) next
  df1 <- if (!is.null(c1)) {
    c1$LR.contribution %>%
      mutate(group = opt$group1)
  } else {
    tibble()
  }
  df2 <- if (!is.null(c2)) {
    c2$LR.contribution %>%
      mutate(group = opt$group2)
  } else {
    tibble()
  }
  dff <- bind_rows(df1, df2)
  if (nrow(dff) == 0) next

  dff <- dff %>%
    arrange(desc(contribution)) %>%
    mutate(
      name = fct_reorder(name, contribution),
      group = factor(group, levels = c(opt$group2, opt$group1)),
      group_disp = display_map(group)
    )

  n_lr <- n_distinct(dff$name)
  h <- height_from_n(n_lr)

  p_contrib <- ggplot(dff, aes(x = contribution, y = name, fill = group)) +
    geom_bar(stat = "identity", width = 0.75) +
    facet_wrap(~group_disp, ncol = 2, scales = "free_y") +
    labs(
      title = paste("Contribution –", pathway),
      x = "Relative contribution", y = NULL
    ) +
    scale_fill_manual(
      values = c("HT" = opt$color2, "CNT" = opt$color1), guide = "none"
    ) +
    theme_classic() +
    theme(
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      strip.text = element_text(size = 12, face = "bold")
    )

  out_file <- file.path(
    opt$outdir,
    paste0(
      "individual_pathways/",
      "contrib_",
      gsub(
        "[^A-Za-z0-9_]+",
        "_",
        pathway
      ),
      ".", opt$format
    )
  )
  ggsave(out_file, p_contrib, width = 12, height = h, dpi = opt$dpi)
}

# ------------------------- Heatmaps (incoming/outgoing) ----------------------
# Use the significant pathway set from df_abs/df_rel to keep report compact
path_keep <- unique(df_abs$name)

# Both groups must have CellChat objects:
g1 <- ccc$cc_list[[opt$group1]]
g2 <- ccc$cc_list[[opt$group2]]
if (length(path_keep) > 0) {
  # outgoing
  ht1_out <- netAnalysis_signalingRole_heatmap(
    g1,
    pattern = "outgoing",
    signaling = path_keep,
    title = opt$groupname1,
    width = 6,
    height = 14
  )
  ht2_out <- netAnalysis_signalingRole_heatmap(
    g2,
    pattern = "outgoing",
    signaling = path_keep,
    title = opt$groupname2,
    width = 6,
    height = 14
  )
  png(
    file.path(
      opt$outdir,
      paste0(
        "heatmap_outgoing_",
        opt$groupname1,
        "_vs_",
        opt$groupname2,
        ".png"
      )
    ),
    width = 3000,
    height = 2500,
    res = 300
  )
  draw(ht2_out + ht1_out, ht_gap = unit(0.5, "cm"))
  dev.off()

  # incoming
  ht1_in <- netAnalysis_signalingRole_heatmap(
    g1,
    pattern = "incoming",
    signaling = path_keep,
    title = opt$groupname1,
    width = 6,
    height = 14
  )
  ht2_in <- netAnalysis_signalingRole_heatmap(
    g2,
    pattern = "incoming",
    signaling = path_keep,
    title = opt$groupname2,
    width = 6,
    height = 14
  )
  png(
    file.path(
      opt$outdir,
      paste0(
        "heatmap_incoming_",
        opt$groupname1,
        "_vs_",
        opt$groupname2,
        ".png"
      )
    ),
    width = 3000,
    height = 2500,
    res = 300
  )
  draw(ht2_in + ht1_in, ht_gap = unit(0.5, "cm"))
  dev.off()
} else {
  message("No significant pathways to display in heatmaps at the chosen alpha.")
}

# ------------------------- Circle plots: selected pathways (fixed) -----------
dir.create(
  file.path(opt$outdir, "circle_plots"),
  showWarnings = FALSE,
  recursive = TRUE
)

# candidate list from CLI
paths_requested <- trimws(unlist(strsplit(opt$paths, ",")))

# helpers
sanitize <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

present_in <- function(obj) {
  # try canonical storage
  pw <- tryCatch(obj@netP$pathways, error = function(e) character())
  if (length(pw) > 0) {
    return(unique(as.character(pw)))
  }
  # fallbacks: names on prob/weight lists
  p1 <- tryCatch(names(obj@net$prob), error = function(e) character())
  p2 <- tryCatch(names(obj@net$weight), error = function(e) character())
  unique(as.character(c(p1, p2)))
}

get_edge_vals <- function(obj, pathway) {
  df <- tryCatch(
    subsetCommunication(obj, signaling = pathway),
    error = function(e) NULL
  )
  if (is.null(df) || nrow(df) == 0) {
    return(numeric(0))
  }
  v <- if ("prob" %in% names(df)) {
    df$prob
  } else if ("weight" %in% names(df)) {
    df$weight
  } else {
    numeric(0)
  }
  as.numeric(v[is.finite(v)])
}

open_device <- function(file, width, height, dpi, fmt) {
  switch(tolower(fmt),
    "png" = png(
      file,
      width = width * dpi,
      height = height * dpi,
      res = dpi,
      bg = "white"
    ),
    "pdf" = pdf(
      file,
      width = width,
      height = height,
      onefile = FALSE
    ),
    "svg" = svg(
      file,
      width = width,
      height = height,
      bg = "white"
    ),
    {
      png(
        file,
        width = width * dpi,
        height = height * dpi,
        res = dpi,
        bg = "white"
      )
    }
  )
}

# determine which requested pathways are present
# in at least one group (case-insensitive)
present_g1 <- present_in(ccc$cc_list[[opt$group1]])
present_g2 <- present_in(ccc$cc_list[[opt$group2]])
present_any <- unique(c(present_g1, present_g2))

match_ci <- function(need, have) have[tolower(have) %in% tolower(need)]
paths_to_plot <- match_ci(paths_requested, present_any)

if (length(paths_to_plot) == 0) {
  message(
    "None of the requested pathways available to plot in either CellChat object"
  )
} else {
  for (pathway in paths_to_plot) {
    # common edge scaling across the two groups via communication table
    vals1 <- get_edge_vals(ccc$cc_list[[opt$group1]], pathway)
    vals2 <- get_edge_vals(ccc$cc_list[[opt$group2]], pathway)
    wmax <- suppressWarnings(max(c(vals1, vals2), na.rm = TRUE))
    if (!is.finite(wmax)) wmax <- NA_real_
    ewmax <- if (is.na(wmax)) NULL else wmax

    for (grp in c(opt$group1, opt$group2)) {
      obj <- ccc$cc_list[[grp]]
      # skip if pathway not present in this group
      if (!(tolower(pathway) %in% tolower(present_in(obj)))) {
        message("Skipping ", pathway, " in ", grp, " (not present).")
        next
      }
      label <- if (grp == opt$group1) opt$groupname1 else opt$groupname2
      outfile <- file.path(
        opt$outdir, "circle_plots",
        paste0(
          "circle_",
          sanitize(pathway),
          "_",
          sanitize(grp),
          ".",
          opt$format
        )
      )
      open_device(
        outfile,
        width = opt$circle_w,
        height = opt$circle_h,
        dpi = opt$dpi,
        fmt = opt$format
      )

      op <- par(no.readonly = TRUE)
      par(xpd = NA, mar = c(4, 4, 4, 4)) # <- key: more room + allow overdraw

      # NOTE: this draws directly and returns NULL
      netVisual_aggregate(
        obj,
        signaling = pathway,
        layout = "circle",
        edge.weight.max = ewmax,
        vertex.size.max = 10,
        edge.width.max = 5,
        vertex.label.cex = 2,
        pt.title = 10
      )
      par(op)
      dev.off()
    }
  }
}

cat("Done\n")
