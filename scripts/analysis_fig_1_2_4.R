# =============================================================================
# Single-cell Telomere Analysis Pipeline
# =============================================================================
# Description:
#   Generates figures for the single-cell telomere manuscript. Analyses are
#   based on the Satpathy et al. single-cell ATAC-seq dataset and TCGA 
#   single-cell atlas dataset.
#   Figures cover:
#     - Figure 1: Telomere read depth, coverage thresholds, and iTRPM by
#                 cell type and patient
#     - Figure 2: Chromatin state correlations, openness metric definition,
#                 and relationship with telomere content
#     - Figure 4: Openness-adjusted telomere content by treatment response
#     - Supplementary figures (PF2S, PF3S–PF5S)
# =============================================================================


# ── Packages ──────────────────────────────────────────────────────────────────

library(readr)
library(forcats)
library(tidyverse)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(tidyr)
library(rstatix)
library(reshape2)
library(scales)
library(viridis)


# ── Color palettes & cell-type ordering ───────────────────────────────────────

# Color map for T cell subtypes (used in T-cell-only plots)
tcell_colors <- c(
  "No T cell"         = "#999999",
  "Treg 1"            = "#c6a4d4",
  "Treg 2"            = "#a37ec9",
  "Treg 3"            = "#8356b4",
  "Treg 4"            = "#5e3fa0",
  "Early TEx"         = "#b3d4fc",
  "Intermediate TEx"  = "#569fd6",
  "Terminal TEx"      = "#1f78b4",
  "Th1"               = "#fdb462",
  "Th17"              = "#f67c36",
  "Tfh 1"             = "#ffd92f",
  "Tfh 2"             = "#e6ac00",
  "Naive CD4 T"       = "#b9e3b1",
  "Activated CD4 T"   = "#71c476",
  "Memory CD4 T"      = "#2e8b57",
  "Naive CD8 T"       = "#f4a6a6",
  "Effector CD8 T"    = "#e15759",
  "Memory CD8 T"      = "#a63636",
  "Other T"           = "#cccccc"
)

tcell_order <- c(
  "No T cell",
  "Treg 1", "Treg 2", "Treg 3", "Treg 4",
  "Early TEx", "Intermediate TEx", "Terminal TEx",
  "Th1", "Th17",
  "Tfh 1", "Tfh 2",
  "Naive CD4 T", "Activated CD4 T", "Memory CD4 T",
  "Naive CD8 T", "Effector CD8 T", "Memory CD8 T",
  "Other T"
)

# Color map for full tumor microenvironment (T cells + all other cell types)
tumor_tcell_colors <- c(
  "Treg 1"           = "#c29dda",
  "Treg 2"           = "#aa78d4",
  "Treg 3"           = "#944ec4",
  "Treg 4"           = "#6c1ab6",
  "Early TEx"        = "#aed6f1",
  "Intermediate TEx" = "#5dade2",
  "Terminal TEx"     = "#2874a6",
  "Th1"              = "#fbb45c",
  "Th17"             = "#e67e22",
  "Tfh 1"            = "#fff176",
  "Tfh 2"            = "#f4d03f",
  "Naive CD4 T"      = "#a9dfbf",
  "Activated CD4 T"  = "#58d68d",
  "Memory CD4 T"     = "#2e8b57",
  "Naive CD8 T"      = "#f5b7b1",
  "Effector CD8 T"   = "#ec7063",
  "Memory CD8 T"     = "#a93226",
  "Other T"          = "#d5d8dc",
  "NK1"              = "#6cc4b9",
  "NK2"              = "#1e968d",
  "B"                = "#63bcd4",
  "Plasma B"         = "#208db5",
  "Myeloid"          = "#8b4513",
  "Endothelial"      = "#6c5ce7",
  "Fibroblasts"      = "#b7b75e",
  "Tumor 1"          = "#f8bbd0",
  "Tumor 2"          = "#f06292",
  "Tumor 3"          = "#d81b60",
  "Tumor 4"          = "#880e4f"
)

tumor_tcell_order <- c(
  "Treg 1", "Treg 2", "Treg 3", "Treg 4",
  "Early TEx", "Intermediate TEx", "Terminal TEx",
  "Th1", "Th17", "Tfh 1", "Tfh 2",
  "Naive CD4 T", "Activated CD4 T", "Memory CD4 T",
  "Naive CD8 T", "Effector CD8 T", "Memory CD8 T",
  "Other T",
  "NK1", "NK2",
  "B", "Plasma B",
  "Myeloid",
  "Endothelial",
  "Fibroblasts",
  "Tumor 1", "Tumor 2", "Tumor 3", "Tumor 4"
)


# ── Shared cell-type group definitions ────────────────────────────────────────

cancer_cell_types  <- c("Tumor 1", "Tumor 2", "Tumor 3", "Tumor 4")

tcell_types <- c(
  "Treg 1", "Treg 2", "Treg 3", "Treg 4",
  "Early TEx", "Intermediate TEx", "Terminal TEx",
  "Th1", "Th17", "Tfh 1", "Tfh 2", "Other T",
  "Memory CD8 T", "Naive CD4 T", "Naive CD8 T",
  "Effector CD8 T", "Activated CD4 T", "Memory CD4 T"
)

immune_cell_types  <- c("B", "NK1", "NK2", "Myeloid", "Plasma B")
stromal_cell_types <- c("Endothelial", "Fibroblasts")

# Per-cell-type color lookup (used in lollipop axis coloring)
cell_colors <- c(
  setNames(rep("#B2182B", length(cancer_cell_types)),  cancer_cell_types),
  setNames(rep("#1B9E77", length(tcell_types)),         tcell_types),
  setNames(rep("#2166AC", length(immune_cell_types)),   immune_cell_types),
  setNames(rep("#969696", length(stromal_cell_types)),  stromal_cell_types)
)

# Broad cell-type color map (used in density/boxplots)
broad_cell_type_colors <- c(
  "Cancer cell"  = "#B2182B",
  "Immune cell"  = "#2166AC",
  "T cell"       = "#1B9E77",
  "Stromal cell" = "#969696"
)

# Response group colors
response_colors <- c("Non-responder" = "#D6604D", "Responder" = "#4393C3")

# Amplification bias reference lines and legend points
amplification_bias_df  <- tibble(y = c(152, 304, 456, 608), f = c("1", "2", "3", "4"))
amplification_legend_df <- data.frame(f = factor(c("1", "2", "3", "4")))

amplification_bias_colors <- c(
  "1" = "#f9a71bff",
  "2" = "#ddc729ff",
  "3" = "#39a14eff",
  "4" = "#016f5aff"
)


# ── Helper: classify cells into broad groups ───────────────────────────────────

assign_broad_cell_type <- function(cell_sub_type) {
  case_when(
    cell_sub_type %in% cancer_cell_types  ~ "Cancer cell",
    cell_sub_type %in% tcell_types        ~ "T cell",
    cell_sub_type %in% immune_cell_types  ~ "Immune cell",
    cell_sub_type %in% stromal_cell_types ~ "Stromal cell",
    TRUE ~ NA_character_
  )
}

# ── Helper: compute openness-adjusted telomere content ─────────────────────────
# Regresses log2(telomere content) against openness and returns
# back-transformed residuals re-centred at the original mean.

adjust_telomere_for_openness <- function(df) {
  df %>%
    mutate(
      log_tel             = log2(tel_content),
      telomere_content_corr = 2^(residuals(lm(log_tel ~ openness)) +
                                   mean(log_tel, na.rm = TRUE))
    )
}

# ── Helper: recode & factor pre/post treatment labels ─────────────────────────

recode_pre_post <- function(df) {
  df %>%
    mutate(
      pre_post = recode(
        pre_post,
        "pre treatment"  = "Pre-treatment",
        "post treatment" = "Post-treatment"
      ),
      pre_post = factor(pre_post, levels = c("Pre-treatment", "Post-treatment"))
    )
}


# ── Load & prepare Satpathy data ──────────────────────────────────────────────

satpathy_df <- read.csv("~/path/to/project/Satpathy_masterdf_Chromatinstate.csv") %>%
  mutate(openness = log2(X13_Heterochrom / Active_Sites))


# =============================================================================
# Figure 1
# =============================================================================

# ── PF1b: Detection probability vs. total reads by amplification bias ─────────

telomere_length_fixed <- 10000
total_reads_vals      <- c(2000, 5000, 10000, 20000, 30000, 50000, 70000, 100000)
rb_vals               <- seq(1, 4, by = 1)
rb_colors             <- c("#F9A71A", "#dcc729ff", "#3aa14dff", "#026F5A")

# Probability of observing at least n_intra_val intratelomeric reads
p_tel_calc <- function(telomere_length, ratio_bias) {
  ((telomere_length * 2 * 23) / (3.2e9)) * ratio_bias
}

generate_rb_data <- function(n_intra_val) {
  expand.grid(Total_Reads = total_reads_vals, RB = rb_vals) %>%
    mutate(
      p           = p_tel_calc(telomere_length_fixed, RB),
      Probability = 1 - pbinom(n_intra_val - 1, size = Total_Reads, prob = p)
    )
}

PF1b <- generate_rb_data(5) %>%
  ggplot(aes(Total_Reads, Probability, color = factor(RB))) +
  geom_vline(xintercept = 10000, color = "#000000", linewidth = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 30000, color = "#000000", linewidth = 0.5, linetype = "dashed") +
  geom_point(size = 3) +
  geom_line(aes(group = RB), linewidth = 1.2) +
  scale_x_log10(
    breaks = total_reads_vals,
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  scale_color_manual(values = rb_colors, name = "Amplification\nbias") +
  labs(x = expression(Log[10] ~ "total reads"), y = "Probability") +
  theme_bw() +
  theme(
    text             = element_text(size = 15),
    legend.position  = c(.2, .8),
    legend.background = element_rect(color = "lightgrey", fill = "white", linewidth = .5)
  )


# ── PF1c: Distribution of total reads per cell with 30k threshold ─────────────

masterdataframe_all_cells <- readRDS("~/path/to/project/masterdataframe_all_cells.rds")

# Pre-compute annotation labels for the two sides of the threshold
label_df <- masterdataframe_all_cells %>%
  summarise(
    left_n  = sum(total_reads <= 30000),
    right_n = sum(total_reads > 30000),
    total_n = n()
  ) %>%
  mutate(
    left_text  = sprintf("%s cells\n(%.1f%%)", left_n,  left_n  / total_n * 100),
    right_text = sprintf("%s cells\n(%.1f%%)", right_n, right_n / total_n * 100)
  )

PF1c <- masterdataframe_all_cells %>%
  ggplot(aes(total_reads)) +
  geom_histogram(bins = 100, fill = "darkgrey") +
  scale_x_log10(labels = label_number(scale_cut = cut_short_scale())) +
  annotation_logticks(side = "b") +
  geom_vline(xintercept = 30000, linetype = "dashed") +
  geom_text(data = label_df, aes(x = 9000,   y = 1000, label = left_text),  inherit.aes = FALSE, size = 5) +
  geom_text(data = label_df, aes(x = 500000, y = 1000, label = right_text), inherit.aes = FALSE, size = 5) +
  labs(x = "Total reads per cell", y = "Number of cells") +
  theme_bw() +
  theme(text = element_text(size = 15))


# ── PF1e: iTRPM by broad cell type ────────────────────────────────────────────

PF1e <- satpathy_df %>%
  mutate(
    new_cell_type = assign_broad_cell_type(cell_sub_type),
    itrpm         = intratel_reads / total_reads * 1e6,
    new_cell_type = fct_reorder(new_cell_type, itrpm, .fun = median, .desc = TRUE)
  ) %>%
  ggplot(aes(new_cell_type, itrpm, fill = new_cell_type)) +
  geom_boxplot() +
  geom_hline(data = amplification_bias_df, aes(yintercept = y, color = f),
             linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  geom_hline(yintercept = 200,  linetype = "dashed") +
  geom_hline(yintercept = 1000, linetype = "dashed") +
  geom_point(data = amplification_legend_df,
             aes(x = -Inf, y = -Inf, color = f),
             shape = 15, size = 6, inherit.aes = FALSE) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  scale_fill_manual(values  = broad_cell_type_colors) +
  scale_color_manual(name   = "Amplification\nbias", values = amplification_bias_colors) +
  labs(x = "Cell type", y = "Intratelomeric reads per million (iTRPM)", color = "Amplification bias") +
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(
    text         = element_text(size = 15),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.background = element_blank()
  )


# ── PF1d: iTRPM by patient and treatment timepoint ────────────────────────────

PF1d <- satpathy_df %>%
  recode_pre_post() %>%
  mutate(
    itrpm           = intratel_reads / total_reads * 1e6,
    patientID_panel = if_else(
      pre_post == "Pre-treatment",
      paste0(patientID, "_pre"),
      paste0(patientID, "_post")
    )
  ) %>%
  group_by(pre_post) %>%
  mutate(patientID_panel = fct_reorder(patientID_panel, itrpm, .fun = median, .desc = TRUE)) %>%
  ungroup() %>%
  ggplot(aes(patientID_panel, itrpm, fill = response)) +
  geom_boxplot() +
  geom_hline(data = amplification_bias_df, aes(yintercept = y, color = f),
             linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  geom_hline(yintercept = 200,  linetype = "dashed") +
  geom_hline(yintercept = 1000, linetype = "dashed") +
  geom_point(data = amplification_legend_df,
             aes(x = -Inf, y = -Inf, color = f),
             shape = 15, size = 6, inherit.aes = FALSE) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  facet_wrap(~pre_post, scales = "free_x", space = "free_x") +
  scale_x_discrete(labels = function(x) sub("(_pre|_post)$", "", x)) +
  scale_fill_manual(values  = response_colors) +
  scale_color_manual(name   = "Amplification\nbias", values = amplification_bias_colors) +
  labs(
    x     = "PatientID",
    y     = "Intratelomeric reads per million (iTRPM)",
    color = "Amplification bias",
    fill  = NULL
  ) +
  theme_bw() +
  theme(
    text         = element_text(size = 15),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.background = element_blank()
  )


# ── PF1f: iTRPM by cell subtype ───────────────────────────────────────────────

PF1f <- satpathy_df %>%
  mutate(
    itrpm        = intratel_reads / total_reads * 1e6,
    cell_sub_type = fct_reorder(cell_sub_type, itrpm, .fun = median, .desc = TRUE)
  ) %>%
  ggplot(aes(cell_sub_type, itrpm, fill = cell_sub_type)) +
  geom_boxplot() +
  geom_hline(data = amplification_bias_df, aes(yintercept = y, color = f),
             linetype = "dashed", linewidth = 0.8, show.legend = FALSE) +
  geom_hline(yintercept = 200,  linetype = "dashed") +
  geom_hline(yintercept = 1000, linetype = "dashed") +
  geom_point(data = amplification_legend_df,
             aes(x = -Inf, y = -Inf, color = f),
             shape = 15, size = 6, inherit.aes = FALSE) +
  scale_y_log10() +
  annotation_logticks(sides = "l") +
  scale_fill_manual(values  = tumor_tcell_colors) +
  scale_color_manual(name   = "Amplification\nbias", values = amplification_bias_colors) +
  labs(x = "Cell sub type", y = "Intratelomeric reads per million (iTRPM)", color = "Amplification bias") +
  guides(fill = "none", color = "none") +
  theme_bw() +
  theme(
    text         = element_text(size = 15),
    axis.text.x  = element_text(angle = 45, hjust = 1),
    strip.background = element_blank()
  )


# ── Assemble & save Figure 1 ──────────────────────────────────────────────────

PF1 <- plot_grid(
  plot_grid(ggplot() + theme_void(), PF1b, PF1c,
            nrow = 1, labels = c("a", "b", "c")),
  plot_grid(PF1d, PF1e, PF1f,
            nrow = 1, rel_widths = c(.4, .12, .48),
            axis = "tb", align = "hv", labels = c("d", "e", "f")),
  nrow = 2, rel_heights = c(.45, .55)
)

ggsave(
  filename = "PF1.pdf",
  plot     = PF1,
  path     = "~/path/to/project/img/",
  width    = 15, height = 11, units = "in", dpi = 600
)


# =============================================================================
# Figure 2
# =============================================================================

# ── PF2a: Correlation matrix of chromatin state read counts ───────────────────

satpathy_joined_all_chromatin_regions <- read.csv(
  "~/path/to/project/Satpathy_joined_all_Chromatin_Regions.csv"
)
corr_plot_chromatin_df <- read.csv(
  "~/path/to/project/Corr_plot_Chromatin_df.csv"
)

corr_mat     <- cor(corr_plot_chromatin_df, use = "pairwise.complete.obs", method = "pearson")
melted_corr  <- melt(corr_mat)

strip_x <- function(x) sub("^X", "", x)

# Color axis labels: red = active states (1–7), blue = heterochromatin (13)
axis_colors_PF2a <- sapply(levels(melted_corr$Var1), function(x) {
  if (grepl("^X[1-7]_", x)) "red"
  else if (grepl("^X13_", x)) "blue"
  else "grey40"
})

PF2a <- melted_corr %>%
  ggplot(aes(Var1, Var2, fill = value, size = abs(value))) +
  geom_point(pch = 21) +
  scale_fill_gradient2(
    name     = "Pearson\nCorrelation\nCoefficient",
    low      = "blue",
    high     = "red",
    mid      = "white",
    midpoint = 0,
    limit    = c(-1, 1)
  ) +
  scale_size_continuous(range = c(1, 8)) +
  guides(
    size = "none",
    fill = guide_colorbar(
      barheight      = unit(15, "lines"),
      frame.colour   = "black", frame.linewidth = 0.25,
      ticks.colour   = "black", ticks.linewidth = 0.25
    )
  ) +
  scale_x_discrete(labels = strip_x) +
  scale_y_discrete(labels = strip_x) +
  theme_bw() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, colour = axis_colors_PF2a),
    axis.text.y  = element_text(colour = axis_colors_PF2a[levels(melted_corr$Var2)]),
    axis.title   = element_blank(),
    text         = element_text(size = 15)
  )


# ── PF2b: Scatter of heterochromatin vs. active-site reads, coloured by openness

PF2b <- satpathy_df %>%
  ggplot(aes(X13_Heterochrom, Active_Sites, col = openness)) +
  geom_point() +
  scale_color_gradient2(low = "#1a9850", mid = "lightgrey", high = "#762a83", midpoint = 0) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Heterochromatin reads", y = "Reads from active sites", col = "Openness") +
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = "bottom") +
  guides(
    col = guide_colorbar(
      barwidth       = unit(15, "lines"),
      frame.colour   = "black", frame.linewidth = 0.25,
      ticks.colour   = "black", ticks.linewidth = 0.25
    )
  )


# ── PF2c: Openness density by broad cell type, split by treatment timepoint ───

PF2c <- satpathy_df %>%
  recode_pre_post() %>%
  mutate(new_cell_type = assign_broad_cell_type(cell_sub_type)) %>%
  ggplot(aes(openness, col = new_cell_type, fill = new_cell_type)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~pre_post) +
  scale_fill_manual(values  = broad_cell_type_colors) +
  scale_color_manual(values = broad_cell_type_colors) +
  labs(x = "Openness", y = "Density", col = "Cell type", fill = "Cell type") +
  theme_bw() +
  theme(
    text             = element_text(size = 15),
    legend.position  = "bottom",
    strip.background = element_blank()
  )


# ── PF2d: Volcano plot — change in openness pre vs. post treatment ────────────

openness_change_df <- satpathy_df %>%
  group_by(cell_sub_type) %>%
  summarise(
    diff = median(openness[pre_post == "post treatment"]) -
      median(openness[pre_post == "pre treatment"]),
    p    = wilcox.test(openness ~ pre_post)$p.value,
    .groups = "drop"
  ) %>%
  mutate(
    p_adj = p.adjust(p, method = "BH"),
    sig   = p_adj < 0.01 & diff > 0.5
  )

PF2d <- openness_change_df %>%
  ggplot(aes(diff, -log10(p), col = sig)) +
  geom_point(size = 3) +
  geom_text_repel(
    data = subset(openness_change_df, sig),
    aes(label = cell_sub_type),
    size = 5, max.overlaps = Inf, box.padding = 0.5, point.padding = 0.3
  ) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  labs(
    x = "Change in Openness (Post – Pre)",
    y = expression(-Log[10] ~ "p-value")
  ) +
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = "none")


# ── PF2e: Lollipop — per-subtype openness correlation ─────────────────────────

openness_corr_df <- satpathy_df %>%
  group_by(cell_sub_type) %>%
  summarise(
    openness_correlation = abs(cor(X13_Heterochrom, Active_Sites,
                                   method = "pearson", use = "complete.obs")),
    .groups = "drop"
  ) %>%
  mutate(cell_sub_type = fct_reorder(cell_sub_type, openness_correlation))

PF2e <- openness_corr_df %>%
  ggplot(aes(openness_correlation, cell_sub_type)) +
  geom_segment(aes(x = 0, xend = openness_correlation,
                   y = cell_sub_type, yend = cell_sub_type), color = "grey") +
  geom_point(size = 3, color = "black") +
  labs(x = "Openness correlation", y = "Cell sub type") +
  theme_bw() +
  theme(
    text         = element_text(size = 15),
    axis.text.y  = element_text(
      color = cell_colors[levels(fct_reorder(openness_corr_df$cell_sub_type,
                                             openness_corr_df$openness_correlation))]
    )
  )


# ── PF2f & PF2g: Heterochromatin/active-site scatter & openness vs. telomere ──
# Restricted to high-coverage cells (≥30k reads)

options(scipen = 10000)

PF2f <- satpathy_df %>%
  filter(total_reads >= 30000) %>%
  arrange(tel_content) %>%
  ggplot(aes(X13_Heterochrom, Active_Sites, col = tel_content)) +
  geom_point() +
  scale_color_viridis_c(
    option = "A",
    limits = c(0, 25000),
    oob    = scales::squish,
    labels = label_number(scale_cut = cut_short_scale())
  ) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Heterochromatin reads", y = "Reads from active sites", col = "Telomere content") +
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = "bottom") +
  guides(
    col = guide_colorbar(
      barwidth       = unit(15, "lines"),
      frame.colour   = "black", frame.linewidth = 0.25,
      ticks.colour   = "black", ticks.linewidth = 0.25
    )
  )

PF2g <- satpathy_df %>%
  filter(total_reads >= 30000) %>%
  ggplot(aes(openness, tel_content)) +
  geom_point() +
  stat_cor(method = "spearman", p.accuracy = 0.0001) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Openness", y = "Telomere content") +
  theme_bw() +
  theme(text = element_text(size = 15))


# ── Assemble & save Figure 2 ──────────────────────────────────────────────────

PF2 <- plot_grid(
  plot_grid(PF2a, PF2b, PF2c, PF2d,
            align = "hv", labels = c("a", "b", "c", "d"),
            rel_widths = c(.55, .45), rel_heights = c(.6, .46)),
  plot_grid(PF2e,
            plot_grid(PF2f, PF2g, ncol = 1, labels = c("f", "g")),
            labels = c("e", "")),
  ncol = 1, rel_heights = c(.55, .45)
)

ggsave(
  filename = "PF2.pdf",
  plot     = PF2,
  path     = "~/path/to/project/img/",
  width    = 12, height = 18, units = "in", dpi = 600
)


# =============================================================================
# Supplementary Figure 2 (PF2S): Openness-adjusted telomere content validation
# =============================================================================

# PF2Sa: Spearman correlation of openness vs. openness-adjusted telomere content
PF2Sa <- satpathy_df %>%
  filter(total_reads >= 30000, tel_content > 0, !is.na(openness)) %>%
  mutate(
    log_tel               = log2(tel_content),
    telomere_content_corr = residuals(lm(log_tel ~ openness)) + mean(log_tel, na.rm = TRUE)
  ) %>%
  ggplot(aes(openness, telomere_content_corr)) +
  geom_point() +
  stat_cor(
    method       = "spearman",
    p.accuracy   = 0.001,
    p.label      = "p.signif",
    label.x.npc  = "left",
    label.y.npc  = "top"
  ) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Openness", y = "Openness-adjusted telomere content") +
  theme_bw() +
  theme(text = element_text(size = 15))

# PF2Sb & PF2Sc: Raw vs. adjusted telomere content distributions by open/closed chromatin
open_closed_df <- satpathy_df %>%
  filter(total_reads >= 30000, tel_content <= 20000, tel_content > 0, !is.na(openness)) %>%
  mutate(
    log_tel               = log2(tel_content),
    telomere_content_corr = 2^(residuals(lm(log_tel ~ openness)) + mean(log_tel, na.rm = TRUE)),
    openness_f            = if_else(openness >= 0, "open", "closed")
  )

open_closed_colors <- c("closed" = "#1a9850", "open" = "#762a83")

PF2Sb <- open_closed_df %>%
  ggplot(aes(tel_content, col = openness_f, fill = openness_f)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = open_closed_colors) +
  scale_fill_manual(values  = open_closed_colors) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Telomere content", y = "Density", fill = NULL, col = NULL) +
  theme_bw() +
  theme(
    text              = element_text(size = 15),
    legend.position   = c(.8, .8),
    legend.background = element_rect(color = "lightgrey", fill = "white", linewidth = .5)
  )

PF2Sc <- open_closed_df %>%
  ggplot(aes(telomere_content_corr, col = openness_f, fill = openness_f)) +
  geom_density(alpha = 0.5) +
  scale_color_manual(values = open_closed_colors) +
  scale_fill_manual(values  = open_closed_colors) +
  scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Openness-adjusted telomere content", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 15), legend.position = "none")

PF2S <- plot_grid(PF2Sa, PF2Sb, PF2Sc, nrow = 1)

ggsave(
  filename = "PF2S.pdf",
  plot     = PF2S,
  path     = "~/path/to/project/img/",
  width    = 14, height = 5, units = "in", dpi = 600
)


# =============================================================================
# Figure 4
# =============================================================================

# Shared pre-processing: filter, recode, and compute adjusted telomere content
satpathy_adj <- satpathy_df %>%
  filter(total_reads >= 30000, tel_content > 0, !is.na(openness)) %>%
  recode_pre_post() %>%
  mutate(new_cell_type = assign_broad_cell_type(cell_sub_type)) %>%
  adjust_telomere_for_openness()


# ── PF4a: Adjusted telomere content by response, broad cell type, and timepoint

PF4a <- satpathy_adj %>%
  ggplot(aes(response, telomere_content_corr, fill = response)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  scale_fill_manual(values = response_colors) +
  facet_grid(pre_post ~ new_cell_type) +
  coord_cartesian(ylim = c(0, 20000)) +
  stat_compare_means(
    method       = "wilcox.test",
    paired       = FALSE,
    label        = "p.signif",
    size         = 5,
    label.y      = 14000,
    label.x      = 1.3,
    comparisons  = list(c("Non-responder", "Responder"))
  ) +
  labs(x = "Response", y = "Openness-adjusted telomere content") +
  theme_bw() +
  theme(
    legend.position  = "none",
    text             = element_text(size = 15),
    strip.background = element_blank(),
    axis.text.x      = element_text(angle = 45, vjust = 1, hjust = 1)
  )


# ── PF4b: Volcano plot — adjusted telomere content, responders vs. non-responders

# Helper: run Wilcoxon tests per cell subtype and attach log2FC
run_wilcoxon <- function(df, label) {
  logfc_df <- df %>%
    group_by(cell_sub_type, response) %>%
    summarise(median_tel = median(tel_content_corr, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = response, values_from = median_tel) %>%
    mutate(log2FC = log2(Responder / `Non-responder`)) %>%
    select(cell_sub_type, log2FC)
  
  df %>%
    group_by(cell_sub_type) %>%
    rstatix::wilcox_test(tel_content_corr ~ response) %>%
    adjust_pvalue(method = "hochberg") %>%
    left_join(logfc_df, by = "cell_sub_type") %>%
    mutate(stage = label)
}

# Prepare data with `tel_content_corr` (naming convention expected by run_wilcoxon)
satpathy_volcano <- satpathy_df %>%
  filter(total_reads >= 30000, tel_content > 0, !is.na(openness)) %>%
  recode_pre_post() %>%
  mutate(
    log_tel          = log2(tel_content),
    tel_content_corr = 2^(residuals(lm(log_tel ~ openness)) + mean(log_tel, na.rm = TRUE))
  )

pre_results <- satpathy_volcano %>%
  filter(pre_post == "Pre-treatment") %>%
  group_by(cell_sub_type) %>%
  filter(n_distinct(response) == 2) %>%
  ungroup() %>%
  run_wilcoxon("Pre-treatment")

post_results <- satpathy_volcano %>%
  filter(pre_post == "Post-treatment") %>%
  group_by(cell_sub_type) %>%
  filter(n_distinct(response) == 2) %>%
  ungroup() %>%
  run_wilcoxon("Post-treatment")

all_results <- bind_rows(pre_results, post_results)

PF4b <- all_results %>%
  mutate(
    stage = factor(stage, levels = c("Pre-treatment", "Post-treatment")),
    sig   = ifelse(p.adj < 0.01 & abs(log2FC) >= 0.5, "significant", "not_significant")
  ) %>%
  ggplot(aes(log2FC, -log10(p), color = sig)) +
  geom_point(size = 2) +
  scale_color_manual(values = c("significant" = "red", "not_significant" = "grey")) +
  geom_label_repel(
    data  = . %>% filter(p.adj < 0.01 & abs(log2FC) >= 0.5),
    aes(label = cell_sub_type),
    color = "black", label.size = NA
  ) +
  coord_cartesian(xlim = c(-4, 4)) +
  facet_wrap(~stage) +
  labs(
    x = expression("Openness-adjusted telomere content " * Log[2] ~ "FC (Responder / Non-responder)"),
    y = expression("-" * Log[10] ~ "p-value")
  ) +
  theme_bw() +
  theme(
    legend.position  = "none",
    text             = element_text(size = 15),
    strip.background = element_blank()
  )


# ── PF4c: UMAP coloured by cell subtype, faceted by response and timepoint ────

PF4c <- satpathy_df %>%
  mutate(
    cell_sub_type = factor(cell_sub_type, levels = tumor_tcell_order)
  ) %>%
  recode_pre_post() %>%
  ggplot(aes(UMAP1, UMAP2, col = cell_sub_type)) +
  geom_point(size = 0.7) +
  scale_color_manual(values = tumor_tcell_colors) +
  facet_grid(pre_post ~ response) +
  coord_equal() +
  labs(col = "") +
  theme_classic() +
  theme(
    legend.position  = "none",
    axis.title       = element_blank(),
    axis.text        = element_blank(),
    axis.ticks       = element_blank(),
    axis.line        = element_blank(),
    text             = element_text(size = 15),
    strip.background = element_blank()
  )


# ── PF4d: Change in cell count vs. change in telomere content (SU008/SU009) ───

PF4d <- satpathy_adj %>%
  filter(patientID %in% c("SU008", "SU009")) %>%
  mutate(cell_sub_type = factor(cell_sub_type, levels = tumor_tcell_order)) %>%
  group_by(patientID, response, cell_sub_type, pre_post) %>%
  summarise(
    cell_count = n(),
    med_tel    = median(telomere_content_corr),
    .groups    = "drop"
  ) %>%
  filter(cell_count >= 10) %>%
  pivot_longer(cols = c(cell_count, med_tel)) %>%
  pivot_wider(names_from = "pre_post", values_from = "value") %>%
  na.omit() %>%
  mutate(diff_value = `Post-treatment` - `Pre-treatment`) %>%
  select(-c(`Post-treatment`, `Pre-treatment`)) %>%
  pivot_wider(names_from = "name", values_from = "diff_value") %>%
  mutate(
    significance = if_else(
      cell_sub_type %in% c("Intermediate TEx", "Terminal TEx"),
      "Significant", "Insignificant        "
    )
  ) %>%
  ggplot(aes(cell_count, med_tel, fill = cell_sub_type, shape = significance)) +
  geom_hline(yintercept = 0, col = "lightgrey") +
  geom_vline(xintercept = 0, col = "lightgrey") +
  geom_point(size = 3) +
  facet_wrap(response ~ patientID, scales = "free", nrow = 2) +
  scale_fill_manual(values = tumor_tcell_colors, drop = FALSE) +
  scale_shape_manual(values = c("Significant" = 24, "Insignificant        " = 21)) +
  labs(
    x = "Change in cell count (Post - Pre)",
    y = "Change in median openness-adjusted telomere content (Post - Pre)"
  ) +
  guides(
    fill  = guide_legend(ncol = 2, override.aes = list(shape = 21, size = 5)),
    shape = guide_legend(ncol = 2, override.aes = list(size = 5))
  ) +
  theme_bw() +
  theme(
    legend.title    = element_blank(),
    legend.spacing.y = unit(0.25, "lines"),
    legend.margin   = margin(0, 0, 0, 0),
    text            = element_text(size = 15),
    strip.background = element_blank()
  )


# ── Assemble & save Figure 4 ──────────────────────────────────────────────────

PF4 <- plot_grid(
  PF4a, PF4b, PF4c, PF4d,
  labels = c("a", "b", "c", "d"),
  rel_heights = c(.45, .55)
)

ggsave(
  filename = "PF4.pdf",
  plot     = PF4,
  path     = "~/path/to/project/img/",
  width    = 15, height = 13, units = "in", dpi = 600
)


# =============================================================================
# TCGA supplementary figures (PF3S, PF4S, PF5S)
# =============================================================================

# ── Load & combine TCGA cohort data ───────────────────────────────────────────

tcga_cohorts <- list(
  BRCA = "BRCA_joined_Heterochrom_Reg_Regions.csv",
  BLCA = "BLCA_joined_Heterochrom_Reg_Regions.csv",
  GBMx = "GBMx_joined_Heterochrom_Reg_Regions.csv",
  KIRP = "KIRP_joined_Heterochrom_Reg_Regions.csv",
  KIRC = "KIRC_joined_Heterochrom_Reg_Regions.csv",
  SKCM = "SKCM_joined_Heterochrom_Reg_Regions.csv",
  COAD = "COAD_joined_Heterochrom_Reg_Regions.csv",
  LUAD = "LUAD_joined_Heterochrom_Reg_Regions.csv"
)

tcga_df <- map_dfr(
  tcga_cohorts,
  ~ read.csv(file.path("~/path/to/project", .x))
) %>%
  mutate(openness = log2(X13_Heterochrom / Active_Sites))

# Shared filter and aesthetics for all three TCGA supplementary figures
tcga_base_theme <- list(
  theme_bw(),
  theme(
    text             = element_text(size = 15),
    legend.position  = "bottom",
    strip.background = element_blank()
  ),
  guides(
    col = guide_colorbar(
      barwidth       = unit(15, "lines"),
      frame.colour   = "black", frame.linewidth = 0.25,
      ticks.colour   = "black", ticks.linewidth = 0.25
    )
  )
)

openness_scatter <- function(df) {
  df %>%
    ggplot(aes(X13_Heterochrom, Active_Sites, col = openness)) +
    geom_point() +
    stat_cor(method = "pearson", p.accuracy = 0.0001) +
    scale_color_gradient2(low = "#1a9850", mid = "lightgrey", high = "#762a83", midpoint = 0) +
    scale_x_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
    labs(x = "Heterochromatin reads", y = "Reads from active sites", col = "Openness")
}

tcga_filtered <- tcga_df %>%
  filter(!is.na(Cohort), total_reads >= 30000)


# PF3S: Openness scatter faceted by TCGA cohort (cancer cells only)
PF3S <- tcga_filtered %>%
  filter(celltype_supplement == "Cancer cell") %>%
  openness_scatter() +
  facet_wrap(~Cohort, nrow = 2) +
  tcga_base_theme

# PF4S: Openness scatter faceted by cell type (all cells)
PF4S <- tcga_filtered %>%
  openness_scatter() +
  facet_wrap(~celltype_supplement, nrow = 2) +
  tcga_base_theme

# PF5S: Openness vs. telomere content, faceted by cohort (cancer cells only)
PF5S <- tcga_filtered %>%
  filter(celltype_supplement == "Cancer cell") %>%
  ggplot(aes(openness, telomere_content)) +
  geom_point() +
  stat_cor(method = "spearman", p.accuracy = 0.0001) +
  scale_y_continuous(labels = label_number(scale_cut = cut_short_scale())) +
  labs(x = "Openness", y = "Telomere content") +
  facet_wrap(~Cohort, scales = "free_y", nrow = 2) +
  theme_bw() +
  theme(text = element_text(size = 15), strip.background = element_blank())


for (fig_name in c("PF3S", "PF4S", "PF5S")) {
  ggsave(
    filename = paste0(fig_name, ".pdf"),
    plot     = get(fig_name),
    path     = "~/path/to/project/img/",
    width    = 12, height = 8, units = "in", dpi = 600
  )
}
