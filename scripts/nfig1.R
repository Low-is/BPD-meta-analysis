# Figure 1: study design and cohorts.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})
OUT <- "./figs"
dir.create(OUT, showWarnings = FALSE)
FONT <- "Liberation Sans"
INK <- "#222222"
MUTED <- "#5f5f5f"
C_BPD <- "#1f3a6e"
C_NO <- "#c8553d"
LIGHTL <- "#b5b3ab"
BOXF <- "#f4f6fa"
BOXC <- "#9fb0cc"
KEYF <- "#e4e9f2"
TESTF <- "#fbefe9"
TESTC <- "#d79a86"
TH <- theme_classic(base_size = 7, base_family = FONT) +
  theme(text = element_text(colour = INK), axis.text = element_text(colour = MUTED, size = 6.2),
    axis.title = element_text(size = 6.8),
    axis.line = element_line(colour = MUTED, linewidth = .3),
    axis.ticks = element_line(colour = MUTED, linewidth = .3),
    plot.title = element_text(size = 7.6, face = "bold", hjust = 0, margin = margin(b = 4)),
    legend.text = element_text(size = 6.2), legend.key.size = unit(7, "pt"),
    strip.background = element_blank(),
    plot.tag = element_text(size = 10, face = "bold", family = FONT))

bx <- function(x, y, w, h, lab, fill = BOXF, col = BOXC, face = "plain", size = 2.2)
list(annotate("rect", xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
  fill = fill, colour = col, linewidth = .35),
  annotate("text", x = x, y = y, label = lab, size = size, lineheight = .95, fontface = face,
    family = FONT, colour = INK))

ar <- function(x1, y1, x2, y2) annotate("segment", x = x1, y = y1, xend = x2, yend = y2, colour = MUTED,
  linewidth = .35,
  arrow = arrow(length = unit(3, "pt"), type = "closed"))

lab_side <- function(y, txt) annotate("text", x = .1, y = y, label = txt, hjust = 0, size = 2.25,
  fontface = "bold", colour = C_BPD, family = FONT)
fA <- ggplot() + coord_cartesian(xlim = c(0, 10), ylim = c(.2, 12.4),
  expand = FALSE) + theme_void(base_family = FONT) +
  bx(5.9, 11.75, 7.8, 1.0,
    "4 independent first-week blood cohorts\n196 infants, 89 with BPD; 2 followed to day 28",
    fill = KEYF, col = C_BPD, face = "bold") +
  lab_side(11.75, "Data") + ar(5.9, 11.23, 5.9, 10.83) +
  bx(5.9, 9.9, 7.8, 0.95,
    "Multicohort meta-analysis\nsignature, pathways, four biological programs") + lab_side(9.9,
      "1") + ar(5.9, 9.41, 5.9, 9.01) +
  bx(3.3, 8.25, 2.45, 1.3, "Held-out\ngenes chosen\nwithout the\ncohort tested", fill = TESTF,
    col = TESTC, size = 2.0) +
  bx(5.9, 8.25, 2.45, 1.3, "GA matching\nand adjustment", fill = TESTF, col = TESTC, size = 2.0) +
  bx(8.5, 8.25, 2.45, 1.3, "Cell\ncomposition", fill = TESTF, col = TESTC, size = 2.0) + lab_side(8.25,
    "2") + ar(5.9, 7.58, 5.9, 7.18) +
  bx(5.9, 6.0, 7.8, 1.05,
    "What do published signatures measure?\n15 signatures, reconstructed and scored only\nin cohorts not used to derive them", size = 2.05) + lab_side(6.0, "3") + ar(5.9, 5.48, 5.9, 5.08) +
  bx(5.9, 4.25, 7.8, 1.25,
    "GA-adjusted discovery (exploratory)\nall time points in each longitudinal cohort; genes from one\ncohort tested in the other; pre-specified confirmatory test", fill = KEYF, col = C_BPD, size = 2.05) +
  lab_side(4.25, "4") + ar(5.9, 3.6, 5.9, 3.2) +
  bx(5.9, 2.5, 7.8, 1.1,
    "Prediction beyond gestational age\nfully out of sample; GA gray zone 27-30 weeks", fill = "#eeeeea",
    col = "#b5b3ab") +
  labs(title = "Study design") + theme(plot.title = element_text(size = 7.6, face = "bold",
    family = FONT, margin = margin(b = 4)),
    plot.tag = element_text(size = 10, face = "bold", family = FONT))
ct <- tribble(~cohort, ~plat, ~tissue, ~bpd, ~no, ~ga, ~role,
  "GSE32472", "Microarray", "Periph. blood", 62, 35, 28.0, "Meta-analysis",
  "GSE220135", "RNA-seq", "Cord blood", 13, 54, 30.0, "Meta-analysis",
  "GSE225881", "Microarray", "Cord blood", 9, 13, 29.1, "Meta-analysis",
  "GSE125873", "RNA-seq", "Periph. blood", 5, 5, NA, "Meta-analysis",
  "GSE302178", "RNA-seq", "Cord blood", 8, 8, NA, "External") |>
  mutate(role = factor(role, levels = c("Meta-analysis", "External")),
    lab = paste0(cohort, "  ", plat, "\n", tissue), lab = factor(lab, levels = rev(lab)),
    long = c("Longitudinal to day 28\n(294 samples)", "Longitudinal to day 28\n(200 samples)", "", "",
      ""))
d <- ct |> pivot_longer(c(bpd, no), names_to = "grp", values_to = "n") |>
  mutate(grp = factor(ifelse(grp == "bpd", "BPD", "No BPD"), levels = c("No BPD", "BPD")))
fB <- ggplot(d, aes(n, lab, fill = grp)) +
  geom_col(width = .68, colour = "white", linewidth = .35) +
  geom_text(aes(label = ifelse(n >= 5, n, "")), position = position_stack(vjust = .5), colour = "white",
    size = 2.05) +
  geom_text(data = ct,
    aes(x = bpd + no + 2, y = lab, label = paste0(ifelse(is.na(ga), "", sprintf("GA %.1f wk", ga)),
      ifelse(long != "", paste0(ifelse(is.na(ga), "", "\n"), long), ""))), inherit.aes = FALSE,
    hjust = 0, size = 1.9, colour = MUTED, lineheight = .9) +
  facet_grid(role ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c("No BPD" = C_NO, "BPD" = C_BPD), name = NULL,
    guide = guide_legend(reverse = TRUE)) +
  scale_x_continuous(expand = expansion(mult = c(0, .75)), breaks = seq(0, 100, 25)) +
  labs(x = "Infants (first-week sample)", y = NULL, title = "Cohorts and their roles") +
  TH + theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
    axis.text.y = element_text(colour = INK, size = 5.9, lineheight = .9),
    strip.text.y = element_text(angle = 0, hjust = 0, size = 6, face = "bold"),
    panel.spacing.y = unit(5, "pt"))
fig <- (fA | fB) + plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "A")
ggsave(file.path(OUT, "NewFigure1.pdf"), fig, width = 7.2, height = 5.0, device = cairo_pdf)
ggsave(file.path(OUT, "NewFigure1.png"), fig, width = 7.2, height = 5.0, dpi = 300, bg = "white")
