# Figure 4: the four programs over the first month (GSE32472 days 4-7/14/28; GSE220135 cord/14/28)
#  Scores from the final core genes (week-1 samples contributed to selection; day 14/28 did not).
#  Expression z-scored within cohort across all time points
source("./00_setup.R")
suppressPackageStartupMessages({
  library(nlme)
  library(ggplot2)
  library(patchwork)
})
load(file.path(OUT, "programs4.RData"))
load(file.path(OUT, "heldout4_relaxed.RData"))
L <- readRDS(file.path(DATA, "longitudinal_GSE32472_GSE220135.rds"))
DISC <- setNames(lapply(P4, function(pg) PR4$core[[which(PR4$program == pg)]] |> transmute(gene, es, p)),
  P4)

psc <- function(core, x) {
  g <- core$gene[core$gene %in% rownames(x)]
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  colMeans(z * sign(core$es[match(g, core$gene)]))
}
SC <- c("Combined score", P4)
ann <- L$annotation |> filter(!is.na(t), !is.na(severity),
  id %in% c(colnames(L$expression$GSE32472), colnames(L$expression$GSE220135)))
LT <- bind_rows(lapply(c("GSE32472", "GSE220135"), function(cc) {
  a <- ann |> filter(cohort == cc)
  x <- L$expression[[cc]][, a$id]
  x <- x[apply(x, 1, sd) > 0, ]
  S <- sapply(P4, function(pg) as.numeric(scale(psc(DISC[[pg]], x))))
  S <- cbind(`Combined score` = rowMeans(S), S)
  hoT <- as.numeric(scale(psc(FOLD[[cc]][["T cell activation"]],
    x)))  # held-out T cell program (cores derived without this cohort)
  data.frame(a, S, `T cell activation (held-out)` = hoT, check.names = FALSE)
}))
LT <- LT |> mutate(grp = case_when(severity == "No BPD" ~ "No BPD", severity == "Mild" ~ "Mild",
  TRUE ~ ifelse(cohort == "GSE220135", "BPD", "Moderate/severe")),
  y = as.integer(grp %in% c("Moderate/severe", "BPD")), tf = factor(t))
TL <- list(GSE32472 = c("D4-7", "D14", "D28"), GSE220135 = c("Cord", "D14", "D28"))
# effect sizes per time point (moderate/severe or BPD vs no BPD)
EFF <- bind_rows(lapply(c("GSE32472", "GSE220135"),
  function(cc) bind_rows(lapply(c(SC, "T cell activation (held-out)"),
    function(s) bind_rows(lapply(0:2, function(tt) {
      d <- LT |> filter(cohort == cc, t == tt, grp != "Mild")
      h <- hedges(d[[s]], d$y)
      data.frame(cohort = cc, score = s, t = tt, time = TL[[cc]][tt + 1], g = h[["g"]],
        p = 2 * pnorm(-abs(h[["g"]] / sqrt(h[["v"]]))))
}))))))
MM <- bind_rows(lapply(c("GSE32472", "GSE220135"),
  function(cc) bind_rows(lapply(c(SC, "T cell activation (held-out)"), function(s) {
    d <- LT |> filter(cohort == cc, grp != "Mild", !is.na(ga))
    d$s <- d[[s]]
    fm <- function(f) lme(f, random = ~ 1 | infant, data = d, method = "ML", na.action = na.omit)
    lr <- function(a, b) anova(a, b)$`p-value`[2]
    m0 <- fm(s ~ tf)
    m1 <- fm(s ~ y + tf)
    m2 <- fm(s ~ y * tf)
    g0 <- fm(s ~ tf + ga)
    g1 <- fm(s ~ y + tf + ga)
    data.frame(cohort = cc, score = s, p_group = lr(m0, m1), p_int = lr(m1, m2), p_group_ga = lr(g0, g1))
}))))
print(EFF |> filter(score %in% c("Combined score", "T cell activation",
  "T cell activation (held-out)")) |> mutate(across(where(is.numeric), ~ round(.x, 3))))
print(MM |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
# trajectory classes (combined score above cohort median at week 1 and day 28)
TR <- LT |> group_by(cohort,
  t) |> mutate(high = as.integer(`Combined score` > median(`Combined score`))) |> ungroup() |> filter(t %in% c(0, 2)) |>
  select(cohort, infant, grp, t, high) |> tidyr::pivot_wider(names_from = t, values_from = high,
    names_prefix = "t") |> filter(!is.na(t0), !is.na(t2)) |>
  mutate(cls = factor(case_when(t0 == 1 & t2 == 1 ~ "Persistent", t0 == 0 & t2 == 1 ~ "Emerging",
    t0 == 1 & t2 == 0 ~ "Resolving", TRUE ~ "Never high"),
    levels = c("Persistent", "Emerging", "Resolving", "Never high")))
FI <- TR |> filter(grp != "Mild") |> group_by(cohort) |> summarise(pers_case = mean(cls[grp != "No BPD"] == "Persistent"), pers_none = mean(cls[grp == "No BPD"] == "Persistent"),
  n_case = sum(grp != "No BPD"), n_none = sum(grp == "No BPD"),
  p = fisher.test(table(grp != "No BPD", cls == "Persistent"))$p.value)
print(FI)
save(LT, EFF, MM, TR, FI, file = file.path(OUT, "fig4.RData"))

# ---- plot ----
FIG <- "./figs"
FONT <- "Liberation Sans"
INK <- "#222222"
MUTED <- "#5f5f5f"
LIGHT <- "#c9c7c0"
C_BPD <- "#1f3a6e"
C_MILD <- "#8aa2cc"
C_NO <- "#c8553d"
update_geom_defaults("text", list(family = FONT, colour = INK, size = 2.1))
TH <- theme_classic(base_size = 7, base_family = FONT) +
  theme(text = element_text(colour = INK), axis.text = element_text(colour = MUTED, size = 6),
    axis.title = element_text(size = 6.6),
    axis.line = element_line(colour = MUTED, linewidth = .3),
    axis.ticks = element_line(colour = MUTED, linewidth = .3),
    plot.title = element_text(size = 7.4, face = "bold", hjust = 0, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 6.1, colour = MUTED, margin = margin(b = 3)),
    plot.caption = element_text(size = 5.8, colour = MUTED, hjust = 0),
    legend.text = element_text(size = 6.1), legend.key.size = unit(7, "pt"),
    strip.background = element_blank(), strip.text = element_text(size = 6.4, face = "bold"),
    plot.margin = margin(4, 6, 4, 4),
    plot.tag = element_text(size = 10, face = "bold", family = FONT))

star <- function(p) ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", "")))
SH <- c("Combined score" = "Combined score",
  "RNA processing and translation" = "RNA processing & translation",
  "T cell activation" = "T cell activation",
  "Cell cycle" = "Cell cycle", "Neutrophil antimicrobial" = "Neutrophil antimicrobial",
  "T cell activation (held-out)" = "T cell activation, held-out")
dA <- EFF |> mutate(col = factor(paste(cohort, time),
  levels = c(paste("GSE32472", TL$GSE32472), paste("GSE220135", TL$GSE220135))),
  score = factor(score, levels = rev(c(SC, "T cell activation (held-out)"))),
  l = paste0(sprintf("%.2f", g), star(p)))
fA <- ggplot(dA, aes(col, score, fill = g)) + geom_tile(colour = "white", linewidth = .8) +
  geom_text(aes(label = l, colour = abs(g) > 1.1,
    fontface = ifelse(score == "Combined score", "bold", "plain")), size = 1.95) +
  geom_vline(xintercept = 3.5, colour = "white", linewidth = 2.2) +
  annotate("text", x = c(2, 5), y = length(levels(dA$score)) + .8,
    label = c("GSE32472: moderate/severe vs none", "GSE220135: BPD vs none"), size = 2,
    fontface = "bold") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colours = c("#a63a2a", "#e8a493", "#f5f5f5", "#9fb3d6", "#1f3a6e"),
    limits = c(-2.2, 2.2), oob = scales::squish, name = "Hedges' g",
    guide = guide_colourbar(barheight = unit(40, "pt"), barwidth = unit(5, "pt"))) +
  scale_x_discrete(labels = function(x) sub("^GSE\\d+ ", "",
    x)) + scale_y_discrete(labels = SH) + coord_cartesian(clip = "off",
      ylim = c(.5, length(levels(dA$score)) + .5)) +
  labs(x = NULL, y = NULL, title = "Program differences over the first month",
    subtitle = "Higher = more BPD-like; day 14 and 28 samples were not used to select genes") +
  TH + theme(axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.y = element_text(colour = INK), plot.margin = margin(14, 4, 4, 4))

traj <- function(cc) {
  d <- LT |> group_by(grp, t) |> filter(cohort == cc) |> summarise(m = mean(`Combined score`),
    se = sd(`Combined score`) / sqrt(n()), n = n_distinct(infant), .groups = "drop")
  lv <- intersect(c("No BPD", "Mild", "Moderate/severe", "BPD"), unique(d$grp))
  d$grp <- factor(d$grp, levels = lv)
  nn <- LT |> filter(cohort == cc) |> distinct(infant, grp) |> count(grp)
  mm <- MM |> filter(cohort == cc, score == "Combined score")
  fmt <- function(p) ifelse(p < .001, formatC(p, format = "e", digits = 1), sprintf("%.2g", p))
  ggplot(d |> mutate(x = t + (as.integer(grp) - 2) * .06),
    aes(x, m, colour = grp)) + geom_hline(yintercept = 0, colour = "#e4e3df", linewidth = .3) +
    geom_errorbar(aes(ymin = m - se, ymax = m + se), width = .06,
      linewidth = .4) + geom_line(linewidth = .7) + geom_point(size = 1.4) +
    scale_colour_manual(values = c("No BPD" = C_NO, "Mild" = C_MILD, "Moderate/severe" = C_BPD,
      BPD = C_BPD), breaks = lv,
      labels = paste0(lv, " (n = ", nn$n[match(lv, nn$grp)], ")"), name = NULL) +
    scale_x_continuous(breaks = 0:2, labels = TL[[cc]]) +
    labs(x = NULL, y = "Combined score (mean, SE)",
      title = paste0(cc, if (cc == "GSE32472") " (microarray)" else " (RNA-seq)"),
      caption = sprintf("Group P = %s; GA-adjusted P = %s\nGroup-by-time P = %s", fmt(mm$p_group),
        fmt(mm$p_group_ga), fmt(mm$p_int))) +
    TH + theme(legend.position = "top", legend.justification = "left", legend.direction = "vertical",
      legend.margin = margin(0, 0, 0, 0))
}
TCOL <- c(Persistent = "#1f3a6e", Emerging = "#d4a017", Resolving = "#e8a493", `Never high` = "#e6e6e6")
dE <- TR |> count(cohort, grp, cls, .drop = FALSE) |> group_by(cohort, grp) |> mutate(N = sum(n),
  f = n / N) |> ungroup() |> filter(N > 0) |>
  mutate(cohort = factor(cohort, levels = c("GSE32472", "GSE220135")),
    grp = factor(grp, levels = c("No BPD", "Mild", "Moderate/severe", "BPD"))) |> arrange(cohort, grp) |>
  mutate(row = paste0(grp, " (n = ", N, ")"))
dE$row <- factor(dE$row, levels = rev(unique(dE$row)))
fE <- ggplot(dE, aes(f, row, fill = cls)) + geom_col(width = .7, colour = "white", linewidth = .4,
  position = position_stack(reverse = TRUE)) +
  geom_text(aes(label = ifelse(f >= .1, paste0(round(100 * f), "%"), ""), colour = cls == "Persistent",
    group = cls), position = position_stack(vjust = .5, reverse = TRUE), size = 1.95) +
  facet_grid(cohort ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK),
    guide = "none") + scale_fill_manual(values = TCOL, name = NULL) +
  scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "Infants sampled in week 1 and at day 28", y = NULL, title = "Trajectory class by BPD status",
    caption = paste0("Persistent class, BPD vs no BPD (Fisher):\n",
      paste0(FI$cohort, " P = ", formatC(FI$p, format = "e", digits = 1), collapse = "; "))) +
  TH + theme(legend.position = "top", legend.justification = "left", strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, size = 6, face = "bold"),
    axis.text.y = element_text(colour = INK), legend.margin = margin(0, 0, 0, 0))
TDEF <- c(Persistent = "high to high", Emerging = "low to high", Resolving = "high to low",
  `Never high` = "low to low")
keyd <- tibble(cls = factor(rep(names(TCOL), each = 2), levels = names(TCOL)), x = rep(0:1, 4),
  y = c(1, 1, 0, 1, 1, 0, 0, 0))
fK <- ggplot(keyd, aes(x, y, colour = cls)) + geom_hline(yintercept = .5, linetype = "22",
  colour = MUTED, linewidth = .25) +
  geom_line(linewidth = .8) + geom_point(size = 1.3) +
  facet_wrap(~ cls, ncol = 1, strip.position = "left",
    labeller = as_labeller(setNames(paste0(names(TDEF), "\n", TDEF), names(TDEF)))) +
  scale_colour_manual(values = replace(TCOL, "Never high", "#9a9a9a"), guide = "none") +
  scale_x_continuous(breaks = 0:1, labels = c("Week 1", "D28"), expand = expansion(add = .25)) +
  scale_y_continuous(breaks = 0:1, labels = c("low", "high"), limits = c(-.3, 1.3)) +
  labs(x = NULL, y = NULL, title = "Trajectory classes", subtitle = "High = above the\ncohort median") +
  TH + theme(strip.placement = "outside",
    strip.text.y.left = element_text(angle = 0, hjust = 1, size = 6, face = "plain"),
    axis.text = element_text(size = 5.4), panel.spacing = unit(2, "pt"))
fig <- free(fA) / (traj("GSE32472") | traj("GSE220135")) / ((fK | fE) + plot_layout(widths = c(1,
  3.2))) + plot_layout(heights = c(1.15, 1.05, .95)) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIG, "Figure4.pdf"), fig, width = 7.2, height = 8, device = cairo_pdf)
ggsave(file.path(FIG, "Figure4.png"), fig, width = 7.2, height = 8, dpi = 300, bg = "white")
