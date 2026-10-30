# Compare program validation, GA matching and cell composition.
# Primary validation and sensitivity analyses
source("./00_setup.R")
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
load(file.path(OUT, "programs4.RData"))
load(file.path(OUT, "heldout4_relaxed.RData"))
FOLD_A <- FOLD
FIG <- "./figs"
DISC <- setNames(lapply(P4, function(pg) PR4$core[[which(PR4$program == pg)]] |> transmute(gene, es, p)),
  P4)
FOLD_B <- lapply(setNames(BLOOD, BLOOD), function(h) {
  tr <- setdiff(BLOOD, h)
  Mt <- meta_all(S_blood |> filter(gse %in% tr), kmin = length(tr) - 1)
  setNames(lapply(P4,
    function(pg) Mt |> filter(gene %in% PR4$genes[[which(PR4$program == pg)]], p < .05,
      ifelse(es > 0, n_pos, n_neg) >= length(tr) - 1) |>
    arrange(p) |> head(30) |> transmute(gene, es, p)), P4)
})

psc <- function(core, x) {
  if (is.null(core) || !nrow(core)) return(rep(NA_real_, ncol(x)))
  g <- core$gene[core$gene %in% rownames(x)]
  if (!length(g)) return(rep(NA_real_, ncol(x)))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  colMeans(z * sign(core$es[match(g, core$gene)]))
}
SC4 <- c(P4, "Combined score")

score <- function(nm, k, type) {
  cs <- switch(type, discovery = DISC, `held-out` = FOLD_A[[k]], `held-out B` = FOLD_B[[k]])
  if (nm == "Combined score") {
    m <- sapply(P4, function(pg) {
      s <- psc(cs[[pg]], E[[k]])
      if (all(is.na(s))) s else as.numeric(scale(s))
    })
    r <- rowMeans(m, na.rm = TRUE)
    r[is.nan(r)] <- NA
    return(r)
  }
  psc(cs[[nm]], E[[k]])
}
TYPES <- c("discovery", "held-out", "held-out B")

# ---- A: AUC ----
AUCT <- bind_rows(lapply(SC4, function(nm) bind_rows(lapply(TYPES, function(tp) {
  r <- t(sapply(BLOOD, function(k) roc_stats(score(nm, k, tp), lab(k))))
  ok <- !is.na(r[, "auc"])
  pa <- pool_auc(r[ok, "auc"], r[ok, "var"])
  data.frame(score = nm, type = tp, auc = pa[1], lo = pa[2], hi = pa[3], k = sum(ok),
    cohorts = paste(round(r[, "auc"], 2), collapse = "/"))
}))))

# ---- B: GA-matched strata, and Model B ----
strat <- function(s, y, ga) {
  w <- floor(ga)
  out <- NULL
  for (st in unique(w)) {
    i <- w == st
    if (sum(y[i] == 1) >= 1 && sum(y[i] == 0) >= 1 && sum(i) >= 3) {
      n1 <- sum(i & y == 1)
      n0 <- sum(i & y == 0)
      out <- rbind(out, c(mean(s[i & y == 1]) - mean(s[i & y == 0]), n1 * n0 / (n1 + n0), n1, n0))
    }
  }
  c(d = sum(out[, 1] * out[, 2]) / sum(out[, 2]), v = sd(s)^2 / sum(out[, 2]), n1 = sum(out[, 3]),
    n0 = sum(out[, 4]))
}

# ---- cell scores (markers as in analysis 3, plus eosinophils) ----
MK <- list(`T cells` = c("CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD6", "LCK", "TRAT1", "CD247"),
  `CD4 T` = c("CD4", "CD40LG", "IL7R", "MAL"),
  `CD8 T` = c("CD8A", "CD8B"),
  `B cells` = c("MS4A1", "CD79A", "CD79B", "CD19", "CD22", "BANK1", "FCRLA", "PAX5", "BLK", "TCL1A"),
  `NK cells` = c("KLRF1", "NCR1", "SH2D1B", "KLRD1", "PRF1"),
  Monocytes = c("CD14", "FCN1", "CSF1R", "VCAN", "LYZ", "CD68"),
  Neutrophils = c("FCGR3B", "CXCR2", "CSF3R", "FPR1", "FPR2", "CMTM2", "MME", "CXCR1", "PROK2"),
  Eosinophils = c("CLC", "IL5RA", "CCR3", "SIGLEC8", "EPX", "PRG3"),
  Erythroid = c("ALAS2", "GYPA", "CA1", "SLC4A1", "AHSP", "HBD", "EPB42", "SPTA1"))
CS <- lapply(setNames(BLOOD, BLOOD), function(k) sapply(MK, function(g) {
  g <- intersect(g, rownames(E[[k]]))
  if (!length(g)) return(rep(NA, ncol(E[[k]])))
  z <- t(scale(t(E[[k]][g, , drop = FALSE])))
  z[is.na(z)] <- 0
  colMeans(z)
}))
COMP <- list(Lymphoid = c("T cells", "B cells", "NK cells"), Myeloid = c("Monocytes", "Neutrophils"),
  Erythroid = "Erythroid")

comp_mat <- function(k) sapply(COMP,
  function(cc) rowMeans(scale(CS[[k]][, cc, drop = FALSE]), na.rm = TRUE))
DEVT <- bind_rows(lapply(SC4, function(nm) bind_rows(lapply(c("discovery", "held-out"), function(tp) {
  r <- t(sapply(GA_COH, function(k) {
    s0 <- score(nm, k, tp)
    if (all(is.na(s0))) return(rep(NA, 10))
    s <- as.numeric(scale(s0))
    y <- lab(k)
    ga <- PH[[k]]$ga_weeks
    X <- comp_mat(k)
    fb <- lm(s ~ y + ga)
    fc <- lm(s ~ y + ga + X)
    c(strat(s, y, ga), coef(fb)[["y"]], vcov(fb)["y", "y"], coef(fc)[["y"]], vcov(fc)["y", "y"],
      summary(lm(s ~ X))$r.squared, 1)
  }))
  ok <- !is.na(r[, 1])
  S <- re_pool(r[ok, 1], r[ok, 2])
  B <- re_pool(r[ok, 5], r[ok, 6])
  C <- re_pool(r[ok, 7], r[ok, 8])
  data.frame(score = nm, type = tp, k = sum(ok), strata = S[["es"]], s_lo = S[["lo"]], s_hi = S[["hi"]],
    s_p = S[["p"]], n1 = sum(r[ok, 3]), n0 = sum(r[ok, 4]),
    B = B[["es"]], B_lo = B[["lo"]], B_hi = B[["hi"]], B_p = B[["p"]], C = C[["es"]], C_lo = C[["lo"]],
    C_hi = C[["hi"]], C_p = C[["p"]], r2 = mean(r[ok, 9]))
}))))
CORR <- bind_rows(lapply(SC4, function(nm) bind_rows(lapply(names(MK), function(cn) {
  zz <- sapply(BLOOD,
    function(k) atanh(cor(score(nm, k, "discovery"), CS[[k]][, cn], use = "complete.obs")))
  w <- sapply(BLOOD, function(k) ncol(E[[k]]) - 3)
  ok <- is.finite(zz)
  data.frame(score = nm, cell = cn, r = tanh(sum(zz[ok] * w[ok]) / sum(w[ok])))
}))))
print(AUCT |> mutate(across(where(is.numeric), ~ round(.x, 2))))
print(DEVT |> mutate(across(where(is.numeric), ~ signif(.x, 2))), width = 250)
print(CORR |> tidyr::pivot_wider(names_from = cell,
  values_from = r) |> mutate(across(where(is.numeric), ~ round(.x, 2))), width = 250)
save(AUCT, DEVT, CORR, FOLD_B, CS, MK, file = file.path(OUT, "fig3.RData"))

# ---- plot ----
FONT <- "Liberation Sans"
INK <- "#222222"
MUTED <- "#5f5f5f"
LIGHT <- "#c9c7c0"
C_BPD <- "#1f3a6e"
update_geom_defaults("text", list(family = FONT, colour = INK, size = 2.1))
TH <- theme_classic(base_size = 7, base_family = FONT) +
  theme(text = element_text(colour = INK), axis.text = element_text(colour = MUTED, size = 6),
    axis.title = element_text(size = 6.6),
    axis.line = element_line(colour = MUTED, linewidth = .3),
    axis.ticks = element_line(colour = MUTED, linewidth = .3),
    plot.title = element_text(size = 7.4, face = "bold", hjust = 0, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 6.1, colour = MUTED, margin = margin(b = 3)),
    legend.text = element_text(size = 6.1), legend.key.size = unit(7, "pt"),
    strip.background = element_blank(), strip.text = element_text(size = 6.4, face = "bold"),
    plot.margin = margin(4, 6, 4, 4), plot.tag = element_text(size = 10, face = "bold", family = FONT))
PCOL <- c("RNA processing and translation" = "#6a3d9a", "T cell activation" = "#1d9fb0",
  "Cell cycle" = "#8a8a82", "Neutrophil antimicrobial" = "#e31a1c", "Combined score" = "#a6611a")
SH <- c("RNA processing and translation" = "RNA processing & translation",
  "T cell activation" = "T cell activation", "Cell cycle" = "Cell cycle",
  "Neutrophil antimicrobial" = "Neutrophil antimicrobial", "Combined score" = "Combined score")
lev <- rev(SC4)
dA <- AUCT |> mutate(score = factor(score, levels = lev))
dA1 <- dA |> filter(type == "discovery")
dA2 <- dA |> filter(type == "held-out")
dA3 <- dA |> filter(type == "held-out B")
fA <- ggplot(mapping = aes(y = score)) +
  geom_vline(xintercept = .5, linetype = "22", colour = MUTED, linewidth = .3) +
  geom_segment(data = dA1 |> left_join(dA2 |> select(score, ho = auc), by = "score"),
    aes(x = ho, xend = auc, yend = score), colour = LIGHT, linewidth = 1.2) +
  geom_point(data = dA3, aes(x = auc), shape = 23, size = 1.5, colour = "#9a9a9a", fill = "white",
    position = position_nudge(y = -.25)) +
  geom_errorbarh(data = dA2, aes(xmin = lo, xmax = hi, colour = score), height = 0, linewidth = .45) +
  geom_point(data = dA1, aes(x = auc), shape = 21, fill = "white", colour = "#7a7a7a", size = 1.7,
    stroke = .6) +
  geom_point(data = dA2, aes(x = auc, colour = score), size = 1.9) +
  geom_text(data = dA2, aes(x = 1.0, label = paste0(k, "/4")), hjust = 1, colour = MUTED, size = 1.9) +
  annotate("text", x = 1.0, y = 5.6, label = "Cohorts", hjust = 1, colour = MUTED, size = 1.85) +
  scale_colour_manual(values = PCOL,
    guide = "none") + scale_y_discrete(labels = SH) + scale_x_continuous(limits = c(.3, 1),
      breaks = seq(.3, 1, .1)) +
  coord_cartesian(clip = "off") +
  labs(x = "Pooled AUC", y = NULL, title = "Held-out validation",
    subtitle = "Open, in-sample; filled, held-out (95% CI);\ndiamond, sensitivity with fixed program gene sets") +
  TH + theme(axis.text.y = element_text(colour = INK), plot.margin = margin(10, 6, 4, 4))
dB <- DEVT |> mutate(score = factor(score, levels = lev),
  type = factor(type, levels = c("discovery", "held-out")))
fB <- ggplot(dB, aes(strata, score, colour = type)) + geom_vline(xintercept = 0, colour = LIGHT,
  linewidth = .3) +
  geom_errorbarh(aes(xmin = s_lo, xmax = s_hi), height = 0, linewidth = .45,
    position = position_dodge(.55)) + geom_point(size = 1.6, position = position_dodge(.55)) +
  scale_colour_manual(values = c(discovery = "#9a9a9a", `held-out` = C_BPD),
    labels = c("In-sample genes", "Held-out genes"), name = NULL, guide = guide_legend(reverse = TRUE)) +
  scale_y_discrete(labels = SH) +
  labs(x = "Difference within GA week (SD, 95% CI)", y = NULL, title = "GA-matched comparison",
    subtitle = sprintf("Up to %d BPD vs %d no BPD in shared GA weeks", dB$n1[1], dB$n0[1])) +
  TH + theme(legend.position = "top", legend.justification = "left", legend.direction = "vertical",
    legend.margin = margin(0, 0, 0, 0), axis.text.y = element_text(colour = INK))
dC <- CORR |> mutate(score = factor(score, levels = lev), cell = factor(cell, levels = names(MK)))
fC <- ggplot(dC, aes(cell, score, fill = r)) + geom_tile(colour = "white", linewidth = .6) +
  geom_text(aes(label = sprintf("%.1f", r), colour = abs(r) > .55), size = 1.9) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = INK), guide = "none") +
  scale_fill_gradientn(colours = c("#1f3a6e", "#9fb3d6", "#f5f5f5", "#e8a493", "#a63a2a"),
    limits = c(-1, 1), name = "Pooled r",
    guide = guide_colourbar(barheight = unit(36, "pt"), barwidth = unit(5, "pt"))) +
  scale_y_discrete(labels = SH) +
  labs(x = NULL, y = NULL, title = "Programs track immune-cell composition",
    subtitle = "Correlation with marker-based cell scores, 4 cohorts") +
  TH + theme(axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text = element_text(colour = INK), axis.text.x = element_text(angle = 35, hjust = 1))
dD <- bind_rows(DEVT |> transmute(score, type, est = B, lo = B_lo, hi = B_hi, model = "GA"),
  DEVT |> transmute(score, type, est = C, lo = C_lo, hi = C_hi, model = "GA + cell composition")) |>
  mutate(score = factor(score, levels = lev),
    type = factor(type, levels = c("discovery", "held-out"),
      labels = c("In-sample genes", "Held-out genes")))
fD <- ggplot(dD, aes(est, score, shape = model)) + geom_vline(xintercept = 0, colour = LIGHT,
  linewidth = .3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .4, colour = MUTED,
    position = position_dodge(.6)) +
  geom_point(aes(colour = score), size = 1.6, fill = "white", stroke = .55,
    position = position_dodge(.6)) +
  facet_wrap(~ type, nrow = 1) + scale_shape_manual(values = c(19, 21),
    name = NULL) + scale_colour_manual(values = PCOL, guide = "none") +
  scale_y_discrete(labels = SH) +
  labs(x = "BPD effect on score (SD, 95% CI)", y = NULL,
    title = "Adjustment for GA and cell composition",
    subtitle = "3 cohorts with GA; lymphoid, myeloid, erythroid axes") +
  TH + theme(legend.position = "top", legend.justification = "left", legend.margin = margin(0, 0, 0, 0),
    axis.text.y = element_text(colour = INK), panel.spacing.x = unit(8, "pt"))
fig <- wrap_plots(A = fA, B = fB, C = fC, D = fD, design = "AB\nCD",
  widths = c(1, 1.1)) + plot_annotation(tag_levels = "A")
ggsave(file.path(FIG, "Figure3.pdf"), fig, width = 7.2, height = 6.4, device = cairo_pdf)
ggsave(file.path(FIG, "Figure3.png"), fig, width = 7.2, height = 6.4, dpi = 300, bg = "white")
