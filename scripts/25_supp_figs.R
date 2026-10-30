# Supplementary figures S1-S7 for the v7 manuscript
src_until <- function(file, from_pat = NULL, until_pat = "^fig <-") {
  l <- readLines(file)
  i0 <- if (is.null(from_pat)) 1 else grep(from_pat, l)[1]
  i1 <- grep(until_pat, l)[1] - 1
  eval(parse(text = l[i0:i1]), envir = globalenv())
}
src_until("./fig2.R")
s1 <- fC
load(file.path(OUT, "fig3.RData"))
SC4 <- c(P4, "Combined score")
src_until("./11_fig3.R", from_pat = "^FONT <-")
s2B <- fB
s2C <- fC
suppressPackageStartupMessages(library(nlme))
source("./figtheme.R")
SD <- file.path(FIG, "supp")
dir.create(SD, showWarnings = FALSE)

sv <- function(p, n, w, h) {
  ggsave(file.path(SD, paste0(n, ".pdf")), p, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(SD, paste0(n, ".png")), p, width = w, height = h, dpi = 250, bg = "white")
}
# S1: core genes per program
sv(s1 + labs(title = "Top 12 core genes of each program"), "FigureS1", 7.4, 3.4)
# S2: GA-matched comparison and composition correlations
sv((s2B | s2C) + plot_layout(widths = c(1, 1.15)) + plot_annotation(tag_levels = "A"), "FigureS2", 7.4,
  3.4)
# S3: programs over the first month (previous Figure 4 build)
file.copy(file.path(FIG, "Figure4.pdf"), file.path(SD, "FigureS3.pdf"), overwrite = TRUE)
file.copy(file.path(FIG, "Figure4.png"), file.path(SD, "FigureS3.png"), overwrite = TRUE)
# S4: signature size
load(file.path(OUT, "fig5.RData"))
ga <- as.numeric(GA_AUC[1:3])
COLT <- c("in-sample" = "#9a9a92", "held-out" = C_BPD)
s4 <- ggplot(TOPN, aes(N, auc, colour = type, fill = type)) + annotate("rect", xmin = 2.4, xmax = 1250,
  ymin = ga[2], ymax = ga[3], fill = C_GA, alpha = .12) +
  geom_hline(yintercept = ga[1], colour = C_GA, linetype = "22",
    linewidth = .5) + geom_hline(yintercept = .5, colour = LIGHT, linewidth = .3) +
  annotate("text", x = 2.6, y = ga[1] + .022, label = sprintf("GA alone %.2f (3 cohorts)", ga[1]),
    colour = C_GA, hjust = 0, size = 2.1) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = .14,
    colour = NA) + geom_line(linewidth = .6) + geom_point(size = 1.3) +
  scale_x_log10(breaks = c(3, 10, 50, 250, 1000)) + scale_y_continuous(limits = c(.45,
    1)) + scale_colour_manual(values = COLT, name = NULL) + scale_fill_manual(values = COLT,
      name = NULL) +
  labs(x = "Number of top genes (log scale)", y = "Pooled AUC (95% CI)",
    title = "Signature size and held-out discrimination",
    subtitle = "Genes ranked in training cohorts only") +
  TH + theme(legend.position = "bottom")
sv(s4, "FigureS4", 3.8, 3.4)
# S5: neutrophil vs T-cell expression
load(file.path(OUT, "v7_partB.RData"))
s5 <- ggplot(D4, aes(T, N, colour = BPD, size = ga)) + geom_point(alpha = .7) + facet_wrap(~ cohort,
  nrow = 1) + scale_size(range = c(2.2, .5), name = "GA (wk)") +
  scale_colour_manual(values = c("No BPD" = C_NO, BPD = C_BPD),
    name = NULL) + labs(x = "T-cell marker expression (z)", y = "Neutrophil marker expression (z)",
      title = "Neutrophil and T-cell expression by cohort",
      subtitle = "Point size, gestational age (smaller = more mature); GSE125873 has no GA") + TH
sv(s5, "FigureS5", 7.4, 2.6)
# S6: protocol v7 pre-specified scores, week 1 per cohort (GA-adjusted) and longitudinal
ga_of <- function(k) PH[[k]]$ga_weeks
NAIVE <- c("CCR7", "TCF7", "LEF1", "BCL11B", "CD27", "CD28", "IL7R", "BACH2", "SATB1")
PAN_T <- c("CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK")
RTE <- c("PTK7", "CR2", "PECAM1", "SOX4", "TOX", "AIF1")
NEU <- c("CSF3R", "CXCR2", "FPR1", "FCGR3B", "MME", "CMTM2")

mz <- function(g, x) {
  g <- intersect(g, rownames(x))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z)))
}

scs <- function(x) {
  n <- mz(NAIVE, x)
  p <- mz(PAN_T, x)
  list(`Naive T` = n, `Pan-T` = p, `Naive within T` = as.numeric(scale(resid(lm(n ~ p)))),
    `Thymic emigrant (literature)` = mz(RTE, x),
    Neutrophil = mz(NEU, x), SPON1 = as.numeric(scale(x["SPON1", ])))
}
W <- bind_rows(lapply(GA_COH, function(k) {
  S <- scs(E[[k]])
  y <- lab(k)
  g <- ga_of(k)
  bind_rows(lapply(names(S), function(v) {
    f <- coef(summary(lm(S[[v]] ~ y + g)))["y", ]
    data.frame(score = v, row = paste0(k, ", week 1"), b = f[1], se = f[2])
  }))
}))
L <- readRDS(file.path(DATA, "longitudinal_GSE32472_GSE220135.rds"))
for (cc in c("GSE32472", "GSE220135")) {
  a <- L$annotation |> filter(cohort == cc, !is.na(t), !is.na(ga), !is.na(severity),
    id %in% colnames(L$expression[[cc]])) |> mutate(y = as.integer(severity != "No BPD"))
  x <- L$expression[[cc]][, a$id]
  x <- x[apply(x, 1, sd) > 0, ]
  S <- scs(x)
  for (v in names(S)) {
    a$s <- S[[v]]
    m <- summary(lme(s ~ y + factor(t) + ga, random = ~ 1 | infant, data = a,
      method = "REML"))$tTable["y", ]
    W <- rbind(W, data.frame(score = v, row = paste0(cc, ", all time points"), b = m[1], se = m[2]))
  }
}
W <- W |> mutate(score = factor(score, levels = names(scs(E$GSE32472))),
  row = factor(row, levels = rev(unique(row))))
s6 <- ggplot(W, aes(b, row)) + geom_vline(xintercept = 0, colour = LIGHT, linewidth = .3) +
  geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0, linewidth = .4,
    colour = MUTED) + geom_point(colour = C_BPD, size = 1.5) +
  facet_wrap(~ score, nrow = 2) + labs(x = "BPD effect after GA (SD, 95% CI); any BPD", y = NULL,
    title = "Pre-specified scores (protocol v7)",
    subtitle = "Week 1: linear models in each GA cohort; all time points: mixed models with infant random intercept") + TH + theme(axis.text.y = element_text(colour = INK))
sv(s6, "FigureS6", 7.4, 4.6)
# S7: exploratory glucocorticoid-response score and time-course adjustment.
load(file.path(OUT, "timecourse.RData"))
g7 <- LT |> group_by(cohort, t, grp = ifelse(y == 1, "BPD", "No BPD")) |> summarise(m = mean(gc),
  se = sd(gc) / sqrt(n())) |>
  mutate(time = factor(case_when(t == 1 ~ "Day 14", t == 2 ~ "Day 28", cohort == "GSE32472" ~ "Day 4-7",
    TRUE ~ "Cord"), levels = c("Cord", "Day 4-7", "Day 14", "Day 28")))
s7a <- ggplot(g7,
  aes(time, m, colour = grp,
    group = grp)) + geom_line(linewidth = .6) + geom_pointrange(aes(ymin = m - se, ymax = m + se),
      size = .25) + facet_wrap(~ cohort, scales = "free_x") +
  scale_colour_manual(values = c("No BPD" = C_NO, BPD = C_BPD), name = NULL) + labs(x = NULL,
    y = "Glucocorticoid-response score (mean, SE)",
    title = "Exploratory glucocorticoid-response score over time") + TH + theme(legend.position = "top")
s7b <- ggplot(ADJ |> filter(score == "comb",
  !grepl("Week 1", model)) |> mutate(model = factor(model,
    levels = rev(c("GA", "GA + GC response", "GA + composition", "GA + GC + composition")))),
  aes(b, model)) + geom_vline(xintercept = 0, colour = LIGHT) + geom_point(colour = C_BPD,
    size = 1.8) + geom_text(aes(label = sprintf("p %.1g", p)), hjust = -.3, size = 2) +
  facet_wrap(~ cohort) + scale_x_continuous(limits = c(-.1,
    1.3)) + labs(x = "BPD effect on combined score, all time points (SD)", y = NULL,
      title = "Adjustment for the proxy and cell composition",
      subtitle = "Held-out cores; mixed models") + TH + theme(axis.text.y = element_text(colour = INK))
sv(s7a / s7b + plot_annotation(tag_levels = "A"), "FigureS7", 7.4, 5.2)
cat("done\n")
