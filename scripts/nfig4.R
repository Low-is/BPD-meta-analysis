# Figure 4: maturity-adjusted discovery and naive T-cell expression.
source("./00_setup.R")
suppressPackageStartupMessages(library(nlme))
source("./figtheme.R")
load(file.path(OUT, "newfigs.RData"))
load(file.path(OUT, "v7_partA.RData"))

ga_of <- function(k) PH[[k]]$ga_weeks
L <- readRDS(file.path(DATA, "longitudinal_GSE32472_GSE220135.rds"))

modscore <- function(M, x) {
  g <- M$gene[M$gene %in% rownames(x)]
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z * sign(M$t[match(g, M$gene)]))))
}

# ---- A: GA-adjusted effects in the two longitudinal cohorts ----
LABG <- c("CCR7", "TCF7", "LEF1", "BCL11B", "CD27", "CD28", "IL7R", "SATB1", "SPON1", "NOG", "CHI3L2",
  "GPA33", "MAN1C1", "CD6")
q <- with(SCAT,
  c(dd = sum(p32 < .05 & p35 < .05 & t32 < 0 & t35 < 0),
    uu = sum(p32 < .05 & p35 < .05 & t32 > 0 & t35 > 0),
    du = sum(p32 < .05 & p35 < .05 & t32 < 0 & t35 > 0),
    ud = sum(p32 < .05 & p35 < .05 & t32 > 0 & t35 < 0)))
dA <- SCAT |> mutate(cls = ifelse(gene %in% LABG, "Naive T-cell genes",
  "Other genes")) |> arrange(cls == "Naive T-cell genes")
lim <- 7.2
fA <- ggplot(dA, aes(t32, t35)) +
  geom_hline(yintercept = 0, colour = LIGHT, linewidth = .3) + geom_vline(xintercept = 0, colour = LIGHT,
    linewidth = .3) +
  geom_point(data = dA |> filter(cls == "Other genes"), colour = "#bdbdb5", size = .35, alpha = .5) +
  geom_point(data = dA |> filter(cls != "Other genes"), colour = C_NAIVE, size = 1.4) +
  annotate("text", x = -lim, y = lim, hjust = 0, vjust = 1, size = 2.0, colour = C_NAIVE,
    lineheight = .95,
    label = paste0("Naive T-cell genes (teal):\n",
      paste(strwrap(paste(LABG, collapse = ", "), 34), collapse = "\n"))) +
  annotate("text", x = -lim, y = -lim, label = paste0(q["dd"], " lower in both"), hjust = 0, vjust = 0,
    size = 2.1, fontface = "bold", colour = C_BPD) +
  annotate("text", x = lim, y = lim, label = paste0(q["uu"], " higher in both"), hjust = 1, vjust = 1,
    size = 2.1, colour = MUTED) +
  annotate("text", x = lim, y = -lim, label = paste0(q["ud"] + q["du"], " discordant"), hjust = 1,
    vjust = 0, size = 2.1, colour = MUTED) +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = "GSE32472: BPD effect after GA (t)", y = "GSE220135: BPD effect after GA (t)",
    title = "After GA, both cohorts agree on naive T-cell genes",
    subtitle = sprintf("%s genes; all time points, infant as blocking factor; counts at p < 0.05 in both", format(nrow(dA), big.mark = ","))) + TH

# ---- B: cross-cohort replication of the 50-gene module ----
RB <- NULL
for (src in names(out)) {
  M <- out[[src]]
  tgt <- setdiff(names(out), src)
  a <- L$annotation |> filter(cohort == tgt, !is.na(t), !is.na(ga), !is.na(severity),
    id %in% colnames(L$expression[[tgt]])) |> mutate(y = as.integer(severity != "No BPD"))
  a$s <- modscore(M, L$expression[[tgt]][, a$id])
  w <- a |> filter(t == 0)
  f <- coef(summary(lm(s ~ y + ga, w)))["y", ]
  m <- summary(lme(s ~ y + factor(t) + ga, random = ~ 1 | infant, data = a, method = "REML"))$tTable["y",
  ]
  RB <- rbind(RB,
    data.frame(src = src, tgt = tgt, when = c("Week 1", "All time points"), b = c(f[1], m[1]),
      se = c(f[2], m[2]), p = c(f[4], m[5])))
  for (k in c("GSE225881", "GSE125873")) {
    s <- modscore(M, E[[k]])
    y <- lab(k)
    if (k == "GSE225881") {
      f <- coef(summary(lm(s ~ y + ga_of(k))))[2, ]
      RB <- rbind(RB, data.frame(src = src, tgt = k, when = "Week 1", b = f[1], se = f[2], p = f[4]))
    }
    else {
      h <- hedges(s, y)
      RB <- rbind(RB,
        data.frame(src = src, tgt = k, when = "Week 1 (no GA)", b = h[["g"]], se = sqrt(h[["v"]]),
          p = 2 * pnorm(-abs(h[["g"]] / sqrt(h[["v"]])))))
    }
  }
}
RB <- RB |> mutate(row = paste0(sub("GSE", "", tgt), ", ", sub("ga\\)", "GA)", tolower(when))),
  src = factor(paste("Genes found in", src),
    levels = paste("Genes found in", c("GSE32472", "GSE220135"))),
  big = tgt %in% c("GSE32472", "GSE220135"), row = factor(row, levels = rev(unique(row))))
fB <- ggplot(RB, aes(b, row, colour = big)) + geom_vline(xintercept = 0, colour = LIGHT,
  linewidth = .3) +
  geom_errorbarh(aes(xmin = b - 1.96 * se, xmax = b + 1.96 * se), height = 0,
    linewidth = .45) + geom_point(size = 1.7) +
  geom_text(aes(x = 2.6, label = ifelse(p < .001, sprintf("p %.0e", p), sprintf("p %.2f", p))),
    hjust = 1, size = 1.95, colour = MUTED) +
  facet_wrap(~ src, ncol = 1,
    scales = "free_y") + scale_colour_manual(values = c(`TRUE` = C_BPD, `FALSE` = C_PUB),
      guide = "none") +
  scale_x_continuous(limits = c(-2, 2.65), breaks = -2:2) +
  labs(x = "BPD effect on 50-gene module, GA-adjusted (SD, 95% CI)", y = NULL,
    title = "Module from one cohort, tested in the others",
    subtitle = "Positive = more BPD-like; tested cohorts listed") +
  TH + theme(axis.text.y = element_text(colour = INK))

# ---- C/D: time course ----
TL <- list(GSE32472 = c("Day 4-7", "Day 14", "Day 28"), GSE220135 = c("Cord", "Day 14", "Day 28"))

tc <- function(cc) {
  d <- TC |> filter(cohort == cc) |> mutate(grp = factor(grp, levels = c("No BPD", "Mild", "BPD")))
  if (cc == "GSE32472") d$grp <- factor(ifelse(d$grp == "BPD", "Moderate/severe", as.character(d$grp)),
    levels = c("No BPD", "Mild", "Moderate/severe"))
  s <- d |> group_by(grp, t) |> summarise(m = mean(s), se = sd(s) / sqrt(n()),
    n = n()) |> mutate(time = factor(TL[[cc]][t + 1], levels = TL[[cc]]))
  pp <- TCP |> filter(cohort == cc)
  ggplot(s, aes(time, m, colour = grp, group = grp)) + geom_hline(yintercept = 0, colour = LIGHT,
    linewidth = .3) +
    geom_line(linewidth = .6,
      position = position_dodge(.15)) + geom_pointrange(aes(ymin = m - se, ymax = m + se), size = .25,
        linewidth = .45, position = position_dodge(.15)) +
    scale_colour_manual(values = c("No BPD" = C_NO, Mild = "#8aa2cc", "Moderate/severe" = C_BPD,
      BPD = C_BPD), name = NULL) +
    labs(x = NULL, y = "Naive T-cell module (mean, SE)",
      title = paste0(cc, if (cc == "GSE32472") ": genes from GSE220135" else ": genes from GSE32472"),
      subtitle = sprintf("BPD vs none after GA, all time points: %.2f SD, p %s", -pp$b,
        format(signif(pp$p, 2)))) +
    coord_cartesian(ylim = c(-1.75, .9)) + TH + theme(legend.position = "top",
      legend.justification = "left", legend.margin = margin(0, 0, -2, 0))
}
fig <- ((fA | fB) + plot_layout(widths = c(1.15,
  1))) / (tc("GSE32472") | tc("GSE220135")) + plot_layout(heights = c(1.25,
    .9)) + plot_annotation(tag_levels = "A")
savefig(fig, "NewFigure4", 7.4, 6.8)
print(RB |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
