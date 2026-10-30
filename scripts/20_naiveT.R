# Explore naive T-cell expression relative to total T-cell marker expression and GA.
source("./00_setup.R")
suppressPackageStartupMessages(library(nlme))
load(file.path(OUT, "bigstory.RData"))
load(file.path(OUT, "fig3.RData"))

ga_of <- function(k) PH[[k]]$ga_weeks
PAN_T <- c("CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK")
NAIVE <- c("CCR7", "SELL", "TCF7", "LEF1", "BCL11B", "CD27", "CD28", "IL7R", "BACH2",
  "SATB1")  # naive / stem-like T
RTE   <- c("PTK7", "CR2", "PECAM1", "SOX4", "TOX", "AIF1", "NOG", "CHI3L2", "SPON1", "GPA33", "MAN1C1",
  "DENND2D")  # recent thymic emigrant / naive-specific (lit.)
EFF   <- c("GZMB", "GZMH", "PRF1", "GNLY", "NKG7", "CCL5", "KLRG1", "CX3CR1", "FGFBP2",
  "GZMA")  # effector / memory cytotoxic

mz <- function(g, x) {
  g <- intersect(g, rownames(x))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z)))
}
cat("RTE genes in replicated GA-adjusted list (", nrow(TOPALL), "genes):", intersect(RTE, TOPALL$gene),
  "| naive:", intersect(NAIVE, TOPALL$gene), "| effector:", intersect(EFF, TOPALL$gene), "\n")
L <- readRDS(file.path(DATA, "longitudinal_GSE32472_GSE220135.rds"))
res <- NULL
for (cc in c("GSE32472", "GSE220135")) {
  a <- L$annotation |> filter(cohort == cc, !is.na(t), !is.na(ga), severity != "Mild",
    id %in% colnames(L$expression[[cc]])) |> mutate(y = as.integer(severity != "No BPD"))
  x <- L$expression[[cc]][, a$id]
  x <- x[apply(x, 1, sd) > 0, ]
  a$panT <- mz(PAN_T, x)
  a$naive <- mz(NAIVE, x)
  a$rte <- mz(RTE, x)
  a$eff <- mz(EFF, x)
  a$neu <- mz(MK$Neutrophils, x)
  a$rte_in_T <- as.numeric(scale(resid(lm(rte ~ panT, a))))
  a$naive_in_T <- as.numeric(scale(resid(lm(naive ~ panT, a))))
  for (v in c("panT", "naive", "rte", "eff", "rte_in_T", "naive_in_T")) for (tt in list(0, 0:2)) {
    d <- a |> filter(t %in% tt)
    d$s <- d[[v]]
    if (length(tt) == 1) {
      f <- lm(s ~ y + ga, d)
      b <- coef(summary(f))["y", c(1, 4)]
      f2 <- lm(s ~ y + ga + neu + panT, d)
      b2 <- if (v %in% c("panT")) c(NA, NA) else coef(summary(f2))["y", c(1, 4)]
    }
    else {
      m <- lme(s ~ y + factor(t) + ga, random = ~ 1 | infant, data = d, method = "REML")
      b <- summary(m)$tTable["y", c(1, 5)]
      m2 <- if (v == "panT") NULL else lme(s ~ y + factor(t) + ga + neu + panT, random = ~ 1 | infant,
        data = d, method = "REML")
      b2 <- if (is.null(m2)) c(NA, NA) else summary(m2)$tTable["y", c(1, 5)]
    }
    res <- rbind(res,
      data.frame(cohort = cc, score = v, samples = ifelse(length(tt) == 1, "week 1", "all"), b_GA = b[1],
        p_GA = b[2], b_GA_cells = b2[1], p_GA_cells = b2[2]))
  }
  cat(cc, "r(rte, GA) in no-BPD week 1:", round(with(a |> filter(t == 0, y == 0), cor(rte, ga)), 2),
    " r(rte_in_T, GA):", round(with(a |> filter(t == 0, y == 0), cor(rte_in_T, ga)), 2), "\n")
}
rownames(res) <- NULL
print(res |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
# 4-cohort week-1 check of the RTE-within-T score (GA-adjusted where GA exists)
cat("\nWeek-1 RTE-within-T score, all cohorts:\n")
for (k in BLOOD) {
  x <- E[[k]]
  s <- as.numeric(scale(resid(lm(mz(RTE, x) ~ mz(PAN_T, x)))))
  y <- lab(k)
  ga_p <- if (k %in% GA_COH) signif(coef(summary(lm(s ~ y + ga_of(k))))["y", 4], 2) else NA
  cat(k, " AUC (lower = BPD):",
    round(as.numeric(auc(roc(y, -s, quiet = TRUE, direction = "<", levels = c(0, 1)))), 2), " GA-adj p",
    ga_p, " RTE genes found:", sum(RTE %in% rownames(x)), "\n")
}
save(res, file = file.path(OUT, "naiveT.RData"))
