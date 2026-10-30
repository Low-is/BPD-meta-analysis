# Evaluate neutrophil-minus-T-cell expression scores.
# Marker score: mean z(neutrophil markers) minus mean z(T-cell markers), without outcome-based selection.
# Program score: signed program cores, with a full-data fallback when the neutrophil core is missing.
source("./00_setup.R")
suppressPackageStartupMessages(library(nlme))
load(file.path(OUT, "fig3.RData"))
load(file.path(OUT, "heldout4_relaxed.RData"))
load(file.path(OUT, "programs4.RData"))
load(file.path(OUT, "fig5.RData"))
set.seed(20261029)

ga_of <- function(k) PH[[k]]$ga_weeks

mz <- function(g, x) {
  g <- intersect(g, rownames(x))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  colMeans(z)
}

NT_mark <- function(x) as.numeric(scale(mz(MK$Neutrophils,
  x) - mz(c(MK$`T cells`, MK$`CD4 T`, MK$`CD8 T`), x)))

psc <- function(core, x) {
  if (is.null(core)) return(rep(NA_real_, ncol(x)))
  g <- core$gene[core$gene %in% rownames(x)]
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  colMeans(z * sign(core$es[match(g, core$gene)]))
}

NT_prog <- function(k, x = E[[k]]) {
  a <- psc(FOLD[[k]][["Neutrophil antimicrobial"]], x)
  b <- psc(FOLD[[k]][["T cell activation"]], x)
  if (all(is.na(a))) a <- psc(PR4$core[[which(PR4$program == "Neutrophil antimicrobial")]],
    x)  # Use the full-data core when the fold-specific neutrophil core is missing.
  as.numeric(scale(as.numeric(scale(a)) + as.numeric(scale(b))))
}  # both cores are signed BPD-positive, so the sum = neutrophil high + T cell low
cat("marker genes found:", sapply(BLOOD, function(k) sum(MK$Neutrophils %in% rownames(E[[k]]))), "/",
  sapply(BLOOD, function(k) sum(c(MK$`T cells`, MK$`CD4 T`, MK$`CD8 T`) %in% rownames(E[[k]]))), "\n")

A <- function(y, s) as.numeric(auc(roc(y, s, quiet = TRUE, direction = "<", levels = c(0, 1))))
SC <- list("NT markers (label-free)" = function(k) NT_mark(E[[k]]),
  "NT programs (full-data fallback)" = NT_prog,
  "Neutrophil markers only" = function(k) as.numeric(scale(mz(MK$Neutrophils, E[[k]]))),
  "T cell markers only (reversed)" = function(k) -as.numeric(scale(mz(c(MK$`T cells`, MK$`CD4 T`,
    MK$`CD8 T`), E[[k]]))),
  "SPON1 (reversed)" = function(k) -as.numeric(scale(E[[k]]["SPON1", ])))
for (nm in names(SC)) {
  r <- t(sapply(BLOOD, function(k) roc_stats(SC[[nm]](k), lab(k))))
  pa <- pool_auc(r[, "auc"], r[, "var"])
  ad <- t(sapply(GA_COH, function(k) {
    s <- SC[[nm]](k)
    f <- lm(s ~ lab(k) + ga_of(k))
    c(coef(f)[2], vcov(f)[2, 2])
  }))
  pg <- re_pool(ad[, 1], ad[, 2])
  ho <- t(sapply(GA_COH, function(h) {
    tr <- setdiff(GA_COH, h)
    dtr <- bind_rows(lapply(tr,
      function(k) data.frame(y = lab(k), ga = ga_of(k), s = SC[[nm]](k), coh = k)))
    f0 <- glm(y ~ coh + ga, binomial, dtr)
    f1 <- glm(y ~ coh + ga + s, binomial, dtr)
    d <- data.frame(y = lab(h), ga = ga_of(h), s = SC[[nm]](h))
    lp0 <- coef(f0)["ga"] * d$ga
    lp1 <- coef(f1)["ga"] * d$ga + coef(f1)["s"] * d$s
    B <- replicate(1000, {
      i <- c(sample(which(d$y == 1), replace = TRUE), sample(which(d$y == 0), replace = TRUE))
      A(d$y[i], lp1[i]) - A(d$y[i], lp0[i])
    })
    c(d = A(d$y, lp1) - A(d$y, lp0), v = var(B))
  }))
  ph <- re_pool(ho[, 1], ho[, 2])
  cat(sprintf("%-32s AUC %.2f (%.2f-%.2f) [%s] | GA-adj effect %.2f SD (p %.2g) | held-out dAUC over GA %+.3f (%.3f to %.3f, p %.2f) [%s]\n", nm, pa[1], pa[2], pa[3],
    paste(round(r[, "auc"], 2), collapse = "/"), pg[["es"]], pg[["p"]], ph[["es"]], ph[["lo"]],
    ph[["hi"]], ph[["p"]], paste(sprintf("%+.2f", ho[, 1]), collapse = "/")))
}
# correlation of the label-free NT score with published signatures and our programs (GSE220135 + GSE32472)
cat("\nr(NT markers, GA):", sapply(GA_COH, function(k) round(cor(NT_mark(E[[k]]), ga_of(k)), 2)), "\n")
# longitudinal: does NT stay high in BPD after GA (mixed model)?
L <- readRDS(file.path(DATA, "longitudinal_GSE32472_GSE220135.rds"))
for (cc in c("GSE32472", "GSE220135")) {
  a <- L$annotation |> filter(cohort == cc, !is.na(t), !is.na(ga), severity != "Mild",
    id %in% colnames(L$expression[[cc]])) |> mutate(y = as.integer(severity != "No BPD"))
  a$s <- NT_mark(L$expression[[cc]][, a$id])
  m <- lme(s ~ y + factor(t) + ga, random = ~ 1 | infant, data = a, method = "REML")
  m2 <- lme(s ~ y * factor(t) + ga, random = ~ 1 | infant, data = a, method = "ML")
  m1 <- update(m2, s ~ y + factor(t) + ga)
  cat(cc, "mixed, GA-adj: b", round(fixef(m)[["y"]], 2), "p", signif(summary(m)$tTable["y", 5], 2),
    "| interaction p", signif(anova(m1, m2)$`p-value`[2], 2),
    "| AUC by time:", sapply(0:2, function(tt) {
      d <- a |> filter(t == tt)
      round(A(d$y, d$s), 2)
    }), "\n")
}
