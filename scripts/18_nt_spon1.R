# Combine the marker-based NT score and SPON1 using equal weights or training-cohort model coefficients.
# SPON1 was selected in the three GA cohorts; its selection is not independent of those cohorts.
source("./00_setup.R")
suppressPackageStartupMessages(library(nlme))
load(file.path(OUT, "fig3.RData"))
set.seed(20261029)

ga_of <- function(k) PH[[k]]$ga_weeks

mz <- function(g, x) {
  g <- intersect(g, rownames(x))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  colMeans(z)
}

NT <- function(x) as.numeric(scale(mz(MK$Neutrophils, x) - mz(c(MK$`T cells`, MK$`CD4 T`, MK$`CD8 T`),
  x)))

SP <- function(x) -as.numeric(scale(x["SPON1", ]))

CMB <- function(x) as.numeric(scale(NT(x) + SP(x)))

A <- function(y, s) as.numeric(auc(roc(y, s, quiet = TRUE, direction = "<", levels = c(0, 1))))
cat("r(NT, -SPON1):", sapply(BLOOD, function(k) round(cor(NT(E[[k]]), SP(E[[k]])), 2)), "\n")
for (nm in c("NT", "SP", "CMB")) {
  f <- get(nm)
  r <- t(sapply(BLOOD, function(k) roc_stats(f(E[[k]]), lab(k))))
  pa <- pool_auc(r[, "auc"], r[, "var"])
  ad <- t(sapply(GA_COH, function(k) {
    s <- f(E[[k]])
    m <- lm(s ~ lab(k) + ga_of(k))
    c(coef(m)[2], vcov(m)[2, 2])
  }))
  pg <- re_pool(ad[, 1], ad[, 2])
  cat(sprintf("%-4s AUC %.2f (%.2f-%.2f) [%s] | GA-adj %.2f SD p %.2g\n", nm, pa[1], pa[2], pa[3],
    paste(round(r[, "auc"], 2), collapse = "/"), pg[["es"]], pg[["p"]]))
}
# held-out gain over GA: (i) equal-weight combined score as one term; (ii) GA + NT + SPON1 with separate weights
for (spec in c("one term: CMB", "two terms: NT + SPON1")) {
  ho <- t(sapply(GA_COH, function(h) {
    tr <- setdiff(GA_COH, h)
    mk <- function(k) data.frame(y = lab(k), ga = ga_of(k), a = NT(E[[k]]), b = SP(E[[k]]),
      c = CMB(E[[k]]), coh = k)
    dtr <- bind_rows(lapply(tr, mk))
    d <- mk(h)
    f0 <- glm(y ~ coh + ga, binomial, dtr)
    lp0 <- coef(f0)["ga"] * d$ga
    if (spec == "one term: CMB") {
      f1 <- glm(y ~ coh + ga + c, binomial, dtr)
      lp1 <- coef(f1)["ga"] * d$ga + coef(f1)["c"] * d$c
    }
    else {
      f1 <- glm(y ~ coh + ga + a + b, binomial, dtr)
      lp1 <- coef(f1)["ga"] * d$ga + coef(f1)["a"] * d$a + coef(f1)["b"] * d$b
    }
    B <- replicate(2000, {
      i <- c(sample(which(d$y == 1), replace = TRUE), sample(which(d$y == 0), replace = TRUE))
      A(d$y[i], lp1[i]) - A(d$y[i], lp0[i])
    })
    c(ga = A(d$y, lp0), full = A(d$y, lp1), d = A(d$y, lp1) - A(d$y, lp0), v = var(B))
  }))
  p <- re_pool(ho[, "d"], ho[, "v"])
  p2 <- re_pool(ho[2:3, "d"], ho[2:3, "v"])
  cat(sprintf("\n%s: held-out AUC GA -> GA+score: %s\n  pooled dAUC %+.3f (%.3f to %.3f, p %.3f) | folds where SPON1 would be selected (225881, 220135): %+.3f (p %.2f)\n", spec,
    paste(sprintf("%.2f->%.2f", ho[, "ga"], ho[, "full"]), collapse = ", "), p[["es"]], p[["lo"]],
    p[["hi"]], p[["p"]], p2[["es"]], p2[["p"]]))
}
# longitudinal
L <- readRDS(file.path(DATA, "longitudinal_GSE32472_GSE220135.rds"))
for (cc in c("GSE32472", "GSE220135")) {
  a <- L$annotation |> filter(cohort == cc, !is.na(t), !is.na(ga), severity != "Mild",
    id %in% colnames(L$expression[[cc]])) |> mutate(y = as.integer(severity != "No BPD"))
  x <- L$expression[[cc]][, a$id]
  a$s <- CMB(x)
  m <- lme(s ~ y + factor(t) + ga, random = ~ 1 | infant, data = a, method = "REML")
  cat(cc, "combined, mixed GA-adj: b", round(fixef(m)[["y"]], 2), "p",
    signif(summary(m)$tTable["y", 5], 2), "| AUC by time:", sapply(0:2, function(tt) {
      d <- a |> filter(t == tt)
      round(A(d$y, d$s), 2)
    }),
    "| p(score | GA) by time:", sapply(0:2, function(tt) {
      d <- a |> filter(t == tt)
      signif(coef(summary(glm(y ~ ga + s, binomial, d)))["s", 4], 2)
    }), "\n")
}
