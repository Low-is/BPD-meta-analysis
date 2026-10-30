# Evaluate SPON1 after its selection in the three GA cohorts; in-sample performance is optimistic.
# Assess leave-one-GA-cohort-out selection, other cohorts, later time points and transferred model weights.
source("./00_setup.R")
suppressPackageStartupMessages(library(nlme))
load(file.path(OUT, "fig3.RData"))
set.seed(20261029)

ga_of <- function(k) PH[[k]]$ga_weeks
G <- "SPON1"
cat("present:", sapply(c(BLOOD, "GSE302178"), function(k) G %in% rownames(E[[k]])), "\n")
S_adj <- lapply(setNames(GA_COH, GA_COH),
  function(k) cohort_stats(k, lab(k), covars = data.frame(ga = ga_of(k))))
# (a) leave-one-GA-cohort-out ranking of SPON1 by GA-adjusted p among genes same-direction in both training cohorts
for (h in GA_COH) {
  M <- meta_all(bind_rows(S_adj[setdiff(GA_COH, h)]),
    kmin = 2) |> filter(n_pos == k | n_neg == k) |> arrange(p)
  r <- which(M$gene == G)
  cat(h, "held out: SPON1 rank", r, "of", nrow(M), " g =", round(M$es[r], 2), " p =", signif(M$p[r], 2),
    " FDR =", signif(M$fdr[r], 2), "\n")
}
M3 <- meta_all(bind_rows(S_adj), kmin = 3) |> filter(gene == G)
print(M3)
# per-cohort unadjusted and adjusted, AUC
for (k in c(BLOOD, "GSE302178")) {
  if (!G %in% rownames(E[[k]])) next
  s <- E[[k]][G, ]
  y <- lab(k)
  r <- roc_stats(s, y)
  ga_p <- if (!all(is.na(ga_of(k)))) signif(coef(summary(glm(y ~ ga_of(k) + s, binomial)))["s", 4],
    2) else NA
  cat(k, " n1/n0", sum(y), sum(1 - y), " AUC", round(r[["auc"]], 2), "(", round(r[["lo"]], 2), "-",
    round(r[["hi"]], 2), ")  p(SPON1 | GA)", ga_p,
    " r neutro", if (!is.null(CS[[k]])) round(cor(s, CS[[k]][, "Neutrophils"]), 2) else NA, " r Tcell",
    if (!is.null(CS[[k]])) round(cor(s, CS[[k]][, "T cells"]), 2) else NA,
    " r GA", if (!all(is.na(ga_of(k)))) round(cor(s, ga_of(k), use = "c"), 2) else NA, "\n")
}
cat("\ncell correlations GSE32472:\n")
print(round(cor(E$GSE32472[G, ], CS$GSE32472), 2))
# (c) held-out gain over GA (weights from the other two GA cohorts; SPON1 z-scored within cohort)
A <- function(y, s) as.numeric(auc(roc(y, s, quiet = TRUE, direction = "<", levels = c(0, 1))))
for (h in GA_COH) {
  tr <- setdiff(GA_COH, h)
  dtr <- bind_rows(lapply(tr,
    function(k) data.frame(y = lab(k), ga = ga_of(k), s = as.numeric(scale(E[[k]][G, ])), coh = k)))
  f0 <- glm(y ~ coh + ga, binomial, dtr)
  f1 <- glm(y ~ coh + ga + s, binomial, dtr)
  d <- data.frame(y = lab(h), ga = ga_of(h), s = as.numeric(scale(E[[h]][G, ])))
  lp0 <- coef(f0)["ga"] * d$ga
  lp1 <- coef(f1)["ga"] * d$ga + coef(f1)["s"] * d$s
  B <- replicate(2000, {
    i <- c(sample(which(d$y == 1), replace = TRUE), sample(which(d$y == 0), replace = TRUE))
    A(d$y[i], lp1[i]) - A(d$y[i], lp0[i])
  })
  cat(h, "held out: AUC GA", round(A(d$y, lp0), 3), " GA+SPON1", round(A(d$y, lp1), 3), " dAUC",
    round(A(d$y, lp1) - A(d$y, lp0), 3), "(", round(quantile(B, .025), 3), "to",
    round(quantile(B, .975), 3), ")\n")
}
# (d) later time points (not used for selection)
L <- readRDS(file.path(DATA, "longitudinal_GSE32472_GSE220135.rds"))
for (cc in c("GSE32472", "GSE220135")) {
  a <- L$annotation |> filter(cohort == cc, !is.na(t), !is.na(ga), severity != "Mild",
    id %in% colnames(L$expression[[cc]])) |> mutate(y = as.integer(severity != "No BPD"))
  a$s <- as.numeric(scale(L$expression[[cc]][G, a$id]))
  for (tt in 0:2) {
    d <- a |> filter(t == tt)
    p <- coef(summary(glm(y ~ ga + s, binomial, d)))["s", 4]
    cat(cc, "t", tt, " n1/n0", sum(d$y), sum(1 - d$y), " AUC SPON1", round(A(d$y, d$s), 2), " AUC GA",
      round(A(d$y, -d$ga), 2), " p(SPON1 | GA)", signif(p, 2), "\n")
  }
  m <- lme(s ~ y + factor(t) + ga, random = ~ 1 | infant, data = a, method = "REML")
  cat(cc, "mixed model, all time points, GA-adjusted: b", round(fixef(m)[["y"]], 2), " p",
    signif(summary(m)$tTable["y", 5], 2), "\n")
}
