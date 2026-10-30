# Neutrophil-T-cell expression balance analyses.
source("./00_setup.R")
suppressPackageStartupMessages({
  library(nlme)
  library(ggplot2)
})
set.seed(20261029)

ga_of <- function(k) PH[[k]]$ga_weeks
NEU <- c("CSF3R", "CXCR2", "FPR1", "FCGR3B", "MME", "CMTM2")
PAN_T <- c("CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK")
PAIRS <- list(c("CSF3R", "CD3E"), c("CXCR2", "CD3D"), c("FPR1", "CD3G"), c("CSF3R", "LCK"))

A <- function(y, s) as.numeric(auc(roc(y, s, quiet = TRUE, direction = "<", levels = c(0, 1))))

lm_mean <- function(g,
  x) colMeans(x[intersect(g, rownames(x)), ,
    drop = FALSE])  # mean log expression (supplied log scale)

mz <- function(g, x) {
  g <- intersect(g, rownames(x))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z)))
}

comp <- function(x) {
  n <- lm_mean(NEU, x)
  t <- lm_mean(PAN_T, x)
  list(N = as.numeric(scale(n)), Tneg = -as.numeric(scale(t)), BAL = as.numeric(scale(n - t)), Nraw = n,
    Traw = t)
}
C <- lapply(setNames(BLOOD, BLOOD), function(k) comp(E[[k]]))

# ---- B1 ----
cat("== B1: balance (week 1) ==\n")
for (k in BLOOD) {
  y <- lab(k)
  b <- C[[k]]$BAL
  ga <- ga_of(k)
  cat(sprintf("%-10s AUC %.2f | raw balance BPD %.2f vs none %.2f", k, A(y, b),
    mean(C[[k]]$Nraw[y == 1] - C[[k]]$Traw[y == 1]), mean(C[[k]]$Nraw[y == 0] - C[[k]]$Traw[y == 0])),
    if (k %in% GA_COH) sprintf("| GA-adj b %+.2f p %.2g | r(GA) no-BPD %+.2f", coef(lm(b ~ y + ga))[2],
      coef(summary(lm(b ~ y + ga)))[2, 4], cor(b[y == 0], ga[y == 0])) else "", "\n")
}

pool_w1 <- function(v) {
  r <- t(sapply(BLOOD, function(k) roc_stats(C[[k]][[v]], lab(k))))
  pool_auc(r[, "auc"], r[, "var"])
}

ga_adj <- function(v) {
  r <- t(sapply(GA_COH, function(k) {
    s <- C[[k]][[v]]
    f <- lm(s ~ lab(k) + ga_of(k))
    c(coef(f)[2], vcov(f)[2, 2])
  }))
  re_pool(r[, 1], r[, 2])
}

# ---- B2: identical held-out infants ----
heldout <- function(make) {
  ho <- t(sapply(GA_COH, function(h) {
    tr <- setdiff(GA_COH, h)
    dtr <- bind_rows(lapply(tr, function(k) data.frame(y = lab(k), ga = ga_of(k), make(k), coh = k)))
    d <- data.frame(y = lab(h), ga = ga_of(h), make(h))
    vars <- setdiff(names(d), c("y", "ga"))
    f0 <- glm(y ~ coh + ga, binomial, dtr)
    f1 <- glm(as.formula(paste("y ~ coh + ga +", paste(vars, collapse = "+"))), binomial, dtr)
    lp0 <- coef(f0)["ga"] * d$ga
    lp1 <- coef(f1)["ga"] * d$ga + as.matrix(d[, vars, drop = FALSE]) %*% coef(f1)[vars]
    B <- replicate(1000, {
      i <- c(sample(which(d$y == 1), replace = TRUE), sample(which(d$y == 0), replace = TRUE))
      A(d$y[i], lp1[i]) - A(d$y[i], lp0[i])
    })
    c(d = A(d$y, lp1) - A(d$y, lp0), v = var(B))
  }))
  p <- re_pool(ho[, "d"], ho[, "v"])
  c(dAUC = p[["es"]], lo = p[["lo"]], hi = p[["hi"]], p = p[["p"]])
}
SPEC <- list("T cell alone" = function(k) data.frame(a = C[[k]]$Tneg),
  "Neutrophil alone" = function(k) data.frame(a = C[[k]]$N),
  "Fixed balance" = function(k) data.frame(a = C[[k]]$BAL),
  "Both, fitted weights" = function(k) data.frame(a = C[[k]]$N, b = C[[k]]$Tneg))
B2 <- bind_rows(lapply(names(SPEC), function(nm) {
  v <- c("T cell alone" = "Tneg", "Neutrophil alone" = "N", "Fixed balance" = "BAL",
    "Both, fitted weights" = NA)[nm]
  pa <- if (!is.na(v)) pool_w1(v) else rep(NA, 4)
  ga <- if (!is.na(v)) ga_adj(v) else c(es = NA, p = NA)
  data.frame(model = nm, AUC_alone_4coh = pa[1], GA_adj_b = ga[["es"]], GA_adj_p = ga[["p"]],
    t(heldout(SPEC[[nm]])))
}))
cat("\n== B2: components vs balance (held-out over GA, 3 GA cohorts) ==\n")
print(B2 |> mutate(across(where(is.numeric), ~ signif(.x, 2))), row.names = FALSE)

# ---- B3: declared pairs; selection within training cohorts ----
pd <- function(pr, x) if (all(pr %in% rownames(x))) as.numeric(x[pr[1], ] - x[pr[2], ]) else rep(NA,
  ncol(x))
B3 <- bind_rows(lapply(BLOOD, function(h) {
  tr <- setdiff(BLOOD, h)
  sc <- sapply(PAIRS, function(pr) {
    r <- t(sapply(tr, function(k) roc_stats(pd(pr, E[[k]]), lab(k))))
    if (all(is.na(r[, "auc"]))) NA else pool_auc(r[, "auc"], r[, "var"])[1]
  })
  best <- PAIRS[[which.max(sc)]]
  s <- pd(best, E[[h]])
  y <- lab(h)
  data.frame(heldout = h, chosen = paste(best, collapse = "-"), train_AUC = max(sc, na.rm = TRUE),
    test_AUC_diff = A(y, s),
    test_AUC_rule = A(y, as.numeric(s > 0) + rnorm(length(s), 0, 1e-9)), frac_rule_true = mean(s > 0))
}))
cat("\n== B3: declared pairs (best chosen in training) ==\n")
print(B3 |> mutate(across(where(is.numeric), ~ signif(.x, 2))), row.names = FALSE)
cat("Pooled held-out AUC, chosen pair difference:",
  round(pool_auc(B3$test_AUC_diff,
    sapply(BLOOD,
      function(h) roc_stats(pd(strsplit(B3$chosen[B3$heldout == h], "-")[[1]], E[[h]]),
        lab(h))[["var"]]))[1:3], 2), "\n")

# ---- B4: 2D view ----
D4 <- bind_rows(lapply(BLOOD,
  function(k) data.frame(cohort = k, N = C[[k]]$N, T = -C[[k]]$Tneg,
    BPD = factor(lab(k), 0:1, c("No BPD", "BPD")), ga = ga_of(k))))
p4 <- ggplot(D4, aes(T, N, colour = BPD, size = ga)) + geom_point(alpha = .7) + facet_wrap(~ cohort,
  nrow = 1) + scale_size(range = c(2.2, .5), name = "GA (wk)") +
  scale_colour_manual(values = c("No BPD" = "#c8553d",
    BPD = "#1f3a6e")) + labs(x = "T-cell expression (z)",
      y = "Neutrophil expression (z)") + theme_classic(base_size = 8)
ggsave("./figs/v7_B4_neutrophil_vs_T.png", p4, width = 8, height = 2.6, dpi = 200, bg = "white")
cat("\n== B4: r(N, T) by cohort:", sapply(BLOOD, function(k) round(cor(C[[k]]$N, -C[[k]]$Tneg), 2)),
  "\n")

# ---- B5: exploratory PCA; updated in 28_axis_pca.R ----
S32 <- S_blood |> filter(gse == "GSE32472")
PUB <- list(Luo_PLoS = c("ACSL1", "GALNT14", "WIPI1", "MAPK14", "PROK2", "CREB5"),
  Fang = c("LPIN1", "ACADSB", "WIPI1", "SLC7A11"),
  Zou = c("ARG1", "CSF3R", "IL1R1", "IL1R2", "MMP9", "RETN", "S100A12", "SOCS3"),
  Zhang = c("MMP9", "SIGLEC5", "DYSF", "MGAM", "S100A12"),
  Guo = c("IL7R", "CXCL10", "DEFA4", "PRTN3", "NCAPG", "CCNB1"), Tu_FP = c("PELI1", "PYGL", "STAT4"),
  Sun = c("CCDC141", "CHI3L2", "RIMKLB", "PDLIM7"), Jia = c("NFATC3", "ERMN", "PLA2G4A"),
  Abilbayeva = c("BTN3A3", "AK5", "NOG", "GCSAM"),
  Authors5 = c("PNPO", "MSANTD2", "CD4", "SNX1", "P2RX7"),
  Th2 = c("UGT2B28", "SYNJ1", "SORD", "PTGDR2", "LGALS7", "FBN1", "TIPARP", "KCNJ8", "TSPAN11"))

sig <- function(g, x) {
  d <- sign(S32$g[match(g, S32$gene)])
  ok <- g %in% rownames(x) & !is.na(d)
  if (sum(ok) == 0) return(rep(NA, ncol(x)))
  z <- t(scale(t(x[g[ok], , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z * d[ok])))
}
B5 <- bind_rows(lapply(setdiff(BLOOD, "GSE32472"), function(k) {
  M <- sapply(PUB, sig, x = E[[k]])
  M <- M[, colSums(is.na(M)) == 0, drop = FALSE]
  pc <- prcomp(M, scale. = TRUE)
  bal <- C[[k]]$BAL
  pc1 <- pc$x[, 1] * sign(cor(pc$x[, 1], bal))
  data.frame(cohort = k, n_sigs = ncol(M), PC1_var = summary(pc)$importance[2, 1],
    r_BAL_PC1 = cor(bal, pc1),
    median_R2_sig_on_BAL = median(apply(M, 2, function(s) summary(lm(s ~ bal))$r.squared)))
}))
cat("\n== B5: published signatures vs balance (test cohorts) ==\n")
print(B5 |> mutate(across(where(is.numeric), ~ signif(.x, 2))), row.names = FALSE)
B5b <- bind_rows(lapply(names(PUB), function(nm) {
  r <- t(sapply(c("GSE225881", "GSE220135"), function(k) {
    s <- sig(PUB[[nm]], E[[k]])
    if (all(is.na(s))) return(c(NA, NA))
    f <- lm(s ~ lab(k) + ga_of(k) + C[[k]]$BAL)
    c(coef(f)[2], vcov(f)[2, 2])
  }))
  r <- r[!is.na(r[, 1]), , drop = FALSE]
  if (!nrow(r)) return(NULL)
  p <- re_pool(r[, 1], r[, 2])
  data.frame(signature = nm, b_after_GA_BAL = p[["es"]], p = p[["p"]])
}))
cat("Signature-BPD association after GA + balance (GSE225881 + GSE220135):\n")
print(B5b |> mutate(across(where(is.numeric), ~ signif(.x, 2))), row.names = FALSE)

# ---- B6: methylation cross-check (GSE220135) ----
L <- readRDS(file.path(DATA, "longitudinal_GSE32472_GSE220135.rds"))
CF <- readRDS(file.path(DATA, "GSE220135_cell_fractions.rds"))
a <- L$annotation |> filter(cohort == "GSE220135", !is.na(t), !is.na(ga), !is.na(severity),
  id %in% colnames(L$expression$GSE220135)) |> mutate(y = as.integer(severity != "No BPD")) |>
  left_join(CF, by = c("id" = "gsm")) |> filter(!is.na(Gran))
x <- L$expression$GSE220135[, a$id]
a$BAL <- as.numeric(scale(lm_mean(NEU, x) - lm_mean(PAN_T, x)))
a$mNLR <- log((a$Gran + 1) / (a$CD4T + a$CD8T + a$B + a$NK + 1))
a$mNTR <- log((a$Gran + 1) / (a$CD4T + a$CD8T + 1))
cat("\n== B6: GSE220135 methylation vs RNA balance (", nrow(a), "samples ) ==\n")
cat("r(BAL, log methylation NLR):", round(cor(a$BAL, a$mNLR), 2), "| r(BAL, log methylation N/T):",
  round(cor(a$BAL, a$mNTR), 2), "| by time:",
  sapply(0:2, function(tt) round(with(a[a$t == tt, ], cor(BAL, mNTR)), 2)), "\n")
for (v in c("mNLR", "mNTR", "BAL")) {
  a$s <- as.numeric(scale(a[[v]]))
  m0 <- lme(s ~ y + factor(t), random = ~ 1 | infant, data = a, method = "REML")
  m1 <- update(m0, s ~ y + factor(t) + ga)
  w <- a |> filter(t == 0)
  cat(sprintf("%-5s mixed: BPD b %+.2f (p %.2g) | +GA %+.2f (p %.2g) | week-1 AUC %.2f\n", v,
    fixef(m0)[["y"]], summary(m0)$tTable["y", 5], fixef(m1)[["y"]], summary(m1)$tTable["y", 5],
    A(w$y, w$s)))
}
save(B2, B3, B5, B5b, D4, file = file.path(OUT, "v7_partB.RData"))
