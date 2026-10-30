# Evaluate top-gene prediction, gain beyond GA and published blood signatures.
# Rank top-N genes in training cohorts using random-effects meta-analysis and directional concordance.
# Evaluate transferred GA and score coefficients with paired bootstrap intervals for the AUC difference.
source("./00_setup.R")
load(file.path(OUT, "heldout4_relaxed.RData"))
load(file.path(OUT, "programs4.RData"))
load(file.path(OUT, "fig3.RData"))
set.seed(20261029)

ga_of <- function(k) PH[[k]]$ga_weeks

sig_score <- function(genes, dirs, x) {
  ok <- genes %in% rownames(x)
  g <- genes[ok]
  if (!length(g)) return(rep(NA_real_, ncol(x)))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z * dirs[ok])))
}

rank_tr <- function(tr) {
  m <- length(tr)
  meta_all(S_blood |> filter(gse %in% tr), kmin = max(1, m - 1)) |> filter(I2 < 75,
    ifelse(es > 0, n_pos, n_neg) >= max(1, m - 1)) |> arrange(p)
}
NGRID <- c(3, 5, 10, 20, 50, 100, 250, 500, 1000)

# ---- A: top-N held-out ----
RK <- lapply(setNames(BLOOD, BLOOD), function(h) rank_tr(setdiff(BLOOD, h)))
RK_all <- rank_tr(BLOOD)
TOPN <- bind_rows(lapply(NGRID, function(N) bind_rows(lapply(c("held-out", "in-sample"), function(tp) {
  r <- t(sapply(BLOOD, function(h) {
    R <- head(if (tp == "held-out") RK[[h]] else RK_all, N)
    roc_stats(sig_score(R$gene, sign(R$es), E[[h]]), lab(h))
  }))
  pa <- pool_auc(r[, "auc"], r[, "var"])
  data.frame(N = N, type = tp, auc = pa[1], lo = pa[2], hi = pa[3],
    cohorts = paste(round(r[, "auc"], 2), collapse = "/"))
}))))
GA_AUC <- {
  r <- t(sapply(GA_COH, function(k) roc_stats(-ga_of(k), lab(k))))
  pa <- pool_auc(r[, "auc"], r[, "var"])
  c(pa[1:3], cohorts = paste(round(r[, "auc"], 2), collapse = "/"))
}
print(TOPN |> mutate(across(where(is.numeric), ~ round(.x, 3))))
print(GA_AUC)

# ---- B: beyond GA, fully held out ----
heldout_pred <- function(h, sfun) {
  tr <- setdiff(GA_COH, h)
  dtr <- bind_rows(lapply(tr, function(k) data.frame(y = lab(k), ga = ga_of(k), s = sfun(k), coh = k)))
  f0 <- glm(y ~ coh + ga, binomial, dtr)
  f1 <- glm(y ~ coh + ga + s, binomial, dtr)
  dh <- data.frame(y = lab(h), ga = ga_of(h), s = sfun(h))
  lp0 <- coef(f0)["ga"] * dh$ga
  lp1 <- coef(f1)["ga"] * dh$ga + coef(f1)["s"] * dh$s
  A <- function(i,
    lp) as.numeric(auc(roc(dh$y[i], lp[i], quiet = TRUE, direction = "<", levels = c(0, 1))))
  ii <- seq_len(nrow(dh))
  B <- replicate(2000, {
    i <- c(sample(which(dh$y == 1), replace = TRUE), sample(which(dh$y == 0), replace = TRUE))
    A(i, lp1) - A(i, lp0)
  })
  c(auc_ga = A(ii, lp0), auc_ga_score = A(ii, lp1), d_auc = A(ii, lp1) - A(ii, lp0),
    d_lo = quantile(B, .025), d_hi = quantile(B, .975), d_var = var(B), n1 = sum(dh$y),
    n0 = sum(1 - dh$y))
}
S_adj <- lapply(setNames(GA_COH, GA_COH),
  function(k) cohort_stats(k, lab(k), covars = data.frame(ga = ga_of(k))))

rank_adj <- function(tr) meta_all(bind_rows(S_adj[tr]),
  kmin = length(tr)) |> filter(n_pos == k | n_neg == k) |> arrange(p)

psc <- function(core, x) {
  if (is.null(core)) return(rep(NA_real_, ncol(x)))
  g <- core$gene[core$gene %in% rownames(x)]
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  as.numeric(scale(colMeans(z * sign(core$es[match(g, core$gene)]))))
}
NTOP <- 50
METH <- list(
  "Top 50 genes" = function(h) {
    R <- head(RK[[h]], NTOP)
    function(k) sig_score(R$gene, sign(R$es), E[[k]])
  },
  "Top 50 genes, GA-adjusted ranking" = function(h) {
    R <- head(rank_adj(setdiff(GA_COH, h)), NTOP)
    function(k) sig_score(R$gene, sign(R$es), E[[k]])
  },
  "T cell activation program" = function(h) function(k) psc(FOLD[[h]][["T cell activation"]], E[[k]]),
  "Combined program score" = function(h) function(k) as.numeric(scale(rowMeans(sapply(P4,
    function(pg) psc(FOLD[[h]][[pg]], E[[k]])), na.rm = TRUE))))
BEY <- bind_rows(lapply(names(METH), function(nm) bind_rows(lapply(GA_COH, function(h) {
  sf <- METH[[nm]](h)
  if (all(is.na(sf(h)))) return(NULL)
  data.frame(method = nm, heldout = h, t(heldout_pred(h, sf)), check.names = FALSE)
}))))
names(BEY) <- sub("\\.2\\.5%|\\.97\\.5%", "", names(BEY))
BEYP <- BEY |> group_by(method) |> summarise(r = list(re_pool(d_auc, d_var)), auc_ga = mean(auc_ga),
  auc_gs = mean(auc_ga_score), k = n()) |>
  mutate(d = sapply(r, `[[`, "es"), lo = sapply(r, `[[`, "lo"), hi = sapply(r, `[[`, "hi"),
    p = sapply(r, `[[`, "p")) |> select(-r)
print(BEY |> mutate(across(where(is.numeric), ~ round(.x, 3))))
print(BEYP |> mutate(across(where(is.numeric), ~ round(.x, 3))))

# ---- C: published-signature audit ----
PUB <- list(
  "Luo 2024, PLoS One" = c("ACSL1", "GALNT14", "WIPI1", "MAPK14", "PROK2", "CREB5"),
  "Luo 2024, Medicine" = c("CYYR1", "GALNT14", "OLAH"),
  "Fang 2024" = c("LPIN1", "ACADSB", "WIPI1", "SLC7A11"),
  "Zou 2024" = c("ARG1", "CSF3R", "IL1R1", "IL1R2", "MMP9", "RETN", "S100A12", "SOCS3"),
  "Zhang 2025" = c("MMP9", "SIGLEC5", "DYSF", "MGAM", "S100A12"),
  "Guo 2025" = c("IL7R", "CXCL10", "DEFA4", "PRTN3", "NCAPG", "CCNB1"),
  "Tu 2025, Front Pediatr" = c("PELI1", "PYGL", "STAT4"),
  "Sun 2026" = c("CCDC141", "CHI3L2", "RIMKLB", "PDLIM7"),
  "Jia 2023" = c("NFATC3", "ERMN", "PLA2G4A", "MTMR9LP", "LOC440700"),
  "Li 2024 (APP)" = "APP",
  "Abilbayeva 2026" = c("BTN3A3", "AK5", "NOG", "GCSAM"),
  "Ann Med 2024 (autophagy)" = c("WIPI1", "TOMM70A", "BAG3", "PRKCQ"),
  "Windhorst 2023" = c("LAT", "LOC100505876", "MRVI1", "KLHDC4", "ATP6V0E2", "FAM102A", "FBXO3", "OGDH",
    "SLC25A3", "GAL", "CAT", "SOAT2",
    "TNFRSF10C", "BCL2A1", "DUSP2", "PDLIM7", "POLD2", "CPLX1", "SCML1", "TMED1", "HEBP1", "RRN3",
    "ZEB1", "SLC2A11"),
  "Authors 2023, 5 genes" = c("PNPO", "MSANTD2", "CD4", "SNX1", "P2RX7"),
  "Authors 2025, Th2 (asthma)" = c("UGT2B28", "SYNJ1", "SORD", "PTGDR2", "LGALS7", "FBN1", "TIPARP",
    "KCNJ8", "TSPAN11"))
TEST <- setdiff(BLOOD, "GSE32472")
S32 <- S_blood |> filter(gse == "GSE32472")
for (k in names(E)) if ("TOMM70" %in% rownames(E[[k]]) && !"TOMM70A" %in% rownames(E[[k]])) rownames(E[[k]])[rownames(E[[k]]) == "TOMM70"] <- "TOMM70A"  # symbol alias
# Windhorst et al. 2023: derived in GSE225881; directions from the authors' PAM class scores (Supplementary Table 1; positive "Score BPD" = higher in BPD)
DIRS <- list("Windhorst 2023" = c(LAT = -1, LOC100505876 = 1, MRVI1 = 1, KLHDC4 = -1, ATP6V0E2 = -1,
  FAM102A = -1, FBXO3 = -1, OGDH = -1, SLC25A3 = -1, GAL = 1, CAT = 1,
  SOAT2 = 1, TNFRSF10C = 1, BCL2A1 = 1, DUSP2 = -1, PDLIM7 = 1, POLD2 = 1, CPLX1 = -1, SCML1 = -1,
  TMED1 = -1, HEBP1 = 1, RRN3 = -1, ZEB1 = -1, SLC2A11 = 1))
TESTS <- list("Windhorst 2023" = setdiff(BLOOD, "GSE225881"))
# Guo et al. 2025 (PLoS One): trained in GSE220135 (validated in GSE32472) -> directions from GSE220135; scored in the other three cohorts
gg <- c("IL7R", "CXCL10", "DEFA4", "PRTN3", "NCAPG", "CCNB1")
S35 <- S_blood |> filter(gse == "GSE220135")
DIRS[["Guo 2025"]] <- setNames(sign(S35$g[match(gg, S35$gene)]), gg)
TESTS[["Guo 2025"]] <- c("GSE32472", "GSE225881", "GSE125873")
R32 <- S32 |> arrange(p) |> filter(is.finite(g))
REF <- list("De novo top 5 (GSE32472)" = head(R32$gene, 5),
  "De novo top 50 (GSE32472)" = head(R32$gene, 50))
ALL <- c(PUB, REF)
grp <- c(rep("Published", 13), rep("Authors' prior", 2), rep("Reference", 2))
CELLS <- c("Neutrophils", "T cells")
AUD <- bind_rows(lapply(seq_along(ALL), function(i) {
  genes <- ALL[[i]]
  nm <- names(ALL)[i]
  dirs <- if (!is.null(DIRS[[nm]])) unname(DIRS[[nm]][genes]) else sign(S32$g[match(genes, S32$gene)])
  dirs[is.na(dirs)] <- 0
  tst <- if (!is.null(TESTS[[nm]])) TESTS[[nm]] else TEST
  per <- bind_rows(lapply(tst, function(k) {
    s <- sig_score(genes, dirs, E[[k]])
    ncov <- sum(genes %in% rownames(E[[k]]) & dirs != 0)
    if (ncov == 0) return(data.frame(cohort = k, auc = NA, var = NA, ncov = 0, g_adj = NA, v_adj = NA,
      rN = NA, rT = NA))
    st <- roc_stats(s, lab(k))
    ga <- if (k %in% GA_COH) ga_of(k) else NULL
    if (!is.null(ga)) {
      f <- lm(s ~ y + ga, data.frame(s = s, y = lab(k), ga = ga))
      ga_e <- c(coef(f)[["y"]], vcov(f)["y", "y"])
    } else ga_e <- c(NA, NA)
    cs <- CS[[k]]
    data.frame(cohort = k, auc = st[["auc"]], var = st[["var"]], ncov = ncov, g_adj = ga_e[1],
      v_adj = ga_e[2],
      rN = cor(s, cs[, "Neutrophils"]), rT = cor(s, cs[, "T cells"]))
  }))
  pa <- pool_auc(per$auc, per$var)
  pg <- if (sum(!is.na(per$g_adj))) re_pool(per$g_adj[!is.na(per$g_adj)],
    per$v_adj[!is.na(per$v_adj)]) else c(es = NA, lo = NA, hi = NA, p = NA)
  data.frame(signature = names(ALL)[i], group = grp[i], n_genes = length(genes),
    dir32 = paste(sum(dirs > 0), sum(dirs < 0), sep = "+/-"),
    derivation_cohort = setdiff(BLOOD, tst), cov = paste(per$ncov, collapse = "/"),
    auc_32 = per$auc[match("GSE32472", per$cohort)], auc_81 = per$auc[match("GSE225881", per$cohort)],
    auc_73 = per$auc[match("GSE125873", per$cohort)],
    auc_35 = per$auc[match("GSE220135", per$cohort)], auc = pa[1], lo = pa[2], hi = pa[3],
    ga_g = pg[["es"]], ga_lo = pg[["lo"]], ga_hi = pg[["hi"]], ga_p = pg[["p"]],
    rN = mean(per$rN, na.rm = TRUE), rT = mean(per$rT, na.rm = TRUE))
}))
print(AUD |> mutate(across(where(is.numeric), ~ round(.x, 2))), width = 250)
save(TOPN, GA_AUC, BEY, BEYP, AUD, RK, RK_all, NTOP, file = file.path(OUT, "fig5.RData"))
