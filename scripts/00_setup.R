# Shared data loading, statistical functions and program scoring.
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(pROC)
  library(limma)
})
select <- dplyr::select
filter <- dplyr::filter
options(dplyr.summarise.inform = FALSE, warn = -1)

DATA <- "./data"
OUT <- "./out"
dir.create(OUT, showWarnings = FALSE)

BLOOD <- c("GSE32472", "GSE225881", "GSE125873",
  "GSE220135")  # 4 independent cohorts (GSE188944 excluded)
GA_COH <- c("GSE32472", "GSE225881", "GSE220135")
D <- lapply(setNames(c(BLOOD, "GSE302178"), c(BLOOD, "GSE302178")),
  function(k) readRDS(file.path(DATA, paste0(k, "_analysis_ready.rds"))))
E <- lapply(D, `[[`, "expression")

lab <- function(k) as.integer(D[[k]]$pheno$bpd)
PH <- lapply(D, `[[`, "pheno")
MODEL5 <- readRDS(file.path(DATA, "model_and_meta_analysis.rds"))
ANCHOR <- MODEL5$anchor_genes  # GO anchors (+ signature genes), fixed before analysis
SIG <- MODEL5$signature_21
PROGS <- names(ANCHOR)
BPD_PROG <- setdiff(PROGS, "Proliferation")

# ---- statistics ----
cohort_stats <- function(k, y, covars = NULL, m = E[[k]]) {
  keep <- !is.na(y)
  m <- m[, keep, drop = FALSE]
  y <- y[keep]
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  J <- 1 - 3 / (4 * (n1 + n0) - 9)
  if (is.null(covars)) {
    v1 <- apply(m[, y == 1, drop = FALSE], 1, var)
    v0 <- apply(m[, y == 0, drop = FALSE], 1, var)
    sp <- sqrt(((n1 - 1) * v1 + (n0 - 1) * v0) / (n1 + n0 - 2))
    g <- (rowMeans(m[, y == 1, drop = FALSE]) - rowMeans(m[, y == 0, drop = FALSE])) / sp * J
  } else {
    X <- model.matrix(~ y + ., data = data.frame(y = y, covars[keep, , drop = FALSE]))
    g <- eBayes(lmFit(m, X))$t[, "y"] * sqrt(1 / n1 + 1 / n0) * J
  }
  v <- (n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0))
  data.frame(gene = rownames(m), gse = k, n1 = n1, n0 = n0, g = g, v = v,
    p = 2 * pnorm(-abs(g / sqrt(v))))
}

meta_all <- function(S, kmin = 2, hk = FALSE) {
  r <- S |> filter(is.finite(g), is.finite(v)) |> group_by(gene) |> filter(n() >= kmin) |>
    summarise(k = n(), w_sum = sum(1 / v), mu_f = sum(g / v) / w_sum, Q = sum((g - mu_f)^2 / v),
      tau2 = pmax(0, (Q - (k - 1)) / (w_sum - sum(1 / v^2) / w_sum)),
      es = sum(g / (v + tau2)) / sum(1 / (v + tau2)), se = sqrt(1 / sum(1 / (v + tau2))),
      hk_q = sum((g - es)^2 / (v + tau2)) / pmax(k - 1, 1),
      n_pos = sum(g > 0), n_neg = sum(g < 0)) |>
    mutate(p = 2 * pnorm(-abs(es / se)),
      se_hk = se * sqrt(pmax(1, hk_q)), p_hk = 2 * pt(-abs(es / se_hk), pmax(k - 1, 1)),
      I2 = ifelse(Q > 0, pmax(0, (Q - (k - 1)) / Q) * 100, 0)) |> select(-w_sum, -mu_f)
  if (hk) r <- r |> mutate(p = p_hk)
  r |> mutate(fdr = p.adjust(p, "BH"))
}
NA_ROC <- c(auc = NA, lo = NA, hi = NA, var = NA, n1 = NA, n0 = NA)

roc_stats <- function(score, y) {
  keep <- !is.na(y) & !is.na(score)
  if (!sum(keep) || sd(score[keep]) == 0) return(NA_ROC)
  r <- roc(y[keep], score[keep], quiet = TRUE, direction = "<", levels = c(0, 1))
  ci <- ci.auc(r, method = "delong")
  c(auc = as.numeric(auc(r)), lo = ci[1], hi = min(1, ci[3]), var = max(pROC::var(r), 1e-4),
    n1 = sum(y[keep] == 1), n0 = sum(y[keep] == 0))
}

pool_auc <- function(a, v) {
  ok <- !is.na(a) & !is.na(v)
  a <- a[ok]
  v <- v[ok]
  w <- 1 / v
  mu <- sum(w * a) / sum(w)
  Q <- sum(w * (a - mu)^2)
  k <- length(a)
  t2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  ws <- 1 / (v + t2)
  est <- sum(ws * a) / sum(ws)
  se <- sqrt(1 / sum(ws))
  c(auc = est, lo = est - 1.96 * se, hi = min(1, est + 1.96 * se),
    I2 = if (Q > 0) max(0, (Q - (k - 1)) / Q) * 100 else 0)
}

re_pool <- function(g, v, hk = FALSE) {
  ok <- is.finite(g) & is.finite(v)
  g <- g[ok]
  v <- v[ok]
  w <- 1 / v
  mu <- sum(w * g) / sum(w)
  Q <- sum(w * (g - mu)^2)
  k <- length(g)
  t2 <- max(0, (Q - (k - 1)) / (sum(w) - sum(w^2) / sum(w)))
  ws <- 1 / (v + t2)
  e <- sum(ws * g) / sum(ws)
  se <- sqrt(1 / sum(ws))
  if (hk && k > 1) {
    q <- sum(ws * (g - e)^2) / (k - 1)
    se2 <- se * sqrt(max(1, q))
    crit <- qt(.975, k - 1)
    return(c(es = e, se = se2, lo = e - crit * se2, hi = e + crit * se2,
      p = 2 * pt(-abs(e / se2), k - 1), k = k))
  }
  c(es = e, se = se, lo = e - 1.96 * se, hi = e + 1.96 * se, p = 2 * pnorm(-abs(e / se)), k = k)
}

hedges <- function(s, y) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  sp <- sqrt(((n1 - 1) * var(s[y == 1]) + (n0 - 1) * var(s[y == 0])) / (n1 + n0 - 2))
  g <- (mean(s[y == 1]) - mean(s[y == 0])) / sp * (1 - 3 / (4 * (n1 + n0) - 9))
  c(g = g, v = (n1 + n0) / (n1 * n0) + g^2 / (2 * (n1 + n0)))
}

# ---- programs ----
S_blood <- bind_rows(lapply(BLOOD, function(k) cohort_stats(k, lab(k))))

core_genes <- function(tr, prog, S = S_blood) {
  m <- length(tr)
  S <- S |> filter(gse %in% tr, gene %in% ANCHOR[[prog]])
  M <- meta_all(S, kmin = m - 1) |> filter(p < .05) |> filter(ifelse(es > 0, n_pos, n_neg) >= m - 1)
  if (!nrow(M)) return(data.frame(gene = character(), dir = numeric(), es = numeric(), p = numeric()))
  ok <- rep(TRUE, nrow(M))
  for (d in tr) {
    L <- meta_all(S |> filter(gse != d, gene %in% M$gene), kmin = 1)
    L <- L[match(M$gene, L$gene), ]
    ok <- ok & !is.na(L$es) & sign(L$es) == sign(M$es)
  }
  M[ok, ] |> arrange(p) |> head(30) |> transmute(gene, dir = sign(es), es, p)
}

prog_score <- function(core, x) {
  g <- core$gene[core$gene %in% rownames(x)]
  if (!length(g)) return(rep(NA_real_, ncol(x)))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  setNames(colMeans(z * core$dir[match(g, core$gene)]), colnames(x))
}

build_model <- function(tr) {
  cores <- lapply(setNames(PROGS, PROGS), function(p) core_genes(tr, p))
  tra <- bind_rows(lapply(PROGS, function(p) bind_rows(lapply(tr, function(k) {
    st <- roc_stats(prog_score(cores[[p]], E[[k]]), lab(k))
    data.frame(program = p, cohort = k, auc = st[["auc"]], var = st[["var"]])
  }))))
  ent <- tra |> filter(!is.na(auc)) |> group_by(program) |> summarise(pa = list(pool_auc(auc, var))) |>
    mutate(auc = sapply(pa, `[`, 1), lo = sapply(pa, `[`, 2)) |> select(-pa)
  list(cores = cores, entry = ent,
    progs = ent |> filter(program %in% BPD_PROG, lo > .5) |> pull(program))
}
SCORES <- c("BPD tracking index", PROGS)

score_fun <- function(model, nm, x) if (nm == "BPD tracking index")
rowMeans(sapply(model$progs, function(p) prog_score(model$cores[[p]], x)),
  na.rm = TRUE) else prog_score(model$cores[[nm]], x)
