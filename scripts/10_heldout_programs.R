# Rebuild signatures, pathways and programs within each leave-one-cohort-out training fold.
# Default signature filters: FDR < 0.05, |g| > 0.5, I2 < 75%, nonzero CI and all-but-one concordance.
# Conserved pathways require FDR < 0.05 and all-but-one directional concordance in training cohorts.
# Cluster leading genes at Jaccard distance 0.9; require >= 3 pathways and >= 8 genes.
# Name clusters by same-direction overlap with final programs at Jaccard >= 0.05; retain up to 30 genes.
source("./00_setup.R")
load(file.path(OUT, "pathways_gsva.RData"))
load(file.path(OUT, "programs4.RData"))
P4dir <- c("RNA processing and translation" = -1, "T cell activation" = -1, "Cell cycle" = 1,
  "Neutrophil antimicrobial" = 1)
FINAL <- setNames(lapply(P4, function(pg) PR4$genes[[which(PR4$program == pg)]]), P4)

derive <- function(tr, relaxed = FALSE) {
  m <- length(tr)
  Ms <- meta_all(S_blood |> filter(gse %in% tr), kmin = m - 1) |> mutate(lo = es - 1.96 * se,
    hi = es + 1.96 * se) |>
    filter(if (relaxed) p < .05 else (fdr < .05 & abs(es) > .5), I2 < 75, lo > 0 | hi < 0,
      ifelse(es > 0, n_pos, n_neg) >= m - 1)
  pw <- PW |> filter(gse %in% tr) |> group_by(set) |> filter(n() >= m - 1) |>
    summarise(r = list(re_pool(g, v)), n_pos = sum(g > 0), n_neg = sum(g < 0)) |>
    mutate(es = sapply(r, `[[`, "es"), p = sapply(r, `[[`, "p"), fdr = p.adjust(p, "BH"),
      same = ifelse(es > 0, n_pos, n_neg)) |>
    filter(if (relaxed) p < .05 else fdr < .05, same >= m - 1) |> mutate(dir = sign(es))
  pw$lead <- lapply(seq_len(nrow(pw)),
    function(i) Ms$gene[Ms$gene %in% GOL[[pw$set[i]]] & sign(Ms$es) == pw$dir[i]])
  pw <- pw |> filter(lengths(lead) >= 3)
  cores <- setNames(vector("list", 4), P4)
  for (d in c(-1, 1)) {
    i <- which(pw$dir == d)
    if (length(i) < 3) next
    L <- pw$lead[i]
    n <- length(L)
    J <- matrix(0, n, n)
    for (a in 1:n) for (b in 1:n) J[a,
      b] <- length(intersect(L[[a]], L[[b]])) / length(union(L[[a]], L[[b]]))
    ct <- cutree(hclust(as.dist(1 - J), "average"), h = .9)
    for (cc in unique(ct)) {
      ii <- i[ct == cc]
      g <- unique(unlist(pw$lead[ii]))
      if (length(ii) < 3 || length(g) < 8) next
      cand <- P4[P4dir == d]
      jj <- sapply(cand, function(pg) length(intersect(g, FINAL[[pg]])) / length(union(g, FINAL[[pg]])))
      if (max(jj) >= .05) {
        pg <- cand[which.max(jj)]
        cores[[pg]] <- unique(c(cores[[pg]], g))
      }
    }
  }
  lapply(cores,
    function(g) if (is.null(g)) NULL else Ms |> filter(gene %in% g) |> arrange(p) |> head(30) |> transmute(gene, es, p))
}
RELAX <- as.logical(Sys.getenv("RELAX", "FALSE"))
FOLD <- lapply(setNames(BLOOD, BLOOD), function(h) derive(setdiff(BLOOD, h), relaxed = RELAX))
save(FOLD, file = file.path(OUT, if (RELAX) "heldout4_relaxed.RData" else "heldout4_strict.RData"))
DISC <- setNames(lapply(P4, function(pg) PR4$core[[which(PR4$program == pg)]] |> transmute(gene, es, p)),
  P4)
for (h in BLOOD) cat(h, "held out -> core genes per program:",
  paste(sapply(FOLD[[h]], function(x) if (is.null(x)) 0 else nrow(x)), collapse = "/"),
  "| overlap with final core:",
  paste(sapply(P4,
    function(pg) if (is.null(FOLD[[h]][[pg]])) 0 else length(intersect(FOLD[[h]][[pg]]$gene,
      DISC[[pg]]$gene))), collapse = "/"), "\n")

psc <- function(core, x) {
  if (is.null(core)) return(rep(NA_real_, ncol(x)))
  g <- core$gene[core$gene %in% rownames(x)]
  if (!length(g)) return(rep(NA_real_, ncol(x)))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  colMeans(z * sign(core$es[match(g, core$gene)]))
}
SC4 <- c(P4, "Combined score")

score <- function(nm, k, type) {
  cs <- if (type == "held-out") FOLD[[k]] else DISC
  if (nm == "Combined score") return(rowMeans(sapply(P4,
    function(pg) as.numeric(scale(psc(cs[[pg]], E[[k]])))), na.rm = TRUE))
  psc(cs[[nm]], E[[k]])
}

# ---- AUC: discovery vs held-out ----
AUCT <- bind_rows(lapply(SC4, function(nm) bind_rows(lapply(c("discovery", "held-out"), function(tp) {
  r <- t(sapply(BLOOD, function(k) roc_stats(score(nm, k, tp), lab(k))))
  pa <- pool_auc(r[, "auc"], r[, "var"])
  data.frame(score = nm, type = tp, auc = pa[1], lo = pa[2], hi = pa[3], I2 = pa[4],
    cohort_auc = paste(round(r[, "auc"], 2), collapse = "/"))
}))))
print(AUCT |> mutate(across(where(is.numeric), ~ round(.x, 2))))
if (!RELAX) quit(save = "no")

# ---- GA-matched strata ----
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
GMT <- bind_rows(lapply(SC4, function(nm) bind_rows(lapply(c("discovery", "held-out"), function(tp) {
  r <- t(sapply(GA_COH, function(k) {
    s <- as.numeric(scale(score(nm, k, tp)))
    y <- lab(k)
    ga <- PH[[k]]$ga_weeks
    fa <- lm(s ~ y)
    fb <- lm(s ~ y + ga)
    c(strat(s, y, ga), coef(fa)[["y"]], vcov(fa)["y", "y"], coef(fb)[["y"]], vcov(fb)["y", "y"])
  }))
  S <- re_pool(r[, "d"], r[, "v"])
  A <- re_pool(r[, 5], r[, 6])
  B <- re_pool(r[, 7], r[, 8])
  data.frame(score = nm, type = tp, strata = S[["es"]], s_lo = S[["lo"]], s_hi = S[["hi"]],
    s_p = S[["p"]], n1 = sum(r[, "n1"]), n0 = sum(r[, "n0"]),
    A = A[["es"]], B = B[["es"]], B_lo = B[["lo"]], B_hi = B[["hi"]], B_p = B[["p"]])
}))))
print(GMT |> mutate(across(where(is.numeric), ~ signif(.x, 3))))
save(FOLD, DISC, AUCT, GMT, P4, SC4,
  file = file.path(OUT, if (RELAX) "heldout4_relaxed.RData" else "heldout4.RData"))
