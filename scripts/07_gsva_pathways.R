# Estimate pathway activity with a custom GSVA implementation in four blood cohorts.
# Use Gaussian-kernel CDFs, log odds, symmetric rank weights 
# Use the sum of maximum positive and negative deviations as the enrichment score.
# Test GO:BP sets of 15-500 measured genes via msigdbr.
# Pool cohort-level Hedges' g with DerSimonian-Laird random effects.
# Define conserved pathways by FDR < 0.05 and concordant direction in >= 3 measured cohorts.
source("./00_setup.R")
suppressPackageStartupMessages({
  library(msigdbr)
  library(data.table)
})
GO <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") |> select(gs_name,
  gene_symbol) |> distinct()
GOL <- split(GO$gene_symbol, GO$gs_name)

gsva_core <- function(x, sets, minsz = 15, maxsz = 500) {
  x <- x[apply(x, 1, sd) > 0, , drop = FALSE]
  p <- nrow(x)
  n <- ncol(x)
  h <- apply(x, 1, sd) / 4
  Z <- matrix(0, p, n)
  for (j in seq_len(n)) Z[,
    j] <- rowMeans(pnorm((x[, j] - x) / h))  # kernel CDF of sample j's value, per gene
  Z <- pmin(pmax(Z, 1e-6), 1 - 1e-6)
  Z <- log(Z / (1 - Z))
  idx <- lapply(sets, function(g) which(rownames(x) %in% g))
  sz <- lengths(idx)
  keep <- sz >= minsz & sz <= maxsz
  idx <- idx[keep]
  mem <- data.table(set = rep(seq_along(idx), lengths(idx)), gene = unlist(idx))
  m <- tabulate(mem$set, length(idx))
  ES <- matrix(NA_real_, length(idx), n, dimnames = list(names(idx), colnames(x)))
  for (j in seq_len(n)) {
    r <- integer(p)
    r[order(Z[, j], decreasing = TRUE)] <- seq_len(p)  # rank 1 = highest
    s <- abs(p / 2 - r)
    d <- copy(mem)
    d[, `:=`(pos = r[gene], w = s[gene])]
    setorder(d, set, pos)
    d[, `:=`(k = seq_len(.N), cw = cumsum(w) / sum(w), m = .N), by = set]
    d[, miss := (pos - k) / (p - m)]
    d[, `:=`(after = cw - miss, before = (cw - w / sum(w)) - miss), by = set]
    e <- d[, .(mx = max(0, after), mn = min(0, before)), by = set]
    ES[e$set, j] <- e$mx + e$mn
  }
  ES
}
GS <- lapply(setNames(BLOOD, BLOOD), function(k) {
  t0 <- Sys.time()
  r <- gsva_core(E[[k]], GOL)
  cat(k, nrow(r), "sets", round(as.numeric(Sys.time() - t0, units = "secs")), "s\n")
  r
})
PW <- bind_rows(lapply(BLOOD, function(k) {
  G <- GS[[k]]
  y <- lab(k)
  bind_rows(lapply(rownames(G), function(s) {
    h <- hedges(G[s, ], y)
    data.frame(set = s, gse = k, g = h[["g"]], v = h[["v"]])
  }))
}))
PM <- PW |> group_by(set) |> filter(n() >= 3) |> summarise(k = n(), r = list(re_pool(g, v)),
  r_hk = list(re_pool(g, v, hk = TRUE)),
  n_pos = sum(g > 0), n_neg = sum(g < 0), I2 = {
    w <- 1 / v
    mu <- sum(w * g) / sum(w)
    Q <- sum(w * (g - mu)^2)
    if (Q > 0) max(0, (Q - (n() - 1)) / Q) * 100 else 0
}) |>
  mutate(es = sapply(r, `[[`, "es"), lo = sapply(r, `[[`, "lo"), hi = sapply(r, `[[`, "hi"),
    p = sapply(r, `[[`, "p"), p_hk = sapply(r_hk, `[[`, "p")) |>
  select(-r, -r_hk) |> mutate(fdr = p.adjust(p, "BH"), same = ifelse(es > 0, n_pos, n_neg),
    conserved = fdr < .05 & same >= 3)
cat("Pathways tested:", nrow(PM), "| conserved:", sum(PM$conserved), "(lower",
  sum(PM$conserved & PM$es < 0), ", higher", sum(PM$conserved & PM$es > 0), ")\n")
print(PM |> filter(conserved) |> arrange(p) |> head(40) |> transmute(set = sub("GOBP_", "", set),
  g = round(es, 2), fdr = signif(fdr, 2), I2 = round(I2), same, k), n = 40)
save(GS, PW, PM, GOL, file = file.path(OUT, "pathways_gsva.RData"))
