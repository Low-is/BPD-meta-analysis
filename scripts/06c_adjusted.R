# Like-for-like: same 3 GA cohorts, same criteria, unadjusted vs GA-adjusted selection
source("./00_setup.R")
S_un  <- S_blood |> filter(gse %in% GA_COH)
S_adj <- bind_rows(lapply(GA_COH,
  function(k) cohort_stats(k, lab(k), covars = data.frame(ga = PH[[k]]$ga_weeks))))

pick <- function(S, ndir) meta_all(S, kmin = 3) |> mutate(lo = es - 1.96 * se, hi = es + 1.96 * se) |>
  filter(fdr < .05, abs(es) > .5, I2 < 75, lo > 0 | hi < 0, ifelse(es > 0, n_pos, n_neg) >= ndir)
for (nd in c(3, 2)) {
  U <- pick(S_un, nd)
  A <- pick(S_adj, nd)
  cat(sprintf("Same direction in >= %d of 3 cohorts: unadjusted %d (up %d, down %d) | GA-adjusted %d (up %d, down %d) | overlap %d\n",
    nd, nrow(U), sum(U$es > 0), sum(U$es < 0), nrow(A), sum(A$es > 0), sum(A$es < 0),
    length(intersect(U$gene, A$gene))))
}
load(file.path(OUT, "signature_rebuild.RData"))
A3 <- pick(S_adj, 3)
cat("Of the 1,151 (4-cohort unadjusted):", sum(SIGR$gene %in% A3$gene),
  "also pass full criteria GA-adjusted (3 cohorts, 3/3)\n")
MA <- meta_all(S_adj, kmin = 3) |> mutate(lo = es - 1.96 * se, hi = es + 1.96 * se)
MU <- meta_all(S_un, kmin = 3) |> mutate(lo = es - 1.96 * se, hi = es + 1.96 * se)
for (nm in c("unadjusted", "GA-adjusted")) {
  X <- if (nm == "unadjusted") MU else MA
  cat(sprintf("%-11s | CI excl 0: %5d | FDR<0.05: %5d | |g|>0.5: %5d | FDR & |g|>0.5: %4d | min FDR %.2g\n", nm,
    sum(X$lo > 0 | X$hi < 0), sum(X$fdr < .05), sum(abs(X$es) > .5 & (X$lo > 0 | X$hi < 0)),
    sum(X$fdr < .05 & abs(X$es) > .5), min(X$fdr)))
}
print(A3 |> select(gene, es, lo, hi, I2, fdr))
cat("Genes tested (3 GA cohorts):", nrow(MA), "| expected nominal p<0.05 by chance:",
  round(.05 * nrow(MA)), "\n")
