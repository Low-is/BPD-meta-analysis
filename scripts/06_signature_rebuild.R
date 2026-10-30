# Rebuild the blood signature in four independent cohorts.
# Require FDR < 0.05, |pooled g| > 0.5, I2 < 75%, a nonzero 95% CI and concordance in >= 3 cohorts.
# Include genes measured in >= 3 cohorts; assess leave-one-cohort-out and GA-adjusted selection.
source("./00_setup.R")
M <- meta_all(S_blood, kmin = 3) |>
  mutate(lo = es - 1.96 * se, hi = es + 1.96 * se, n_same = pmax(n_pos, n_neg),
    dir_ok = ifelse(es > 0, n_pos, n_neg) >= 3,
    ci_ok = lo > 0 | hi < 0)

crit <- function(M, gthr = .5) M |> filter(fdr < .05, abs(es) > gthr, I2 < 75, ci_ok, dir_ok)
SIGR <- crit(M)
cat("Genes tested (measured in >=3 cohorts):", nrow(M), "\n")
cat("Pass all criteria (|g| > 0.5):", nrow(SIGR), " up:", sum(SIGR$es > 0), " down:", sum(SIGR$es < 0),
  "\n")
for (t in c(.6, .7, .8, .9, 1.0)) cat(sprintf("  |g| > %.1f: %d genes (up %d, down %d)\n", t,
  nrow(crit(M, t)), sum(crit(M, t)$es > 0), sum(crit(M, t)$es < 0)))
# leave-one-cohort-out survival
loo <- sapply(BLOOD, function(d) {
  Ml <- meta_all(S_blood |> filter(gse != d), kmin = 2) |> mutate(lo = es - 1.96 * se,
    hi = es + 1.96 * se)
  Ml <- Ml[match(SIGR$gene, Ml$gene), ]
  !is.na(Ml$es) & sign(Ml$es) == sign(SIGR$es) & (Ml$lo > 0 | Ml$hi < 0)
})
SIGR$loo_all <- rowSums(loo) == length(BLOOD)
cat("Survive omission of every cohort (same sign, CI excludes 0):", sum(SIGR$loo_all), "\n")
# GA-adjusted (3 cohorts)
S_adj <- bind_rows(lapply(GA_COH,
  function(k) cohort_stats(k, lab(k), covars = data.frame(ga = PH[[k]]$ga_weeks))))
MA <- meta_all(S_adj, kmin = 2) |> mutate(lo = es - 1.96 * se, hi = es + 1.96 * se)
MA <- MA[match(SIGR$gene, MA$gene), ]
SIGR$ga_es <- MA$es
SIGR$ga_sig <- !is.na(MA$es) & sign(MA$es) == sign(SIGR$es) & (MA$lo > 0 | MA$hi < 0)
cat("Same sign and CI excludes 0 after GA adjustment:", sum(SIGR$ga_sig, na.rm = TRUE),
  " | median |g| unadj", round(median(abs(SIGR$es)), 2),
  "-> GA-adj", round(median(abs(SIGR$ga_es), na.rm = TRUE), 2), "\n")
# original 21
o21 <- M[match(SIG$gene, M$gene), ] |> mutate(pass = gene %in% SIGR$gene)
cat("Original 21 genes passing rebuilt criteria:", sum(o21$pass), ":",
  paste(o21$gene[o21$pass], collapse = ", "), "\n")
# top genes
print(SIGR |> arrange(p) |> head(25) |> transmute(gene, g = round(es, 2), lo = round(lo, 2),
  hi = round(hi, 2), I2 = round(I2), k, n_same, fdr = signif(fdr, 2), loo_all, ga_g = round(ga_es, 2),
  ga_sig))
save(M, SIGR, file = file.path(OUT, "signature_rebuild.RData"))
