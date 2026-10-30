# Apply signature-specific direction overrides and test-cohort exclusions.
source("./source_section.R")
source_section("./23_newfigs_data.R", "Signature scoring setup", "Signature association analysis")

r2 <- function(s, Z) summary(lm(s ~ ., data.frame(s = s, Z)))$r.squared
PCA <- bind_rows(lapply(BLOOD, function(k) {
  nms <- names(PUB)[sapply(names(PUB),
    function(nm) k %in% (if (!is.null(TESTS[[nm]])) TESTS[[nm]] else TEST))]
  M <- sapply(nms, function(nm) sig(PUB[[nm]], E[[k]], nm))
  M <- M[, colSums(is.na(M)) == 0, drop = FALSE]
  if (ncol(M) < 3) return(data.frame(cohort = k, n_samples = nrow(M), n_sigs = ncol(M)))
  pc <- prcomp(M, scale. = TRUE)
  imp <- summary(pc)$importance
  b <- BAL(E[[k]])
  N <- CS[[k]][, "Neutrophils"]
  Tc <- CS[[k]][, "T cells"]
  p1 <- pc$x[, 1] * sign(cor(pc$x[, 1], b))
  p2 <- pc$x[, 2] * sign(cor(pc$x[, 2], b))
  y <- lab(k)
  ga <- if (k %in% GA_COH) ga_of(k) else NULL
  data.frame(cohort = k, n_samples = nrow(M), n_sigs = ncol(M), PC1_var = imp[2, 1], PC2_var = imp[2, 2],
    PC1_2_cum = imp[3, 2],
    r_PC1_bal = cor(p1, b), r_PC2_bal = cor(p2, b), r_PC1_neu = cor(p1, N), r_PC1_T = cor(p1, Tc),
    r_PC2_neu = cor(p2, N), r_PC2_T = cor(p2, Tc),
    R2_bal_on_PC12 = r2(b, data.frame(p1, p2)),
    meanR2_sigs_on_bal = mean(apply(M, 2, r2, Z = data.frame(b))),
    medR2_sigs_on_bal = median(apply(M, 2, r2, Z = data.frame(b))),
    meanR2_sigs_on_neu_T = mean(apply(M, 2, r2, Z = data.frame(N, Tc))),
    meanR2_sigs_on_GA = if (!is.null(ga)) mean(apply(M, 2, r2, Z = data.frame(ga))) else NA,
    r_PC1_GA = if (!is.null(ga)) cor(p1, ga) else NA, AUC_PC1 = A(y, p1))
}))
print(PCA |> mutate(across(where(is.numeric), ~ round(.x, 3))), width = 250)
# per-signature R2 on balance and on (neutrophil, T) markers, pooled (mean over test cohorts)
SIGR2 <- bind_rows(lapply(names(PUB), function(nm) {
  tst <- if (!is.null(TESTS[[nm]])) TESTS[[nm]] else TEST
  bind_rows(lapply(tst, function(k) {
    s <- sig(PUB[[nm]], E[[k]], nm)
    if (all(is.na(s))) return(NULL)
    data.frame(signature = nm, cohort = k, R2_bal = r2(s, data.frame(BAL(E[[k]]))),
      R2_neuT = r2(s, data.frame(CS[[k]][, c("Neutrophils", "T cells")])))
  }))
})) |>
  group_by(signature) |> summarise(across(c(R2_bal, R2_neuT), mean), n_cohorts = n())
print(SIGR2 |> mutate(across(where(is.numeric), ~ round(.x, 2))), n = 30)
save(PCA, SIGR2, file = file.path(OUT, "axis_pca.RData"))
# export for Supplementary Table S14 (b, c)
dir.create(file.path(OUT, "supp_tables"), showWarnings = FALSE)
write.csv(PCA |> filter(!is.na(PC1_var)) |> mutate(across(where(is.numeric), ~ round(.x, 3))),
  file.path(OUT, "supp_tables", "S14b_signature_PCA_and_balance.csv"), row.names = FALSE)
write.csv(SIGR2 |> mutate(across(where(is.numeric), ~ round(.x, 3))),
  file.path(OUT, "supp_tables", "S14c_signature_variance_explained.csv"), row.names = FALSE)
