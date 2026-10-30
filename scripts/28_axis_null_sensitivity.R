# Assess random gene-set nulls and sensitivity to shared genes.
# Null A uses random signs; null B uses derivation-cohort signs; each uses 1,000 draws.
# Remove balance genes and genes shared by >= 2 signatures for the sensitivity analysis.
source("./source_section.R")
source_section("./23_newfigs_data.R", "Signature scoring setup", "Signature association analysis")
set.seed(20261030)

derv <- function(nm) if (nm == "Guo 2025") "GSE220135" else if (nm == "Windhorst 2023") "GSE225881" else "GSE32472"

tests_of <- function(nm) if (!is.null(TESTS[[nm]])) TESTS[[nm]] else TEST

pcstats <- function(M, b) {
  M <- M[, colSums(is.na(M)) == 0 & apply(M, 2, sd) > 0, drop = FALSE]
  pc <- prcomp(M, scale. = TRUE)
  c(n_sigs = ncol(M), PC1_var = summary(pc)$importance[2, 1], PC2_var = summary(pc)$importance[2, 2],
    r_PC1_bal = abs(cor(pc$x[, 1], b)),
    r_PC2_bal = abs(cor(pc$x[, 2], b)),
    R2_bal_on_PC1_PC2 = summary(lm(b ~ pc$x[, 1] + pc$x[, 2]))$r.squared,
    mean_abs_r_bal = mean(abs(cor(M, b))))
}
PCT <- c("GSE225881", "GSE125873", "GSE220135")

# ---- observed (identical to 28_axis_pca.R) and the signature set / sizes per cohort ----
OBS <- list()
SIZES <- list()
for (k in PCT) {
  nms <- names(PUB)[sapply(names(PUB), function(nm) k %in% tests_of(nm))]
  M <- sapply(nms, function(nm) sig(PUB[[nm]], E[[k]], nm))
  keep <- colSums(is.na(M)) == 0
  nms <- nms[keep]
  OBS[[k]] <- pcstats(M[, keep, drop = FALSE], BAL(E[[k]]))
  SIZES[[k]] <- sapply(nms, function(nm) {
    g <- PUB[[nm]]
    d <- if (!is.null(DIRS[[nm]])) unname(DIRS[[nm]][g]) else sign(S32$g[match(g, S32$gene)])
    sum(g %in% rownames(E[[k]]) & !is.na(d) & d != 0)
  })
}
print(do.call(rbind, OBS))

# ---- C1 nulls ----
U <- Reduce(intersect, lapply(BLOOD, function(k) rownames(E[[k]])))
GD <- lapply(setNames(c("GSE32472", "GSE220135", "GSE225881"), c("GSE32472", "GSE220135", "GSE225881")),
  function(k) {
    s <- S_blood |> filter(gse == k)
    setNames(sign(s$g), s$gene)
})
U <- U[U %in% names(GD$GSE32472) & U %in% names(GD$GSE220135) & U %in% names(GD$GSE225881)]
Z <- lapply(setNames(PCT, PCT), function(k) {
  z <- t(scale(t(E[[k]][U, ])))
  z[!is.finite(z)] <- 0
  z
})
cat("Universe:", length(U), "genes\n")
B <- 1000
NULLS <- bind_rows(lapply(c("A", "B"),
  function(type) bind_rows(lapply(seq_len(B), function(it) bind_rows(lapply(PCT, function(k) {
    sz <- SIZES[[k]]
    M <- sapply(names(sz), function(nm) {
      g <- sample(U, sz[[nm]])
      d <- if (type == "A") sample(c(-1, 1), length(g), replace = TRUE) else GD[[derv(nm)]][g]
      as.numeric(scale(colMeans(Z[[k]][g, , drop = FALSE] * d)))
    })
    data.frame(type = type, it = it, cohort = k, t(pcstats(M, BAL(E[[k]]))))
}))))))
C1 <- bind_rows(lapply(c("A", "B"), function(type) bind_rows(lapply(PCT, function(k) {
  nl <- NULLS |> filter(type == !!type, cohort == k)
  o <- OBS[[k]]
  data.frame(null = c(A = "A: random genes, random signs",
    B = "B: random genes, derivation-cohort signs")[type], cohort = k,
    obs_PC1_var = o[["PC1_var"]], null_PC1_var_median = median(nl$PC1_var),
    null_PC1_var_p95 = quantile(nl$PC1_var, .95),
    p_PC1_var = (1 + sum(nl$PC1_var >= o[["PC1_var"]])) / (B + 1),
    obs_r_PC1_bal = o[["r_PC1_bal"]], null_r_PC1_bal_median = median(nl$r_PC1_bal),
    null_r_PC1_bal_p95 = quantile(nl$r_PC1_bal, .95),
    p_r_PC1_bal = (1 + sum(nl$r_PC1_bal >= o[["r_PC1_bal"]])) / (B + 1),
    obs_mean_abs_r_bal = o[["mean_abs_r_bal"]], null_mean_abs_r_bal_median = median(nl$mean_abs_r_bal),
    null_mean_abs_r_bal_p95 = quantile(nl$mean_abs_r_bal, .95),
    p_mean_abs_r_bal = (1 + sum(nl$mean_abs_r_bal >= o[["mean_abs_r_bal"]])) / (B + 1))
}))))
rownames(C1) <- NULL
print(C1 |> mutate(across(where(is.numeric), ~ round(.x, 3))), width = 250)

# ---- C2 sensitivity ----
tab <- table(unlist(lapply(PUB, unique)))
SHARED <- names(tab)[tab >= 2]
DROP <- union(SHARED, c(NEU, PAN_T))
cat("Removed genes:", paste(sort(intersect(DROP, unlist(PUB))), collapse = ", "), "\n")
source_section("./23_newfigs_data.R", "Signature association analysis", "Methylation comparison")
PUB0 <- PUB
SIGMAP0 <- SIGMAP
PUB <- lapply(PUB0, function(g) setdiff(g, DROP))
PUB <- PUB[sapply(PUB, length) > 0]
cat("Signatures changed:",
  paste(names(PUB0)[sapply(names(PUB0), function(n) !identical(PUB0[[n]], PUB[[n]]))], collapse = "; "),
  "\n")
source_section("./23_newfigs_data.R", "Signature association analysis", "Methylation comparison")
SIGMAP2 <- SIGMAP
SIGMAP <- SIGMAP0
C2PCA <- bind_rows(lapply(PCT, function(k) {
  nms <- names(PUB)[sapply(names(PUB), function(nm) k %in% tests_of(nm))]
  M <- sapply(nms, function(nm) sig(PUB[[nm]], E[[k]], nm))
  data.frame(cohort = k, t(pcstats(M, BAL(E[[k]]))))
}))

axis_of <- function(n, t) ifelse(pmax(abs(n), abs(t)) < .5, "Weak / neither",
  ifelse(abs(n) >= abs(t), "Neutrophil-associated", "T-cell-associated"))
C2SIG <- SIGMAP0 |> transmute(signature, n_genes_before = sapply(PUB0[signature], length),
  axis_before = axis_of(r_neu, r_T), r_neu_before = r_neu, r_T_before = r_T,
  b_GAbal_before = b_GAbal, p_GAbal_before = p_GAbal) |>
  left_join(SIGMAP2 |> transmute(signature, n_genes_after = sapply(PUB[signature], length),
    axis_after = axis_of(r_neu, r_T), r_neu_after = r_neu, r_T_after = r_T,
    b_GA_after = b_GA, p_GA_after = p_GA, b_GAbal_after = b_GAbal, p_GAbal_after = p_GAbal),
    by = "signature")
print(C2PCA |> mutate(across(where(is.numeric), ~ round(.x, 3))))
print(C2SIG |> mutate(across(where(is.numeric), ~ round(.x, 3))), width = 250)
save(OBS, NULLS, C1, C2PCA, C2SIG, SHARED, DROP, U, file = file.path(OUT, "axis_null_sensitivity.RData"))
dir.create(file.path(OUT, "supp_tables"), showWarnings = FALSE)
write.csv(C1 |> mutate(across(where(is.numeric), ~ round(.x, 4))),
  file.path(OUT, "supp_tables", "S14e_random_gene_set_nulls.csv"), row.names = FALSE)
write.csv(C2PCA |> mutate(across(where(is.numeric), ~ round(.x, 3))),
  file.path(OUT, "supp_tables", "S14f_shared_gene_sensitivity_PCA.csv"), row.names = FALSE)
write.csv(C2SIG |> mutate(across(where(is.numeric), ~ round(.x, 4))),
  file.path(OUT, "supp_tables", "S14g_shared_gene_sensitivity_signatures.csv"), row.names = FALSE)
