# Numbers used in the manuscript text + external cohort (GSE302178) check of the v7 scores
source("./00_setup.R")
for (f in c("signature_rebuild", "programs4", "fig3", "fig4", "fig5", "v7_partA", "v7_partB", "newfigs",
  "timecourse")) load(file.path(OUT, paste0(f, ".RData")))

A <- function(y, s) as.numeric(auc(roc(y, s, quiet = TRUE, direction = "<", levels = c(0, 1))))
cat("== signature ==\n")
print(dim(SIGR))
print(table(sign(SIGR$es)))
cat("== programs ==\n")
print(PR4[, c("program", "n_path")])
print(sapply(PR4$core, nrow))
cat("== AUCT ==\n")
print(AUCT |> mutate(across(where(is.numeric), ~ round(.x, 2))))
cat("== DEVT ==\n")
print(DEVT |> select(score, type, k, strata, s_lo, s_hi, s_p, B, B_p, C,
  C_p) |> mutate(across(where(is.numeric), ~ signif(.x, 2))))

psc <- function(core, x) {
  g <- core$gene[core$gene %in% rownames(x)]
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  colMeans(z * sign(core$es[match(g, core$gene)]))
}
cat("== in-sample pooled g of program scores ==\n")
for (i in seq_len(nrow(PR4))) {
  r <- t(sapply(BLOOD, function(k) hedges(psc(PR4$core[[i]], E[[k]]), lab(k))))
  p <- re_pool(r[, 1], r[, 2])
  cat(PR4$program[i], round(p, 3), "\n")
}
cat("== longitudinal programs ==\n")
print(MM |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
print(FI)
cat("== GA AUC ==\n")
print(GA_AUC)
cat("== cell CORR (pooled r) ==\n")
print(CORR |> tidyr::pivot_wider(names_from = cell,
  values_from = r) |> mutate(across(where(is.numeric), ~ round(.x, 2))))

# ---- external GSE302178 ----
x <- E$GSE302178
y <- lab("GSE302178")
cat("\n== GSE302178: n", length(y), "BPD", sum(y), " genes", nrow(x), "\n")
NAIVE <- c("CCR7", "TCF7", "LEF1", "BCL11B", "CD27", "CD28", "IL7R", "BACH2", "SATB1")
PAN_T <- c("CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK")
NEU <- c("CSF3R", "CXCR2", "FPR1", "FCGR3B", "MME", "CMTM2")

mz <- function(g, x) {
  g <- intersect(g, rownames(x))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z)))
}

lmn <- function(g, x) colMeans(x[intersect(g, rownames(x)), , drop = FALSE])

modscore <- function(M, x) {
  g <- M$gene[M$gene %in% rownames(x)]
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z * sign(M$t[match(g, M$gene)]))))
}
ext <- c(naive_low = A(y, -mz(NAIVE, x)), panT_low = A(y, -mz(PAN_T, x)),
  balance = A(y, as.numeric(lmn(NEU, x) - lmn(PAN_T, x))),
  module_from32472 = A(y, modscore(out$GSE32472, x)),
  module_from220135 = A(y, modscore(out$GSE220135, x)),
  SPON1_low = if ("SPON1" %in% rownames(x)) A(y, -x["SPON1", ]) else NA)
cov <- c(naive = sum(NAIVE %in% rownames(x)), panT = sum(PAN_T %in% rownames(x)),
  neu = sum(NEU %in% rownames(x)), mod32 = sum(out$GSE32472$gene %in% rownames(x)),
  mod35 = sum(out$GSE220135$gene %in% rownames(x)))
print(round(ext, 2))
print(cov)
EXT <- list(auc = ext, cov = cov)
save(EXT, file = file.path(OUT, "external_v7.RData"))
