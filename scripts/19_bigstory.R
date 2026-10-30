# Assess GA distributions, maturity associations, GA-restricted discrimination and longitudinal discovery.
source("./00_setup.R")
load(file.path(OUT, "fig3.RData"))
load(file.path(OUT, "fig5.RData"))
load(file.path(OUT, "heldout4_relaxed.RData"))

ga_of <- function(k) PH[[k]]$ga_weeks

A <- function(y, s) as.numeric(auc(roc(y, s, quiet = TRUE, direction = "<", levels = c(0, 1))))
cat("== 1) GA by group ==\n")
for (k in GA_COH) {
  y <- lab(k)
  g <- ga_of(k)
  cat(k, " BPD median", median(g[y == 1]), "[", range(g[y == 1]), "]  no BPD median", median(g[y == 0]),
    "[", range(g[y == 0]), "]  AUC GA", round(A(y, -g), 2),
    "| overlap zone n:", sum(g >= max(min(g[y == 0]), 0) & g <= max(g[y == 1])), "\n")
}

mz <- function(g, x) {
  g <- intersect(g, rownames(x))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  colMeans(z)
}

NT <- function(x) as.numeric(scale(mz(MK$Neutrophils, x) - mz(c(MK$`T cells`, MK$`CD4 T`, MK$`CD8 T`),
  x)))

sig_score <- function(genes, dirs, x) {
  ok <- genes %in% rownames(x)
  g <- genes[ok]
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z * dirs[ok])))
}
S32 <- S_blood |> filter(gse == "GSE32472")
PUB <- list(Zou = c("ARG1", "CSF3R", "IL1R1", "IL1R2", "MMP9", "RETN", "S100A12", "SOCS3"),
  Zhang = c("MMP9", "SIGLEC5", "DYSF", "MGAM", "S100A12"),
  Guo = c("IL7R", "CXCL10", "DEFA4", "PRTN3", "NCAPG", "CCNB1"),
  Sun = c("CCDC141", "CHI3L2", "RIMKLB", "PDLIM7"), `Tu FP` = c("PELI1", "PYGL", "STAT4"),
  `Luo PLoS` = c("ACSL1", "GALNT14", "WIPI1", "MAPK14", "PROK2", "CREB5"),
  `Authors 5` = c("PNPO", "MSANTD2", "CD4", "SNX1", "P2RX7"))
SC <- c(lapply(PUB, function(g) {
  d <- sign(S32$g[match(g, S32$gene)])
  d[is.na(d)] <- 0
  function(x) sig_score(g, d, x)
}),
  list(`Top 50 (all 4)` = function(x) {
    R <- head(RK_all, 50)
    sig_score(R$gene, sign(R$es), x)
  }, `NT score` = NT, SPON1 = function(x) -as.numeric(scale(x["SPON1", ]))))
cat("\n== 2) maturity test: r(score, GA) among NO-BPD infants  |  BPD AUC ==\n")
MAT <- bind_rows(lapply(names(SC), function(nm) {
  r <- sapply(GA_COH, function(k) {
    y <- lab(k)
    s <- SC[[nm]](E[[k]])
    c(r0 = cor(s[y == 0], ga_of(k)[y == 0]), auc = A(y, s))
  })
  data.frame(score = nm, r_GA_noBPD = paste(sprintf("%+.2f", r["r0", ]), collapse = "/"),
    mean_r = mean(r["r0", ]), auc_BPD = paste(sprintf("%.2f", r["auc", ]), collapse = "/"))
}))
print(MAT |> arrange(mean_r) |> mutate(mean_r = round(mean_r, 2)), row.names = FALSE)
cat("\n== 3) gray zone: GA 27-30 wk (both groups present) ==\n")
for (k in GA_COH) {
  y <- lab(k)
  g <- ga_of(k)
  i <- g >= 27 & g < 31
  if (sum(y[i]) < 3 || sum(1 - y[i]) < 3) next
  cat(k, " n", sum(i), "(", sum(y[i]), "BPD)  AUC GA", round(A(y[i], -g[i]), 2), " NT",
    round(A(y[i], NT(E[[k]])[i]), 2), " SPON1", round(A(y[i], SC$SPON1(E[[k]])[i]), 2),
    " Top50", round(A(y[i], SC$`Top 50 (all 4)`(E[[k]])[i]), 2), "\n")
}
cat("\n== 4) GA-adjusted discovery: week 1 only vs all longitudinal samples ==\n")
L <- readRDS(file.path(DATA, "longitudinal_GSE32472_GSE220135.rds"))

fitc <- function(cc, all = TRUE) {
  a <- L$annotation |> filter(cohort == cc, !is.na(t), !is.na(ga), severity != "Mild",
    id %in% colnames(L$expression[[cc]])) |> mutate(y = as.integer(severity != "No BPD"))
  if (!all) a <- a |> filter(t == 0)
  x <- L$expression[[cc]][, a$id]
  x <- x[apply(x, 1, sd) > 0, ]
  X <- if (all) model.matrix(~ y + ga + factor(t), a) else model.matrix(~ y + ga, a)
  if (all) {
    dc <- duplicateCorrelation(x, X, block = a$infant)
    f <- lmFit(x, X, block = a$infant, correlation = dc$consensus)
  } else f <- lmFit(x, X)
  tt <- topTable(eBayes(f), coef = "y", number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), t = tt$t, p = tt$P.Value, fdr = tt$adj.P.Val)
}
for (all in c(FALSE, TRUE)) {
  a <- fitc("GSE32472", all)
  b <- fitc("GSE220135", all)
  m <- inner_join(a, b, by = "gene", suffix = c("32", "35"))
  rep <- m |> filter(p32 < .05, p35 < .05, sign(t32) == sign(t35))
  chance <- nrow(m) * .05 * .05 / 2 * 2
  z <- (m$t32 + m$t35) / sqrt(2)
  pz <- 2 * pnorm(-abs(z))
  fz <- p.adjust(pz, "BH")
  cat(ifelse(all, "ALL samples", "Week 1 only"), ": genes", nrow(m), "| FDR<.05 in 32472:",
    sum(a$fdr < .05), " in 220135:", sum(b$fdr < .05),
    "| p<.05 both, same sign:", nrow(rep), "(chance ~", round(nrow(m) * .05 * .05 / 2),
    ") | combined z FDR<.05:", sum(fz < .05), "\n")
  if (all) {
    top <- m |> mutate(z = z, fdr = fz) |> filter(fdr < .05) |> arrange(-abs(z))
    print(head(top, 25))
    assign("TOPALL", top)
  }
}
save(MAT, TOPALL, file = file.path(OUT, "bigstory.RData"))
