# ---- Signature scoring setup ----
# Prepare data for the revised main Figures 3-5.
source("./00_setup.R")
suppressPackageStartupMessages(library(nlme))
load(file.path(OUT, "fig5.RData"))
load(file.path(OUT, "v7_partA.RData"))
load(file.path(OUT, "v7_partB.RData"))
load(file.path(OUT, "fig3.RData"))
set.seed(20261029)

ga_of <- function(k) PH[[k]]$ga_weeks

A <- function(y, s) as.numeric(auc(roc(y, s, quiet = TRUE, direction = "<", levels = c(0, 1))))
NEU <- c("CSF3R", "CXCR2", "FPR1", "FCGR3B", "MME", "CMTM2")
PAN_T <- c("CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK")

lmn <- function(g, x) colMeans(x[intersect(g, rownames(x)), , drop = FALSE])

BAL <- function(x) as.numeric(scale(lmn(NEU, x) - lmn(PAN_T, x)))

mz <- function(g, x) {
  g <- intersect(g, rownames(x))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z)))
}
S32 <- S_blood |> filter(gse == "GSE32472")
PUB <- list("Luo 2024, PLoS One" = c("ACSL1", "GALNT14", "WIPI1", "MAPK14", "PROK2", "CREB5"),
  "Luo 2024, Medicine" = c("CYYR1", "GALNT14", "OLAH"),
  "Fang 2024" = c("LPIN1", "ACADSB", "WIPI1", "SLC7A11"),
  "Zou 2024" = c("ARG1", "CSF3R", "IL1R1", "IL1R2", "MMP9", "RETN", "S100A12", "SOCS3"),
  "Zhang 2025" = c("MMP9", "SIGLEC5", "DYSF", "MGAM", "S100A12"),
  "Guo 2025" = c("IL7R", "CXCL10", "DEFA4", "PRTN3", "NCAPG", "CCNB1"),
  "Tu 2025, Front Pediatr" = c("PELI1", "PYGL", "STAT4"),
  "Sun 2026" = c("CCDC141", "CHI3L2", "RIMKLB", "PDLIM7"),
  "Jia 2023" = c("NFATC3", "ERMN", "PLA2G4A", "MTMR9LP", "LOC440700"), "Li 2024 (APP)" = "APP",
  "Abilbayeva 2026" = c("BTN3A3", "AK5", "NOG", "GCSAM"),
  "Ann Med 2024 (autophagy)" = c("WIPI1", "TOMM70A", "BAG3", "PRKCQ"),
  "Windhorst 2023" = c("LAT", "LOC100505876", "MRVI1", "KLHDC4", "ATP6V0E2", "FAM102A", "FBXO3", "OGDH",
    "SLC25A3", "GAL", "CAT", "SOAT2",
    "TNFRSF10C", "BCL2A1", "DUSP2", "PDLIM7", "POLD2", "CPLX1", "SCML1", "TMED1", "HEBP1", "RRN3",
    "ZEB1", "SLC2A11"),
  "Authors 2023, 5 genes" = c("PNPO", "MSANTD2", "CD4", "SNX1", "P2RX7"),
  "Authors 2025, Th2" = c("UGT2B28", "SYNJ1", "SORD", "PTGDR2", "LGALS7", "FBN1", "TIPARP", "KCNJ8",
    "TSPAN11"))
for (k in names(E)) if ("TOMM70" %in% rownames(E[[k]]) && !"TOMM70A" %in% rownames(E[[k]])) rownames(E[[k]])[rownames(E[[k]]) == "TOMM70"] <- "TOMM70A"
DIRS <- list("Windhorst 2023" = c(LAT = -1, LOC100505876 = 1, MRVI1 = 1, KLHDC4 = -1, ATP6V0E2 = -1,
  FAM102A = -1, FBXO3 = -1, OGDH = -1, SLC25A3 = -1, GAL = 1, CAT = 1,
  SOAT2 = 1, TNFRSF10C = 1, BCL2A1 = 1, DUSP2 = -1, PDLIM7 = 1, POLD2 = 1, CPLX1 = -1, SCML1 = -1,
  TMED1 = -1, HEBP1 = 1, RRN3 = -1, ZEB1 = -1, SLC2A11 = 1))
TESTS <- list("Windhorst 2023" = c("GSE32472", "GSE125873", "GSE220135"))
# Guo et al. 2025 (PLoS One): trained in GSE220135 (validated in GSE32472) -> directions from GSE220135; scored in the other three cohorts
gg <- c("IL7R", "CXCL10", "DEFA4", "PRTN3", "NCAPG", "CCNB1")
S35 <- S_blood |> filter(gse == "GSE220135")
DIRS[["Guo 2025"]] <- setNames(sign(S35$g[match(gg, S35$gene)]), gg)
TESTS[["Guo 2025"]] <- c("GSE32472", "GSE225881", "GSE125873")

sig <- function(g, x, nm = "") {
  d <- if (!is.null(DIRS[[nm]])) unname(DIRS[[nm]][g]) else sign(S32$g[match(g, S32$gene)])
  ok <- g %in% rownames(x) & !is.na(d)
  if (!sum(ok)) return(rep(NA, ncol(x)))
  z <- t(scale(t(x[g[ok], , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z * d[ok])))
}
TEST <- c("GSE225881", "GSE125873", "GSE220135")

# ---- Signature association analysis ----
fz <- function(r, n) {
  ok <- !is.na(r)
  if (!any(ok)) return(NA)
  tanh(sum(atanh(r[ok]) * (n[ok] - 3)) / sum(n[ok] - 3))
}
SIGMAP <- bind_rows(lapply(names(PUB), function(nm) {
  tst <- if (!is.null(TESTS[[nm]])) TESTS[[nm]] else TEST
  r <- t(sapply(tst, function(k) {
    s <- sig(PUB[[nm]], E[[k]], nm)
    if (all(is.na(s))) return(c(NA, NA, NA, NA, ncol(E[[k]])))
    y <- lab(k)
    c(cor(s, BAL(E[[k]])), cor(s, CS[[k]][, "Neutrophils"]), cor(s, CS[[k]][, "T cells"]),
      if (k %in% GA_COH) cor(s[y == 0], ga_of(k)[y == 0]) else NA, ncol(E[[k]]))
  }))
  ga_after <- t(sapply(intersect(tst, GA_COH), function(k) {
    s <- sig(PUB[[nm]], E[[k]], nm)
    if (all(is.na(s))) return(rep(NA, 4))
    b <- BAL(E[[k]])
    f1 <- lm(s ~ lab(k) + ga_of(k))
    f2 <- lm(s ~ lab(k) + ga_of(k) + b)
    c(coef(f1)[2], vcov(f1)[2, 2], coef(f2)[2], vcov(f2)[2, 2])
  }))
  ok <- !is.na(ga_after[, 1])
  p1 <- if (any(ok)) re_pool(ga_after[ok, 1], ga_after[ok, 2]) else c(es = NA, lo = NA, hi = NA, p = NA)
  p2 <- if (any(ok)) re_pool(ga_after[ok, 3], ga_after[ok, 4]) else c(es = NA, lo = NA, hi = NA, p = NA)
  data.frame(signature = nm, r_bal = fz(r[, 1], r[, 5]), r_neu = fz(r[, 2], r[, 5]),
    r_T = fz(r[, 3], r[, 5]), r_GA0 = mean(r[, 4], na.rm = TRUE),
    b_GA = p1[["es"]], lo_GA = p1[["lo"]], hi_GA = p1[["hi"]], p_GA = p1[["p"]], b_GAbal = p2[["es"]],
    lo_GAbal = p2[["lo"]], hi_GAbal = p2[["hi"]], p_GAbal = p2[["p"]])
}))

# ---- Methylation comparison ----
L <- readRDS(file.path(DATA, "longitudinal_GSE32472_GSE220135.rds"))
CF <- readRDS(file.path(DATA, "GSE220135_cell_fractions.rds"))
METH <- L$annotation |> filter(cohort == "GSE220135", !is.na(t), !is.na(ga), !is.na(severity),
  id %in% colnames(L$expression$GSE220135)) |> mutate(y = as.integer(severity != "No BPD")) |>
  left_join(CF, by = c("id" = "gsm")) |> filter(!is.na(Gran))
METH$BAL <- BAL(L$expression$GSE220135[, METH$id])
METH$mNLR <- log((METH$Gran + 1) / (METH$CD4T + METH$CD8T + METH$B + METH$NK + 1))

# ---- Fig 4: maturity-matched discovery ----
disc <- function(cc) {
  a <- L$annotation |> filter(cohort == cc, !is.na(t), !is.na(ga), !is.na(severity),
    id %in% colnames(L$expression[[cc]]))
  if (cc == "GSE32472") a <- a |> filter(severity != "Mild")
  a <- a |> mutate(y = as.integer(severity != "No BPD"))
  x <- L$expression[[cc]][, a$id]
  x <- x[apply(x, 1, sd) > 0, ]
  X <- model.matrix(~ y + ga + factor(t), a)
  dc <- duplicateCorrelation(x, X, block = a$infant)
  f <- eBayes(lmFit(x, X, block = a$infant, correlation = dc$consensus))
  tt <- topTable(f, coef = "y", number = Inf, sort.by = "none")
  data.frame(gene = rownames(tt), t = tt$t, p = tt$P.Value)
}
D32 <- disc("GSE32472")
D35 <- disc("GSE220135")
SCAT <- inner_join(D32, D35, by = "gene", suffix = c("32", "35"))
NAIVE <- c("CCR7", "TCF7", "LEF1", "BCL11B", "CD27", "CD28", "IL7R", "BACH2", "SATB1", "SPON1", "NOG",
  "CHI3L2", "GPA33", "MAN1C1", "CD6")
SCAT$naive <- SCAT$gene %in% NAIVE

modscore <- function(M, x) {
  g <- M$gene[M$gene %in% rownames(x)]
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z * sign(M$t[match(g, M$gene)]))))
}
# maturity check: r(score, GA) among no-BPD infants at week 1 (GA cohorts; module scored with the module from the OTHER longitudinal cohort)
mod_for <- function(k) if (k == "GSE32472") out$GSE220135 else out$GSE32472
MATR <- bind_rows(
  bind_rows(lapply(GA_COH, function(k) {
    y <- lab(k)
    g <- ga_of(k)
    R <- head(RK_all, 50)
    s <- {
      ok <- R$gene %in% rownames(E[[k]])
      z <- t(scale(t(E[[k]][R$gene[ok], ])))
      z[!is.finite(z)] <- 0
      as.numeric(colMeans(z * sign(R$es[ok])))
    }
    data.frame(score = "Top 50 genes, not GA-adjusted", cohort = k, r = cor(s[y == 0], g[y == 0]))
  })),
  bind_rows(lapply(GA_COH, function(k) {
    y <- lab(k)
    g <- ga_of(k)
    data.frame(score = "Neutrophil-T balance", cohort = k, r = cor(BAL(E[[k]])[y == 0], g[y == 0]))
  })),
  bind_rows(lapply(GA_COH, function(k) {
    y <- lab(k)
    g <- ga_of(k)
    data.frame(score = "Naive T module (other cohort)", cohort = k,
      r = cor(modscore(mod_for(k), E[[k]])[y == 0], g[y == 0]))
  })),
  bind_rows(lapply(names(PUB), function(nm) bind_rows(lapply(GA_COH, function(k) {
    y <- lab(k)
    s <- sig(PUB[[nm]], E[[k]])
    if (all(is.na(s))) return(NULL)
    data.frame(score = "Published signatures", sig = nm, cohort = k,
      r = cor(s[y == 0], ga_of(k)[y == 0]))
})))))
# time course of the cross-cohort module
TC <- bind_rows(lapply(c("GSE32472", "GSE220135"), function(cc) {
  a <- L$annotation |> filter(cohort == cc, !is.na(t), !is.na(ga), !is.na(severity),
    id %in% colnames(L$expression[[cc]])) |>
    mutate(y = as.integer(severity != "No BPD"),
      grp = ifelse(severity == "No BPD", "No BPD", ifelse(severity == "Mild", "Mild", "BPD")))
  if (cc == "GSE220135") a$grp[a$grp == "Mild"] <- "BPD"
  a$s <- -modscore(mod_for(cc), L$expression[[cc]][, a$id])
  a$cohort <- cc
  a
}))  # sign flipped: higher = more naive T expression
TCP <- bind_rows(lapply(c("GSE32472", "GSE220135"), function(cc) {
  a <- TC |> filter(cohort == cc)
  if (cc == "GSE32472") a <- a |> filter(grp != "Mild")
  m <- lme(s ~ y + factor(t) + ga, random = ~ 1 | infant, data = a, method = "REML")
  data.frame(cohort = cc, b = fixef(m)[["y"]], p = summary(m)$tTable["y", 5])
}))

# ---- Fig 5: prediction ----
wk1 <- function(cc, M) {
  a <- L$annotation |> filter(cohort == cc, t == 0, !is.na(ga), !is.na(severity),
    id %in% colnames(L$expression[[cc]])) |> mutate(y = as.integer(severity != "No BPD"))
  a$s <- modscore(M, L$expression[[cc]][, a$id])
  a$bal <- BAL(L$expression[[cc]][, a$id])
  a
}
XPRED <- list()
ROCD <- NULL
GZ <- NULL
for (src in names(out)) {
  tgt <- setdiff(names(out), src)
  dtr <- wk1(src, out[[src]])
  d <- wk1(tgt, out[[src]])
  f0 <- glm(y ~ ga, binomial, dtr)
  f1 <- glm(y ~ ga + s, binomial, dtr)
  lp0 <- predict(f0, d)
  lp1 <- predict(f1, d)
  B <- replicate(2000, {
    i <- c(sample(which(d$y == 1), replace = TRUE), sample(which(d$y == 0), replace = TRUE))
    A(d$y[i], lp1[i]) - A(d$y[i], lp0[i])
  })
  XPRED[[tgt]] <- c(d = A(d$y, lp1) - A(d$y, lp0), lo = quantile(B, .025), hi = quantile(B, .975),
    auc0 = A(d$y, lp0), auc1 = A(d$y, lp1))
  for (nm in c("GA alone", "GA + naive T module")) {
    lp <- if (nm == "GA alone") lp0 else lp1
    r <- roc(d$y, lp, quiet = TRUE, direction = "<", levels = c(0, 1))
    ROCD <- rbind(ROCD,
      data.frame(cohort = tgt, model = nm, fpr = 1 - r$specificities, tpr = r$sensitivities,
        auc = as.numeric(auc(r))))
  }
  gz <- d$ga >= 27 & d$ga < 31
  GZ <- rbind(GZ,
    data.frame(cohort = tgt, n = sum(gz), n_bpd = sum(d$y[gz]),
      score = c("GA", "Neutrophil-T balance", "Naive T module"),
      auc = c(A(d$y[gz], -d$ga[gz]), A(d$y[gz], d$bal[gz]), A(d$y[gz], d$s[gz]))))
}
save(AUD, SIGMAP, METH, SCAT, MATR, XC, TC, TCP, XPRED, ROCD, GZ, BEYP, B2, H5, TOPN, GA_AUC,
  file = file.path(OUT, "newfigs.RData"))
print(SIGMAP |> mutate(across(where(is.numeric), ~ round(.x, 2))))
print(TCP)
print(XPRED)
print(GZ |> mutate(auc = round(auc, 2)))
cat("Quadrant counts (p<.05 both):",
  with(SCAT,
    c(down_down = sum(p32 < .05 & p35 < .05 & t32 < 0 & t35 < 0),
      up_up = sum(p32 < .05 & p35 < .05 & t32 > 0 & t35 > 0),
      discordant = sum(p32 < .05 & p35 < .05 & sign(t32) != sign(t35)))), " n genes", nrow(SCAT), "\n")
print(MATR |> group_by(score) |> summarise(r = mean(r), n = n()))
