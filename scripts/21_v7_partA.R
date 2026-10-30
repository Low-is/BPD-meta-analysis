# Naive T-cell analyses, protocol v7 Part A.
source("./00_setup.R")
suppressPackageStartupMessages(library(nlme))
set.seed(20261029)

ga_of <- function(k) PH[[k]]$ga_weeks
NAIVE <- c("CCR7", "TCF7", "LEF1", "BCL11B", "CD27", "CD28", "IL7R", "BACH2", "SATB1")
PAN_T <- c("CD3D", "CD3E", "CD3G", "CD2", "CD247", "LCK")
RTE   <- c("PTK7", "CR2", "PECAM1", "SOX4", "TOX", "AIF1")
NEU   <- c("CSF3R", "CXCR2", "FPR1", "FCGR3B", "MME", "CMTM2")

mz <- function(g, x) {
  g <- intersect(g, rownames(x))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z)))
}

scores <- function(x) {
  n <- mz(NAIVE, x)
  p <- mz(PAN_T, x)
  list(NAIVE = n, PAN_T = p, NAIVE_in_T = as.numeric(scale(resid(lm(n ~ p)))), RTE_lit = mz(RTE, x),
    NEU = mz(NEU, x), SPON1 = as.numeric(scale(x["SPON1", ])))
}

A <- function(y, s) as.numeric(auc(roc(y, s, quiet = TRUE, direction = "<", levels = c(0, 1))))
cat("Gene coverage (NAIVE/PAN_T/RTE/NEU):\n")
print(sapply(BLOOD,
  function(k) sapply(list(NAIVE, PAN_T, RTE, NEU), function(g) sum(g %in% rownames(E[[k]])))))
SCW <- lapply(setNames(BLOOD, BLOOD), function(k) scores(E[[k]]))
out <- list()

# ---- H1 and H4 (week 1): pooled GA-adjusted BPD coefficient; H2: r with GA among no-BPD ----
WK1 <- bind_rows(lapply(names(SCW[[1]]), function(v) {
  r <- t(sapply(GA_COH, function(k) {
    s <- SCW[[k]][[v]]
    y <- lab(k)
    g <- ga_of(k)
    f <- lm(s ~ y + g)
    c(b = coef(f)[["y"]], var = vcov(f)["y", "y"], r0 = cor(s[y == 0], g[y == 0]), n0 = sum(y == 0),
      auc_low = A(y, -s))
  }))
  pb <- re_pool(r[, "b"], r[, "var"])
  zr <- atanh(r[, "r0"])
  pr <- re_pool(zr, 1 / (r[, "n0"] - 3))
  data.frame(score = v, b_GA = pb[["es"]], lo = pb[["lo"]], hi = pb[["hi"]], p = pb[["p"]],
    neg_cohorts = sum(r[, "b"] < 0),
    per_cohort_b = paste(sprintf("%+.2f", r[, "b"]), collapse = "/"), r_GA_noBPD = tanh(pr[["es"]]),
    r_lo = tanh(pr[["lo"]]), r_hi = tanh(pr[["hi"]]),
    auc_lowerBPD = paste(sprintf("%.2f", r[, "auc_low"]), collapse = "/"),
    auc_125873 = A(lab("GSE125873"), -SCW$GSE125873[[v]]))
}))
cat("\n== H1/H4/H2: week 1, GA-adjusted, any BPD (3 GA cohorts) ==\n")
print(WK1 |> mutate(across(where(is.numeric), ~ signif(.x, 2))), row.names = FALSE)

# ---- H3/H4 longitudinal ----
L <- readRDS(file.path(DATA, "longitudinal_GSE32472_GSE220135.rds"))
LONG <- NULL
for (cc in c("GSE32472", "GSE220135")) for (def in c("any", "modsev")) {
  if (cc == "GSE220135" && def == "modsev") next
  a <- L$annotation |> filter(cohort == cc, !is.na(t), !is.na(ga), !is.na(severity),
    id %in% colnames(L$expression[[cc]]))
  if (def == "modsev") a <- a |> filter(severity != "Mild")
  a <- a |> mutate(y = as.integer(severity != "No BPD"))
  x <- L$expression[[cc]][, a$id]
  x <- x[apply(x, 1, sd) > 0, ]
  S <- scores(x)
  for (v in names(S)) {
    a$s <- S[[v]]
    m <- lme(s ~ y + factor(t) + ga, random = ~ 1 | infant, data = a, method = "REML")
    LONG <- rbind(LONG,
      data.frame(cohort = cc, outcome = def, score = v, b = fixef(m)[["y"]],
        p = summary(m)$tTable["y", 5],
        auc_by_time = paste(sapply(0:2, function(tt) {
          i <- a$t == tt
          sprintf("%.2f", A(a$y[i], -a$s[i]))
        }), collapse = "/")))
  }
}
cat("\n== H3/H4: all time points, GA-adjusted mixed models ==\n")
print(LONG |> mutate(across(where(is.numeric), ~ signif(.x, 2))), row.names = FALSE)

# ---- secondary: GSE32472 moderate/severe vs none at week 1 ----
w <- L$annotation |> filter(cohort == "GSE32472", t == 0, !is.na(ga),
  severity %in% c("No BPD", "Moderate/severe"),
  id %in% colnames(L$expression$GSE32472)) |> mutate(y = as.integer(severity != "No BPD"))
Sw <- scores(L$expression$GSE32472[, w$id])
cat("\nGSE32472 week 1, mod/sev vs none (GA-adj):",
  paste(sapply(names(Sw),
    function(v) sprintf("%s %+.2f (p %.2g)", v, coef(summary(lm(Sw[[v]] ~ w$y + w$ga)))[2, 1],
      coef(summary(lm(Sw[[v]] ~ w$y + w$ga)))[2, 4])), collapse = "; "), "\n")

# ---- H5 prediction: held-out dAUC over GA, and gray zone ----
H5 <- bind_rows(lapply(c("NAIVE", "PAN_T", "NAIVE_in_T", "NEU", "SPON1"), function(v) {
  sg <- if (v == "NEU") 1 else -1
  ho <- t(sapply(GA_COH, function(h) {
    tr <- setdiff(GA_COH, h)
    dtr <- bind_rows(lapply(tr,
      function(k) data.frame(y = lab(k), ga = ga_of(k), s = sg * SCW[[k]][[v]], coh = k)))
    d <- data.frame(y = lab(h), ga = ga_of(h), s = sg * SCW[[h]][[v]])
    f0 <- glm(y ~ coh + ga, binomial, dtr)
    f1 <- glm(y ~ coh + ga + s, binomial, dtr)
    lp0 <- coef(f0)["ga"] * d$ga
    lp1 <- coef(f1)["ga"] * d$ga + coef(f1)["s"] * d$s
    B <- replicate(2000, {
      i <- c(sample(which(d$y == 1), replace = TRUE), sample(which(d$y == 0), replace = TRUE))
      A(d$y[i], lp1[i]) - A(d$y[i], lp0[i])
    })
    gz <- d$ga >= 27 & d$ga < 31
    c(d = A(d$y, lp1) - A(d$y, lp0), v = var(B),
      gz_ga = if (length(unique(d$y[gz])) == 2) A(d$y[gz], -d$ga[gz]) else NA,
      gz_s = if (length(unique(d$y[gz])) == 2) A(d$y[gz], d$s[gz]) else NA, gz_n1 = sum(d$y[gz]))
  }))
  p <- re_pool(ho[, "d"], ho[, "v"])
  data.frame(score = v, dAUC = p[["es"]], lo = p[["lo"]], hi = p[["hi"]], p = p[["p"]],
    per_cohort = paste(sprintf("%+.3f", ho[, "d"]), collapse = "/"),
    grayzone_AUC_GA = paste(sprintf("%.2f", ho[, "gz_ga"]), collapse = "/"),
    grayzone_AUC_score = paste(sprintf("%.2f", ho[, "gz_s"]), collapse = "/"),
    gz_BPD_n = paste(ho[, "gz_n1"], collapse = "/"))
}))
cat("\n== H5: held-out gain over GA; gray zone 27-30 wk ==\n")
print(H5 |> mutate(across(where(is.numeric), ~ signif(.x, 2))), row.names = FALSE)

# ---- cross-cohort discovery -> test of the data-driven module ----
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
  tt <- topTable(f, coef = "y", number = Inf)
  tt$gene <- rownames(tt)
  tt
}

modscore <- function(M, x) {
  g <- M$gene[M$gene %in% rownames(x)]
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[!is.finite(z)] <- 0
  as.numeric(scale(colMeans(z * sign(M$t[match(g, M$gene)]))))
}
XC <- NULL
for (src in c("GSE32472", "GSE220135")) {
  DD <- disc(src)
  tgt <- setdiff(c("GSE32472", "GSE220135"), src)
  common <- Reduce(intersect, list(rownames(L$expression$GSE32472), rownames(L$expression$GSE220135)))
  M <- DD |> filter(gene %in% common) |> arrange(P.Value) |> head(50)
  cat("\nModule from", src, "(top 50, GA-adjusted, all samples):",
    paste(head(M$gene, 30), collapse = ", "), "... | lower in BPD:", sum(M$t < 0), "/50\n")
  a <- L$annotation |> filter(cohort == tgt, !is.na(t), !is.na(ga), !is.na(severity),
    id %in% colnames(L$expression[[tgt]])) |> mutate(y = as.integer(severity != "No BPD"))
  a$s <- modscore(M, L$expression[[tgt]][, a$id])
  m <- lme(s ~ y + factor(t) + ga, random = ~ 1 | infant, data = a, method = "REML")
  w1 <- a |> filter(t == 0)
  XC <- rbind(XC,
    data.frame(module_from = src, tested_in = tgt, week1_b_GA = coef(lm(s ~ y + ga, w1))[["y"]],
      week1_p = coef(summary(lm(s ~ y + ga, w1)))["y", 4], week1_AUC = A(w1$y, w1$s),
      all_b_GA = fixef(m)[["y"]], all_p = summary(m)$tTable["y", 5],
      r_GA_noBPD_wk1 = with(w1 |> filter(y == 0), cor(s, ga))))
  for (k in c("GSE225881", "GSE125873")) {
    dk <- data.frame(s = modscore(M, E[[k]]), y = lab(k), ga = if (k == "GSE225881") ga_of(k) else NA)
    cat("  ", k, "module genes found:", sum(M$gene %in% rownames(E[[k]])), " n =", nrow(dk), "\n")
    fk <- if (k == "GSE225881") coef(summary(lm(s ~ y + ga, dk)))["y", ] else NULL
    XC <- rbind(XC,
      data.frame(module_from = src, tested_in = k,
        week1_b_GA = if (!is.null(fk)) fk[1] else hedges(dk$s, dk$y)[["g"]],
        week1_p = if (!is.null(fk)) fk[4] else NA, week1_AUC = A(dk$y, dk$s), all_b_GA = NA, all_p = NA,
        r_GA_noBPD_wk1 = if (k == "GSE225881") with(dk[dk$y == 0, ], cor(s, ga)) else NA))
  }
  out[[src]] <- M
}
cat("\n== Cross-cohort module test (score signed so that higher = BPD-like) ==\n")
print(XC |> mutate(across(where(is.numeric), ~ signif(.x, 2))), row.names = FALSE)
save(WK1, LONG, H5, XC, out, file = file.path(OUT, "v7_partA.RData"))
