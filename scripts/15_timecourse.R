# Evaluate longitudinal prediction beyond GA and adjustment for expression-based covariates.
# Compare moderate/severe versus no BPD in GSE32472 and BPD versus no BPD in GSE220135.
# Use training-fold program cores for held-out time-point and day-28-minus-baseline scores.
# Evaluate within-cohort likelihood-ratio tests, leave-one-out AUC gain and cross-cohort prediction.
# Evaluate glucocorticoid-response gene set.
source("./00_setup.R")
suppressPackageStartupMessages(library(nlme))
load(file.path(OUT, "programs4.RData"))
load(file.path(OUT, "heldout4_relaxed.RData"))
load(file.path(OUT, "fig3.RData"))
set.seed(20261029)
L <- readRDS(file.path(DATA, "longitudinal_GSE32472_GSE220135.rds"))
DISC <- setNames(lapply(P4, function(pg) PR4$core[[which(PR4$program == pg)]] |> transmute(gene, es, p)),
  P4)

psc <- function(core, x) {
  if (is.null(core)) return(rep(NA_real_, ncol(x)))
  g <- core$gene[core$gene %in% rownames(x)]
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  colMeans(z * sign(core$es[match(g, core$gene)]))
}

msc <- function(g, x) {
  g <- intersect(g, rownames(x))
  z <- t(scale(t(x[g, , drop = FALSE])))
  z[is.na(z)] <- 0
  colMeans(z)
}
GC <- c("FKBP5", "TSC22D3", "ZBTB16", "DDIT4", "PER1", "KLF9", "IRS2", "CD163", "VSIG4", "ADAMTS2")
ann <- L$annotation |> filter(!is.na(t), !is.na(severity), !is.na(ga),
  id %in% c(colnames(L$expression$GSE32472), colnames(L$expression$GSE220135)))
LT <- bind_rows(lapply(c("GSE32472", "GSE220135"), function(cc) {
  a <- ann |> filter(cohort == cc)
  x <- L$expression[[cc]][, a$id]
  x <- x[apply(x, 1, sd) > 0, ]
  ho <- sapply(P4, function(pg) as.numeric(scale(psc(FOLD[[cc]][[pg]], x))))
  data.frame(a, comb = as.numeric(scale(rowMeans(ho, na.rm = TRUE))), tcell = ho[, "T cell activation"],
    comb_in = as.numeric(scale(rowMeans(sapply(P4,
      function(pg) as.numeric(scale(psc(DISC[[pg]], x))))))),
    gc = as.numeric(scale(msc(GC, x))), gc_n = sum(GC %in% rownames(x)),
    lym = as.numeric(scale(msc(c(MK$`T cells`, MK$`B cells`, MK$`NK cells`), x))),
    mye = as.numeric(scale(msc(c(MK$Monocytes, MK$Neutrophils), x))))
})) |>
  filter(severity != "Mild") |> mutate(y = as.integer(severity != "No BPD"), tf = factor(t))
cat("GC genes found:", unique(LT$gc_n), "\n")
print(LT |> count(cohort, t, y))

A <- function(y,
  s) if (length(unique(y)) < 2) NA else as.numeric(auc(roc(y, s, quiet = TRUE, direction = "<",
    levels = c(0, 1))))
TL <- list(GSE32472 = c("Day 4-7", "Day 14", "Day 28"), GSE220135 = c("Cord", "Day 14", "Day 28"))

wide <- function(cc, v) LT |> filter(cohort == cc) |> select(infant, y, ga, t,
  all_of(v)) |> tidyr::pivot_wider(names_from = t, values_from = all_of(v), names_prefix = "t")

get_d <- function(cc, tt, v) {
  if (tt == "change") {
    w <- wide(cc, v) |> filter(!is.na(t0), !is.na(t2))
    data.frame(y = w$y, ga = w$ga, s = w$t2 - w$t0)
  }
  else LT |> filter(cohort == cc, t == tt) |> transmute(y, ga, s = .data[[v]])
}

# ---- (1) prediction ----
PRED <- bind_rows(lapply(c("comb", "tcell", "gc"),
  function(v) bind_rows(lapply(c("GSE32472", "GSE220135"),
    function(cc) bind_rows(lapply(list(0, 1, 2, "change"), function(tt) {
      d <- get_d(cc, tt, v)
      oth <- setdiff(c("GSE32472", "GSE220135"), cc)
      dt <- get_d(oth, tt, v)
      f0 <- glm(y ~ ga, binomial, d)
      f1 <- glm(y ~ ga + s, binomial, d)
      lrp <- anova(f0, f1, test = "LRT")$`Pr(>Chi)`[2]
      loo <- t(sapply(seq_len(nrow(d)),
        function(i) c(predict(glm(y ~ ga, binomial, d[-i, ]), d[i, ]),
          predict(glm(y ~ ga + s, binomial, d[-i, ]), d[i, ]))))
      g0 <- glm(y ~ ga, binomial, dt)
      g1 <- glm(y ~ ga + s, binomial, dt)  # fitted in the other cohort
      lp0 <- coef(g0)["ga"] * d$ga
      lp1 <- coef(g1)["ga"] * d$ga + coef(g1)["s"] * d$s
      B <- replicate(1000, {
        i <- c(sample(which(d$y == 1), replace = TRUE), sample(which(d$y == 0), replace = TRUE))
        A(d$y[i], lp1[i]) - A(d$y[i], lp0[i])
      })
      data.frame(score = v, cohort = cc,
        time = if (tt == "change") "Change, day 28 - week 1" else TL[[cc]][tt + 1], n1 = sum(d$y),
        n0 = sum(1 - d$y),
        auc_ga = A(d$y, -d$ga), auc_s = A(d$y, d$s), or_s = exp(coef(f1)[["s"]]), p_lrt = lrp,
        loo_d = A(d$y, loo[, 2]) - A(d$y, loo[, 1]),
        ext_d = A(d$y, lp1) - A(d$y, lp0), ext_lo = quantile(B, .025, na.rm = TRUE),
        ext_hi = quantile(B, .975, na.rm = TRUE))
}))))))
rownames(PRED) <- NULL
print(PRED |> mutate(across(where(is.numeric), ~ signif(.x, 2))), width = 250)

# ---- (2) treatment / composition ----
GCT <- LT |> group_by(cohort, t) |> summarise(gc_bpd = mean(gc[y == 1]), gc_none = mean(gc[y == 0]),
  p = t.test(gc ~ y)$p.value)
print(GCT |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
ADJ <- bind_rows(lapply(c("comb", "tcell"),
  function(v) bind_rows(lapply(c("GSE32472", "GSE220135"), function(cc) {
    d <- LT |> filter(cohort == cc)
    d$s <- d[[v]]
    fm <- function(f) lme(f, random = ~ 1 | infant, data = d, method = "ML")
    est <- function(f) {
      m <- lme(f, random = ~ 1 | infant, data = d, method = "REML")
      c(b = fixef(m)[["y"]], p = summary(m)$tTable["y", "p-value"])
    }
    rbind(data.frame(score = v, cohort = cc, model = "GA", t(est(s ~ y + tf + ga))),
      data.frame(score = v, cohort = cc, model = "GA + GC response", t(est(s ~ y + tf + ga + gc))),
      data.frame(score = v, cohort = cc, model = "GA + composition",
        t(est(s ~ y + tf + ga + lym + mye))),
      data.frame(score = v, cohort = cc, model = "GA + GC + composition",
        t(est(s ~ y + tf + ga + gc + lym + mye))),
      data.frame(score = v, cohort = cc, model = "Week 1 only, GA",
        t(coef(summary(lm(s ~ y + ga, d |> filter(t == 0))))["y", c(1, 4)] |> setNames(c("b", "p")))),
      data.frame(score = v, cohort = cc, model = "Week 1 only, GA + composition",
        t(coef(summary(lm(s ~ y + ga + lym + mye, d |> filter(t == 0))))["y",
          c(1, 4)] |> setNames(c("b", "p")))))
}))))
print(ADJ |> mutate(across(where(is.numeric), ~ signif(.x, 2))))
# Compare marker scores with the supplied GSE220135 deconvolved cell fractions.
CF <- readRDS(file.path(DATA, "GSE220135_cell_fractions.rds"))
d <- LT |> filter(cohort == "GSE220135") |> left_join(CF, by = c("id" = "gsm"))
m <- lme(comb ~ y + tf + ga + Gran + CD4T + CD8T + B + nRBC, random = ~ 1 | infant,
  data = d |> filter(!is.na(Gran)), method = "REML")
CFADJ <- summary(m)$tTable["y", ]
print(CFADJ)
save(LT, PRED, GCT, ADJ, CFADJ, file = file.path(OUT, "timecourse.RData"))
