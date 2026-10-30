# New Figure 5: prediction beyond gestational age
source("./00_setup.R")
source("./figtheme.R")
load(file.path(OUT, "newfigs.RData"))

ga_of <- function(k) PH[[k]]$ga_weeks

A <- function(y, s) as.numeric(auc(roc(y, s, quiet = TRUE, direction = "<", levels = c(0, 1))))

# ---- A: why GA is strong in these cohorts ----
dA <- bind_rows(lapply(GA_COH,
  function(k) data.frame(cohort = factor(k, levels = GA_COH), ga = ga_of(k),
    grp = factor(ifelse(lab(k) == 1, "BPD", "No BPD"), levels = c("No BPD", "BPD")))))
aA <- dA |> group_by(cohort) |> summarise(auc = A(as.integer(grp == "BPD"), -ga))
fA <- ggplot(dA, aes(grp, ga, colour = grp)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 27, ymax = 31, fill = "#f3eee4") +
  geom_jitter(width = .18, height = .08, size = .9, alpha = .75) +
  stat_summary(fun = median, geom = "crossbar", width = .5, colour = INK, linewidth = .3) +
  geom_text(data = aA, aes(x = 1.5, y = 35.3, label = sprintf("GA AUC %.2f", auc)), inherit.aes = FALSE,
    size = 2.1, colour = C_GA) +
  facet_wrap(~ cohort, nrow = 1) + scale_colour_manual(values = c("No BPD" = C_NO, BPD = C_BPD),
    guide = "none") +
  scale_y_continuous(limits = c(21.5, 35.8), breaks = seq(22, 34, 2)) +
  labs(x = NULL, y = "Gestational age (weeks)", title = "GA separates the groups by design",
    subtitle = "Shaded: 27-30 week gray zone; bars, medians") + TH

# ---- B: gain over GA, all scores, fully held out ----
gx <- function(v, nm) c(es = unname(v[1]), lo = unname(v[2]), hi = unname(v[3]))
rows <- list(
  list("Top 50 genes", "Meta-analysis genes", BEYP |> filter(method == "Top 50 genes")),
  list("Top 50 genes, GA-adjusted ranking", "Meta-analysis genes",
    BEYP |> filter(method == "Top 50 genes, GA-adjusted ranking")),
  list("T cell activation program", "Meta-analysis genes",
    BEYP |> filter(method == "T cell activation program")),
  list("Combined program score", "Meta-analysis genes",
    BEYP |> filter(method == "Combined program score")),
  list("Neutrophil-T balance", "Fixed markers",
    B2 |> filter(model == "Fixed balance") |> transmute(d = dAUC, lo, hi, p)),
  list("Naive T-cell markers", "Fixed markers",
    H5 |> filter(score == "NAIVE") |> transmute(d = dAUC, lo, hi, p)),
  list("SPON1 *", "Fixed markers", H5 |> filter(score == "SPON1") |> transmute(d = dAUC, lo, hi, p)),
  list("Naive T module: GSE32472 to GSE220135", "Cross-cohort module",
    data.frame(d = XPRED$GSE220135[["d"]], lo = XPRED$GSE220135[[2]], hi = XPRED$GSE220135[[3]],
      p = NA)),
  list("Naive T module: GSE220135 to GSE32472", "Cross-cohort module",
    data.frame(d = XPRED$GSE32472[["d"]], lo = XPRED$GSE32472[[2]], hi = XPRED$GSE32472[[3]], p = NA)))
dB <- bind_rows(lapply(rows,
  function(r) data.frame(score = r[[1]], family = r[[2]], d = r[[3]]$d, lo = r[[3]]$lo,
    hi = r[[3]]$hi))) |>
  mutate(score = factor(score, levels = rev(score)),
    family = factor(family, levels = c("Meta-analysis genes", "Fixed markers", "Cross-cohort module")))
fB <- ggplot(dB, aes(d, score, colour = family)) + geom_vline(xintercept = 0, colour = LIGHT,
  linewidth = .3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = .5) + geom_point(size = 1.8) +
  geom_text(aes(x = .27, label = sprintf("%+.3f", d)), hjust = 1, size = 2, colour = INK) +
  scale_colour_manual(values = c("Meta-analysis genes" = C_PUB, "Fixed markers" = C_BAL,
    "Cross-cohort module" = C_NAIVE), name = NULL) +
  scale_x_continuous(limits = c(-.13, .28), breaks = seq(-.1, .2, .1)) +
  labs(x = "Change in held-out AUC: GA + score vs GA alone (95% CI)", y = NULL,
    title = "Small gains over GA",
    subtitle = "Weights fitted in training cohorts, applied unchanged",
    caption = "* SPON1 was selected using these cohorts.\nNaive T module, GSE32472 to GSE220135: lower bound at zero.") +
  TH + theme(axis.text.y = element_text(colour = INK), legend.position = "bottom",
    legend.margin = margin(-4, 0, 0, 0))

# ---- C: ROC in GSE220135 ----
dC <- ROCD |> filter(cohort == "GSE220135") |> arrange(model, fpr, tpr)
aC <- dC |> group_by(model) |> summarise(auc = auc[1]) |> mutate(lab = sprintf("%s: AUC %.2f", model,
  auc))
fC <- ggplot(dC, aes(fpr, tpr, colour = model)) + geom_abline(linetype = "22", colour = LIGHT,
  linewidth = .3) + geom_step(linewidth = .65, direction = "vh") +
  scale_colour_manual(values = c("GA alone" = C_GA, "GA + naive T module" = C_NAIVE),
    labels = setNames(aC$lab, aC$model), name = NULL) +
  coord_equal() + labs(x = "1 - specificity", y = "Sensitivity", title = "GSE220135, week 1",
    subtitle = "Module and weights from GSE32472 only") +
  TH + theme(legend.position = c(.98, .03), legend.justification = c(1, 0),
    legend.background = element_blank())

# ---- D: gray zone ----
dD <- GZ |> mutate(score = factor(score, levels = c("GA", "Neutrophil-T balance", "Naive T module")),
  lab = sprintf("%s\n%d BPD of %d", cohort, n_bpd, n))
fD <- ggplot(dD, aes(score, auc, fill = score)) + geom_hline(yintercept = .5, linetype = "22",
  colour = LIGHT, linewidth = .3) +
  geom_col(width = .7) + geom_text(aes(label = sprintf("%.2f", auc)), vjust = -.4,
    size = 2.1) + facet_wrap(~ lab, nrow = 1) +
  scale_fill_manual(values = c(GA = C_GA, "Neutrophil-T balance" = C_BAL, "Naive T module" = C_NAIVE),
    guide = "none") +
  scale_y_continuous(limits = c(0, 1.02),
    breaks = seq(0, 1, .25)) + scale_x_discrete(labels = c("GA", "Neutrophil-T\nbalance",
      "Naive T\nmodule")) +
  labs(x = NULL, y = "AUC within 27-30 weeks", title = "In the GA gray zone",
    subtitle = "Module from the other cohort; small numbers") + TH
fig <- ((fA | fB) + plot_layout(widths = c(1, 1.25))) / ((fC | fD) + plot_layout(widths = c(1,
  1.25))) + plot_layout(heights = c(1, 1)) + plot_annotation(tag_levels = "A")
savefig(fig, "NewFigure5", 7.4, 6.6)
print(dB)
