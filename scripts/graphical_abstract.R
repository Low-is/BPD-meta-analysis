# Draw a four-panel graphical abstract at 1200 x 1200 pixels and 300 dpi.
source("./00_setup.R")
source("./figtheme.R")
load(file.path(OUT, "newfigs.RData"))
C_NEUg <- "#c98a1b"
C_Tg <- "#1d9fb0"
C_WK <- "#a7a59e"

txt <- function(lines, sizes, faces, cols, ys, x = .5, hj = .5) {
  p <- ggplot() + xlim(0, 1) + ylim(0, 1) + theme_void()
  for (i in seq_along(lines)) p <- p + annotate("text", x = x, y = ys[i], label = lines[i],
    size = sizes[i], fontface = faces[i], colour = cols[i], hjust = hj, family = FONT, lineheight = .95)
  p
}
GTH <- theme_classic(base_size = 6,
  base_family = FONT) + theme(axis.text = element_text(size = 4.4, colour = MUTED),
    axis.title = element_text(size = 4.8, colour = INK),
    plot.title = element_text(size = 6.2, face = "bold", hjust = 0, margin = margin(b = 1)),
    plot.subtitle = element_text(size = 4.6, colour = MUTED, hjust = 0, margin = margin(b = 2)),
    axis.line = element_line(linewidth = .25, colour = MUTED),
    axis.ticks = element_line(linewidth = .25, colour = MUTED), axis.ticks.length = unit(1.5, "pt"),
    plot.margin = margin(4, 6, 4, 4), legend.position = "none",
    plot.background = element_rect(fill = "#fafaf7", colour = "#e3e1da", linewidth = .4))

axis_of <- function(n, t) ifelse(pmax(abs(n), abs(t)) < .5, "Weak",
  ifelse(abs(n) >= abs(t), "Neutrophil", "T cell"))
AXC <- c(Neutrophil = C_NEUg, "T cell" = C_Tg, Weak = C_WK)
head <- txt(c("What do blood signatures of bronchopulmonary dysplasia measure?",
  "4 independent first-week blood cohorts  |  196 preterm infants  |  15 published gene signatures,",
  "each reconstructed and tested only in cohorts not used to derive it"),
  c(2.7, 1.95, 1.95), c("bold", "plain", "plain"), c(INK, MUTED, MUTED), c(.78, .40, .14))
# 1: external performance
G <- SIGMAP |> mutate(axis = axis_of(r_neu, r_T)) |> select(signature, axis)
rn <- c("Authors 2025, Th2 (asthma)" = "Authors 2025, Th2")
d1 <- AUD |> filter(group != "Reference") |> mutate(signature = ifelse(signature %in% names(rn),
  rn[signature], signature)) |> left_join(G, by = "signature") |>
  arrange(auc) |> mutate(i = row_number())
ref <- AUD |> filter(grepl("top 5 ", signature))
p1 <- ggplot(d1, aes(auc, i, colour = axis)) + geom_vline(xintercept = .5, colour = MUTED,
  linewidth = .25, linetype = 2) + geom_vline(xintercept = ref$auc, colour = INK, linewidth = .3,
    linetype = 3) +
  annotate("text", x = ref$auc + .015, y = 2, label = "5-gene\nreference", size = 1.35, hjust = 0,
    colour = INK, lineheight = .85) +
  geom_errorbarh(aes(xmin = lo, xmax = pmin(hi, 1)), height = 0, linewidth = .3,
    alpha = .7) + geom_point(size = .9) +
  scale_colour_manual(values = AXC) + scale_x_continuous(limits = c(.05, 1),
    breaks = c(.1, .3, .5, .7, .9)) + scale_y_continuous(breaks = NULL) +
  labs(x = "Area under the curve in new cohorts", y = NULL, title = "1  Limited external performance",
    subtitle = "Pooled area under the curve 0.30-0.70") +
  GTH + theme(axis.line.y = element_blank())
# 2: immune axes
d2 <- SIGMAP |> filter(!is.na(r_neu)) |> mutate(axis = axis_of(r_neu, r_T))
p2 <- ggplot(d2, aes(r_neu, r_T, colour = axis)) + geom_hline(yintercept = 0, colour = LIGHT,
  linewidth = .2) + geom_vline(xintercept = 0, colour = LIGHT, linewidth = .2) +
  geom_point(size = 1.15) + scale_colour_manual(values = AXC) +
  annotate("text", x = .80, y = .14, label = "neutrophil-\nassociated", size = 1.5, colour = C_NEUg,
    lineheight = .85) +
  annotate("text", x = .40, y = -.90, label = "T-cell-\nassociated", size = 1.5, colour = C_Tg,
    lineheight = .85) +
  scale_x_continuous(limits = c(-.25, 1.02),
    breaks = c(0, .5, 1)) + scale_y_continuous(limits = c(-1, .4), breaks = c(-1, -.5, 0)) +
  labs(x = "Correlation with neutrophil genes", y = "Correlation with T-cell genes",
    title = "2  Converge on immune composition",
    subtitle = "Scores track neutrophil and T-cell genes") + GTH
# 3: gestational age
d3 <- bind_rows(lapply(GA_COH,
  function(k) data.frame(cohort = factor(k, levels = c("GSE32472", "GSE225881", "GSE220135")),
    ga = PH[[k]]$ga_weeks, grp = factor(ifelse(lab(k) == 1, "Yes", "No"), levels = c("No", "Yes")))))
set.seed(1)
p3 <- ggplot(d3, aes(grp, ga, colour = grp)) + geom_jitter(width = .22, height = .1, size = .35,
  alpha = .7) +
  stat_summary(fun = median, geom = "crossbar", width = .55, linewidth = .25, colour = INK) +
  facet_wrap(~ cohort, nrow = 1,
    labeller = as_labeller(c(GSE32472 = "Cohort 1", GSE225881 = "Cohort 2",
      GSE220135 = "Cohort 3"))) + scale_colour_manual(values = c(No = C_NO, Yes = C_BPD)) +
  labs(x = "Developed bronchopulmonary dysplasia", y = "Gestational age (weeks)",
    title = "3  Gestational age dominates", subtitle = "Gene scores add little beyond gestational age") +
  GTH + theme(strip.background = element_blank(),
    strip.text = element_text(size = 4.4, colour = MUTED, margin = margin(0, 0, 1, 0)),
    axis.text.x = element_text(size = 4))
# 4: exploratory naive T-cell expression
d4 <- TC |> filter(cohort == "GSE220135") |> group_by(t,
  grp = ifelse(y == 1, "BPD", "No BPD")) |> summarise(m = mean(s), se = sd(s) / sqrt(n()),
    .groups = "drop") |>
  mutate(time = factor(c("Cord", "Day 14", "Day 28")[t + 1], levels = c("Cord", "Day 14", "Day 28")))
p4 <- ggplot(d4,
  aes(time, m, colour = grp,
    group = grp)) + geom_line(linewidth = .45) + geom_pointrange(aes(ymin = m - se, ymax = m + se),
      size = .1, linewidth = .3) +
  scale_colour_manual(values = c("No BPD" = C_NO, BPD = C_BPD)) +
  annotate("text", x = 3.1, y = .82, label = "No disease", size = 1.5, colour = C_NO,
    hjust = 1) + annotate("text", x = 3.1, y = -1.05, label = "Disease", size = 1.5, colour = C_BPD,
      hjust = 1) +
  annotate("label", x = 1.05, y = -1.95, label = "Pre-specified test: not met (p = 0.54)", size = 1.35,
    hjust = 0, colour = INK, fill = "white", label.size = .15, label.padding = unit(1.2, "pt")) +
  coord_cartesian(ylim = c(-2.05, 1)) +
  labs(x = NULL, y = "Naive T-cell gene score", title = "4  Exploratory: lower naive T-cell genes",
    subtitle = "Adjusted for gestational age; seen in 2 cohorts") + GTH
foot <- txt(c("Blood signatures of bronchopulmonary dysplasia largely reflect maturity and immune cells",
  "Next: maturity-matched cohorts with thymic output, blood counts, and cytometry"),
  c(1.95, 1.85), c("bold", "plain"), c(INK, MUTED), c(.70, .25))
fig <- head / ((p1 | p2) / (p3 | p4)) / foot + plot_layout(heights = c(.55, 3.3, .42))
ggsave(file.path(FIG, "Graphical_abstract.png"), fig, width = 4, height = 4, dpi = 300, bg = "white")
ggsave(file.path(FIG, "Graphical_abstract.pdf"), fig, width = 4, height = 4, device = cairo_pdf)
