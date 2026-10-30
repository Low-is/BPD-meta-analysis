# Figure 3: external signature performance and associations with immune expression axes.
source("./00_setup.R")
source("./figtheme.R")
load(file.path(OUT, "newfigs.RData"))
load(file.path(OUT, "axis_pca.RData"))
C_NEUg <- "#c98a1b"
C_Tg <- "#1d9fb0"
C_WK <- "#a7a59e"
AX <- c("Neutrophil-associated" = C_NEUg, "T-cell-associated" = C_Tg, "Weak / neither" = C_WK,
  "Reference list" = INK)
rn <- c("Authors 2025, Th2 (asthma)" = "Authors 2025, Th2")

short_of <- function(s) dplyr::case_when(s == "Tu 2025, Front Pediatr" ~ "Tu 2025",
  s == "Luo 2024, Medicine" ~ "Luo 2024b",
  s == "Luo 2024, PLoS One" ~ "Luo 2024a", grepl("Authors 2023", s) ~ "Moreira 5-gene",
  grepl("Th2", s) ~ "Moreira Th2", grepl("APP", s) ~ "Li 2024",
  grepl("autophagy", s) ~ "Ann Med 2024", TRUE ~ s)
# axis assignment: larger absolute correlation with neutrophil or T-cell marker scores; weak if both |r| < 0.5
G <- SIGMAP |> mutate(short = short_of(signature),
  axis = ifelse(pmax(abs(r_neu), abs(r_T), na.rm = TRUE) < .5, "Weak / neither",
    ifelse(abs(r_neu) >= abs(r_T), "Neutrophil-associated", "T-cell-associated")))

# ---- A: external performance ----
dA <- AUD |> mutate(signature = ifelse(signature %in% names(rn), rn[signature], signature),
  short = ifelse(group == "Reference", sub("De novo ", "Reference: ", signature),
    short_of(signature))) |>
  left_join(G |> select(signature, axis),
    by = "signature") |> mutate(axis = ifelse(group == "Reference", "Reference list", axis),
      recon = dplyr::case_when(group == "Reference" ~ "Reference gene list (this study)",
        signature == "Windhorst 2023" ~ "Genes + directions from published classifier",
        TRUE ~ "Genes + directions from derivation data"),
      lab = sprintf("%s%s [%d]", short,
        ifelse(signature %in% c("Tu 2025, Front Pediatr", "Guo 2025"), " #", ""), n_genes),
      lab = factor(lab, levels = lab[order(group != "Reference", auc)]))
fA <- ggplot(dA, aes(auc, lab, colour = axis)) + annotate("rect", xmin = -Inf, xmax = Inf, ymin = .4,
  ymax = 2.5, fill = "#f1f0ec") +
  geom_vline(xintercept = .5, colour = LIGHT, linewidth = .3, linetype = 2) +
  geom_errorbarh(aes(xmin = lo, xmax = pmin(hi, 1)), height = 0,
    linewidth = .45) + geom_point(aes(shape = recon), size = 1.7, fill = "white", stroke = .6) +
  scale_colour_manual(values = AX, guide = "none") +
  scale_shape_manual(values = c("Genes + directions from derivation data" = 16,
    "Genes + directions from published classifier" = 23, "Reference gene list (this study)" = 15),
    name = NULL) +
  scale_x_continuous(limits = c(.05, 1), breaks = seq(.1, .9, .2)) +
  labs(x = "Pooled external AUC (95% CI)", y = NULL, title = "Limited external performance",
    subtitle = "Reconstructed gene scores (not original models),\nscored only in cohorts not used to derive them",
    caption = "# validated by its authors in a cohort used here for testing") +
  TH + theme(axis.text.y = element_text(colour = INK, size = 6), legend.position = "bottom",
    legend.direction = "vertical", legend.margin = margin(-2, 0, 0, 0)) +
  guides(shape = guide_legend(override.aes = list(colour = MUTED)))

# ---- B: immune axes ----
LP <- tribble(~short, ~nx, ~ny, ~hj,
  "Zhang 2025", -.035, .02, 1, "Zou 2024", -.035, 0, 1, "Luo 2024a", -.035, -.01, 1, "Tu 2025", -.035,
  -.02, 1,
  "Sun 2026", -.03, .05, .5, "Guo 2025", .035, -.03, 0, "Fang 2024", .035, .02, 0, "Moreira 5-gene",
  -.02, .055, .5, "Abilbayeva 2026", .035, 0, 0,
  "Luo 2024b", -.035, 0, 1, "Moreira Th2", 0, .05, .5, "Ann Med 2024", -.01, -.055, .5, "Windhorst 2023",
  .035, -.01, 0, "Jia 2023", 0, -.055, .5, "Li 2024", .035, 0, 0)
dB <- G |> filter(!is.na(r_neu)) |> left_join(LP, by = "short") |> mutate(nx = coalesce(nx, 0),
  ny = coalesce(ny, .05), hj = coalesce(hj, .5),
  kept = dplyr::case_when(!is.na(p_GAbal) & p_GAbal < .05 & b_GAbal > 0 ~ "yes, same direction",
    !is.na(p_GAbal) & p_GAbal < .05 ~ "yes, reversed", TRUE ~ "no"))
pcs <- PCA |> filter(!is.na(PC1_var))
pctxt <- sprintf("First two principal\ncomponents: %.0f-%.0f%%\nof variance; explain\n%.0f-%.0f%% of the\nneutrophil-T balance\n(%d test cohorts)",
  100 * min(pcs$PC1_2_cum), 100 * max(pcs$PC1_2_cum), 100 * min(pcs$R2_bal_on_PC12),
  100 * max(pcs$R2_bal_on_PC12), nrow(pcs))
fB <- ggplot(dB, aes(r_neu, r_T)) +
  geom_hline(yintercept = 0, colour = LIGHT, linewidth = .3) + geom_vline(xintercept = 0, colour = LIGHT,
    linewidth = .3) +
  geom_point(aes(fill = axis, shape = kept), size = 2.5, colour = INK, stroke = .35) +
  geom_text(aes(label = short, x = r_neu + nx, y = r_T + ny, hjust = hj), size = 1.95, colour = INK) +
  annotate("text", x = -.43, y = .44, label = pctxt, hjust = 0, vjust = 1, size = 1.95, colour = MUTED,
    lineheight = .95) +
  scale_fill_manual(values = AX, name = "Immune axis") +
  scale_shape_manual(values = c("no" = 21, "yes, same direction" = 24, "yes, reversed" = 25),
    name = "BPD association after\nGA + balance (p < 0.05)") +
  scale_x_continuous(limits = c(-.45, 1.05),
    breaks = seq(-.2, 1, .2)) + scale_y_continuous(limits = c(-1, .45), breaks = seq(-1, .4, .2)) +
  labs(x = "Correlation with neutrophil marker score", y = "Correlation with T-cell marker score",
    title = "Signatures align with two immune-composition axes",
    subtitle = "Pooled across test cohorts; axis = larger |r| (weak, both |r| < 0.5)") +
  TH + theme(legend.position = "right",
    legend.box = "vertical") + guides(fill = guide_legend(override.aes = list(shape = 21, size = 2.4)),
      shape = guide_legend(override.aes = list(size = 2.2), keyheight = unit(10, "pt")))

# ---- C: adjustment ----
dC <- dB |> filter(!is.na(b_GA)) |> mutate(axis = factor(axis, levels = names(AX)[1:3]),
  short = reorder(short, b_GAbal))
fC <- ggplot(dC, aes(y = short)) + geom_vline(xintercept = 0, colour = LIGHT, linewidth = .3) +
  geom_segment(aes(x = b_GA, xend = b_GAbal, yend = short), colour = LIGHT, linewidth = .9) +
  geom_point(aes(x = b_GA), shape = 21, fill = "white", colour = MUTED, size = 1.6) +
  geom_errorbarh(aes(xmin = lo_GAbal, xmax = hi_GAbal, colour = axis), height = 0, linewidth = .45) +
  geom_point(aes(x = b_GAbal, colour = axis),
    size = 1.8) + geom_text(aes(x = 1.72, label = star(p_GAbal)), colour = INK, size = 2.4) +
  facet_grid(axis ~ ., scales = "free_y", space = "free_y") + scale_colour_manual(values = AX,
    guide = "none") + scale_x_continuous(limits = c(-1.8, 1.8)) +
  labs(x = "BPD effect on signature score (SD; positive, more BPD-like)", y = NULL,
    title = "Few associations survive adjustment",
    subtitle = "Open, adjusted for GA; filled, GA + balance (95% CI);\n3 test cohorts reporting GA; *p < 0.05, **p < 0.01") +
  TH + theme(axis.text.y = element_text(colour = INK, size = 6),
    strip.text.y = element_text(angle = 0, size = 5.6, hjust = 0, face = "plain", colour = MUTED))

# ---- D: methylation cross-check ----
dD <- METH |> mutate(grp = factor(ifelse(y == 1, "BPD", "No BPD"), levels = c("No BPD", "BPD")))
r <- cor(dD$BAL, dD$mNLR)
fD <- ggplot(dD, aes(mNLR, BAL, colour = grp)) + geom_point(size = 1,
  alpha = .75) + geom_smooth(aes(group = 1), method = "lm", se = FALSE, colour = INK, linewidth = .4,
    formula = y ~ x) +
  annotate("text", x = min(dD$mNLR), y = max(dD$BAL),
    label = sprintf("r = %.2f\n%d samples, %d infants", r, nrow(dD), length(unique(dD$infant))),
    hjust = 0, vjust = 1, size = 2.2) +
  scale_colour_manual(values = c("No BPD" = C_NO, BPD = C_BPD), name = NULL) +
  labs(x = "Methylation-inferred neutrophil-to-lymphocyte ratio (log)",
    y = "RNA neutrophil-T-cell balance (z)", title = "RNA balance agrees with an independent estimate",
    subtitle = "GSE220135, cord blood to day 28; inferred proportions, not cell counts") +
  TH + theme(legend.position = c(.98, .04), legend.justification = c(1, 0),
    legend.background = element_blank())
fig <- ((fA | fB) + plot_layout(widths = c(1, 1.3))) / ((fC | fD) + plot_layout(widths = c(1.05,
  1))) + plot_layout(heights = c(1.15, 1)) + plot_annotation(tag_levels = "A")
savefig(fig, "NewFigure3", 7.4, 7.6)
print(table(G$axis))
print(G |> select(short, axis, r_neu, r_T) |> arrange(axis))
