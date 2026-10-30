# Figure 2: the four-cohort multicohort signature and programs, in-sample vs held-out, and adjustment
src_until <- function(file, from_pat = NULL, until_pat = "^fig <-") {
  l <- readLines(file)
  i0 <- if (is.null(from_pat)) 1 else grep(from_pat, l)[1]
  i1 <- grep(until_pat, l)[1] - 1
  eval(parse(text = l[i0:i1]), envir = globalenv())
}
src_until("./fig2.R")
g2A <- fA
g2B <- fB
load(file.path(OUT, "fig3.RData"))
SC4 <- c(P4, "Combined score")
src_until("./11_fig3.R", from_pat = "^FONT <-")
g3A <- fA
g3D <- fD
source("./figtheme.R")
g2A <- g2A + labs(title = "Conserved pathways group into four programs")
g2B <- g2B + labs(title = "Program scores by cohort", subtitle = "Hedges' g, in-sample")
g3A <- g3A + labs(title = "Held-out validation",
  subtitle = "Open, in-sample; filled, held-out (95% CI);\ndiamond, fixed gene sets")
g3D <- g3D + labs(title = "Most of the signal is GA and cell composition") +
  scale_shape_manual(values = c("GA" = 21, "GA + cell composition" = 19),
    name = NULL)  # open, GA; filled, GA + composition
fig <- ((g2A | g2B) + plot_layout(widths = c(1.5, 1))) / ((g3A | g3D) + plot_layout(widths = c(1,
  1.35))) + plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "A")
savefig(fig, "NewFigure2", 7.4, 6.6)
