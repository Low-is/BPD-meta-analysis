# Shared figure style: navy for BPD, terracotta for no BPD and red/blue for effect direction.
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})
FIG <- "./figs"
FONT <- "Liberation Sans"
INK <- "#222222"
MUTED <- "#5f5f5f"
LIGHT <- "#c9c7c0"
C_BPD <- "#1f3a6e"
C_NO <- "#c8553d"
C_GA <- "#b8860b"
C_NAIVE <- "#1d9fb0"
C_BAL <- "#7a5195"
C_PUB <- "#8a8a82"
DIV <- c("#a63a2a", "#e8a493", "#f5f5f5", "#9fb3d6", "#1f3a6e")
update_geom_defaults("text", list(family = FONT, colour = INK, size = 2.2))
TH <- theme_classic(base_size = 7.5, base_family = FONT) +
  theme(text = element_text(colour = INK), axis.text = element_text(colour = MUTED, size = 6.4),
    axis.title = element_text(size = 6.9),
    axis.line = element_line(colour = MUTED, linewidth = .3),
    axis.ticks = element_line(colour = MUTED, linewidth = .3),
    plot.title = element_text(size = 7.8, face = "bold", hjust = 0, margin = margin(b = 3)),
    plot.subtitle = element_text(size = 6.4, colour = MUTED, margin = margin(b = 3)),
    plot.caption = element_text(size = 6, colour = MUTED, hjust = 0),
    legend.text = element_text(size = 6.4), legend.title = element_text(size = 6.4),
    legend.key.size = unit(7, "pt"),
    strip.background = element_blank(), strip.text = element_text(size = 6.6, face = "bold", hjust = 0),
    plot.margin = margin(4, 6, 4, 4), plot.tag = element_text(size = 10, face = "bold", family = FONT))

star <- function(p) ifelse(is.na(p), "",
  ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))

savefig <- function(fig, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), fig, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(FIG, paste0(name, ".png")), fig, width = w, height = h, dpi = 300, bg = "white")
}
