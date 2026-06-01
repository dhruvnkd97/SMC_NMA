library(here)
library(dplyr)
library(forestplot)
library(grid)
options(scipen = 999, OutDec = "·")

my_file <- here("data", "smc_nma_sev malaria_n events_forest.csv")

# ---------- 1) Dataset ----------
nma_forest_main <- read.csv(my_file)
nma_forest_main <- nma_forest_main[order(nma_forest_main$est, na.last = TRUE), ] # Order by effect size

# ---------- 2) Label p-values column ----------
fmt_p <- function(p) ifelse(is.na(p), "",
                            ifelse(p < 0.001, "p<0.001", sprintf("p=%.2f", p)))

# ---------- 3) Weighting box size ----------
rng <- range(nma_forest_main$weight) #range of weights

boxsz <- if (diff(rng) == 0) {
  rep(0.35, nrow(nma_forest_main))                             # Constant box size if all weights equal
} else {
  0.15 + 0.45 * (nma_forest_main$weight - rng[1]) / diff(rng)  # Linear rescale all ranges from 0.15–0.60
}

# ---------- 4) Labels (left/right of the plot) ----------
fmt3 <- function(x) formatC(x, format = "f", digits = 2, decimal.mark = getOption("OutDec"))

labeltext <- rbind(c("Antimalarials vs Placebo/No drug", 
                     "OR [95% CrI]"),
                   cbind(Label = nma_forest_main$label,
                         `OR [95% CrI]` = paste0(fmt3(nma_forest_main$est), " (", 
                                                 fmt3(nma_forest_main$lo), ", ", 
                                                 fmt3(nma_forest_main$hi), ")"
                   )))

mean  <- c(NA, nma_forest_main$est)
lower <- c(NA, nma_forest_main$lo)
upper <- c(NA, nma_forest_main$hi)


is_sum <- c(TRUE, rep(FALSE, nrow(nma_forest_main)))

# ---------- 5) Axes ----------
# 1) Ticks on X-axis
ticks_all <- c(0.1, 0.2, 0.5, 1, 2)

# 2) Compute data range for x-axis
xmin <- min(nma_forest_main$lo)
xmax <- max(nma_forest_main$hi)

clip_left  <- max(0, xmin)
first_right_tick <- min(ticks_all[ticks_all > 0]) 
clip_right <- max(xmax, first_right_tick)

ticks <- ticks_all[ticks_all >= (clip_left - 1) & ticks_all <= (clip_right + 2)]

# ---------- 6) Draw the forest plot ----------

forest_severe_malaria <- forestplot(
  labeltext = labeltext,
  mean = mean, lower = lower, upper = upper,
  is.summary = is_sum,
  graph.pos = 2,                        # put the plot between column 1 & 2
  xlog = TRUE, 
  zero = 1,
  clip = c(0.1, 2), 
  xticks = c(0.1, 0.2, 0.5, 1, 2),
  col = fpColors(box = "grey45", line = "black", summary = "black"),
  txt_gp = fpTxtGp(
    label   = gpar(cex = 1.00, fontfamily = "Arial"),
    ticks   = gpar(cex = 0.90, fontfamily = "Arial"),
    xlab    = gpar(cex = 1.00, fontfamily = "Arial"),
    summary = gpar(fontface = "bold", fontfamily = "Arial")),
  align = c("l", "l", "r"),               # align each label text column
  hrzl_lines = list("2" = gpar(lwd = 1, col = "black")),  # rule under header
  colgap = unit(6, "mm"),
)

show(forest_severe_malaria)

caption <- "Common effect model\nWithin-design heterogeneity: Q = 0·92; p-value = 0·6321; I² = 0%"
grid.newpage()

pushViewport(viewport(
  layout = grid.layout(
    nrow = 2, ncol = 1,
    heights = unit.c(unit(0.86, "npc"), unit(0.14, "npc"))
  )
))

pushViewport(viewport(layout.pos.row = 1))
show(forest_severe_malaria)
upViewport()

pushViewport(viewport(layout.pos.row = 2))
grid.text(
  caption,
  x = unit(4.9, "mm"),
  y = unit(1, "npc"),
  just = c("left", "top"),
  gp = gpar(fontsize = 10, fontfamily = "Arial", lineheight = 1.1)
)
upViewport(2)