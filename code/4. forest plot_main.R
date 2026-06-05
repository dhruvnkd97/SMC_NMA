# Forest plot (main) for NMA


# Load packages ####
library(dplyr)
library(forestplot)
library(grid)
library(svglite)



# Load and clean data ####
nma_forest_main <- read.csv("plots/forest_main.csv")

    # --- Order by effect size
    nma_forest_main <- nma_forest_main[order(nma_forest_main$est, na.last = TRUE), ] 


    # --- Label p-values column 
    fmt_p <- function(p) ifelse(is.na(p), "",
                                ifelse(p < 0.001, "p<0.001", sprintf("p=%.3f", p)))


    # --- Weights for box sizes
    rng <- range(nma_forest_main$weight) #range of weights

    boxsz <- if (diff(rng) == 0) {
      rep(0.35, nrow(nma_forest_main))  
    } else {
      0.15 + 0.45 * (nma_forest_main$weight - rng[1]) / diff(rng)  # Rescale all ranges from 0.15–0.60
    }

    # --- Plot labels
    labeltext <- rbind(c("Antimalarials vs Placebo/No drug", 
                         "IRR [95% CrI]"),
                cbind(Label = nma_forest_main$label,
                `IRR [95% CrI]` = sprintf("%.2f (%.2f, %.2f)", 
                                nma_forest_main$est, 
                                nma_forest_main$lo, 
                                nma_forest_main$hi)
                                ))

    mean  <- c(NA, nma_forest_main$est)
    lower <- c(NA, nma_forest_main$lo)
    upper <- c(NA, nma_forest_main$hi)


    is_sum <- c(TRUE, rep(FALSE, nrow(nma_forest_main)))


    # --- Axes 
    # X-axis (log-scale)
       ticks_all <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10)
    
        # Data range for x-axis
        xmin <- min(nma_forest_main$lo)
        xmax <- max(nma_forest_main$hi)
    
        clip_left  <- max(0.05, xmin / 1.5)
        first_right_tick <- min(ticks_all[ticks_all > 1]) 
        clip_right <- max(xmax * 1.5, first_right_tick)
    
        ticks <- ticks_all[ticks_all >= clip_left & ticks_all <= clip_right]


# Forest plot ####
  svglite("forest_main.svg", width = 10, height = 5.625) 
  
  forestplot(
    labeltext = labeltext,
    mean  = mean, lower = lower, upper = upper,
    is.summary = is_sum,
    graph.pos = 2,                        
    xlog = TRUE, 
    zero = 1, 
    clip = c(clip_left, clip_right), 
    xticks = ticks,
    col = fpColors(box = "grey45", line = "black", summary = "black"),
    txt_gp = fpTxtGp(
      label   = gpar(cex = 1.00, fontfamily = "Arial"),
      ticks   = gpar(cex = 0.90, fontfamily = "Arial"),
      xlab    = gpar(cex = 1.00, fontfamily = "Arial"),
      summary = gpar(fontface = "bold", fontfamily = "Arial")),
    align = c("l","l","r"),  # align each label text column
    hrzl_lines = list("2" = gpar(lwd = 1, col = "black")),  # rule under header
    colgap = unit(6, "mm"),
    xlab = "Incidence rate ratio",
  )
  
  
  dev.off()


