# Heat plot for treatment ranking


# Load packages ####
library(ggplot2)
library(tidyr)
library(dplyr)
library(viridis)
library(RColorBrewer)

# Lead and clean data
  
  # --- Treatment labels
  treatments_heatplot <- c( 
                "ASAQ", 
                "ASAQ (bimonthly)", 
                "DHAPQ", 
                "Placebo/No drug",
                "Seasonal RTS,S", 
                "SP (bimonthly)", 
                "SPAQ",
                "SPAQ + Seasonal RTS,S",
                "SPAS/SP (monthly)",
                "SPPQ"
                )


  # --- Data
    mat_heatplot <- matrix(c(
      0.00000, 0.00000, 0.01450, 0.06825, 0.12425, 0.33675, 0.45550, 0.00075, 0.00000, 0.00000,
      0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00025, 0.00000, 0.05175, 0.94750, 0.00050,
      0.00000, 0.38600, 0.61050, 0.00325, 0.00025, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
      0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00125, 0.99875,
      0.00000, 0.00050, 0.09650, 0.36500, 0.32775, 0.13725, 0.07300, 0.00000, 0.00000, 0.00000,
      0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.01725, 0.93775, 0.04500, 0.00000,
      0.00000, 0.00000, 0.10125, 0.47475, 0.36075, 0.05650, 0.00675, 0.00000, 0.00000, 0.00000,
      0.90375, 0.09625, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00000,
      0.00000, 0.00025, 0.01375, 0.05850, 0.14775, 0.43975, 0.33975, 0.00000, 0.00025, 0.00000,
      0.09625, 0.51700, 0.16350, 0.03025, 0.03925, 0.02950, 0.10775, 0.00975, 0.00600, 0.00075
    ), nrow = length(treatments_heatplot), byrow = TRUE)
    
    colnames(mat_heatplot) <- paste0("Rank", 1:10)
    rownames(mat_heatplot) <- treatments_heatplot

  # ---  Ordering treatments by Rank1..Rank10 l (highest first) 
    order_idx <- do.call(order, 
                         as.data.frame(-mat_heatplot[, paste0("Rank", 1:10)]))
    
    ordered_treatments <- rownames(mat_heatplot)[order_idx]
    
    plot_levels <- rev(ordered_treatments)  # reverse so highest appears at top in ggplot

  # --- Build long data frame
    df <- as.data.frame(mat_heatplot) %>%
      mutate(Treatment = rownames(.)) %>%
      pivot_longer(cols = starts_with("Rank"), 
                   names_to = "Rank", 
                   values_to = "Probability") %>%
      mutate(Treatment = factor(Treatment, 
                                levels = plot_levels),
             Rank = factor(Rank, 
                           levels = paste0("Rank", 1:10)))
    
  # --- Labels for probabilities
    fmt_label <- function(x) {
      case_when(
        x == 0 ~ "0",
        x < 0.001 ~ formatC(x, format = "e", digits = 2),
        TRUE ~ sprintf("%.3f", x)
      )
    }
    df$Label <- fmt_label(df$Probability)
    
    df <- df %>% mutate(p_vis = sqrt(Probability))

  # --- Warm colours
    warm_cols <- c("#ffffff", 
                   "#fff7ec", 
                   "#fee6ce", 
                   "#fdd49e", 
                   "#fdbb84", 
                   "#fc8d59", 
                   "#ef6548", 
                   "#d7301f")

  # --- Legend tick locations in original probability scale (choose small ones to reveal detail)
    legend_breaks_orig <- c(0, 1e-5, 1e-4, 1e-3, 0.01, 0.05, 0.1, 0.5, 1)
    legend_breaks_trans <- sqrt(legend_breaks_orig) 
    
    legend_labels <- sapply(legend_breaks_orig, function(x) {
      if (x == 0) "0" else if (x < 0.001) formatC(x, format="e", digits=2) else sprintf("%.3f", x)
    })
    
# Forest plot ####  
    svglite("plots/rankheatplot_uncomp malaria.svg", width = 10, height = 5.625) 
    
    p_vis <- ggplot(df, aes(x = Rank, y = Treatment, fill = p_vis)) +
      geom_tile(color = "#f0f0f0", 
                linewidth = 0.25) +
      geom_text(data = df %>% filter(p_vis > 0), aes(label = Label), color = "#222222", size = 3) +
      scale_fill_gradientn(colors = warm_cols, limits = c(0, 1), na.value = "white") +
      scale_x_discrete(labels = as.character(1:10)) +
      theme_minimal(base_size = 12) +
      theme(
        panel.background = element_rect(fill = "white", colour = NA),
        plot.background  = element_rect(fill = "white", colour = NA),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, face = "plain", size = 11),
        axis.text.y = element_text(face = "plain", size = 11),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        axis.ticks = element_blank(),
        legend.position = "none"
      )
    
    print(p_vis)
    
    dev.off()




















