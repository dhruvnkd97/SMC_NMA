# Forest plot (split) for NMA


# Load packages ####
library(dplyr)
library(forestplot)
library(grid)
library(svglite)
library(tools)

# Load and clean data ####
forest_split_plot <-  read.csv("forest plot_split_dp.csv")

    # --- Label p-values column 
    fmt_p <- function(p) ifelse(is.na(p), "",
                                ifelse(p < 0.001, "p<0.001", 
                                       sprintf("p=%.3f", p))) #Labelling p values



    forest_split_plot <- forest_split_plot %>%
      mutate(
        label_source = toTitleCase(as.character(source)),
        label_ci     = sprintf("%.2f (%.2f, %.2f)", est, lo, hi),
        label_p      = ifelse(source == "network", fmt_p(p_incoherence), "")
      )




    # --- Order rows: Direct → Indirect → Network
    src_order <- c("direct","indirect","network")
    forest_split_plot <- forest_split_plot %>% 
      mutate(source = factor(source, levels = src_order))  


      # --- Insert bold comparison headers and mark Network as summary (diamond) ----------
      rows <- list()
    
      for (cmp in unique(forest_split_plot$comparison)) {
        
        g <- forest_split_plot %>% filter(comparison == cmp) %>% arrange(source)
        
        # Header row (bold, no estimate)
        header <- tibble(
          comparison = cmp,
          source = factor(NA, levels = src_order),
          est = NA_real_, lo = NA_real_, hi = NA_real_,
          label_source = "", label_ci = "", label_p = "",
          is_summary = TRUE,  # bold header
          boxsize = NA_real_
        )
        
        # Body rows (network row drawn as diamond/summary)
        g <- g %>%
          mutate(is_summary = (source == "network"))
        rows[[length(rows)+1]] <- bind_rows(header, g)
      }
      
      tbl <- bind_rows(rows)



      # --- Plot labels
        labeltext <- cbind(
          Comparison = ifelse(is.na(tbl$source), tbl$comparison, ""), 
          Source     = tbl$label_source,
          `IRR [95% CrI]` = tbl$label_ci,
          `p (incoherence)`     = tbl$label_p
        )

      # ---- Data for forest plot ----------
      mean  <- tbl$est
      lower <- tbl$lo
      upper <- tbl$hi
      is.sum <- tbl$is_summary


      # --- Axes settings 
      xmin <- min(lower[!is.na(lower)], na.rm = TRUE)
      xmax <- max(upper[!is.na(upper)], na.rm = TRUE)
      clip <- c(max(0.05, xmin/1.5), xmax*1.5)
      ticks_all <- c(0.1, 0.2, 0.5, 1, 2, 5, 10)
      ticks <- ticks_all[ticks_all >= clip[1] & ticks_all <= clip[2]]


      tbl$row_class <- ifelse(is.na(tbl$source), "header", as.character(tbl$source))
      is.sum <- is.na(tbl$source)

      # Functions to draw values
      fun_list <- list()
      for (i in 1:nrow(tbl)) {
        if (tbl$row_class[i] == "header") {
          # Headers get an empty function
          fun_list[[i]] <- function(mean, lower, upper, size, ...) { NULL }
        } else if (tbl$row_class[i] == "direct") {
          fun_list[[i]] <- function(mean, lower, upper, size, ...) {
            fpDrawNormalCI(mean, lower, upper, size, gp = gpar(fill = "grey65", col = "grey20"), ...)
          }
        } else if (tbl$row_class[i] == "indirect") {
          fun_list[[i]] <- function(mean, lower, upper, size, ...) {
            fpDrawNormalCI(mean, lower, upper, size, gp = gpar(fill = "white", col = "#1f78b4", lwd = 1.5), ...)
          }
        } else if (tbl$row_class[i] == "network") {
          fun_list[[i]] <- function(mean, lower, upper, size, ...) {
            fpDrawNormalCI(mean, lower, upper, size, gp = gpar(fill = "black", col = "black"), ...)
          }
        }
      }
      
      
   #Forest plot ####   
      svglite("forest_split_dp.svg", width = 11, height = 6.5)
      
      styles_matrix <- matrix("plain", nrow = nrow(tbl), ncol = ncol(labeltext))
      styles_matrix[tbl$row_class == "header" | tbl$row_class == "network", ] <- "bold"
      
      forestplot(
        labeltext  = labeltext,
        mean       = mean, 
        lower      = lower, 
        upper      = upper,
        is.summary = is.sum,
        fn.ci_geom = fun_list,
        xlog       = TRUE, 
        zero       = 1, 
     #   clip       = clip_limits,
        xticks     = ticks,
        graph.pos  = 3,
        lwd.ci     = 1.5, 
        lwd.zero   = 1, 
        vertices   = TRUE,  
        txt_gp     = fpTxtGp(
          label   = gpar(cex = 0.90, fontfamily = "Arial"),
          ticks   = gpar(cex = 0.85, fontfamily = "Arial"),
          xlab    = gpar(cex = 0.90, fontfamily = "Arial"),
          summary = gpar(fontface = "bold", fontfamily = "Arial")         
        ),
        lineheight = unit(0.65, "cm"),               
        colgap     = unit(6, "mm"),               
        xlab       = "Incidence Rate Ratio",
        title      = "Incoherence analysis\nDHAPQ vs Placebo/No SMC",
        align      = c("l", "l", "r", "r")
      )
      
      dev.off()
























