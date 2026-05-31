#SMC Network meta-analysis - tolerability outcomes


# Load packages ####
library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(meta)
library(netmeta)
library(dmetar)
library(gemtc)
library(rjags)



# Load data ####
     #replace 'filename.csv' with relevant files from data/tolerability
      tol_pairwise_input <- read_csv("data/tolerability/filename.csv") 
      
# Pairwise meta-analysis ####

      # To estimate contrast-based treatment effects
      for_nma <- pairwise(
        treat = list(trt1, trt2, trt3, trt4),
        event = list(n1, n2, n3, n4),
        n = list(t1, t2, t3, t4),
        data = tol_pairwise_input, 
        studlab = studyid,
        sm = "RD", # "RR" for relative effects
      )    
      
      
      tol_data <- for_nma %>% 
        select(studlab, treat1, treat2, TE, seTE) %>% 
        rename(
          studyid = studlab,
          logTE = TE)

      
      # Total number of events
      a <- sum(tol_pairwise_input$t1, na.rm = TRUE)
      b <- sum(tol_pairwise_input$t2, na.rm = TRUE)
      c <- sum(tol_pairwise_input$t3, na.rm = TRUE)
      d <- sum(tol_pairwise_input$t4, na.rm = TRUE)
      
      a+b+c+d # Total events


# Frequentist NMA ####
      
  # --- NMA model
      smc_tol_network <- netmeta(
        TE = logTE,
        seTE = seTE,
        treat1 = treat1,
        treat2 = treat2,
        studlab = studyid,
        data = tol_data,
        sm = "RD",  # "RR" for relative effects
        common = TRUE,
        #random = TRUE,
        reference.group = "placebo",
        details.chkmultiarm = TRUE,
        tol.multiarm = 0.5,
        tol.multiarm.se = 1,
        sep.trts = " vs "
      )
      
      summary(smc_tol_network)
      
      
  # --- Plots
    # Network graph
      # Treatment labels
      smc_tol_network$trts
      
      long_labels <- c( #Select treatment represented in the network for each tolerability outcome
        "ASAQ",
        "DHAPQ",
        "Placebo",
        "SP + 1 AS",
        "SP + 3 AS",
        "SPAQ",
        "SPPQ")
      
      netgraph(smc_tol_network, 
               labels = long_labels)
      
    # Direct evidence plot 
      tol_devidence <- direct.evidence.plot(smc_tol_network)
      plot(tol_devidence)
      
    # Forest plot of relative efficacy 
      
      forest(
        smc_tol_network,
        reference.group = "placebo",
        labels = long_labels,
        fontsize = 11,
        pooled = "common",
        sortvar = -TE,
        smlab = "Antimalarials vs placebo\n('OUTCOME', n='TOTAL EVENTS')", #Insert outcome and total events
        #xlim = c(0.05, 10),
        cex = 0.5,
        rightlabs = c("Risk\n difference", "95% CI"), #Risk ratio for relative effects
        overall.hetstat = TRUE
      )
      
      
      
  # --- Inconsistency analysis: Node splitting method
      
      #Node-splitting
      sm_netsplit <- netsplit(smc_sm_network)
      
      #Inconsistency forest plot
      netsplit_placebo <- netsplit(
        smc_sm_network,
        method = "SIDDE",
        reference.group = "placebo/nodrug",
        show = "all",
        baseline.reference = TRUE,
        common = TRUE,
        backtransf = TRUE,
        ci = TRUE)
      
      
      forest(netsplit_placebo)    
      
      
# Bayesian NMA ####

 # --- Shape to GEMTC format

  tol_data_long <- tol_data %>% 
    pivot_longer( #Reshape from wide to long
      cols = c("treat1", "treat2", "logTE", "seTE"), # Select columns to reshape
      names_to = c(".value"), # Use column name parts as output column names 
      names_pattern = "(..)"  # Group column names by first two columns (trt1 & trt2 == grouped into treatments)
    ) %>% 
    rename( #Rename columns
      study = studyid,
      diff = lo,
      std.err = se,
      treatment = tr
    ) %>% 
    arrange(study)



 # --- Multi-arm cleaning 
      
      #Drop redundant multi-arm comparisons (not needed by GEMTC)
    sm_data_long$row_id <- seq_len(nrow(sm_data_long))
    drop_rows <- c(1, 2, 
                   21, 22, 
                   23, 24, 
                   25, 26)
    
          sm_data_long <- sm_data_long %>%
            filter(!(row_id %in% drop_rows)) %>%
            select(-row_id)
    
    sm_data_long <- sm_data_long %>%      # Remove duplicated baseline rows
      group_by(study, treatment) %>%
      filter(!(is.na(diff) & duplicated(treatment))) %>%
      ungroup()
  
    # Standard error in baseline arm for multi-arm studies (required by GEMTC)
    baseline_SE <- tibble(     #Calculated as SE(λ) = 1/sqrt(y); where y = No. events (counts) in arm 
      study = c("Chandramohan2021",
                "Kweku2008"),
      baseline_se = c(0.3015113,
                      0.2294157)     
    )         
      
      #Merge baseline standard errors to main dataset
      sm_data_long <- sm_data_long %>%
        left_join(baseline_SE, by = "study") %>%
        mutate(
          std.err = if_else(is.na(diff) & !is.na(baseline_se), baseline_se, std.err)
        ) %>%
        select(-baseline_se)
  
      
  # --- Treatment labels for plots
      sm_data_long <- sm_data_long %>%
          mutate(treatment = recode(treatment,
                               "placebo/nodrug" = "placebo_nodrug",
                               "spaq+rtss"      = "spaq_rtss"))

      treat.codes <- c(
        "spaq"           = "SPAQ",
        "placebo_nodrug" = "Placebo/No drug control",
        "rtss"           = "RTS,S vaccine",
        "spaq_rtss"      = "SPAQ + RTS,S Vaccine",
        "asaq"           = "ASAQ",
        "sp"             = "SP (bimonthly)",
        "asaq_bi"        = "ASAQ (bimonthly)") %>%
        data.frame(description = ., stringsAsFactors = FALSE) %>%
        rownames_to_column("id")




  # --- NMA 
  
    # Network
    sm_network <- mtc.network(
    data.re    = sm_data_long,
    treatments = treat.codes
  )

    summary(sm_network, use.description = TRUE)

    # NMA model
    sm_model <- mtc.model(
      sm_network,
      likelihood  = "normal",
      link        = "identity",
      linearModel = "fixed",
      n.chain     = 4
    )

    mcmc  <- mtc.run(sm_model, 
                     n.adapt = 5000, 
                     n.iter = 10000, 
                     thin = 10)

    # Model diagnostics
    plot(mcmc)
    gelman.plot(mcmc)
    gelman.diag(mcmc)$mpsrf
    
    summary(mcmc)



 # --- Network graph

  plot(
    sm_network,
    use.description = TRUE,
    vertex.color = "yellow",
    vertex.label.color = "black",
    vertex.label.family = "Arial",
    vertex.label.dist = 2.1,
    vertex.label.cex = 0.8
  )


 # --- Relative effects + forest plot

    #Summary of relative effects vs placebo
    results_vs_placebo <- relative.effect(mcmc, t1 = "placebo_nodrug")
    summary(results_vs_placebo)
    
    
    # Extract results
    output <- summary(results_vs_placebo) 
    
    irr_re_plac <- as.data.frame(output[["summaries"]][["statistics"]]) #Extract logIRRs
    cri_re_plac <- as.data.frame(output[["summaries"]][["quantiles"]]) #Extract credible intervals
    
    re_plac <- cbind(irr_re_plac, cri_re_plac)
    re_plac$comparison <- rownames(re_plac)
    
    
    # Forest plot
    forest(results_vs_placebo, use.description = TRUE)


# --- Treatment ranking
  rank_sm <- rank.probability(mcmc_ce_main, preferredDirection = -1)
  plot(rank_sm, beside = TRUE)



















