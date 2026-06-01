# Contrast-based Bayesian NMA in Stan (rstan)

library(here)
library(dplyr)
library(rstan)
library(stringr)
library(flextable)
library(tibble)
library(scales)
library(bayesplot)
library(ggplot2)
library(truncnorm)
library(tidyr)

set.seed(123)
options(scipen = 999, OutDec = ".", mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

# Paths
data_path <- here("data", "smc_nma_severe malaria_num events.csv")
stan_path <- here("code", "2. severe malaria_smc_nma_num events.stan")
out_dir   <- here("plots")
data_dir  <- here("data")

# Reference treatment (must match CSV)
REF <- "placebo/nodrug"

# ---- Load data ----
smc_data <- read.csv(data_path, stringsAsFactors = FALSE)

# Find columns containing specified prefix
find_cols <- function(df, prefix) {
  cols <- grep(paste0("^", prefix, "[0-9]+$"), names(df), value = TRUE)
  if (!length(cols)) stop("No columns found for prefix: ", prefix)
  cols[order(as.integer(sub(paste0("^", prefix), "", cols)))]
}

event_cols <- find_cols(smc_data, "n")
total_cols <- find_cols(smc_data, "t")
trt_cols   <- find_cols(smc_data, "trt")

# Standardise placebo labels so all aliases map to a single reference treatment
canon_trt <- function(x) {
  x <- trimws(tolower(as.character(x)))
  x[x == ""] <- NA_character_
  x[x %in% c("placebo", "no drug", "nodrug", "placebo/nodrug", "placebo / nodrug")] <- REF
  x
}

# Convert numeric arm columns and standardise treatment labels
smc_data <- smc_data %>%
  mutate(across(all_of(c(event_cols, total_cols)), as.integer),
         across(all_of(trt_cols), canon_trt))

# ---- Index studies ----
studies <- smc_data %>%
  distinct(studyid) %>%
  arrange(studyid) %>%
  mutate(s_id = row_number())

S <- nrow(studies)
s_index <- setNames(studies$s_id, studies$studyid)

smc_data <- smc_data %>%
  mutate(s_id = unname(s_index[studyid])) %>%
  arrange(s_id)

# ---- Prepare treatment list ----
treatments <- sort(unique(unlist(smc_data[, trt_cols], use.names = FALSE)))
treatments <- treatments[!is.na(treatments) & nzchar(treatments)]
treatments <- c(setdiff(treatments, REF), REF)
Tn <- length(treatments)
bt <- Tn # Store baseline treatment

# Build integer matrices for Stan: event counts, totals, and treatment IDs
n <- as.matrix(smc_data[, event_cols, drop = FALSE])
t <- as.matrix(smc_data[, total_cols, drop = FALSE])
trt_raw <- as.matrix(smc_data[, trt_cols, drop = FALSE])
trt_raw <- apply(trt_raw, 2, canon_trt)
trt <- apply(trt_raw, 2, function(x) match(x, treatments))

storage.mode(n) <- storage.mode(t) <- storage.mode(trt) <- "integer"
n[is.na(n)] <- t[is.na(t)] <- trt[is.na(trt)] <- 0L

# Number of arms per study
A <- as.integer(rowSums(!is.na(trt_raw)))

stan_data <- list(T = as.integer(Tn), S = as.integer(S), A = A,
                  n = n, t = t, trt = trt, bt = as.integer(bt))

# ---- Compile and sample from the Bayesian NMA model ----
mod <- stan_model(file = stan_path)

fit <- sampling(object = mod, data = stan_data, chains = 4,
                iter = 2000, warmup = 1000, seed = 123,
                control = list(adapt_delta = 0.99, max_treedepth = 15))

# ---- Convert posterior draws into pairwise OR summaries ----
d_draws <- as.matrix(fit, pars = "d")
OR_draws <- as.matrix(fit, pars = "OR")
LOR_draws <- as.matrix(fit, pars = "LOR")

# Traceplots for treatment effects (d) to assess MCMC convergence
mcmc_trace(d_draws, pars = paste0("d[", 1:(Tn - 1), "]")) +
  xlab("Iteration") + ylab("Log Odds Ratio relative to baseline")

direct_cols <- c(21, 22, 23)

# Compare posterior with prior
for (i in direct_cols) {
  post_plot <- ggplot(LOR_draws, aes(x = LOR_draws[, i])) +
    geom_density(fill = "skyblue", alpha = 0.5) +
    stat_function(fun = dtruncnorm, n = 500, args = list(a = 0, mean = 0, sd = 0.5), linetype = "dashed") +
    geom_vline(xintercept = 0) + theme_minimal() +
    labs(title = "Posterior (shaded) vs Prior (dashed) Density of LOR", x = "Log Odds Ratio", y = "Density") +
    xlim(0, 2) + scale_y_continuous(limits = c(0, 2.5), expand = c(0, 0))
  show(post_plot)
}

# Extract median and 95% credible interval
qs <- function(x) c(med = median(x), l95 = unname(quantile(x, 0.025)),
                    u95 = unname(quantile(x, 0.975)))

# Apply function above to posterior OR draws
summary_raw <- data.frame(t(apply(OR_draws, 2, qs)))

# Extract treatment indices from OR parameter names
parse_idx <- function(nm) {
  t(vapply(nm,
           function(x) as.integer(unlist(regmatches(x, gregexpr("\\d+", x)))),
           integer(2)))
}

idx <- parse_idx(colnames(OR_draws))

# Populate pairwise OR matrices (median, 95% credible interval, standard error)
or_med      <- matrix(NA_real_, Tn, Tn, dimnames = list(treatments, treatments))
or_lo       <- matrix(NA_real_, Tn, Tn, dimnames = list(treatments, treatments))
or_hi       <- matrix(NA_real_, Tn, Tn, dimnames = list(treatments, treatments))
or_se       <- matrix(NA_real_, Tn, Tn, dimnames = list(treatments, treatments))
or_combined <- matrix(NA_character_, Tn, Tn, dimnames = list(treatments, treatments))

or_med[cbind(idx[, 1], idx[, 2])]      <- summary_raw$med
or_lo[cbind(idx[, 1], idx[, 2])]       <- summary_raw$l95
or_hi[cbind(idx[, 1], idx[, 2])]       <- summary_raw$u95
or_se[cbind(idx[, 1], idx[, 2])]       <- (summary_raw$u95 - summary_raw$l95) / 3.92

# Format median and 95% credible interval for table display
fmt_ci <- function(med, lo, hi) sprintf("%.2f (%.2f; %.2f)", med, lo, hi)
or_combined[cbind(idx[, 1], idx[, 2])] <- fmt_ci(summary_raw$med, summary_raw$l95, summary_raw$u95)

# Reorder treatments by the number of pairwise "wins" (OR < 1)
count <- rowSums(or_med < 1, na.rm = TRUE)
indices <- order(count)

# Reorder matrices so treatments are displayed from worst to best
reorder_square <- function(x, idx) x[idx, idx, drop = FALSE]

or_combined <- reorder_square(or_combined, indices)
or_med      <- reorder_square(or_med, indices)
or_lo       <- reorder_square(or_lo, indices)
or_hi       <- reorder_square(or_hi, indices)
or_se       <- reorder_square(or_se, indices)

diag(or_combined) <- treatments[indices]
or_combined[upper.tri(or_combined)] <- "-"

write.csv(or_combined, file.path(out_dir, "league_table_sev malaria_n events.csv"))

# ---- Create flextable for pairwise OR results ----
OR_flextable <- or_combined
treatment_names <- c("SPAQ + Seasonal RTS,S", "SPAQ", "Seasonal RTS,S", "SP (bimonthly)", "Placebo/No drug")
rownames(OR_flextable) <- colnames(OR_flextable) <- treatment_names[indices]
diag(OR_flextable) <- NA

OR_flextable <- flextable(rownames_to_column(as.data.frame(OR_flextable), "Treatment A (row) vs B (column)")) %>%
  theme_box() %>%
  bold(j = 1, bold = TRUE) %>%
  bg(j = 1, bg = "grey80") %>%
  bg(bg = "grey80", part = "header")

for (k in seq_len(Tn)) {
  OR_flextable <- OR_flextable %>%
    bg(i = k, j = k + 1, bg = "grey80") %>%
    compose(i = k, j = k + 1, as_paragraph("NA"), part = "body") %>%
    italic(i = k, j = k + 1, italic = TRUE, part = "body")
}

OR_flextable

# ---- Extract comparisons versus baseline for forest plot ----
ref_idx <- match(REF, rownames(or_med))
other_idx <- setdiff(seq_len(Tn), ref_idx)

severe_malaria_forest <- data.frame(label = treatment_names[other_idx],
                                    est = or_med[ref_idx, other_idx],
                                    lo = or_lo[ref_idx, other_idx],
                                    hi = or_hi[ref_idx, other_idx],
                                    weight = 1 / (((or_hi[ref_idx, other_idx] - or_lo[ref_idx, other_idx]) / 3.92)^2),
                                    type = "study")

write.csv(severe_malaria_forest, file.path(data_dir, "smc_nma_sev malaria_n events_forest.csv"),
          row.names = FALSE)

# --- Compute treatment ranks from posterior draws of d ----
d_rank <- as.matrix(fit, pars = paste0("d[", seq_len(Tn), "]"))
rank_mat <- t(apply(d_rank, 1, rank, ties.method = "first"))

rank_probs <- sapply(seq_len(Tn), function(r) colMeans(rank_mat == r))
mean_rank <- colMeans(rank_mat)

rank_order <- indices
rank_probs <- rank_probs[rank_order, , drop = FALSE]
rownames(rank_probs) <- treatment_names
colnames(rank_probs) <- as.character(1:Tn)
mean_rank <- mean_rank[rank_order]

rank_tbl <- as.data.frame(round(rank_probs, 3)) |>
  rownames_to_column("Treatment")

my_rankogram <- flextable(rank_tbl) |>
  theme_box() |>
  bold(j = 1, bold = TRUE) |>
  colformat_num(j = 2:(Tn + 1), digits = 3, decimal.mark = "·") |>
  add_header_row(values = c("Treatment", "Rank Probability"), colwidths = c(1, Tn)) |>
  align(align = "center", part = "header") |>
  merge_at(i = 1:2, j = 1, part = "header")

colourer <- col_numeric(palette = c("transparent", "red"), domain = c(0, 1))  
my_rankogram <- bg(my_rankogram, bg = colourer, j = 2:(Tn + 1), part = "body")
my_rankogram

# Derive SUCRA (surface under cumulative ranking curve)
cum_probs <- t(apply(rank_probs, 1, cumsum))
sucra <- rowSums(cum_probs[, -ncol(cum_probs), drop = FALSE]) / (Tn - 1)

sucra_df <- data.frame(
  Treatment = rownames(rank_probs),
  SUCRA = round(sucra, 2),
  `Mean rank` = round(mean_rank, 2),
  check.names = FALSE
)

my_sucra <- flextable(sucra_df) |>
  theme_box() |>
  bold(j = 1, bold = TRUE)

my_sucra

# Plot rank probability and cumulative probability curves for each treatment
rank_df <- as.data.frame(rank_probs) |>
  mutate(treatment = rownames(rank_probs)) |>
  pivot_longer(-treatment, names_to = "rank", values_to = "prob") |>
  mutate(rank = as.integer(gsub("\\D", "", rank)))

cum_df <- rank_df |>
  arrange(treatment, rank) |>
  group_by(treatment) |>
  mutate(cum_prob = cumsum(prob)) |>
  ungroup()

plot_prob_curve <- function(my_df, y_col, y_title) {
  ggplot(my_df, aes(x = rank, y = .data[[y_col]])) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ factor(treatment, levels = treatment_names)) +
    scale_x_continuous(breaks = 1:Tn) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Rank", y = y_title) +
    theme_minimal(base_size = 12)
}
  
plot_prob_curve(rank_df, "prob", "Probability")
plot_prob_curve(cum_df, "cum_prob", "Cumulative probability")

# ---- Session info ----
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "sessionInfo_sev malaria_n events.txt"))
