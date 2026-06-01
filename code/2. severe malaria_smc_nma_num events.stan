// ========================================================================================================
// Contrast-based Bayesian NMA
// - Compares multiple treatments using a binomial likelihood for the number of events
// - Models log-odds of events using a normal distribution
// ========================================================================================================

data {
  int<lower=1> T;                        // number of treatments (including baseline)
  int<lower=1> S;                        // number of studies
  array[S] int<lower=1> A;               // number of arms in each study
  array[S, max(A)] int<lower=0> n;       // number of events (padded to max arms)
  array[S, max(A)] int<lower=0> t;       // total participants (padded to max arms)
  array[S, max(A)] int<lower=0> trt;     // treatment IDs (padded to max arms)
  int<lower=1> bt;                       // baseline treatment
}

parameters {
  vector[T - 1] d_raw;                   // all non-baseline treatment effects
  vector[S] alpha;                       // log-odds of baseline treatment
}

transformed parameters {
  vector[T] d = rep_vector(0, T);        // baseline treatment effect is 0
  matrix[S, max(A)] logit_p = rep_matrix(0, S, max(A));

  {
    int j = 1;
    for (treat in 1:T) {
      if (treat != bt) {
        d[treat] = d_raw[j];
        j += 1;
      }
    }
  }

  for (i in 1:S) {
    for (k in 1:A[i]) {
      
      // linear predictor for arm-specific event probability
      logit_p[i, k] = alpha[i] + d[trt[i, k]];
    }
  }
}

model {
  
  // Weakly informative priors
  // Stabilise estimation while allowing data to dominate
  // Baseline event probability ~0.05
  d_raw ~ normal(0, 0.5);
  alpha ~ normal(-2, 0.5);

  for (i in 1:S) {
    for (k in 1:A[i]) {
      
      // Binomial likelihood using log-odds ratio
      n[i, k] ~ binomial_logit(t[i, k], logit_p[i, k]);
    }
  }
}

generated quantities {
  matrix[T, T] OR;
  matrix[T, T] LOR;
  
  // Generate log-odds ratio between all pairs of non-baseline treatments
  for (i in 1:T) {
    for (j in 1:T) {
      LOR[i, j] = d[j] - d[i];
      OR[i, j]  = exp(LOR[i, j]);
    }
  }
}
