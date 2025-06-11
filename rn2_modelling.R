# ---------------------------------------------------------------------------------------
# BAYESIAN MODELING Research Note
# 2: Modelling & Diagnostics (with Radical-Reform)
# ---------------------------------------------------------------------------------------

library(tidyverse)
library(cmdstanr)
library(stringr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(ggdist)
library(forcats)
library(dplyr)

# ----------------------------------------
# LOAD DATA
# ----------------------------------------

setwd("/Users/henrybaker/repositories/bayesian_modeling/research_note")

stan_data <- readRDS("data/stan_data_fit_fast.rds")

if (!dir.exists("plots_qoi")) {
  dir.create("plots_qoi")
}

# ----------------------------------------------------
# SANITY CHECKS ON `stan_data` BEFORE SENDING TO STAN
# ----------------------------------------------------

# 1. Basic existence and required fields
if (!is.list(stan_data)) {
  stop("`stan_data` must be a list.")
}
required_names <- c(
  "N","P","X","R","region_id","Q","party_id",
  "J_opt","N_opt","i_opt","j_opt","y_opt","lower_opt","upper_opt",
  "J_env","N_env","i_env","j_env","y_env","lower_env","upper_env",
  "J_rad","N_rad","i_rad","j_rad","y_rad","lower_rad","upper_rad"
)
missing_names <- setdiff(required_names, names(stan_data))
if (length(missing_names) > 0) {
  stop(paste0("`stan_data` is missing these elements: ", paste(missing_names, collapse = ", ")))
}

# 2. Check N and P
if (!is.numeric(stan_data$N) || length(stan_data$N) != 1 || stan_data$N < 1) {
  stop("`stan_data$N` must be a single positive integer.")
}
if (!is.numeric(stan_data$P) || length(stan_data$P) != 1 || stan_data$P < 1) {
  stop("`stan_data$P` must be a single positive integer.")
}
if (!is.matrix(stan_data$X)) {
  stop("`stan_data$X` must be a matrix.")
}
if (nrow(stan_data$X) != stan_data$N || ncol(stan_data$X) != stan_data$P) {
  stop(paste0("`X` has dimensions ", 
              paste(dim(stan_data$X), collapse = " x "), 
              " but should be ", stan_data$N, " x ", stan_data$P, "."))
}

# 3. Check region_id and party_id
if (!is.integer(stan_data$region_id) && !is.numeric(stan_data$region_id)) {
  stop("`stan_data$region_id` must be an integer or numeric vector.")
}
if (length(stan_data$region_id) != stan_data$N ||
    any(stan_data$region_id < 1) || any(stan_data$region_id > stan_data$R)) {
  stop("`region_id` must be length N and each entry in [1..R].")
}
if (!is.integer(stan_data$party_id) && !is.numeric(stan_data$party_id)) {
  stop("`stan_data$party_id` must be an integer or numeric vector.")
}
if (length(stan_data$party_id) != stan_data$N ||
    any(stan_data$party_id < 1) || any(stan_data$party_id > stan_data$Q)) {
  stop("`party_id` must be length N and each entry in [1..Q].")
}

# 4. “opt” block: lengths, bounds, and no NAs
if (!is.numeric(stan_data$J_opt) || length(stan_data$J_opt) != 1 || stan_data$J_opt < 1) {
  stop("`stan_data$J_opt` must be a single positive integer.")
}
if (!is.numeric(stan_data$N_opt) || length(stan_data$N_opt) != 1 || stan_data$N_opt < 1) {
  stop("`stan_data$N_opt` must be a single positive integer.")
}
if (length(stan_data$i_opt) != stan_data$N_opt ||
    length(stan_data$j_opt) != stan_data$N_opt ||
    length(stan_data$y_opt) != stan_data$N_opt) {
  stop("Lengths of `i_opt`, `j_opt`, `y_opt` must all equal N_opt.")
}
if (any(is.na(stan_data$i_opt)) || any(is.na(stan_data$j_opt)) || any(is.na(stan_data$y_opt))) {
  stop("`i_opt`, `j_opt`, and `y_opt` must not contain NA.")
}
if (any(stan_data$i_opt < 1) || any(stan_data$i_opt > stan_data$N)) {
  stop("All entries of `i_opt` must be in [1..N].")
}
if (any(stan_data$j_opt < 1) || any(stan_data$j_opt > stan_data$J_opt)) {
  stop("All entries of `j_opt` must be in [1..J_opt].")
}
if (length(stan_data$lower_opt) != stan_data$J_opt ||
    length(stan_data$upper_opt) != stan_data$J_opt) {
  stop("`lower_opt` and `upper_opt` must each be length J_opt.")
}
if (any(is.na(stan_data$lower_opt)) || any(is.na(stan_data$upper_opt))) {
  stop("`lower_opt` and `upper_opt` must not contain NA.")
}

# 5. “env” block: lengths, bounds, and no NAs
if (!is.numeric(stan_data$J_env) || length(stan_data$J_env) != 1 || stan_data$J_env < 1) {
  stop("`stan_data$J_env` must be a single positive integer.")
}
if (!is.numeric(stan_data$N_env) || length(stan_data$N_env) != 1 || stan_data$N_env < 1) {
  stop("`stan_data$N_env` must be a single positive integer.")
}
if (length(stan_data$i_env) != stan_data$N_env ||
    length(stan_data$j_env) != stan_data$N_env ||
    length(stan_data$y_env) != stan_data$N_env) {
  stop("Lengths of `i_env`, `j_env`, `y_env` must all equal N_env.")
}
if (any(is.na(stan_data$i_env)) || any(is.na(stan_data$j_env)) || any(is.na(stan_data$y_env))) {
  stop("`i_env`, `j_env`, and `y_env` must not contain NA.")
}
if (any(stan_data$i_env < 1) || any(stan_data$i_env > stan_data$N)) {
  stop("All entries of `i_env` must be in [1..N].")
}
if (any(stan_data$j_env < 1) || any(stan_data$j_env > stan_data$J_env)) {
  stop("All entries of `j_env` must be in [1..J_env].")
}
if (length(stan_data$lower_env) != stan_data$J_env ||
    length(stan_data$upper_env) != stan_data$J_env) {
  stop("`lower_env` and `upper_env` must each be length J_env.")
}
if (any(is.na(stan_data$lower_env)) || any(is.na(stan_data$upper_env))) {
  stop("`lower_env` and `upper_env` must not contain NA.")
}

# 6. “rad” block: lengths, bounds, and no NAs
if (!is.numeric(stan_data$J_rad) || length(stan_data$J_rad) != 1 || stan_data$J_rad < 1) {
  stop("`stan_data$J_rad` must be a single positive integer.")
}
if (!is.numeric(stan_data$N_rad) || length(stan_data$N_rad) != 1 || stan_data$N_rad < 1) {
  stop("`stan_data$N_rad` must be a single positive integer.")
}
if (length(stan_data$i_rad) != stan_data$N_rad ||
    length(stan_data$j_rad) != stan_data$N_rad ||
    length(stan_data$y_rad) != stan_data$N_rad) {
  stop("Lengths of `i_rad`, `j_rad`, `y_rad` must all equal N_rad.")
}
if (any(is.na(stan_data$i_rad)) || any(is.na(stan_data$j_rad)) || any(is.na(stan_data$y_rad))) {
  stop("`i_rad`, `j_rad`, and `y_rad` must not contain NA.")
}
if (any(stan_data$i_rad < 1) || any(stan_data$i_rad > stan_data$N)) {
  stop("All entries of `i_rad` must be in [1..N].")
}
if (any(stan_data$j_rad < 1) || any(stan_data$j_rad > stan_data$J_rad)) {
  stop("All entries of `j_rad` must be in [1..J_rad].")
}
if (length(stan_data$lower_rad) != stan_data$J_rad ||
    length(stan_data$upper_rad) != stan_data$J_rad) {
  stop("`lower_rad` and `upper_rad` must each be length J_rad.")
}
if (any(is.na(stan_data$lower_rad)) || any(is.na(stan_data$upper_rad))) {
  stop("`lower_rad` and `upper_rad` must not contain NA.")
}

# 7. Check that R and Q are positive integers
if (!is.numeric(stan_data$R) || length(stan_data$R) != 1 || stan_data$R < 1) {
  stop("`stan_data$R` must be a single positive integer.")
}
if (!is.numeric(stan_data$Q) || length(stan_data$Q) != 1 || stan_data$Q < 1) {
  stop("`stan_data$Q` must be a single positive integer.")
}

# 8. Everything passed—signal success
message("✅ All sanity checks passed! `stan_data` is consistent.")

# ----------------------------------------
# FULL Stan MODEL CODE 
# ----------------------------------------

# NB: I DROPPED THE TRUNCATION
# IF HAVE TIME: TRY ORDINAL APPROPRIATE MODEL 
# IE NOT WITH NORMAL PRIOR

stan_code_2 <- "
data {
  int<lower=1> N;                // # of individuals
  int<lower=1> P;                // # of numeric covariates
  matrix[N, P] X;                // standardized covariates

  int<lower=1> R;                           // # of regions
  array[N] int<lower=1, upper=R> region_id; // region index

  int<lower=1> Q;                          // # of parties
  array[N] int<lower=1, upper=Q> party_id; // party index

  int<lower=1> J_opt;                          // # of optimism items
  int<lower=1> N_opt;                          // total optimism responses
  array[N_opt] int<lower=1, upper=N>    i_opt;  // respondent index
  array[N_opt] int<lower=1, upper=J_opt> j_opt; // item index
  array[N_opt] real                     y_opt; // observed (standardized)

  int<lower=1> J_env;                          // # of environment items
  int<lower=1> N_env;                          // total environment responses
  array[N_env] int<lower=1, upper=N>    i_env;  // respondent index
  array[N_env] int<lower=1, upper=J_env> j_env; // item index
  array[N_env] real                     y_env; // observed (standardized)

  int<lower=1> J_rad;                          // # of radical-reform items
  int<lower=1> N_rad;                          // total radical-reform responses
  array[N_rad] int<lower=1, upper=N>    i_rad;  // respondent index
  array[N_rad] int<lower=1, upper=J_rad> j_rad; // item index
  array[N_rad] real                     y_rad; // observed (standardized)
}

parameters {
  // A) Latent-factor hierarchy (3-dimensional now: φ, θ, ψ)
  cholesky_factor_corr[3] Lcorr_eta;  // corr(φ, θ, ψ)
  vector<lower=0>[3]      tau_eta;    // half-Normal(0, 0.3)
  matrix[3, N]            z_eta;      // non-centered

  // 2) Region intercepts (3-dimensional)
  cholesky_factor_corr[3] Lcorr_alpha; // corr across (α₁, α₂, α₃)
  matrix[3, R]            z_alpha;     // non-centered
  vector<lower=0>[3]      sigma_alpha; // SD ≥ 0

  // 3) Party intercepts (3-dimensional)
  cholesky_factor_corr[3] Lcorr_delta; // corr across (δ₁, δ₂, δ₃)
  matrix[3, Q]            z_delta;     // non-centered
  vector<lower=0>[3]      sigma_delta; // SD ≥ 0

  // 4) Covariate slopes (3 × P)
  matrix[3, P]            B; // Normal(0, 0.5)

  // B) Measurement: optimism items
  vector[J_opt]           beta_opt;    // intercepts
  vector<lower=0>[J_opt]  lambda_opt;  // loadings ≥ 0
  vector<lower=0>[J_opt]  sigma_opt;   // residual SD ≥ 0

  // C) Measurement: environment items
  vector[J_env]           beta_env;    // intercepts
  vector<lower=0>[J_env]  lambda_env;  // loadings ≥ 0
  vector<lower=0>[J_env]  sigma_env;   // residual SD ≥ 0

  // D) Measurement: radical-reform items
  vector[J_rad]           beta_rad;    // intercepts
  vector<lower=0>[J_rad]  lambda_rad;  // loadings ≥ 0
  vector<lower=0>[J_rad]  sigma_rad;   // residual SD ≥ 0
}

transformed parameters {
  // Expand region & party covariance matrices
  cov_matrix[3] Sigma_alpha =
    diag_pre_multiply(sigma_alpha, Lcorr_alpha)
    * diag_pre_multiply(sigma_alpha, Lcorr_alpha)';
  cov_matrix[3] Sigma_delta =
    diag_pre_multiply(sigma_delta, Lcorr_delta)
    * diag_pre_multiply(sigma_delta, Lcorr_delta)';

  // Vectors for each latent dimension
  vector[N] phi;
  vector[N] theta;
  vector[N] psi;

  for (i in 1:N) {
    // region + party intercept contributions (3-vector)
    vector[3] mu_eta_i =
      diag_pre_multiply(sigma_alpha, Lcorr_alpha) * z_alpha[, region_id[i]] +
      diag_pre_multiply(sigma_delta, Lcorr_delta) * z_delta[, party_id[i]] +
      B * to_vector(X[i]);

    // latent noise (3-vector)
    vector[3] noise =
      diag_pre_multiply(tau_eta, Lcorr_eta) * z_eta[, i];

    vector[3] eta_i = mu_eta_i + noise;
    phi[i]   = eta_i[1];
    theta[i] = eta_i[2];
    psi[i]   = eta_i[3];
  }
}

model {
  // 1) Priors on τ_eta (latent SDs)
  // Optimism (φ) and Environment (θ): existing tuned priors
  tau_eta[1] ~ normal(0, 0.3);
  tau_eta[2] ~ normal(0, 0.3);
  
  // Radical Reform (ψ): tighter prior
  tau_eta[3] ~ normal(0, 0.1);  

  Lcorr_eta        ~ lkj_corr_cholesky(2.0);
  to_vector(z_eta) ~ normal(0, 1);

  // 2) Region intercept priors
  Lcorr_alpha      ~ lkj_corr_cholesky(2.0);
  sigma_alpha      ~ normal(0, 0.1) T[0, ];
  to_vector(z_alpha) ~ normal(0, 1);

  // 3) Party intercept priors
  Lcorr_delta      ~ lkj_corr_cholesky(2.0);
  sigma_delta      ~ normal(0, 0.1) T[0, ];
  to_vector(z_delta) ~ normal(0, 1);

  // 4) Covariate slopes
  to_vector(B)     ~ normal(0, 0.5);

  // 5) Measurement: optimism
  beta_opt        ~ normal(0, 0.5);
  lambda_opt      ~ lognormal(log(1), 0.2);
  sigma_opt       ~ normal(1, 0.2);

  // 6) Measurement: environment
  beta_env        ~ normal(0, 0.5);
  lambda_env      ~ lognormal(log(1), 0.2);
  sigma_env       ~ normal(1, 0.2);

  // 7) Measurement: radical-reform
  beta_rad        ~ normal(0, 0.5);
  lambda_rad      ~ lognormal(log(1), 0.2);
  sigma_rad       ~ normal(1, 0.2);

  // 8) Likelihood: y_opt
  for (n in 1:N_opt) {
    int ii = i_opt[n];
    int jj = j_opt[n];
    real mu_opt = beta_opt[jj] + lambda_opt[jj] * phi[ii];
    y_opt[n] ~ normal(mu_opt, sigma_opt[jj]);
  }

  // 9) Likelihood: y_env
  for (n in 1:N_env) {
    int ii = i_env[n];
    int jj = j_env[n];
    real mu_env = beta_env[jj] + lambda_env[jj] * theta[ii];
    y_env[n] ~ normal(mu_env, sigma_env[jj]);
  }

  // 10) Likelihood: y_rad
  for (n in 1:N_rad) {
    int ii = i_rad[n];
    int jj = j_rad[n];
    real mu_rad = beta_rad[jj] + lambda_rad[jj] * psi[ii];
    y_rad[n] ~ normal(mu_rad, sigma_rad[jj]);
  }
}

generated quantities {
  // A) Posterior-predictive y's for each block
  vector[N_opt] y_opt_sim;
  for (n in 1:N_opt) {
    int ii = i_opt[n];
    int jj = j_opt[n];
    real mu_opt = beta_opt[jj] + lambda_opt[jj] * phi[ii];
    y_opt_sim[n] = normal_rng(mu_opt, sigma_opt[jj]);
  }

  vector[N_env] y_env_sim;
  for (n in 1:N_env) {
    int ii = i_env[n];
    int jj = j_env[n];
    real mu_env = beta_env[jj] + lambda_env[jj] * theta[ii];
    y_env_sim[n] = normal_rng(mu_env, sigma_env[jj]);
  }

  vector[N_rad] y_rad_sim;
  for (n in 1:N_rad) {
    int ii = i_rad[n];
    int jj = j_rad[n];
    real mu_rad = beta_rad[jj] + lambda_rad[jj] * psi[ii];
    y_rad_sim[n] = normal_rng(mu_rad, sigma_rad[jj]);
  }

  // B) R² computations for each latent block
  // 1) R²_opt
  vector[N_opt] yhat_opt_resp;       
  vector[N_opt] var_resid_opt_resp;  

  for (n in 1:N_opt) {
    int ii = i_opt[n];
    int jj = j_opt[n];
    real mu_opt_n = beta_opt[jj] + lambda_opt[jj] * phi[ii];
    yhat_opt_resp[n]      = mu_opt_n;
    var_resid_opt_resp[n] = square(sigma_opt[jj]);
  }

  real Var_pred_opt = variance(yhat_opt_resp);
  real E_resid_opt  = mean(var_resid_opt_resp);
  real R2_opt       = Var_pred_opt / (Var_pred_opt + E_resid_opt);

  // 2) R²_env
  vector[N_env] yhat_env_resp;
  vector[N_env] var_resid_env_resp;

  for (n in 1:N_env) {
    int ii = i_env[n];
    int jj = j_env[n];
    real mu_env_n = beta_env[jj] + lambda_env[jj] * theta[ii];
    yhat_env_resp[n]      = mu_env_n;
    var_resid_env_resp[n] = square(sigma_env[jj]);
  }

  real Var_pred_env = variance(yhat_env_resp);
  real E_resid_env  = mean(var_resid_env_resp);
  real R2_env       = Var_pred_env / (Var_pred_env + E_resid_env);

  // 3) R²_rad
  vector[N_rad] yhat_rad_resp;
  vector[N_rad] var_resid_rad_resp;

  for (n in 1:N_rad) {
    int ii = i_rad[n];
    int jj = j_rad[n];
    real mu_rad_n = beta_rad[jj] + lambda_rad[jj] * psi[ii];
    yhat_rad_resp[n]      = mu_rad_n;
    var_resid_rad_resp[n] = square(sigma_rad[jj]);
  }

  real Var_pred_rad = variance(yhat_rad_resp);
  real E_resid_rad  = mean(var_resid_rad_resp);
  real R2_rad       = Var_pred_rad / (Var_pred_rad + E_resid_rad);

  // === (B) Generate full 3×R alpha_region and 3×Q delta_party ===
  // NB: in the end I did post-estimation in R, so commending this out
    
  // matrix[3, R] alpha_region;
  // matrix[3, Q] delta_party;
  
  // {
  //  matrix[3, 3] AlphaFactor = diag_pre_multiply(sigma_alpha, Lcorr_alpha);
  //  matrix[3, 3] DeltaFactor = diag_pre_multiply(sigma_delta, Lcorr_delta);
    
  //  for (r in 1:R) {
  //    vector[3] tmp_alpha = AlphaFactor * z_alpha[, r];
  //    for (d in 1:3) {
  //      alpha_region[d, r] = tmp_alpha[d];
  //    }
  //  }
    
  //  for (q in 1:Q) {
  //    vector[3] tmp_delta = DeltaFactor * z_delta[, q];
  //    for (d in 1:3) {
  //      delta_party[d, q] = tmp_delta[d];
  //    }
  //  }
  // }
}
"

if (!dir.exists("model_2")) {
  dir.create("model_2", recursive = TRUE)
}

# Write Stan model file
write_lines(stan_code_2, "model_2/model2.stan")

# Compile the model 
model_full_2 <- cmdstan_model("model_2/model2.stan")

# ----------------------------------------
# PRIOR-PREDICTIVE SAMPLING
# ----------------------------------------

fit_prior <- model_full_2$sample(
  data                 = stan_data,
  chains               = 4,
  parallel_chains      = 4,
  iter_warmup          = 2000,
  iter_sampling        = 1000,
  refresh              = 50,
  fixed_param          = TRUE, 
  save_cmdstan_config  = TRUE
)

fit_prior$save_object("model_2/fit_prior_predictive.rds")

#fit_prior <- readRDS("model_2/fit_prior_predictive.rds")

# ----------------------------------------
# FULL MODEL SAMPLING
# ----------------------------------------

fit_full <- model_full_2$sample(
  data                 = stan_data,
  chains               = 4,
  parallel_chains      = 4,
  iter_warmup          = 2000,
  iter_sampling        = 1000,
  refresh              = 50,
  save_cmdstan_config  = TRUE,
  init                 = 0.5,
)

fit_full$save_object("model_2/fit_full.rds")
saveRDS(fit_full, "model_2/fit_full_inference.rds")

#fit_full <- readRDS("model_2/fit_full.rds")
