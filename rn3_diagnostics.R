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

# ----------------------------------------------------
# EXTRACT STAN DATA FOR DOWNSTREAM PROCESSING
# ----------------------------------------------------

# Extract the “opt” data from stan_data:
y_opt   <- stan_data$y_opt      # length N_opt
i_opt   <- stan_data$i_opt      # respondent-indices (1…N) for each of the N_opt rows
j_opt   <- stan_data$j_opt      # item-indices      (1…J_opt) for each of the N_opt rows
N_opt   <- stan_data$N_opt      # total number of opt-responses
J_opt   <- stan_data$J_opt      # number of “opt” items
lower_opt <- stan_data$lower_opt  # vector of length J_opt
upper_opt <- stan_data$upper_opt

# Extract the “env” data:
y_env   <- stan_data$y_env      # length N_env
i_env   <- stan_data$i_env      # respondent-indices (1…N) for each of the N_env rows
j_env   <- stan_data$j_env      # item-indices      (1…J_env) for each of the N_env rows
N_env   <- stan_data$N_env      # total number of env-responses
J_env   <- stan_data$J_env      # number of “env” items
lower_env <- stan_data$lower_env
upper_env <- stan_data$upper_env

# Extract the “rad” (Radical-Reform) data:
y_rad   <- stan_data$y_rad      # length N_rad
i_rad   <- stan_data$i_rad      # respondent-indices (1…N) for each of the N_rad rows
j_rad   <- stan_data$j_rad      # item-indices      (1…J_rad) for each of the N_rad rows
N_rad   <- stan_data$N_rad      # total number of rad-responses
J_rad   <- stan_data$J_rad      # number of “radical-reform” items
lower_rad <- stan_data$lower_rad
upper_rad <- stan_data$upper_rad

# Covariate & grouping indices:
N       <- stan_data$N                # total # of individuals
P       <- stan_data$P                # # of numeric covariates
X       <- stan_data$X                # design matrix (N × P)

R       <- stan_data$R                # # of regions
region_id <- stan_data$region_id      # region index per individual
Q       <- stan_data$Q                # # of parties
party_id  <- stan_data$party_id       # party index per individual

# Build long-format data frames with consistent column names:
opt_long_df <- tibble(
  respondent = i_opt,
  item       = j_opt,
  response   = y_opt
  #lower      = lower_opt[j_opt],
  #upper      = upper_opt[j_opt]
)
env_long_df <- tibble(
  respondent = i_env,
  item       = j_env,
  response   = y_env
  #lower      = lower_env[j_env],
  #upper      = upper_env[j_env]
)
rad_long_df <- tibble(
  respondent = i_rad,
  item       = j_rad,
  response   = y_rad
  #lower      = lower_rad[j_rad],
  #upper      = upper_rad[j_rad]
)

# ----------------------------------------
# SAMPLING
# ----------------------------------------

fit_prior <- readRDS("model_2/fit_prior_predictive.rds")

fit_full <- readRDS("model_2/fit_full.rds")

# ----------------------------------------
# PRIOR-PREDICTIVE CHECKS
# ----------------------------------------

# Create plots directory if it doesn't exist
if (!dir.exists("diagnostic_plots")) {
  dir.create("diagnostic_plots", recursive = TRUE)
}

# Extract prior-predictive simulated values:
y_opt_sim_pr     <- fit_prior$draws(variables = "y_opt_sim")
y_opt_sim_mat_pr <- as_draws_matrix(y_opt_sim_pr)

y_env_sim_pr     <- fit_prior$draws(variables = "y_env_sim")
y_env_sim_mat_pr <- as_draws_matrix(y_env_sim_pr)

y_rad_sim_pr     <- fit_prior$draws(variables = "y_rad_sim")
y_rad_sim_mat_pr <- as_draws_matrix(y_rad_sim_pr)

# Pull observed responses from long data frames:
obs_opt <- opt_long_df$response
obs_env <- env_long_df$response
obs_rad <- rad_long_df$response

# (1) Density overlay for optimism items
set.seed(123)
n_plot_draws_pr     <- 200
draws_to_plot_opt_pr <- sample(nrow(y_opt_sim_mat_pr), size = n_plot_draws_pr)

sim_opt_pr_df <- tibble(
  draw      = rep(1:n_plot_draws_pr, each = ncol(y_opt_sim_mat_pr)),
  obs_index = rep(1:ncol(y_opt_sim_mat_pr), times = n_plot_draws_pr),
  y_opt_sim = as.vector(y_opt_sim_mat_pr[draws_to_plot_opt_pr, ])
)

ppc_prior_opt <- ggplot() +
  geom_density(
    data  = sim_opt_pr_df,
    aes(x = y_opt_sim, group = draw),
    color = "skyblue",
    alpha = 0.1,
    size  = 0.3
  ) +
  geom_density(
    data = tibble(response = obs_opt),
    aes(x = response),
    color = "black",
    size  = 1
  ) +
  labs(
    title    = "Prior Predictive Check: Optimism (Item Densities)",
    subtitle = "Blue = prior draws; black = observed",
    x        = "Standardised response",
    y        = "Density"
  ) +
  theme_minimal()

print(ppc_prior_opt)
ggsave("diagnostic_plots/prior_predictive_optimism_density.png", ppc_prior_opt,
       width = 6, height = 4, dpi = 300)

# (2) Density overlay for environment items
set.seed(456)
draws_to_plot_env_pr <- sample(nrow(y_env_sim_mat_pr), size = n_plot_draws_pr)

sim_env_pr_df <- tibble(
  draw      = rep(1:n_plot_draws_pr, each = ncol(y_env_sim_mat_pr)),
  obs_index = rep(1:ncol(y_env_sim_mat_pr), times = n_plot_draws_pr),
  y_env_sim = as.vector(y_env_sim_mat_pr[draws_to_plot_env_pr, ])
)

ppc_prior_env <- ggplot() +
  geom_density(
    data  = sim_env_pr_df,
    aes(x = y_env_sim, group = draw),
    color = "lightgreen",
    alpha = 0.1,
    size  = 0.3
  ) +
  geom_density(
    data = tibble(response = obs_env),
    aes(x = response),
    color = "black",
    size  = 1
  ) +
  labs(
    title    = "Prior Predictive Check: Environment (Item Densities)",
    subtitle = "Green = prior draws; black = observed",
    x        = "Standardised response",
    y        = "Density"
  ) +
  theme_minimal()

print(ppc_prior_env)
ggsave("diagnostic_plots/prior_predictive_radical_density.png", ppc_prior_env,
       width = 6, height = 4, dpi = 300)

# (3) Density overlay for radical-reform items
set.seed(789)
draws_to_plot_rad_pr <- sample(nrow(y_rad_sim_mat_pr), size = n_plot_draws_pr)

sim_rad_pr_df <- tibble(
  draw      = rep(1:n_plot_draws_pr, each = ncol(y_rad_sim_mat_pr)),
  obs_index = rep(1:ncol(y_rad_sim_mat_pr), times = n_plot_draws_pr),
  y_rad_sim = as.vector(y_rad_sim_mat_pr[draws_to_plot_rad_pr, ])
)

ppc_prior_rad <- ggplot() +
  geom_density(
    data  = sim_rad_pr_df,
    aes(x = y_rad_sim, group = draw),
    color = "salmon",
    alpha = 0.1,
    size  = 0.3
  ) +
  geom_density(
    data = tibble(response = obs_rad),
    aes(x = response),
    color = "black",
    size  = 1
  ) +
  labs(
    title    = "Prior Predictive Check: Radical-Reform (Item Densities)",
    subtitle = "Red = prior draws; black = observed",
    x        = "Standardised response",
    y        = "Density"
  ) +
  theme_minimal()

print(ppc_prior_rad)
ggsave("diagnostic_plots/prior_predictive_radical_density.png", ppc_prior_rad,
       width = 6, height = 4, dpi = 300)

# (4) Person-level prior-predictive: optimism (φ)
person_index_opt <- opt_long_df$respondent

sim_person_opt_list <- lapply(draws_to_plot_opt_pr, function(d_idx) {
  one_draw <- y_opt_sim_mat_pr[d_idx, ]
  tibble(
    respondent = person_index_opt,
    y_pr       = one_draw
  ) %>%
    group_by(respondent) %>%
    summarize(sim_phi_bar = mean(y_pr), .groups = "drop") %>%
    mutate(draw = d_idx)
})

sim_person_opt_df <- bind_rows(sim_person_opt_list)

person_obs_opt <- opt_long_df %>%
  group_by(respondent) %>%
  summarize(obs_phi_bar = mean(response), .groups = "drop")

person_opt_plot <- ggplot() +
  geom_point(
    data = sim_person_opt_df,
    aes(x = respondent, y = sim_phi_bar, group = draw),
    color = "steelblue",
    alpha = 0.05,
    size  = 0.5
  ) +
  geom_point(
    data = person_obs_opt,
    aes(x = respondent, y = obs_phi_bar),
    color = "black",
    size  = 1
  ) +
  labs(
    title    = "Prior Predictive Check: Optimism (Person Means)",
    subtitle = "Black = observed φ̄ᵢ; Blue = prior-predicted draws",
    x        = "Person index",
    y        = "Mean response (φ̄ᵢ)"
  ) +
  theme_minimal()

print(person_opt_plot)
ggsave("diagnostic_plots/prior_predictive_optimism_person_means.png", person_opt_plot,
       width = 6, height = 4, dpi = 300)

# (5) Person-level prior-predictive: environment (θ)
person_index_env <- env_long_df$respondent

sim_person_env_list <- lapply(draws_to_plot_env_pr, function(d_idx) {
  one_draw <- y_env_sim_mat_pr[d_idx, ]
  tibble(
    respondent = person_index_env,
    y_pr       = one_draw
  ) %>%
    group_by(respondent) %>%
    summarize(sim_theta_bar = mean(y_pr), .groups = "drop") %>%
    mutate(draw = d_idx)
})

sim_person_env_df <- bind_rows(sim_person_env_list)

person_obs_env <- env_long_df %>%
  group_by(respondent) %>%
  summarize(obs_theta_bar = mean(response), .groups = "drop")

person_env_plot <- ggplot() +
  geom_point(
    data  = sim_person_env_df,
    aes(x = respondent, y = sim_theta_bar, group = draw),
    color = "lightgreen",
    alpha = 0.05,
    size  = 0.5
  ) +
  geom_point(
    data = person_obs_env,
    aes(x = respondent, y = obs_theta_bar),
    color = "black",
    size  = 1
  ) +
  labs(
    title    = "Prior Predictive Check: Environmentalism (Person Means)",
    subtitle = "Black = observed θ̄ᵢ; Green = prior-predicted draws",
    x        = "Person index",
    y        = "Mean response (θ̄ᵢ)"
  ) +
  theme_minimal()

print(person_env_plot)
ggsave("diagnostic_plots/prior_predictive_environment_person_means.png", person_env_plot,
       width = 6, height = 4, dpi = 300)

# (6) Person-level prior-predictive: radical (ψ)
person_index_rad <- rad_long_df$respondent

sim_person_rad_list <- lapply(draws_to_plot_rad_pr, function(d_idx) {
  one_draw <- y_rad_sim_mat_pr[d_idx, ]
  tibble(
    respondent = person_index_rad,
    y_pr       = one_draw
  ) %>%
    group_by(respondent) %>%
    summarize(sim_psi_bar = mean(y_pr), .groups = "drop") %>%
    mutate(draw = d_idx)
})

sim_person_rad_df <- bind_rows(sim_person_rad_list)

person_obs_rad <- rad_long_df %>%
  group_by(respondent) %>%
  summarize(obs_psi_bar = mean(response), .groups = "drop")

person_rad_plot <- ggplot() +
  geom_point(
    data  = sim_person_rad_df,
    aes(x = respondent, y = sim_psi_bar, group = draw),
    color = "salmon",
    alpha = 0.05,
    size  = 0.5
  ) +
  geom_point(
    data = person_obs_rad,
    aes(x = respondent, y = obs_psi_bar),
    color = "black",
    size  = 1
  ) +
  labs(
    title    = "Prior Predictive Check: Radical-Reform Person Means",
    subtitle = "Black = observed ψ̄ᵢ; red = prior-predicted draws",
    x        = "Person index",
    y        = "Mean response (ψ̄ᵢ)"
  ) +
  theme_minimal()

print(person_rad_plot)
ggsave("diagnostic_plots/prior_predictive_radical_person_means.png", person_rad_plot,
       width = 6, height = 4, dpi = 300)

# (7) Prior vs. Posterior SD ratio for all parameters
prior_draws  <- as_draws_df(fit_prior$draws())
prior_summary <- prior_draws %>%
  summarise_draws(sd) %>%
  select(variable, prior_sd = sd)

post_draws   <- as_draws_df(fit_full$draws())
post_summary <- post_draws %>%
  summarise_draws(sd) %>%
  select(variable, post_sd = sd)

joined <- prior_summary %>%
  inner_join(post_summary, by = "variable") %>%
  mutate(ratio = post_sd / prior_sd) %>%
  arrange(ratio)

cat("\nTop 20 parameters with smallest posterior/prior SD ratio:\n")
joined %>%
  head(20) %>%
  print(n = Inf)

# ----------------------------------------
# CONVERGENCE DIAGNOSTICS
# ----------------------------------------

# 1) Rhat & ESS
summ_full <- fit_full$summary()

bad_rhat     <- summ_full %>% filter(rhat > 1.01)   %>% arrange(desc(rhat))
bad_ess_bulk <- summ_full %>% filter(ess_bulk < 200) %>% arrange(ess_bulk)
bad_ess_tail <- summ_full %>% filter(ess_tail < 200) %>% arrange(ess_tail)

cat("\nParameters with Rhat > 1.01:\n")
print(bad_rhat)

cat("\nParameters with ESS_bulk < 200:\n")
print(bad_ess_bulk)

cat("\nParameters with ESS_tail < 200:\n")
print(bad_ess_tail)

# 2) Sampler diagnostics: divergences, max-treedepth, BFMI
sampler_diag_df <- as_tibble(fit_full$diagnostic_summary())

bad_diverg   <- sampler_diag_df %>% filter(num_divergent > 0)
bad_treedep  <- sampler_diag_df %>% filter(num_max_treedepth > 0)
bad_bfmi     <- sampler_diag_df %>% filter(ebfmi < 0.2)  # low BFMI threshold

cat("\nChains with divergent transitions:\n")
print(bad_diverg)

cat("\nChains that hit max treedepth:\n")
print(bad_treedep)

cat("\nChains with low BFMI:\n")
print(bad_bfmi)

# 3) Rhat histogram
rhat_df <- tibble(rhat = summ_full$rhat)

rhat_hist <- ggplot(rhat_df, aes(x = rhat)) +
  geom_histogram(
    breaks = seq(0.99, 1.05, by = 0.001),
    fill   = "steelblue",
    color  = "white"
  ) +
  geom_vline(
    xintercept = 1.01,
    linetype   = "dotted",
    color      = "red",
    size       = 1
  ) +
  labs(
    title = "Histogram of Rhat for All Parameters",
    x     = "Rhat",
    y     = "Count"
  ) +
  theme_minimal()

print(rhat_hist)
ggsave("diagnostic_plots/rhat_histogram.png", rhat_hist,
       width = 8, height = 6, dpi = 150)

# 4) Traceplots and density plots for selected parameters (including radical-reform)
draws_arr <- as_draws_array(fit_full$draws())

pars_to_plot <- c(
  "tau_eta[1]", "tau_eta[2]", "tau_eta[3]",
  "sigma_alpha[1]", "sigma_alpha[2]", "sigma_alpha[3]",
  "lambda_opt[1]",
  "lambda_env[1]",
  "lambda_rad[1]"
)

traceplot_sel <- mcmc_trace(
  draws_arr,
  pars       = pars_to_plot,
  facet_args = list(ncol = 1, strip.position = "right")
) +
  ggtitle("Traceplots for Selected Parameters (incl. Radical)") +
  theme_minimal()

print(traceplot_sel)
ggsave("diagnostic_plots/traceplots_selected_params.png", traceplot_sel,
       width = 6, height = 10, dpi = 300)

mcmc_dens <- mcmc_dens(
  draws_arr,
  pars       = pars_to_plot,
  facet_args = list(ncol = 1, strip.position = "right")
) +
  ggtitle("Posterior Densities for Selected Parameters") +
  theme_minimal()

print(mcmc_dens)
ggsave("diagnostic_plots/densities_selected_params.png", mcmc_dens,
       width = 6, height = 10, dpi = 300)

# ----------------------------------------
# POSTERIOR-PREDICTIVE CHECKS
# ----------------------------------------

draws_full_df    <- as_draws_df(fit_full$draws())
y_opt_sim_draws  <- draws_full_df %>% select(starts_with("y_opt_sim["))
y_env_sim_draws  <- draws_full_df %>% select(starts_with("y_env_sim["))
y_rad_sim_draws  <- draws_full_df %>% select(starts_with("y_rad_sim["))

# (1) Density overlay: optimism
set.seed(123)
n_plot_draws_pp     <- 200
draws_to_plot_opt_pp <- sample(nrow(y_opt_sim_draws), size = n_plot_draws_pp)

sim_opt_df <- y_opt_sim_draws %>%
  slice(draws_to_plot_opt_pp) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "obs_index",
    values_to = "y_opt_sim"
  ) %>%
  mutate(draw = rep(1:length(draws_to_plot_opt_pp),
                    each = ncol(y_opt_sim_draws)))

ppc_plot_opt <- ggplot() +
  geom_density(
    data  = sim_opt_df,
    aes(x = y_opt_sim, group = draw),
    color = "steelblue",
    alpha = 0.1,
    size  = 0.3
  ) +
  geom_density(
    data = tibble(response = obs_opt),
    aes(x = response),
    color = "black",
    size  = 1
  ) +
  labs(
    title    = "Posterior Predictive Check: Optimism (Item Densities)",
    subtitle = "Blue = simulated; black = observed",
    x        = "Standardised response",
    y        = "Density"
  ) +
  theme_minimal()

print(ppc_plot_opt)
ggsave("diagnostic_plots/posterior_predictive_optimism_density.png", ppc_plot_opt,
       width = 6, height = 4, dpi = 300)

# (2) Density overlay: environment
set.seed(456)
draws_to_plot_env_pp <- sample(nrow(y_env_sim_draws), size = n_plot_draws_pp)

sim_env_df <- y_env_sim_draws %>%
  slice(draws_to_plot_env_pp) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "obs_index",
    values_to = "y_env_sim"
  ) %>%
  mutate(draw = rep(1:length(draws_to_plot_env_pp),
                    each = ncol(y_env_sim_draws)))

ppc_plot_env <- ggplot() +
  geom_density(
    data  = sim_env_df,
    aes(x = y_env_sim, group = draw),
    color = "lightgreen",
    alpha = 0.1,
    size  = 0.3
  ) +
  geom_density(
    data = tibble(response = obs_env),
    aes(x = response),
    color = "black",
    size  = 1
  ) +
  labs(
    title    = "Posterior Predictive Check: Environmentalism (Item Densities)",
    subtitle = "Green = simulated; black = observed",
    x        = "Standardised response",
    y        = "Density"
  ) +
  theme_minimal()

print(ppc_plot_env)
ggsave("diagnostic_plots/posterior_predictive_environment_density.png", ppc_plot_env,
       width = 6, height = 4, dpi = 300)

# (3) Density overlay: radical-reform
set.seed(789)
draws_to_plot_rad_pp <- sample(nrow(y_rad_sim_draws), size = n_plot_draws_pp)

sim_rad_df <- y_rad_sim_draws %>%
  slice(draws_to_plot_rad_pp) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "obs_index",
    values_to = "y_rad_sim"
  ) %>%
  mutate(draw = rep(1:length(draws_to_plot_rad_pp),
                    each = ncol(y_rad_sim_draws)))

ppc_plot_rad <- ggplot() +
  geom_density(
    data  = sim_rad_df,
    aes(x = y_rad_sim, group = draw),
    color = "salmon",
    alpha = 0.1,
    size  = 0.3
  ) +
  geom_density(
    data = tibble(response = obs_rad),
    aes(x = response),
    color = "black",
    size  = 1
  ) +
  labs(
    title    = "Posterior Predictive Check: Radical-Reform (Item Densities)",
    subtitle = "red = simulated; black = observed",
    x        = "Standardised response",
    y        = "Density"
  ) +
  theme_minimal()

print(ppc_plot_rad)
ggsave("diagnostic_plots/posterior_predictive_radical_density.png", ppc_plot_rad,
       width = 6, height = 4, dpi = 300)

# (4) Person-level posterior-predictive: optimism (φ)
person_index_opt <- opt_long_df$respondent
mat_opt_pp       <- as.matrix(y_opt_sim_draws)
stopifnot(ncol(mat_opt_pp) == nrow(opt_long_df))

sim_person_opt_list <- lapply(draws_to_plot_opt_pp, function(d_idx) {
  one_draw <- mat_opt_pp[d_idx, ]
  tibble(
    respondent = person_index_opt,
    y_pp       = one_draw
  ) %>%
    group_by(respondent) %>%
    summarize(sim_phi_bar = mean(y_pp), .groups = "drop") %>%
    mutate(draw = d_idx)
})

sim_person_opt_df <- bind_rows(sim_person_opt_list)

person_obs_opt <- opt_long_df %>%
  group_by(respondent) %>%
  summarize(obs_phi_bar = mean(response), .groups = "drop")

person_pp_opt_plot <- ggplot() +
  geom_point(
    data  = sim_person_opt_df,
    aes(x = respondent, y = sim_phi_bar, group = draw),
    color = "steelblue",
    alpha = 0.05,
    size  = 0.5
  ) +
  geom_point(
    data = person_obs_opt,
    aes(x = respondent, y = obs_phi_bar),
    color = "black",
    size  = 1
  ) +
  labs(
    title    = "Posterior Predictive Check: Optimism Person Means",
    subtitle = "Black = observed φ̄ᵢ; Blue = posterior-predicted draws",
    x        = "Person index",
    y        = "Mean response (φ̄ᵢ)"
  ) +
  theme_minimal()

print(person_pp_opt_plot)
ggsave("diagnostic_plots/posterior_predictive_optimism_person_means.png", person_pp_opt_plot,
       width = 6, height = 4, dpi = 300)

# (5) Person-level posterior-predictive: environment (θ)
person_index_env <- env_long_df$respondent
mat_env_pp       <- as.matrix(y_env_sim_draws)
stopifnot(ncol(mat_env_pp) == nrow(env_long_df))

sim_person_env_list <- lapply(draws_to_plot_env_pp, function(d_idx) {
  one_draw <- mat_env_pp[d_idx, ]
  tibble(
    respondent = person_index_env,
    y_pp       = one_draw
  ) %>%
    group_by(respondent) %>%
    summarize(sim_theta_bar = mean(y_pp), .groups = "drop") %>%
    mutate(draw = d_idx)
})

sim_person_env_df <- bind_rows(sim_person_env_list)

person_obs_env <- env_long_df %>%
  group_by(respondent) %>%
  summarize(obs_theta_bar = mean(response), .groups = "drop")

person_pp_env_plot <- ggplot() +
  geom_point(
    data  = sim_person_env_df,
    aes(x = respondent, y = sim_theta_bar, group = draw),
    color = "lightgreen",
    alpha = 0.05,
    size  = 0.5
  ) +
  geom_point(
    data = person_obs_env,
    aes(x = respondent, y = obs_theta_bar),
    color = "black",
    size  = 1
  ) +
  labs(
    title    = "Posterior Predictive Check: Environment (Person Means)",
    subtitle = "Black = observed θ̄ᵢ; Green = posterior-predicted draws",
    x        = "Person index",
    y        = "Mean response (θ̄ᵢ)"
  ) +
  theme_minimal()

print(person_pp_env_plot)
ggsave("diagnostic_plots/posterior_predictive_environment_person_means.png", person_pp_env_plot,
       width = 6, height = 4, dpi = 300)

# (6) Person-level posterior-predictive: radical (ψ)
person_index_rad <- rad_long_df$respondent
mat_rad_pp       <- as.matrix(y_rad_sim_draws)
stopifnot(ncol(mat_rad_pp) == nrow(rad_long_df))

sim_person_rad_list <- lapply(draws_to_plot_rad_pp, function(d_idx) {
  one_draw <- mat_rad_pp[d_idx, ]
  tibble(
    respondent = person_index_rad,
    y_pp       = one_draw
  ) %>%
    group_by(respondent) %>%
    summarize(sim_psi_bar = mean(y_pp), .groups = "drop") %>%
    mutate(draw = d_idx)
})

sim_person_rad_df <- bind_rows(sim_person_rad_list)

person_obs_rad <- rad_long_df %>%
  group_by(respondent) %>%
  summarize(obs_psi_bar = mean(response), .groups = "drop")

person_pp_rad_plot <- ggplot() +
  geom_point(
    data  = sim_person_rad_df,
    aes(x = respondent, y = sim_psi_bar, group = draw),
    color = "salmon",
    alpha = 0.05,
    size  = 0.5
  ) +
  geom_point(
    data = person_obs_rad,
    aes(x = respondent, y = obs_psi_bar),
    color = "black",
    size  = 1
  ) +
  labs(
    title    = "Posterior Predictive Check: Radical-Reform (Person Means)",
    subtitle = "Black = observed ψ̄ᵢ; red = posterior-predicted draws",
    x        = "Person index",
    y        = "Mean response (ψ̄ᵢ)"
  ) +
  theme_minimal()

print(person_pp_rad_plot)
ggsave("diagnostic_plots/posterior_predictive_radical_person_means.png", person_pp_rad_plot,
       width = 6, height = 4, dpi = 300)

# (7) Item-level posterior-predictive check: environment means
env_obs_df <- env_long_df %>%
  mutate(row_index = row_number())

mat_env_sim <- as.matrix(y_env_sim_draws)
ppc_mean_env_per_index <- colMeans(mat_env_sim)

item_summary_env <- env_obs_df %>%
  mutate(ppc_mean = ppc_mean_env_per_index[row_index]) %>%
  group_by(item) %>%
  summarize(
    obs_mean   = mean(response),
    ppc_mean   = mean(ppc_mean),
    obs_sd     = sd(response),
    ppc_sd_est = sd(ppc_mean),
    .groups    = "drop"
  )

print(item_summary_env)

item_env_plot <- ggplot(item_summary_env, aes(x = obs_mean, y = ppc_mean)) +
  geom_point(size = 3, color = "lightgreen") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Item-Level Check: Environment Means",
    x     = "Observed mean (response)",
    y     = "Posterior-predictive mean"
  ) +
  theme_minimal()

print(item_env_plot)
ggsave("diagnostic_plots/item_level_means_environment.png", item_env_plot,
       width = 6, height = 4, dpi = 300)

# (NEW) Item-level posterior-predictive check: optimism means
opt_obs_df <- opt_long_df %>%
  mutate(row_index = row_number())

mat_opt_sim <- as.matrix(y_opt_sim_draws)
ppc_mean_opt_per_index <- colMeans(mat_opt_sim)

item_summary_opt <- opt_obs_df %>%
  mutate(ppc_mean = ppc_mean_opt_per_index[row_index]) %>%
  group_by(item) %>%
  summarize(
    obs_mean   = mean(response),
    ppc_mean   = mean(ppc_mean),
    obs_sd     = sd(response),
    ppc_sd_est = sd(ppc_mean),
    .groups    = "drop"
  )

print(item_summary_opt)

item_opt_plot <- ggplot(item_summary_opt, aes(x = obs_mean, y = ppc_mean)) +
  geom_point(size = 3, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Item-Level Check: Optimism Means",
    x     = "Observed mean (response)",
    y     = "Posterior-predictive mean"
  ) +
  theme_minimal()

print(item_opt_plot)
ggsave("diagnostic_plots/item_level_means_optimism.png", item_opt_plot,
       width = 6, height = 4, dpi = 300)


# (8) Item-level posterior-predictive check: radical means
rad_obs_df <- rad_long_df %>%
  mutate(row_index = row_number())

mat_rad_sim <- as.matrix(y_rad_sim_draws)
ppc_mean_rad_per_index <- colMeans(mat_rad_sim)

item_summary_rad <- rad_obs_df %>%
  mutate(ppc_mean = ppc_mean_rad_per_index[row_index]) %>%
  group_by(item) %>%
  summarize(
    obs_mean   = mean(response),
    ppc_mean   = mean(ppc_mean),
    obs_sd     = sd(response),
    ppc_sd_est = sd(ppc_mean),
    .groups    = "drop"
  )

print(item_summary_rad)

item_rad_plot <- ggplot(item_summary_rad, aes(x = obs_mean, y = ppc_mean)) +
  geom_point(size = 3, color = "salmon") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  labs(
    title = "Item-Level Check: Radical-Reform Means",
    x     = "Observed mean (response)",
    y     = "Posterior-predictive mean"
  ) +
  theme_minimal()

print(item_rad_plot)
ggsave("diagnostic_plots/item_level_radicalreform.png", item_rad_plot,
       width = 6, height = 4, dpi = 300)

# (9) Category-frequency check for a single optimism item (item = 1)
opt_item1_obs <- opt_long_df %>%
  filter(item == 1) %>%
  pull(response)

item1_row_idx <- which(opt_long_df$item == 1)
sim_item1_matrix <- mat_opt_pp[, item1_row_idx]

round_to_likert <- function(x) {
  x_rounded <- round(x)
  pmin(pmax(x_rounded, -2), 2)
}

obs_cats <- round_to_likert(opt_item1_obs)
obs_freq <- table(factor(obs_cats, levels = c(-2, -1, 0, 1, 2)))
obs_cat_df <- tibble(
  cat        = factor(c(-2, -1, 0, 1, 2), levels = c(-2, -1, 0, 1, 2)),
  proportion = as.numeric(obs_freq) / sum(obs_freq)
)

sim_cat_props <- apply(sim_item1_matrix, 1, function(vec_sim) {
  sim_cats <- round_to_likert(vec_sim)
  tab <- table(factor(sim_cats, levels = c(-2, -1, 0, 1, 2)))
  as.numeric(tab / sum(tab))
})

sim_cat_df <- as_tibble(t(sim_cat_props)) %>%
  set_names(c("-2", "-1", "0", "1", "2")) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(
    cols       = c("-2", "-1", "0", "1", "2"),
    names_to   = "cat",
    values_to  = "proportion"
  )

sim_cat_summary <- sim_cat_df %>%
  group_by(cat) %>%
  summarize(
    p50 = median(proportion),
    p10 = quantile(proportion, 0.10),
    p90 = quantile(proportion, 0.90),
    .groups = "drop"
  )

cat_freq_plot <- ggplot() +
  geom_col(
    data  = obs_cat_df,
    aes(x = cat, y = proportion),
    fill  = "gray80"
  ) +
  geom_errorbar(
    data  = sim_cat_summary,
    aes(x = cat, ymin = p10, ymax = p90),
    width = 0.2,
    color = "steelblue",
    size  = 1
  ) +
  geom_point(
    data  = sim_cat_summary,
    aes(x = cat, y = p50),
    color = "steelblue",
    size  = 2
  ) +
  labs(
    title    = "Category Frequencies for Optimism Item 1",
    subtitle = "Gray bars = observed; Blue = posterior predictive",
    x        = "Rounded Likert Category",
    y        = "Proportion"
  ) +
  theme_minimal()

print(cat_freq_plot)
ggsave("diagnostic_plots/category_frequencies_optimism_item1.png", cat_freq_plot,
       width = 6, height = 4, dpi = 300)

# (10) Category-frequency check for a single radical item (item = 1)
rad_item1_obs <- rad_long_df %>%
  filter(item == 1) %>%
  pull(response)

item_rad1_row_idx <- which(rad_long_df$item == 1)
sim_rad1_matrix <- mat_rad_pp[, item_rad1_row_idx]

obs_cats_rad <- round_to_likert(rad_item1_obs)
obs_freq_rad <- table(factor(obs_cats_rad, levels = c(-2, -1, 0, 1, 2)))
obs_cat_rad_df <- tibble(
  cat        = factor(c(-2, -1, 0, 1, 2), levels = c(-2, -1, 0, 1, 2)),
  proportion = as.numeric(obs_freq_rad) / sum(obs_freq_rad)
)

sim_cat_props_rad <- apply(sim_rad1_matrix, 1, function(vec_sim) {
  sim_cats <- round_to_likert(vec_sim)
  tab <- table(factor(sim_cats, levels = c(-2, -1, 0, 1, 2)))
  as.numeric(tab / sum(tab))
})

sim_cat_df_rad <- as_tibble(t(sim_cat_props_rad)) %>%
  set_names(c("-2", "-1", "0", "1", "2")) %>%
  mutate(draw = row_number()) %>%
  pivot_longer(
    cols       = c("-2", "-1", "0", "1", "2"),
    names_to   = "cat",
    values_to  = "proportion"
  )

sim_cat_summary_rad <- sim_cat_df_rad %>%
  group_by(cat) %>%
  summarize(
    p50 = median(proportion),
    p10 = quantile(proportion, 0.10),
    p90 = quantile(proportion, 0.90),
    .groups = "drop"
  )

cat_freq_rad_plot <- ggplot() +
  geom_col(
    data  = obs_cat_rad_df,
    aes(x = cat, y = proportion),
    fill  = "gray80"
  ) +
  geom_errorbar(
    data  = sim_cat_summary_rad,
    aes(x = cat, ymin = p10, ymax = p90),
    width = 0.2,
    color = "salmon",
    size  = 1
  ) +
  geom_point(
    data  = sim_cat_summary_rad,
    aes(x = cat, y = p50),
    color = "salmon",
    size  = 2
  ) +
  labs(
    title    = "Category Frequencies for Radical Item 1",
    subtitle = "Gray bars = observed; Red = posterior predictive",
    x        = "Rounded Likert Category",
    y        = "Proportion"
  ) +
  theme_minimal()

print(cat_freq_rad_plot)
ggsave("diagnostic_plots/category_frequencies_radical_item1.png", cat_freq_rad_plot,
       width = 6, height = 4, dpi = 300)

# ---------------------------------------------------------------------------------------
# 7. POSTERIOR R² ESTIMATES (including Radical-Reform)
# ---------------------------------------------------------------------------------------
# Extract posterior draws of R² from the full model:
r2_draws_df <- fit_full$draws(
  variables = c("R2_opt", "R2_env", "R2_rad"),
  format    = "draws_df"
)

# Quick numerical summary:
r2_summary <- r2_draws_df %>%
  pivot_longer(cols = everything(), names_to = "block", values_to = "R2") %>%
  group_by(block) %>%
  summarize(
    mean_R2   = mean(R2),
    median_R2 = median(R2),
    lower95   = quantile(R2, 0.025),
    upper95   = quantile(R2, 0.975),
    .groups   = "drop"
  )

print(r2_summary)



r2_long <- r2_draws_df %>%
  pivot_longer(
    cols      = starts_with("R2_"),
    names_to  = "block",
    values_to = "R2"
  )

# Plot posterior density of each R²
r2_density_plot <- ggplot(r2_long, aes(x = R2, fill = block)) +
                      geom_density(alpha = 0.4) +
                      labs(
                        title = "Posterior Densities of R²\n(Optimism vs. Environment vs. Radical)",
                        x     = expression(R^2),
                        y     = "Density"
                      ) +
                      theme_minimal()

print(r2_density_plot)
ggsave("diagnostic_plots/r2_posterior_density.png", r2_density_plot,
       width = 6, height = 4, dpi = 300)

cat("\n=========================================\n")
cat("Modeling & Diagnostics Completed Successfully\n")
cat("All diagnostic_plots saved under 'diagnostic_plots/' directory.\n")
cat("=========================================\n\n")
