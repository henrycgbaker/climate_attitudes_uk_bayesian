# ---------------------------------------------------------------------------------------
# DOWNSTREAM ANALYSIS: QUANTITIES OF INTEREST & STRATEGIC VISUALIZATIONS
# ---------------------------------------------------------------------------------------

# 0. Setup -------------------------------------------------------------------------------
library(tidyverse)
library(cmdstanr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(ggdist)      
library(forcats)
library(cowplot)     
library(corrr)        
library(ggraph)       
library(igraph)      
library(scales)       
library(plotly)
library(ggrepel)
library(ggridges)      

setwd("/Users/henrybaker/repositories/bayesian_modeling/research_note")

if (!dir.exists("qoi_outputs")) dir.create("qoi_outputs")
if (!dir.exists("qoi_outputs/plots")) dir.create("qoi_outputs/plots")

# ---------------------------------------------------------------------------------------
# 1. Load Data & Model Fit --------------------------------------------------------------
# ---------------------------------------------------------------------------------------

stan_data <- readRDS("data/stan_data_fit_fast.rds")
fit_full   <- readRDS("model_2_fit_fast/fit_full.rds")

# 1c. Define region and party labels (must match the wrangling script)
region_levels <- c(
  "Yorkshire & the Humber", "West Midlands", "Scotland", "Wales",
  "North West", "Eastern", "South West", "East Midlands", "London", "South East"
)

party_levels <- c(
  "Green", "Labour", "Plaid Cymru", "Scottish National Party (SNP)",
  "Liberal Democrat", "Conservative", "Reform UK",
  "Another party", "Don’t know", "Won’t vote"
)

# ---------------------------------------------------------------------------------------
# 2. Extract Posterior Draws ------------------------------------------------------------
# ---------------------------------------------------------------------------------------

draws_df <- as_draws_df(fit_full$draws())
n_draws  <- nrow(draws_df)

# ---------------------------------------------------------------------------------------
# 3. Posterior Slope Coefficients -------------------------------------------------------
# ---------------------------------------------------------------------------------------

# 3a. Extract all slope coefficients B[1..3,1..P]
P <- stan_data$P
slope_names <- paste0("B[", rep(1:3, each = P), ",", rep(1:P, times = 3), "]")

# 3b. Summarize posterior means (skip CI columns)
slopes_df <- draws_df %>%
  select(all_of(slope_names)) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "param",
    values_to = "value"
  ) %>%
  group_by(param) %>%
  summarize(
    mean    = mean(value),
    .groups = "drop"
  ) %>%
  arrange(param)

# 3c. Save slopes summary (means only) to CSV and print first rows
write_csv(slopes_df, "qoi_outputs/slopes_summary.csv")
cat("\n=== Posterior Slopes Summary (Means Only) ===\n")
print(n = 50, slopes_df, 50)

# ---------------------------------------------------------------------------------------
# 3d. Covariate Effects: Overlapping Ridgeline Plot (replacing half-eye) ----------------
# ---------------------------------------------------------------------------------------

# 3d.i. Define covariate names in order used in stan_data$X
cov_names <- c(
  "gender",
  "age_65plus", "age_55_64", "age_45_54", "age_35_44", "age_25_34",
  "edu_L4plus", "edu_L3", "edu_L2", "edu_L1", "edu_appr", "edu_other",
  "insecurity"
)

# 3d.ii. Annotate each draw with latent dimension and covariate index (long format)
slopes_draws_long <- draws_df %>%
  select(all_of(slope_names)) %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "param",
    values_to = "value"
  ) %>%
  mutate(
    latent = case_when(
      str_detect(param, "^B\\[1,") ~ "φ (Optimism)",
      str_detect(param, "^B\\[2,") ~ "θ (Environment)",
      str_detect(param, "^B\\[3,") ~ "ψ (Radical-Reform)",
      TRUE                         ~ NA_character_
    ),
    cov_index = as.integer(str_remove_all(param, "B\\[[123],|\\]"))
  ) %>%
  mutate(
    covariate = factor(cov_names[cov_index], levels = cov_names)
  ) %>%
  select(latent, covariate, cov_index, value)

# 3d.iii. Compute a small vertical “nudge” so that φ, θ, ψ densities at the same covariate
#         do not sit exactly on top of one another
slopes_draws_long <- slopes_draws_long %>%
  mutate(
    y_position = case_when(
      latent == "φ (Optimism)"      ~ cov_index - 0.15,
      latent == "θ (Environment)"   ~ cov_index,
      latent == "ψ (Radical-Reform)"~ cov_index + 0.15,
      TRUE                           ~ cov_index
    )
  )

# 3d.iv. Rebuild the ridgeline with a larger vertical scale—use default ggplot discrete colors
covariate_ridge_overlap <- ggplot(
  slopes_draws_long,
  aes(
    x     = value,
    y     = y_position,
    fill  = latent,
    color = latent
  )
) +
  geom_density_ridges(
    aes(
      height = ..density..,
      group  = interaction(latent, covariate)
    ),
    stat            = "density",
    scale           = 3,       # ↑ make ridges taller
    rel_min_height  = 0.01,
    alpha           = 0.6,
    trim            = FALSE
  ) +
  # dashed vertical line at x = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray40") +
  # label each covariate at integer positions
  scale_y_continuous(
    breaks = seq_len(length(levels(slopes_draws_long$covariate))),
    labels = levels(slopes_draws_long$covariate)
  ) +
  # Use default discrete scales so that we match Script 1’s original coloring
  scale_fill_discrete(name = "Latent Trait") +
  scale_color_discrete(name = "Latent Trait") +
  coord_cartesian(expand = FALSE) +
  theme_minimal(base_size = 12) +
  labs(
    title    = "Individual-level Covariate Effects on Latent Traits",
    x        = "Posterior Slope Value",
    y        = "Covariate"
  ) +
  theme(
    axis.title.y       = element_text(size = 11),
    axis.text.y        = element_text(size = 9),
    legend.position    = "bottom",
    legend.key.size    = unit(0.4, "cm"),
    panel.grid.major.y = element_blank()
  )

# 3d.v. Save the new ridgeline plot
ggsave(
  filename = "qoi_outputs/plots/covariate_effects_overlapping_ridgelines.png",
  plot     = covariate_ridge_overlap,
  width    = 9,
  height   = 15,  
  dpi      = 300
)

print(covariate_ridge_overlap)

# ---------------------------------------------------------------------------------------
# 4. Posterior Latent Scores (φᵢ, θᵢ, ψᵢ) ---------------------------------------------------

# 4a. Extract phi[i], theta[i], psi[i] from draws
N <- stan_data$N
phi_names   <- paste0("phi[", 1:N, "]")
theta_names <- paste0("theta[", 1:N, "]")
psi_names   <- paste0("psi[", 1:N, "]")

phi_draws   <- draws_df %>% select(all_of(phi_names))
theta_draws <- draws_df %>% select(all_of(theta_names))
psi_draws   <- draws_df %>% select(all_of(psi_names))

# 4b. Compute posterior means for each respondent
phi_means   <- colMeans(phi_draws)
theta_means <- colMeans(theta_draws)
psi_means   <- colMeans(psi_draws)

latent_scores_df <- tibble(
  respondent = 1:N,
  phi_hat    = phi_means,
  theta_hat  = theta_means,
  psi_hat    = psi_means
)

write_csv(latent_scores_df, "qoi_outputs/posterior_latent_scores.csv")
cat("\n=== Posterior Latent Scores (first 20 respondents) ===\n")
print(n = 20, latent_scores_df)

# ---------------------------------------------------------------------------------------
# 5. Average Marginal Effects (AME) --------------------------------------------------------

# 5a. Since the latent-regression is linear, AME_p^(ℓ) = β_p^(ℓ)
#     We already extracted means in “slopes_df”. Save AME table as is.
ame_df <- slopes_draws_long %>%
  group_by(latent, covariate) %>%
  summarize(mean = mean(value), .groups = "drop") %>%
  arrange(latent, covariate)

write_csv(ame_df, "qoi_outputs/ame_summary.csv")
cat("\n=== Average Marginal Effects Summary (first 10 rows) ===\n")
print(head(ame_df, 10))

# ---------------------------------------------------------------------------------------
# 6. Bayesian R² ----------------------------------------------------------------------------

# 6a. Extract R² draws (R2_opt, R2_env, R2_rad)
r2_draws_df <- draws_df %>%
  select(starts_with("R2_opt"), starts_with("R2_env"), starts_with("R2_rad")) %>%
  rename_with(~ c("R2_opt", "R2_env", "R2_rad"))

# 6b. Summarize posterior mean for each R² (skip CI calc)
r2_summary <- r2_draws_df %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "block",
    values_to = "R2"
  ) %>%
  group_by(block) %>%
  summarize(
    mean_R2   = mean(R2),
    .groups   = "drop"
  )

write_csv(r2_summary, "qoi_outputs/r2_summary.csv")
cat("\n=== Bayesian R² Summary (Means Only) ===\n")
print(r2_summary)

# 6c. Plot posterior densities of R² (φ, θ, ψ) — unchanged
r2_density <- r2_draws_df %>%
  pivot_longer(cols = everything(), names_to = "block", values_to = "R2") %>%
  ggplot(aes(x = R2, fill = block)) +
  geom_density(alpha = 0.4) +
  labs(
    title = "Posterior Densities of Bayesian R²",
    x     = expression(R^2),
    y     = "Density",
    fill  = "Latent Block"
  ) +
  theme_minimal() +
  scale_fill_brewer(palette = "Set2")

ggsave(
  filename = "qoi_outputs/plots/r2_posteriors.png",
  plot     = r2_density,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# ---------------------------------------------------------------------------------------
# 7. Group-Level Intercepts: Regions & Parties ------------------------------------------------

# 7a. Extract draws of z_alpha[1,1], …, z_alpha[3,R], sigma_alpha[1..3], Lcorr_alpha (lower-triangular entries)
z_alpha_names <- expand.grid(l = 1:3, r = 1:stan_data$R) %>%
  mutate(name = paste0("z_alpha[", l, ",", r, "]")) %>%
  pull(name)

sigma_alpha_names <- paste0("sigma_alpha[", 1:3, "]")

Lcorr_alpha_names <- expand.grid(i = 1:3, j = 1:3) %>%
  filter(i >= j) %>%
  mutate(name = paste0("Lcorr_alpha[", i, ",", j, "]")) %>%
  pull(name)

# Convert draws to array form for alpha reconstruction
draws_array <- as_draws_array(fit_full$draws())

# Extract z_alpha and sigma_alpha as matrices
z_alpha_mat     <- draws_df %>% select(all_of(z_alpha_names)) %>% as.matrix()       # n_draws × (3*R)
sigma_alpha_mat <- draws_df %>% select(all_of(sigma_alpha_names)) %>% as.matrix()   # n_draws × 3

# Build a list of Lcorr_alpha matrices, one 3×3 per draw
Lcorr_alpha_list <- vector("list", length = n_draws)
for (s in seq_len(n_draws)) {
  L <- matrix(0, nrow = 3, ncol = 3)
  for (nm in Lcorr_alpha_names) {
    coords <- str_match(nm, "^Lcorr_alpha\\[(\\d+),(\\d+)\\]$")[, 2:3] %>% as.integer()
    L[coords[1], coords[2]] <- draws_df[[nm]][s]
  }
  Lcorr_alpha_list[[s]] <- L
}

# Reconstruct α_draws: array of dims (n_draws, R, 3)
alpha_draws <- array(
  NA_real_,
  dim = c(n_draws, stan_data$R, 3),
  dimnames = list(NULL, NULL, c("phi", "theta", "psi"))
)

for (s in seq_len(n_draws)) {
  D_alpha <- diag(sigma_alpha_mat[s, ])
  L_alpha <- Lcorr_alpha_list[[s]]
  M_alpha <- D_alpha %*% L_alpha   # 3×3
  z_s     <- matrix(z_alpha_mat[s, ], nrow = 3, ncol = stan_data$R, byrow = FALSE)
  alpha_draws[s, , ] <- t(M_alpha %*% z_s)  # R×3
}

# 7b. Similarly reconstruct δ (party intercepts)
z_delta_names     <- expand.grid(l = 1:3, q = 1:stan_data$Q) %>%
  mutate(name = paste0("z_delta[", l, ",", q, "]")) %>%
  pull(name)

sigma_delta_names <- paste0("sigma_delta[", 1:3, "]")

Lcorr_delta_names <- expand.grid(i = 1:3, j = 1:3) %>%
  filter(i >= j) %>%
  mutate(name = paste0("Lcorr_delta[", i, ",", j, "]")) %>%
  pull(name)

z_delta_mat     <- draws_df %>% select(all_of(z_delta_names)) %>% as.matrix()
sigma_delta_mat <- draws_df %>% select(all_of(sigma_delta_names)) %>% as.matrix()

Lcorr_delta_list <- vector("list", length = n_draws)
for (s in seq_len(n_draws)) {
  Ld <- matrix(0, nrow = 3, ncol = 3)
  for (nm in Lcorr_delta_names) {
    coords <- str_match(nm, "^Lcorr_delta\\[(\\d+),(\\d+)\\]$")[, 2:3] %>% as.integer()
    Ld[coords[1], coords[2]] <- draws_df[[nm]][s]
  }
  Lcorr_delta_list[[s]] <- Ld
}

delta_draws <- array(
  NA_real_,
  dim = c(n_draws, stan_data$Q, 3),
  dimnames = list(NULL, NULL, c("phi", "theta", "psi"))
)

for (s in seq_len(n_draws)) {
  D_delta <- diag(sigma_delta_mat[s, ])
  L_delta <- Lcorr_delta_list[[s]]
  M_delta <- D_delta %*% L_delta
  z_s     <- matrix(z_delta_mat[s, ], nrow = 3, ncol = stan_data$Q, byrow = FALSE)
  delta_draws[s, , ] <- t(M_delta %*% z_s)  # Q×3
}

# 7c. Compute posterior means for α and δ (for CSV outputs)
alpha_means <- apply(alpha_draws, c(2, 3), mean)  # R×3 matrix
delta_means <- apply(delta_draws, c(2, 3), mean)  # Q×3 matrix

# 7d. Build data frames with means
region_intercepts_df <- tibble(
  region_id   = 1:stan_data$R,
  region_name = region_levels,
  phi_alpha   = alpha_means[, 1],
  theta_alpha = alpha_means[, 2],
  psi_alpha   = alpha_means[, 3]
)

party_intercepts_df <- tibble(
  party_id    = 1:stan_data$Q,
  party_name  = party_levels,
  phi_delta   = delta_means[, 1],
  theta_delta = delta_means[, 2],
  psi_delta   = delta_means[, 3]
)

# 7e. Instead of computing 95% CIs, we will prepare full-draw data frames for plotting

# 7e.i. Region intercepts: long format of all draws for each region × latent
region_draws_long <- map_dfr(
  seq_len(stan_data$R),
  function(r) {
    tibble(
      region_id = r,
      draw      = 1:n_draws,
      phi       = alpha_draws[, r, 1],
      theta     = alpha_draws[, r, 2],
      psi       = alpha_draws[, r, 3]
    )
  }
) %>%
  pivot_longer(
    cols      = c(phi, theta, psi),
    names_to  = "latent",
    values_to = "value"
  ) %>%
  mutate(
    region_name = region_levels[region_id],
    latent = case_when(
      latent == "phi"   ~ "φ (Optimism)",
      latent == "theta" ~ "θ (Environment)",
      latent == "psi"   ~ "ψ (Radical-Reform)"
    )
  ) %>%
  select(region_name, latent, value)

# 7e.ii. Party intercepts: long format of all draws for each party × latent
party_draws_long <- map_dfr(
  seq_len(stan_data$Q),
  function(q) {
    tibble(
      party_id = q,
      draw     = 1:n_draws,
      phi      = delta_draws[, q, 1],
      theta    = delta_draws[, q, 2],
      psi      = delta_draws[, q, 3]
    )
  }
) %>%
  pivot_longer(
    cols      = c(phi, theta, psi),
    names_to  = "latent",
    values_to = "value"
  ) %>%
  mutate(
    party_name = party_levels[party_id],
    latent = case_when(
      latent == "phi"   ~ "φ (Optimism)",
      latent == "theta" ~ "θ (Environment)",
      latent == "psi"   ~ "ψ (Radical-Reform)"
    )
  ) %>%
  select(party_name, latent, value)

# 7f. Save raw intercept means to CSV (no CIs included)
write_csv(region_intercepts_df, "qoi_outputs/region_intercepts.csv")
write_csv(party_intercepts_df, "qoi_outputs/party_intercepts.csv")

cat("\n=== Region Intercepts (means only, first 6) ===\n")
print(head(region_intercepts_df, 6))
cat("\n=== Party Intercepts (means only, first 6) ===\n")
print(head(party_intercepts_df, 6))

# ---------------------------------------------------------------------------------------
# 7g. Group-Level Effect Plots: Distributions (regions and parties, 3 traits) -------------
# ---------------------------------------------------------------------------------------

# 7g.i. Region-level distribution plots (half-eye; unchanged)
region_plot_dist <- ggplot(
  region_draws_long,
  aes(
    x = value,
    y = fct_rev(region_name),
    fill = latent,
    colour = latent
  )
) +
  stat_halfeye(
    position   = position_nudge(y = 0.2),
    slab_alpha = 0.6,
    slab_size  = 0.5,
    adjust     = 0.8
  ) +
  facet_wrap(~ latent, scales = "free_x", ncol = 1) +
  theme_minimal() +
  labs(
    title    = "Region-Level Intercepts: Posterior Distributions",
    subtitle = "Latent traits: φ (Optimism), θ (Environment), ψ (Radical-Reform)",
    x        = "Intercept Value",
    y        = "Region"
  ) +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = "qoi_outputs/plots/region_group_effects_distributions.png",
  plot     = region_plot_dist,
  width    = 10,
  height   = 12,
  dpi      = 300
)

# 7g.ii. Party-level distribution plots (half-eye; unchanged)
party_plot_dist <- ggplot(
  party_draws_long,
  aes(
    x = value,
    y = fct_rev(party_name),
    fill = latent,
    colour = latent
  )
) +
  stat_halfeye(
    position   = position_nudge(y = 0.2),
    slab_alpha = 0.6,
    slab_size  = 0.5,
    adjust     = 0.8
  ) +
  facet_wrap(~ latent, scales = "free_x", ncol = 1) +
  theme_minimal() +
  labs(
    title    = "Party-Level Intercepts: Posterior Distributions",
    subtitle = "Latent traits: φ (Optimism), θ (Environment), ψ (Radical-Reform)",
    x        = "Intercept Value",
    y        = "Party"
  ) +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = "qoi_outputs/plots/party_group_effects_distributions.png",
  plot     = party_plot_dist,
  width    = 10,
  height   = 12,
  dpi      = 300
)

# ---------------------------------------------------------------------------------------
# 8. Residual Correlations Among Latent Traits ------------------------------------------------

# 8a. Extract Lcorr_eta draws (“Lcorr_eta[i,j]” for i>=j)
Lcorr_eta_names <- expand.grid(i = 1:3, j = 1:3) %>%
  filter(i >= j) %>%
  mutate(name = paste0("Lcorr_eta[", i, ",", j, "]")) %>%
  pull(name)

Lcorr_eta_list <- vector("list", length = n_draws)
for (s in seq_len(n_draws)) {
  Le <- matrix(0, nrow = 3, ncol = 3)
  for (nm in Lcorr_eta_names) {
    coords <- str_match(nm, "^Lcorr_eta\\[(\\d+),(\\d+)\\]$")[, 2:3] %>% as.integer()
    Le[coords[1], coords[2]] <- draws_df[[nm]][s]
  }
  Lcorr_eta_list[[s]] <- Le
}

# 8b. Compute correlation matrix R_eta = Lcorr_eta %*% t(Lcorr_eta) for each draw
R_eta_list <- lapply(Lcorr_eta_list, function(Le) {
  Ce <- Le %*% t(Le)  # this is a correlation matrix
  Ce
})

# 8c. Compute posterior mean correlation matrix (for later network plot)
R_eta_array <- simplify2array(R_eta_list)  # dims: 3×3×n_draws
R_eta_mean  <- apply(R_eta_array, c(1, 2), mean)

# 8d. Build data frame of residual correlations (means)
resid_corr_df <- tibble(
  pair        = c("φ–θ", "φ–ψ", "θ–ψ"),
  correlation = c(
    R_eta_mean[1, 2],
    R_eta_mean[1, 3],
    R_eta_mean[2, 3]
  )
)

write_csv(resid_corr_df, "qoi_outputs/residual_correlations.csv")
cat("\n=== Residual Correlations Among Latents (means) ===\n")
print(resid_corr_df)

# ---------------------------------------------------------------------------------------
# 8e. Posterior Density Plots of Latent Correlations (unchanged, already densities)
# ---------------------------------------------------------------------------------------

# 1. Compute point estimates (posterior means) for each pair:
point_estimates <- corr_long %>%
  group_by(pair) %>%
  summarize(est = mean(corr)) %>%
  ungroup()

# 2. Plot all three densities + overlay a vertical line at each point estimate:
corr_density_with_point <- ggplot(corr_long, aes(x = corr, fill = pair, color = pair)) +
  geom_density(alpha = 0.4, size = 1) +
  # Add a vertical dashed line at each posterior mean:
  geom_vline(
    data = point_estimates,
    aes(xintercept = est, color = pair),
    linetype = "dashed",
    size = 1
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Latent Correlations: Posterior Densities",
    x     = "Correlation",
    y     = "Density",
    fill  = "Pair",
    color = "Pair"
  ) +
  scale_fill_manual(
    values = c(
      "φ–θ" = "#666666",
      "φ–ψ" = "#6a3d9a",
      "θ–ψ" = "#ff7f00"
    )
  ) +
  scale_color_manual(
    values = c(
      "φ–θ" = "#666666",
      "φ–ψ" = "#6a3d9a",
      "θ–ψ" = "#ff7f00"
    )
  ) +
  theme(
    legend.position = "top",
    legend.title     = element_text(size = 12),
    legend.text      = element_text(size = 11)
  )

# 3. Save or display:
ggsave(
  filename = "qoi_outputs/plots/latent_corr_with_point_estimate.png",
  plot     = corr_density_with_point,
  width    = 8,
  height   = 6,
  dpi      = 300
)

# ---------------------------------------------------------------------------------------
# 9. Measurement Reliability (McDonald’s ω) & Communalities -----------------------------------

# 9a. Extract λ and σ for each block
lambda_opt_names <- paste0("lambda_opt[", 1:stan_data$J_opt, "]")
sigma_opt_names  <- paste0("sigma_opt[", 1:stan_data$J_opt, "]")

lambda_env_names <- paste0("lambda_env[", 1:stan_data$J_env, "]")
sigma_env_names  <- paste0("sigma_env[", 1:stan_data$J_env, "]")

lambda_rad_names <- paste0("lambda_rad[", 1:stan_data$J_rad, "]")
sigma_rad_names  <- paste0("sigma_rad[", 1:stan_data$J_rad, "]")

lambda_opt_mat <- draws_df %>% select(all_of(lambda_opt_names)) %>% as.matrix()
sigma_opt_mat  <- draws_df %>% select(all_of(sigma_opt_names)) %>% as.matrix()

lambda_env_mat <- draws_df %>% select(all_of(lambda_env_names)) %>% as.matrix()
sigma_env_mat  <- draws_df %>% select(all_of(sigma_env_names)) %>% as.matrix()

lambda_rad_mat <- draws_df %>% select(all_of(lambda_rad_names)) %>% as.matrix()
sigma_rad_mat  <- draws_df %>% select(all_of(sigma_rad_names)) %>% as.matrix()

# 9b. Compute ω for each draw & each block
omega_opt_draws <- apply(
  cbind(lambda_opt_mat, sigma_opt_mat),
  1,
  function(row) {
    J <- stan_data$J_opt
    λ <- row[1:J]
    σ <- row[(J+1):(2*J)]
    num <- sum(λ)^2
    den <- num + sum(σ^2)
    num / den
  }
)

omega_env_draws <- apply(
  cbind(lambda_env_mat, sigma_env_mat),
  1,
  function(row) {
    J <- stan_data$J_env
    λ <- row[1:J]
    σ <- row[(J+1):(2*J)]
    num <- sum(λ)^2
    den <- num + sum(σ^2)
    num / den
  }
)

omega_rad_draws <- apply(
  cbind(lambda_rad_mat, sigma_rad_mat),
  1,
  function(row) {
    J <- stan_data$J_rad
    λ <- row[1:J]
    σ <- row[(J+1):(2*J)]
    num <- sum(λ)^2
    den <- num + sum(σ^2)
    num / den
  }
)

omega_df <- tibble(
  draw       = 1:n_draws,
  omega_opt  = omega_opt_draws,
  omega_env  = omega_env_draws,
  omega_rad  = omega_rad_draws
)

# 9b.ii. Summarize posterior means (skip CIs)
omega_summary <- omega_df %>%
  pivot_longer(
    cols      = c(omega_opt, omega_env, omega_rad),
    names_to  = "block",
    values_to = "omega"
  ) %>%
  group_by(block) %>%
  summarize(
    mean_omega  = mean(omega),
    .groups     = "drop"
  )

write_csv(omega_summary, "qoi_outputs/reliability_omega.csv")
cat("\n=== McDonald’s ω (Reliability Means Only) ===\n")
print(omega_summary)

# 9c. Communalities & uniquenesses for each item: compute full-draw distributions
summarize_items_distributions <- function(lambda_mat, sigma_mat, block_name) {
  J <- ncol(lambda_mat)
  map_dfr(1:J, function(j) {
    λ_draws    <- lambda_mat[, j]
    σ_draws    <- sigma_mat[, j]
    commun     <- λ_draws^2
    uniqueness <- σ_draws^2 / (λ_draws^2 + σ_draws^2)
    tibble(
      block         = block_name,
      item          = j,
      commun_draw   = list(commun),
      uniq_draw     = list(uniqueness)
    )
  })
}

opt_communi_dist <- summarize_items_distributions(lambda_opt_mat, sigma_opt_mat, "Optimism")
env_communi_dist <- summarize_items_distributions(lambda_env_mat, sigma_env_mat, "Environment")
rad_communi_dist <- summarize_items_distributions(lambda_rad_mat, sigma_rad_mat, "Radical-Reform")

comm_uni_draws <- bind_rows(opt_communi_dist, env_communi_dist, rad_communi_dist)

# For reference, also compute summary means (skip CIs)
comm_uni_summary <- comm_uni_draws %>%
  unnest(c(commun_draw, uniq_draw)) %>%
  group_by(block, item) %>%
  summarize(
    comm_mean = mean(commun_draw),
    uniq_mean = mean(uniq_draw),
    .groups   = "drop"
  )

write_csv(comm_uni_summary, "qoi_outputs/communalities_uniquenesses.csv")
cat("\n=== Communalities & Uniquenesses (Means Only, first 6 rows) ===\n")
print(head(comm_uni_summary, 6))

# ---------------------------------------------------------------------------------------
# 9d. Item Discrimination (λ) Distributions (unchanged; half-eye) ------------------------
# ---------------------------------------------------------------------------------------

# 9d.i. Build a long data frame of all λ draws for each item × latent
lambda_draws_long <- bind_rows(
  map_dfr(1:stan_data$J_opt, function(j) {
    tibble(
      item   = j,
      latent = "φ (Optimism)",
      value  = lambda_opt_mat[, j]
    )
  }),
  map_dfr(1:stan_data$J_env, function(j) {
    tibble(
      item   = j,
      latent = "θ (Environment)",
      value  = lambda_env_mat[, j]
    )
  }),
  map_dfr(1:stan_data$J_rad, function(j) {
    tibble(
      item   = j,
      latent = "ψ (Radical-Reform)",
      value  = lambda_rad_mat[, j]
    )
  })
)

# 9d.ii. Plot full posterior distributions of λ for each item × latent
lambda_dist_plot <- ggplot(
  lambda_draws_long,
  aes(
    x      = value,
    y      = factor(item),
    fill   = latent,
    colour = latent        
  )
) +
  stat_halfeye(
    position      = position_nudge(y = 0.2),
    slab_alpha    = 0.6,
    slab_size     = 0.5,
    adjust        = 0.8,
  ) +
  facet_wrap(~ latent, scales = "free_x", ncol = 1) +
  theme_minimal() +
  labs(
    title    = "Item Discrimination (λ): Posterior Distributions",
    subtitle = "Latent traits: φ (Optimism), θ (Environment), ψ (Radical-Reform)",
    x        = "Discrimination (λ)",
    y        = "Item Index"
  ) +
  theme(
    legend.position = "none"
  )

ggsave(
  filename = "qoi_outputs/plots/item_discrimination_lambda_distributions.png",
  plot     = lambda_dist_plot,
  width    = 8,
  height   = 10,
  dpi      = 300
)

# ---------------------------------------------------------------------------------------
# 10. Visualization: Region-Level Portraits (distributions reflected in error bars removed)
# ---------------------------------------------------------------------------------------

# 10a. 3D Scatterplot of (φ, θ, ψ) for regions (posterior means) — unchanged
region_plot_df <- region_intercepts_df %>%
  select(region_name, phi_alpha, theta_alpha, psi_alpha)

# Compute number of observations per region 
region_counts <- tibble(region_id = stan_data$region_id) %>%
  count(region_id, name = "n_obs")

region_plot_df <- region_intercepts_df %>%
  left_join(region_counts, by = "region_id")

# 10b. Pairwise scatterplots with means only (omit error bars)

# 1) φ vs ψ (x = ψ, y = φ), color = θ
p_psi_phi_col_env <- ggplot(region_plot_df, aes(
  x     = psi_alpha,
  y     = phi_alpha,
  color = theta_alpha,
  size  = n_obs
)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = region_name), size = 3, max.overlaps = 15) +
  scale_color_gradient(
    low   = "grey80",
    high  = "darkgreen",
    name  = "Standardised θ\n(Environment)"
  ) +
  scale_size_continuous(
    range = c(3, 12),
    name  = "Respondents\n(n_obs)"
  ) +
  labs(
    title = "Regions: ψ (Radical-Reform) vs φ (Optimism), colored by θ (Environment)",
    x     = "Standardised ψ (intercept)",
    y     = "Standardised φ (intercept)"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = "qoi_outputs/plots/region_psi_vs_phi_col_env.png",
  plot     = p_psi_phi_col_env,
  width    = 8,
  height   = 6,
  dpi      = 300
)

# 2) φ vs θ (x = θ, y = φ), color = ψ
p_theta_phi_col_rad <- ggplot(region_plot_df, aes(
  x     = theta_alpha,
  y     = phi_alpha,
  color = psi_alpha,
  size  = n_obs
)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = region_name), size = 3, max.overlaps = 15) +
  scale_color_gradient(
    low   = "grey80",
    high  = "red",
    name  = "Standardised ψ\n(Radical-Reform)"
  ) +
  scale_size_continuous(
    range = c(3, 12),
    name  = "Respondents\n(n_obs)"
  ) +
  labs(
    title = "Regions: θ (Environment) vs φ (Optimism), colored by ψ (Radical-Reform)",
    x     = "Standardised θ (intercept)",
    y     = "Standardised φ (intercept)"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = "qoi_outputs/plots/region_theta_vs_phi_col_rad.png",
  plot     = p_theta_phi_col_rad,
  width    = 8,
  height   = 6,
  dpi      = 300
)

# 3) θ vs ψ (x = ψ, y = θ), color = φ
p_psi_theta_col_opt <- ggplot(region_plot_df, aes(
  x     = psi_alpha,
  y     = theta_alpha,
  color = phi_alpha,
  size  = n_obs
)) +
  geom_point(alpha = 0.8) +
  geom_text_repel(aes(label = region_name), size = 3, max.overlaps = 15) +
  scale_color_gradient(
    low   = "grey80",
    high  = "blue",
    name  = "Standardised φ\n(Optimism)"
  ) +
  scale_size_continuous(
    range = c(3, 12),
    name  = "Respondents\n(n_obs)"
  ) +
  labs(
    title = "Regions: ψ (Radical-Reform) vs θ (Environment), colored by φ (Optimism)",
    x     = "Standardised ψ (intercept)",
    y     = "Standardised θ (intercept)"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = "qoi_outputs/plots/region_psi_vs_theta_col_opt.png",
  plot     = p_psi_theta_col_opt,
  width    = 8,
  height   = 6,
  dpi      = 300
)

# 3D scatter, sizing by actual number of respondents (n_obs) — unchanged
fig_3d <- plot_ly(
  data   = region_plot_df,
  x      = ~phi_alpha,
  y      = ~theta_alpha,
  z      = ~psi_alpha,
  size   = ~n_obs,
  text   = ~region_name,
  mode   = "markers",
  type   = "scatter3d",
  marker = list(
    sizemode = "diameter",
    sizeref  = 2 * max(region_plot_df$n_obs) / (20^2), 
    opacity  = 0.8
  )
) %>%
  layout(
    title = "3D Region Portraits: Standardised φ vs θ vs ψ",
    scene = list(
      xaxis = list(title = "Standardised φ (Optimism)"),
      yaxis = list(title = "Standardised θ (Environment)"),
      zaxis = list(title = "Standardised ψ (Radical-Reform)")
    )
  )

htmlwidgets::saveWidget(
  fig_3d,
  file = "qoi_outputs/plots/region_portraits_3d.html",
  selfcontained = TRUE
)

# ---------------------------------------------------------------------------------------
# 11. Visualization: Party-Level Radar Chart ------------------------------------------------

# 11a. Prepare data for radar plot: one row per party, columns φ, θ, ψ (means only)
party_radar_df <- party_intercepts_df %>%
  select(party_name, phi_delta, theta_delta, psi_delta) %>%
  column_to_rownames(var = "party_name")

# 11b. Helper to create radar chart for a single party — unchanged
radar_plot_list <- map(
  rownames(party_radar_df),
  function(p_name) {
    df <- tibble(
      trait = c("Optimism", "Environment", "Radical-Reform"),
      value = as.numeric(party_radar_df[p_name, ])
    )
    df <- rbind(df, df[1, ])  # close the polygon
    ggplot(df, aes(x = trait, y = value, group = 1)) +
      geom_polygon(fill = "#66c2a5", alpha = 0.4) +
      geom_line(color = "#1b9e77", size = 1) +
      geom_point(size = 2) +
      ylim(min(df$value) - 0.5, max(df$value) + 0.5) +
      labs(
        title = paste0("Party Portrait: ", p_name),
        y     = "Intercept Value"
      ) +
      theme_minimal() +
      theme(
        axis.title.x = element_blank(),
        axis.text.x  = element_text(size = 10),
        plot.title   = element_text(size = 12, face = "bold")
      )
  }
)

# 11c. Arrange radar plots in a grid — unchanged
n_parties <- nrow(party_radar_df)
ncols     <- 2
nrows     <- ceiling(n_parties / ncols)
radar_grid <- plot_grid(plotlist = radar_plot_list, ncol = ncols)

ggsave(
  filename = "qoi_outputs/plots/party_radar_charts.png",
  plot     = radar_grid,
  width    = 12,
  height   = 6 * nrows / 2,
  dpi      = 300
)

# ---------------------------------------------------------------------------------------
# 12. Visualization: Covariate Effects Heatmap (posterior means only) -------------------------

ame_wide <- slopes_draws_long %>%
  select(latent, covariate, value) %>%
  group_by(latent, covariate) %>%
  summarize(mean = mean(value), .groups = "drop") %>%
  pivot_wider(
    names_from  = latent,
    values_from = mean
  )

ame_long <- ame_wide %>%
  pivot_longer(
    cols      = -covariate,            # picks up the three latent columns
    names_to  = "latent",
    values_to = "slope"
  )

heatmap_cov <- ggplot(ame_long, aes(x = latent, y = covariate, fill = slope)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low      = "red",
    mid      = "white",
    high     = "blue",
    midpoint = 0,
    name     = "Slope"
  ) +
  labs(
    title    = "Heatmap of Posterior Mean Slopes (AME)",
    subtitle = "Latent traits: φ (Optimism), θ (Environment), ψ (Radical-Reform)",
    x        = "Latent Trait",
    y        = "Covariate"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

print(heatmap_cov)

ggsave(
  filename = "qoi_outputs/plots/covariate_effects_heatmap.png",
  plot     = heatmap_cov,
  width    = 8,
  height   = 6,
  dpi      = 300
)

# ---------------------------------------------------------------------------------------
# 13. Visualization: Latent Correlation Network ------------------------------------------------

# 13a. Build a 3×3 correlation matrix of latent correlations (R_eta_mean)
corr_mat <- R_eta_mean
colnames(corr_mat) <- c("Optimism", "Environment", "Radical-Reform")
rownames(corr_mat) <- c("Optimism", "Environment", "Radical-Reform")

# 13b. Convert to long format for plotting
corr_long_net <- as_tibble(corr_mat, rownames = "trait1") %>%
  pivot_longer(-trait1, names_to = "trait2", values_to = "corr")

# 13c. Create correlation network graph — unchanged
edges <- corr_long_net %>%
  filter(trait1 != trait2, abs(corr) > 0.2)

graph <- graph_from_data_frame(edges, vertices = tibble(name = rownames(corr_mat)))

corr_network_plot <- ggraph(graph, layout = "circle") +
  geom_edge_link(aes(width = abs(corr), color = corr), alpha = 0.8) +
  geom_node_point(size = 10, color = "#66a61e") +
  geom_node_text(aes(label = name), repel = TRUE) +
  scale_edge_color_gradient2(
    low      = "red",
    mid      = "grey80",
    high     = "blue",
    midpoint = 0,
    name     = "r"
  ) +
  scale_edge_width(range = c(0.5, 2), guide = "none") +
  labs(
    title = "Residual Correlation Network Among Latent Traits"
  ) +
  theme_void()

ggsave(
  filename = "qoi_outputs/plots/latent_correlation_network.png",
  plot     = corr_network_plot,
  width    = 6,
  height   = 6,
  dpi      = 300
)

# ---------------------------------------------------------------------------------------
# 14. Visualization: Predicted Latent Distributions for Key Demographics ----------------------

# 14a. Define a function to simulate latent distribution for a hypothetical demographic profile
simulate_latent_profile <- function(profile_vec, region = NULL, party = NULL) {
  # profile_vec: named vector length P
  # region, party: integer IDs (1..R, 1..Q). If NULL, set intercepts = 0.
  
  # Extract posterior draws of B (3×P), tau_eta, and Lcorr_eta_list
  B_draws <- array(NA, dim = c(n_draws, 3, P))
  for (p in 1:P) {
    for (ℓ in 1:3) {
      nm <- paste0("B[", ℓ, ",", p, "]")
      B_draws[, ℓ, p] <- draws_df[[nm]]
    }
  }
  tau_eta_mat  <- draws_df %>% select(starts_with("tau_eta")) %>% as.matrix()  # n_draws × 3
  Lcorr_eta_ls <- Lcorr_eta_list  # from section 8
  
  sim_results <- tibble(
    draw      = 1:n_draws,
    phi_sim   = NA_real_,
    theta_sim = NA_real_,
    psi_sim   = NA_real_
  )
  
  for (s in 1:n_draws) {
    # Region & party intercepts
    alpha_vec <- if (!is.null(region)) alpha_draws[s, region, ] else rep(0, 3)
    delta_vec <- if (!is.null(party))  delta_draws[s, party, ]  else rep(0, 3)
    
    # Covariate linear part
    beta_mat <- B_draws[s, , ]               # 3×P
    x_vec    <- profile_vec                  # length P
    cov_part <- beta_mat %*% x_vec           # 3-vector
    
    # Latent noise
    z_vec    <- rnorm(3)
    M_eta    <- diag(tau_eta_mat[s, ]) %*% Lcorr_eta_ls[[s]]  # 3×3
    noise    <- M_eta %*% z_vec                                # 3-vector
    
    eta      <- alpha_vec + delta_vec + cov_part + noise
    sim_results$phi_sim[s]   <- eta[1]
    sim_results$theta_sim[s] <- eta[2]
    sim_results$psi_sim[s]   <- eta[3]
  }
  
  return(sim_results)
}

# ── 1. Define 10 descriptive profiles ───────────────────────────────────────────────────

make_X <- function(
    gender = 0,
    age_bracket = c("18-24","25-34","35-44","45-54","55-64","65+"),
    edu_level   = c("No qualifications","Level 1","Level 2","Level 3","Level 4+","Apprenticeship","Other"),
    insecurity  = 0
) {
  v <- numeric(13)
  names(v) <- c(
    "gender",
    "age_65plus","age_55_64","age_45_54","age_35_44","age_25_34",
    "edu_L4plus","edu_L3","edu_L2","edu_L1","edu_appr","edu_other",
    "insecurity"
  )
  v["gender"] <- gender
  age_map <- list(
    "65+"   = "age_65plus",
    "55-64" = "age_55_64",
    "45-54" = "age_45_54",
    "35-44" = "age_35_44",
    "25-34" = "age_25_34",
    "18-24" = NULL
  )
  if (!age_bracket %in% names(age_map)) stop("Invalid age_bracket")
  if (!is.null(age_map[[age_bracket]])) {
    v[ age_map[[age_bracket]] ] <- 1
  }
  edu_map <- list(
    "Level 4+"       = "edu_L4plus",
    "Level 3"        = "edu_L3",
    "Level 2"        = "edu_L2",
    "Level 1"        = "edu_L1",
    "Apprenticeship" = "edu_appr",
    "Other"          = "edu_other",
    "No qualifications" = NULL
  )
  if (!edu_level %in% names(edu_map)) stop("Invalid edu_level")
  if (!is.null(edu_map[[edu_level]])) {
    v[ edu_map[[edu_level]] ] <- 1
  }
  v["insecurity"] <- insecurity
  return(v)
}

profiles <- tibble(
  profile_name = c(
    "Elderly Male, Reform UK, No Qualifications, High Insecurity, South East",
    "Young Female, Green Party, Univ Grad, Low Insecurity, London",
    "Mid-Age Male, Conservative, Level 1 Edu, Moderate Insecurity, East Midlands",
    "Mid-Age Female, Labour, Univ Grad, Low Insecurity, North West",
    "Senior Female, SNP, No Qualifications, High Insecurity, Scotland",
    "Senior Male, Lib Dem, Univ Grad, Low Insecurity, South West",
    "Young Male, Conservative, Level 2 Edu, High Insecurity, North West",
    "Young Female, Reform UK, Level 3 Edu, Low Insecurity, Yorkshire & the Humber",
    "Mid-Age Female, Plaid Cymru, Apprenticeship, Moderate Insecurity, Wales",
    "Senior Male, Don’t Know, Other Qualifications, Moderate Insecurity, Eastern"
  ),
  gender = c(0, 1, 0, 1, 1, 0, 0, 1, 1, 0),
  age_bracket = c(
    "65+",   # profile 1
    "18-24", # profile 2
    "35-44", # profile 3
    "35-44", # profile 4
    "55-64", # profile 5
    "55-64", # profile 6
    "18-24", # profile 7
    "18-24", # profile 8
    "45-54", # profile 9
    "65+"    # profile 10
  ),
  edu_level = c(
    "No qualifications",  # profile 1
    "Level 4+",           # profile 2
    "Level 1",            # profile 3
    "Level 4+",           # profile 4
    "No qualifications",  # profile 5
    "Level 4+",           # profile 6
    "Level 2",            # profile 7
    "Level 3",            # profile 8
    "Apprenticeship",     # profile 9
    "Other"               # profile 10
  ),
  insecurity = c(
    1.5,   # very high insecurity
    -1.2,  # very low insecurity
    0.5,   # moderate insecurity
    -0.8,  # low insecurity
    1.3,   # high insecurity
    -1.0,  # low insecurity
    1.4,   # very high insecurity
    -1.1,  # very low insecurity
    0.6,   # moderate insecurity
    0.4    # moderate insecurity
  ),
  region_id = c(10, 9, 8, 5, 3, 7, 5, 1, 4, 6),
  party_id  = c(7, 1, 6, 2, 4, 5, 6, 7, 3, 9)
)

# Build the X-vector for each profile
profiles <- profiles %>%
  rowwise() %>%
  mutate(
    X_vec = list(make_X(
      gender      = gender,
      age_bracket = age_bracket,
      edu_level   = edu_level,
      insecurity  = insecurity
    ))
  ) %>%
  ungroup()

# ── 2. Simulate posterior draws for each profile ──────────────────────────────────────────
profiles <- profiles %>%
  rowwise() %>%
  mutate(
    sim_df = list(simulate_latent_profile(X_vec, region_id, party_id))
  ) %>%
  ungroup()

# ── 3. Prepare combined data for plotting ────────────────────────────────────────────────

# 3a. φ-draws combined
phi_combined <- profiles %>%
  select(profile_name, sim_df) %>%
  unnest(sim_df) %>%
  select(profile_name, phi_sim) %>%
  rename(value = phi_sim) %>%
  mutate(latent = "φ (Optimism)")

# 3b. θ-draws combined
theta_combined <- profiles %>%
  select(profile_name, sim_df) %>%
  unnest(sim_df) %>%
  select(profile_name, theta_sim) %>%
  rename(value = theta_sim) %>%
  mutate(latent = "θ (Environment)")

# 3c. ψ-draws combined
psi_combined <- profiles %>%
  select(profile_name, sim_df) %>%
  unnest(sim_df) %>%
  select(profile_name, psi_sim) %>%
  rename(value = psi_sim) %>%
  mutate(latent = "ψ (Radical-Reform)")

# 3d. Bind all three
all_combined <- bind_rows(phi_combined, theta_combined, psi_combined)

# Factorize ‘profile_name’ so that ordering is the same as in `profiles`
all_combined$profile_name <- factor(
  all_combined$profile_name,
  levels = profiles$profile_name
)

# ── 4. Plot overlapping densities: one facet per latent dimension ─────────────────────────

density_plot <- ggplot(all_combined, aes(x = value, color = profile_name, fill = profile_name)) +
  geom_density(alpha = 0.3, size = 0.5) +
  facet_wrap(~ latent, ncol = 1, scales = "free") +
  scale_color_brewer(palette = "Dark2", name = "Profile") +
  scale_fill_brewer(palette = "Dark2", name = "Profile") +
  labs(
    title = "Predicted Latent Distributions Across Ten Profiles",
    x     = "Latent Value",
    y     = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    strip.text      = element_text(face = "bold", size = 12),
    axis.text.x     = element_text(size = 10),
    axis.text.y     = element_text(size = 10)
  )

ggsave(
  filename = "qoi_outputs/plots/multiple_profiles_latent_densities.png",
  plot     = density_plot,
  width    = 8,
  height   = 14,
  dpi      = 300
)

# ── 5. Print a summary of profiles and remind where the plot is saved ───────────────────

cat("\n============================================\n")
cat("Generated overlapping latent-density plots for 10 profiles.\n")
cat("Each facet shows one latent (φ, θ, or ψ) with all profiles overlaid.\n")
cat("See: qoi_outputs/plots/multiple_profiles_latent_densities.png\n")
cat("Profiles:\n")
print(profiles %>% select(profile_name, gender, age_bracket, edu_level, insecurity, region_id, party_id))
cat("============================================\n\n")

# ---------------------------------------------------------------------------------------
# SEPARATE OVERLAID LATENT DENSITY PLOTS FOR EACH PROFILE (unchanged) -------------------
# ---------------------------------------------------------------------------------------

# 1. Create a directory to hold per-profile plots
if (!dir.exists("qoi_outputs/plots/profiles_separate")) {
  dir.create("qoi_outputs/plots/profiles_separate", recursive = TRUE)
}

# 2. Define colors for each latent
latent_colors <- c(
  "φ (Optimism)"       = "lightblue",
  "θ (Environment)"    = "darkgreen",
  "ψ (Radical-Reform)" = "salmon"
)

# 3. Loop over each profile and generate one plot with three overlaid densities
all_vals <- bind_rows(
  lapply(profiles$sim_df, function(df) {
    tibble(
      value = c(df$phi_sim, df$theta_sim, df$psi_sim)
    )
  })
)
global_min <- min(all_vals$value)
global_max <- max(all_vals$value)
y_max <- 1.1

for (i in seq_len(nrow(profiles))) {
  prof_name <- profiles$profile_name[i]
  sim_df    <- profiles$sim_df[[i]]
  
  sim_long <- sim_df %>%
    select(phi_sim, theta_sim, psi_sim) %>%
    pivot_longer(
      cols      = everything(),
      names_to  = "latent",
      values_to = "value"
    ) %>%
    mutate(
      latent = case_when(
        latent == "phi_sim"   ~ "φ (Optimism)",
        latent == "theta_sim" ~ "θ (Environment)",
        latent == "psi_sim"   ~ "ψ (Radical-Reform)",
        TRUE                  ~ latent
      )
    )
  
  p <- ggplot(sim_long, aes(x = value, color = latent, fill = latent)) +
    geom_density(alpha = 0.3, size = 0.5) +
    scale_color_manual(values = latent_colors, name = "Latent") +
    scale_fill_manual(values = latent_colors, name = "Latent") +
    xlim(global_min, global_max) +
    ylim(0, y_max) +
    labs(
      title = paste0("Latent Distributions for Profile:\n", prof_name),
      x     = "Latent Value",
      y     = "Density"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      plot.title      = element_text(face = "bold", hjust = 0.5)
    )
  
  filename <- prof_name %>%
    str_replace_all("[^A-Za-z0-9_]+", "_") %>%
    str_replace_all("_+", "_") %>%
    str_replace_all("^_|_$", "") %>%
    paste0(".png")
  
  ggsave(
    file   = file.path("qoi_outputs/plots/profiles_separate", filename),
    plot   = p,
    width  = 6,
    height = 4,
    dpi    = 300
  )
}

cat("\n============================================\n")
cat("Saved separate overlaid latent-density plots for each profile under:\n")
cat("  qoi_outputs/plots/profiles_separate/\n")
cat("============================================\n\n")

# ---------------------------------------------------------------------------------------
# 15. Save Key Data Frames ------------------------------------------------------------------

write_csv(region_intercepts_df, "qoi_outputs/region_intercepts.csv")
write_csv(party_intercepts_df,  "qoi_outputs/party_intercepts.csv")
write_csv(latent_scores_df,     "qoi_outputs/latent_scores.csv")
write_csv(ame_df,               "qoi_outputs/ame_slopes.csv")
write_csv(r2_summary,           "qoi_outputs/r2_summary.csv")
write_csv(omega_summary,        "qoi_outputs/reliabilities.csv")
write_csv(comm_uni_summary,     "qoi_outputs/communalities_uniquenesses.csv")
write_csv(resid_corr_df,        "qoi_outputs/residual_correlations.csv")
write_csv(slopes_df,            "qoi_outputs/slopes_summary.csv")

# ---------------------------------------------------------------------------------------
# 16. Completion Message -------------------------------------------------------------------

cat("\n=========================================\n")
cat("Quantities of Interest & Visualizations Generated Successfully.\n")
cat("Covariate-effects plot now uses overlapping ridgelines (with Script 1’s original colors).\n")
cat("Browse the “qoi_outputs/” folder for CSV summaries and “qoi_outputs/plots/” for figures.\n")
cat("=========================================\n\n")
