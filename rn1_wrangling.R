# ---------------------------------------------------------------------------------------
# BAYESIAN MODELING Research Note
# Henry Baker
# ---------------------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(cmdstanr)
library(stringr)
library(posterior)
library(bayesplot)
library(ggplot2)
library(ggdist)
library(forcats)

# ── Helpers ────────────────────────────────────────────────────────────────────────────
recode_map <- function(x, map) map[x]
z_score    <- function(x) (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)

#### 0. Data Loading =============================================================================

setwd("/Users/henrybaker/repositories/bayesian_modeling/research_note")

data <- read_excel("data/LfG Data.xlsx") 
data <- data %>% mutate(obs_id = row_number())

# fit_fast - for model dev:
# set.seed(42)
#sample_ids <- sample(data$obs_id, size = 300)
#data <- data %>% filter(obs_id %in% sample_ids)

# ── 1. Economic Optimism (ϕᵢ) ────────────────────────────────────────────────────────

# Identify the six optimism columns (Q43–Q48)
opt_cols <- grep(
  "^To what extent do you feel optimistic or pessimistic about the following",
  names(data),
  value = TRUE
)

# Recode text → numeric on {−2,−1,0,1,2}
opt_map <- c(
  "Very pessimistic"                   = -2L,
  "Fairly pessimistic"                 = -1L,
  "Neither optimistic nor pessimistic" =  0L,
  "Fairly optimistic"                  =  1L,
  "Very optimistic"                    =  2L
)

# Step 1: Create raw numeric columns
data <- data %>%
  mutate(across(all_of(opt_cols),
                ~ recode_map(., opt_map),
                .names = "num_{.col}"
  ))

# Step 2: Z‐score each of those six
num_opt_cols <- paste0("num_", opt_cols)
data <- data %>%
  mutate(across(all_of(num_opt_cols),
                ~ z_score(.),
                .names = "z_{.col}"
  ))

# Now gather into long format
z_opt_cols <- paste0("z_", num_opt_cols)

opt_long_df <- data %>%
  select(obs_id, all_of(z_opt_cols)) %>%
  pivot_longer(
    -obs_id,
    names_to  = "item_opt",
    values_to = "y_opt"
  ) %>%
  filter(!is.na(y_opt)) %>%
  mutate(
    j_opt = as.integer(factor(item_opt, levels = z_opt_cols))
  ) %>%
  select(obs_id, j_opt, y_opt)

# opt_cols is a character vector of your six “Q43…Q48” column names
# The z‐scored versions live in columns named "z_num_<Qxx>"

# Compute lower and upper bounds for each of the six optimism items:
lower_opt <- sapply(opt_cols, function(col) {
  zname <- paste0("z_num_", col)
  min(data[[zname]], na.rm = TRUE)
})
upper_opt <- sapply(opt_cols, function(col) {
  zname <- paste0("z_num_", col)
  max(data[[zname]], na.rm = TRUE)
})

# Note: j_opt runs from 1..6, corresponding to the six optimism items

# ── 2. Environmental Priority (θᵢ) ─────────────────────────────────────────────────

#  2a) Forced choice Q121
env_forced_col <- grep("\\.\\.\\.121$", names(data), value = TRUE)

#  2b) Four Likert items Q122–Q125
env_likert_cols <- grep(
  "^To what extent do you agree or disagree with the following statements",
  names(data),
  value = TRUE
)

# Recode Likert {Strongly disagree → -2, …, Strongly agree → +2}
env_map <- c(
  "Strongly disagree"           = -2L,
  "Somewhat disagree"           = -1L,
  "Neither agree nor disagree"  =  0L,
  "Somewhat agree"              =  1L,
  "Strongly agree"              =  2L
)

# Which of those four needs reversing?  (the “cost‐first” statement is the 2nd of the four)
anti_col <- env_likert_cols[2]

data <- data %>%
  # 1) Recode all four to {−2,…,+2}
  mutate(across(all_of(env_likert_cols),
                ~ recode_map(., env_map),
                .names = "num_{.col}"
  )) %>%
  # 2) Reverse‐score the anti‐environment item
  mutate(
    !!paste0("num_", anti_col) := -get(paste0("num_", anti_col))
  ) %>%
  # 3) Recode the forced‐choice into {0,1}
  mutate(
    num_forced = case_when(
      .data[[env_forced_col]] ==
        "We should prioritise protecting the natural environment, even if that sometimes prevents action to reduce the cost of living" ~ 1L,
      .data[[env_forced_col]] ==
        "We should prioritise action to reduce the cost of living, even if that sometimes comes at the cost of the environment"      ~ 0L,
      TRUE ~ NA_integer_
    )
  ) %>%
  # 4) Z‐score all five (forced + four recoded Likert)
  mutate(across(
    c("num_forced", paste0("num_", env_likert_cols)),
    ~ z_score(.),
    .names = "z_{.col}"
  ))

# Build long format for environment
num_env_cols <- c("num_forced", paste0("num_", env_likert_cols))
z_env_cols   <- paste0("z_", num_env_cols)

env_long_df <- data %>%
  select(obs_id, all_of(z_env_cols)) %>%
  pivot_longer(
    -obs_id,
    names_to  = "item_env",
    values_to = "y_env"
  ) %>%
  filter(!is.na(y_env)) %>%
  mutate(
    j_env = as.integer(factor(item_env, levels = z_env_cols))
  ) %>%
  select(obs_id, j_env, y_env)

# Now j_env runs from 1..5

# opt_cols is a character vector of your six “Q43…Q48” column names
# The z‐scored versions live in columns named "z_num_<Qxx>"

# Compute lower and upper bounds for each of the environmentalism items:
lower_env <- sapply(z_env_cols, function(col) {
  min(data[[col]], na.rm = TRUE)
})
upper_env <- sapply(z_env_cols, function(col) {
  max(data[[col]], na.rm = TRUE)
})

# ── 3. Radical‐Reform (ψᵢ) ────────────────────────────────────────────────────────

# The eight forced‐choice questions Q72–Q75, Q78–Q79, Q82–Q83 each yield one binary indicator.
# We assume the raw text in each of these columns exactly matches the “radical” vs “status‐quo” phrasing.

# 3a) Identify the eight “For the following pairs…” columns (Q72, Q73, Q74, Q75, Q78, Q79, Q82, Q83).
#     Each of these column names ends with “…72”, “…73”, etc., so we match those suffixes:
rad_cols <- colnames(data)[c(72, 73, 74, 75, 78, 79, 82, 83)]
length(rad_cols)    # should be 8
rad_cols

# 3b) Create a named vector that maps each column to “1” if the respondent chose the radical statement,
#     or “0” if they chose the status‐quo statement.  We need to inspect the exact text in data to fill in both sides.
#     Below, replace the right‐hand strings with the exact phrasing from your Excel file.

# First, the original mapping (keys = "Q72","Q73",…,"Q83"):
radical_map_list <- list(
  Q72 = c(
    "Westminster politicians have moved too slowly on initiatives to help the economy" = 1L,
    "Westminster politicians have generally moved at a good pace to get things done to help the economy" = 0L
  ),
  Q73 = c(
    "I am looking for politicians who show that they can get things done" = 1L,
    "I am looking for politicians who take their time and work for long term goals even if that means there is less getting done in the short term" = 0L
  ),
  Q74 = c(
    "I would be more likely to support politicians who show they are prepared to take radical action to improve everyday people’s lives" = 1L,
    "I would be more likely to support politicians who do not want to change the system and keep things broadly the same with no surprises" = 0L
  ),
  Q75 = c(
    "Britain has had a rough few years but is broadly doing fine and just needs to have a steady few years" = 0L,
    "Britain is broken and needs radical action to fix it" = 1L
  ),
  Q78 = c(
    "The UK is not taking the right steps to achieve more economic growth and raise living standards" = 1L,
    "The UK is taking the right steps to achieve more economic growth and raise living standards" = 0L
  ),
  Q79 = c(
    "The UK is broadly speaking on the wrong track and needs radical reform to deliver a good quality of life for its citizens" = 1L,
    "The UK is broadly speaking headed in the right direction and should continue as normal to deliver a good quality of life for its citizens" = 0L
  ),
  Q82 = c(
    "It is acceptable for unelected bodies and not elected representatives to make the final decisions about things that affect the UK" = 1L,
    "Ultimately final decisions about things that affect the UK should always lie with an elected representative" = 0L
  ),
  Q83 = c(
    "Politicians have allowed unelected bodies and officials to make too many decisions to shield them from accountability" = 1L,
    "Politicians have taken too many decisions and should trust experts and unelected officials to do what is best for the country" = 0L
  )
)

# We know rad_cols[1] corresponds to Q72, rad_cols[2]→Q73, …, rad_cols[8]→Q83.
# So we can assign, in order, the radical_map_list elements to those eight full strings:
qs <- c("Q72","Q73","Q74","Q75","Q78","Q79","Q82","Q83")
radical_map_list2 <- setNames(radical_map_list[qs], rad_cols)

# Now 'radical_map_list2' has names exactly = the eight full column names. For example:
names(radical_map_list2)

# 3c) Recode each of the eight raw columns into {0,1} according to radical_map_list,
#     then standardize (z‐score) each resulting numeric column.

for (col in rad_cols) {
  recode_vec <- radical_map_list2[[col]]
  data <- data %>%
    mutate(
      !!paste0("num_", col) := recode_map(.data[[col]], recode_vec)
    )
}

num_rad_cols <- paste0("num_", rad_cols)

data[num_rad_cols]

# 3d) Compute the z‐score for each “num_<rad_col>”
data <- data %>%
  mutate(across(
    all_of(num_rad_cols),
    ~ z_score(.),
    .names = "z_{.col}"
  ))

# Verify they now exist:
z_rad_cols <- paste0("z_num_", rad_cols)

data[z_rad_cols]

# 3e) Assemble into long format, analogous to opt_long_df and env_long_df

rad_long_df <- data %>%
  select(obs_id, all_of(z_rad_cols)) %>%
  pivot_longer(
    -obs_id,
    names_to  = "item_rad",
    values_to = "y_rad"
  ) %>%
  filter(!is.na(y_rad)) %>%
  mutate(
    # j_rad runs from 1..8, corresponding to the eight items
    j_rad = as.integer(factor(item_rad, levels = z_rad_cols))
  ) %>%
  select(obs_id, j_rad, y_rad)

# At this point:
# - `rad_long_df$i` will be created later (after dropping NAs) in the same way as opt_long_df2 and env_long_df2.
# - j_rad ∈ {1,…,8} and y_rad is the standardized response.


# Compute lower and upper bounds for each of the environmentalism items:
lower_rad <- sapply(z_rad_cols, function(col) {
  min(data[[col]], na.rm = TRUE)
})
upper_rad <- sapply(z_rad_cols, function(col) {
  max(data[[col]], na.rm = TRUE)
})

# ── 4. Material‐Insecurity (Mᵢ) ────────────────────────────────────────────────────────

# Q90–Q97 frequency items
insec_cols <- grep(
  "^For each of the following, say how they apply or do not apply to you:",
  names(data),
  value = TRUE
)

insec_map <- c(
  "Never"        = 0L,
  "Occasionally" = 1L,
  "Fairly often" = 2L,
  "Very often"   = 3L
)

# 1) Recode each to 0–3
data <- data %>%
  mutate(across(
    all_of(insec_cols),
    ~ recode_map(., insec_map),
    .names = "num_{.col}"
  )) %>%
  # 2) Build composite (row mean) and then standardize
  mutate(
    M_full     = rowMeans(select(., paste0("num_", insec_cols)), na.rm = TRUE),
    insecurity = z_score(M_full)
  )

# ── 4. Demographics & Covariates (Xᵢ) ────────────────────────────────────────────────

#  4a) Combine “vote tomorrow” / “vote recent” into one
vote_tomorrow <- grep(
  "Which party would you vote for if there were a general election tomorrow",
  names(data), value = TRUE, ignore.case = TRUE
)
vote_recent <- grep(
  "Which party did you vote for in the most recent general election",
  names(data), value = TRUE, ignore.case = TRUE
)

data <- data %>%
  mutate(
    vote_tomorrow = na_if(.data[[vote_tomorrow]], ""),
    vote_recent   = na_if(.data[[vote_recent]], ""),
    vote          = coalesce(vote_tomorrow, vote_recent)
  )

#  4b) Define party categories (keep “Another party,” “Don’t know,” “Won’t vote” separate)
party_levels <- c(
  "Green", "Labour", "Plaid Cymru",
  "Scottish National Party (SNP)",
  "Liberal Democrat", "Conservative", "Reform UK",
  "Another party", "Dont know", "Won't vote"
)
data <- data %>%
  mutate(
    party = factor(vote, levels = party_levels),
    party_id = as.integer(party) 
    # this maps each level above to 1..10
  )

#  4c) Region → integer 1..10
region_col <- grep("^Region$", names(data), value = TRUE, ignore.case = TRUE)
# Make sure the factor levels of region exactly match your 10 region names
region_levels <- c(
  "Yorkshire & the Humber", "West Midlands", "Scotland", "Wales",
  "North West", "Eastern", "South West", "East Midlands", "London", "South East"
)
data <- data %>%
  mutate(
    region = factor(.data[[region_col]], levels = region_levels),
    region_id = as.integer(region) 
    # now region_id ∈ 1..10
  )

#  4d) Gender → binary dummy (0 = Male, 1 = Female)
gender_col <- grep("^Gender$", names(data), value = TRUE, ignore.case = TRUE)
data <- data %>%
  mutate(
    gender = case_when(
      .data[[gender_col]] == "Male"   ~ 0L,
      .data[[gender_col]] == "Female" ~ 1L,
      TRUE                            ~ NA_integer_
    )
  )

#  4e) Age bracket (6 levels → 5 dummies, reference = “18–24”)
age_bracket_cols_both <- grep("^Age.*\\d+$", names(data), value = TRUE)
age_bracket_col <- age_bracket_cols_both[2]
# Re‐label bracket strings if needed to match exactly; assume they are:
# "65+" (code 1), "55-64" (2), "45-54" (3), "25-34" (4), "35-44" (5), "18-24" (6)
# We create five dummy columns, leaving “18–24” as all zeros.
data <- data %>%
  mutate(
    age_raw = factor(.data[[age_bracket_col]],
                     levels = c("65+", "55-64", "45-54", "25-34", "35-44", "18-24"))
  ) %>%
  mutate(
    age_65plus  = as.integer(age_raw == "65+"),
    age_55_64   = as.integer(age_raw == "55-64"),
    age_45_54   = as.integer(age_raw == "45-54"),
    age_25_34   = as.integer(age_raw == "25-34"),
    age_35_44   = as.integer(age_raw == "35-44")
    # if age_raw == "18-24", all five = 0
  )

#  4f) Education level (7 categories → 6 dummies, reference = “No qualifications”)
edu_col <- grep("educational.*qualification", names(data), value = TRUE, ignore.case = TRUE)
# Assume straightforward matching of the seven strings exactly as in your data
edu_levels <- c(
  "Level 4 qualifications or above, for example a bachelor’s degree or above",
  "Level 2 qualifications, for example 5 or more GCSE passes (formerly O levels)",
  "Level 1 and entry level qualifications, for example 1 to 4 GCSE passes (formerly O levels)",
  "Level 3 qualifications, for example 2 or more A levels",
  "Apprenticeship",
  "No qualifications",
  "Other qualifications"
)
data <- data %>%
  mutate(
    edu_raw = factor(.data[[edu_col]], levels = edu_levels)
  ) %>%
  mutate(
    edu_L4plus  = as.integer(edu_raw == edu_levels[1]),
    edu_L2      = as.integer(edu_raw == edu_levels[2]),
    edu_L1      = as.integer(edu_raw == edu_levels[3]),
    edu_L3      = as.integer(edu_raw == edu_levels[4]),
    edu_appr    = as.integer(edu_raw == edu_levels[5]),
    edu_other   = as.integer(edu_raw == edu_levels[7])
    # If edu_raw == edu_levels[6] (“No qualifications”) → all six = 0
  )

#  4g) Insecurity is already computed & standardized (in earlier step)
#      Check if any NAs remain in insecurity
sum(is.na(data$insecurity)) # is 0

# ── 5. Assemble Covariate DataFrame (cov_df) ─────────────────────────────────────────

cov_df <- data %>%
  select(
    obs_id,
    party_id,
    region_id,
    gender,
    age_65plus, age_55_64, age_45_54, age_35_44, age_25_34,
    edu_L4plus, edu_L3, edu_L2, edu_L1, edu_appr, edu_other,
    insecurity
  )

# Drop rows where any covariate is NA
cov_df <- cov_df %>% drop_na()

# After dropping NAs, we'll need to re‐align opt_long_df and env_long_df to the remaining obs_id

# ── 6. Re‐index respondents and merge with long responses ──────────────────────────────

# 6a) Determine which obs_id remain
valid_ids <- cov_df$obs_id

# 6b) Re‐label them to 1..N
cov_df <- cov_df %>%
  arrange(obs_id) %>%
  mutate(i = row_number())  # new index i ∈ 1..N

# 6c) Build a lookup table: old obs_id → new i
id_map <- tibble(old = cov_df$obs_id, new = cov_df$i)

# 6d) Filter and remap opt_long_df
opt_long_df2 <- opt_long_df %>%
  filter(obs_id %in% valid_ids) %>%
  left_join(id_map, by = c("obs_id" = "old")) %>%
  rename(i = new) %>%
  select(i, j_opt, y_opt)

# 6e) Filter and remap env_long_df
env_long_df2 <- env_long_df %>%
  filter(obs_id %in% valid_ids) %>%
  left_join(id_map, by = c("obs_id" = "old")) %>%
  rename(i = new) %>%
  select(i, j_env, y_env)

# ── 6f) Filter and remap rad_long_df (after you have `id_map` defined) ────────────────
rad_long_df2 <- rad_long_df %>%
  filter(obs_id %in% valid_ids) %>%
  left_join(id_map, by = c("obs_id" = "old")) %>%
  rename(i = new) %>%
  select(i, j_rad, y_rad)

# Now j_rad ∈ 1:8 and i ∈ 1:N_final

# Check how many subjects remain
N_final <- nrow(cov_df)
N_final # is 2933

# ── 7. Extract region_id and party_id vectors (length = N_final) ─────────────────────────

# cov_df is already ordered by obs_id ascending, with new index i = 1..N_final
resp_region <- cov_df$region_id     # length N_final - is also 2933
resp_party  <- cov_df$party_id      # length N_final - is also 2933

# ── 8. Build design matrix X (N_final × 13) in the correct column order ────────────────

# Columns of cov_df in order:
#  gender,
#  age_65plus, age_55_64, age_45_54, age_35_44, age_25_34,
#  edu_L4plus, edu_L3, edu_L2, edu_L1, edu_appr, edu_other,
#  insecurity

Xmat <- cov_df %>%
  select(
    gender,
    age_65plus, age_55_64, age_45_54, age_35_44, age_25_34,
    edu_L4plus, edu_L3, edu_L2, edu_L1, edu_appr, edu_other,
    insecurity
  ) %>%
  as.matrix()

P_cov <- ncol(Xmat)  # should be 13 (confirmed)
P_cov

# ── 9. Final sanity checks ─────────────────────────────────────────────────────────────

# Check that j_opt ranges 1..6 and j_env ranges 1..5
range(opt_long_df2$j_opt) # confirmed: 1-6
range(env_long_df2$j_env) # confirmed: 1-5

# Check that i indices in opt_long_df2 and env_long_df2 both range 1..N_final
range(opt_long_df2$i) 
range(env_long_df2$i) 
range(env_long_df2$i) == c(1, N_final) # confirmed
range(env_long_df2$i) == c(1, N_final) # confirmed

# ── 10. Prepare the stan_data list ─────────────────────────────────────────────────────

stan_data <- list(
  N         = N_final,
  P         = P_cov,
  X         = Xmat,
  R         = length(unique(resp_region)),
  region_id = resp_region,
  Q         = max(resp_party),
  party_id  = resp_party,
  
  J_opt     = length(unique(opt_long_df2$j_opt)),
  N_opt     = nrow(opt_long_df2),
  i_opt     = opt_long_df2$i,
  j_opt     = opt_long_df2$j_opt,
  y_opt     = opt_long_df2$y_opt,
  lower_opt = as.vector(lower_opt), 
  upper_opt = as.vector(upper_opt),  
  
  J_env     = length(unique(env_long_df2$j_env)),
  N_env     = nrow(env_long_df2),
  i_env     = env_long_df2$i,
  j_env     = env_long_df2$j_env,
  y_env     = env_long_df2$y_env,
  lower_env = as.vector(lower_env), 
  upper_env = as.vector(upper_env), 
  
  J_rad     = length(unique(rad_long_df2$j_rad)),
  N_rad     = nrow(rad_long_df2),
  i_rad     = rad_long_df2$i,
  j_rad     = rad_long_df2$j_rad,
  y_rad     = rad_long_df2$y_rad,
  lower_rad = as.vector(lower_rad),   
  upper_rad = as.vector(upper_rad)
)

# save data for Stan
saveRDS(stan_data, "data/stan_data_full.rds")
