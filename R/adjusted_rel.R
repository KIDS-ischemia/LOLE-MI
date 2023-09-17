#### Setup ####

# Packages
library(tidyverse)

# Modeling pkgs
library(survival)
library(rstpm2)
library(cuRe) # LOLE functions

# Read data
source("R/import_data.R")
source("R/helper_functions.R")

# Cases only
study_population <- study_population %>%
  filter(c_case == 1)

#### Adjusted LOLE/MRL ####
##### Setup for relative analyses #####
# Read custom ratetable objects
custom_ratetable_adjusted <- readRDS(file = "data/adjusted_ratetable_list.RDS")

# Get hazards from ratetable object
# yearly scale
study_population$background_hazard_yearly_adjusted_all <- general.haz(
  time = study_population$c_timetodeath_days,
  rmap = list(
    age = study_population$rt_agedays,
    sex = study_population$rt_sex,
    year = study_population$rt_year
  ),
  data = study_population,
  ratetable = custom_ratetable_adjusted$time_split_means_all,
  scale = 365.24
)
study_population$background_hazard_yearly_adjusted_NORMAL_EF <- general.haz(
  time = study_population$c_timetodeath_days,
  rmap = list(
    age = study_population$rt_agedays,
    sex = study_population$rt_sex,
    year = study_population$rt_year
  ),
  data = study_population,
  ratetable = custom_ratetable_adjusted$time_split_means_normal_ef,
  scale = 365.24
)
study_population$background_hazard_yearly_adjusted_IMPAIRED_EF <- general.haz(
  time = study_population$c_timetodeath_days,
  rmap = list(
    age = study_population$rt_agedays,
    sex = study_population$rt_sex,
    year = study_population$rt_year
  ),
  data = study_population,
  ratetable = custom_ratetable_adjusted$time_split_means_impaired_ef,
  scale = 365.24
)

##### Relative models #####
# Stratifying by age, sex, index date
rel_mod_all_adjusted <-
  study_population %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + nsx(d_indexyear, df = 2) + rt_sex,
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 6,
    bhazard = .$background_hazard_yearly_adjusted_all
  )

# stratifying also by LVEF
rel_mod_normal_ef_adjusted <-
  study_population %>%
  filter(normal_ef == ">=50%") %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + rt_sex + nsx(d_indexyear, df = 2),
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 4,
    bhazard = .$background_hazard_yearly_adjusted_NORMAL_EF
  )

rel_mod_impaired_ef_adjusted <-
  study_population %>%
  filter(normal_ef == "<50%") %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + rt_sex + nsx(d_indexyear, df = 2),
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 5,
    bhazard = .$background_hazard_yearly_adjusted_IMPAIRED_EF
  )


# Prepare newdata grid for predictions
newdata <- expand_grid(
  c_timetodeath_years = 0:300 / 10,
  c_age = c(50, 65, 80),
  rt_sex = factor(c("male", "female")),
  d_indexyear = c(2000, 2004, 2008, 2016)
)
# Predict [relative] survival
rel_data_all_adjusted <-
  newdata %>%
  mutate(
    Estimate =
      predict(rel_mod_all_adjusted, newdata = newdata)
  ) %>%
  mutate(normal_ef = "All")

rel_data_all_normal_lvef_adjusted <-
  newdata %>%
  mutate(
    Estimate =
      predict(rel_mod_normal_ef_adjusted, newdata = newdata)
  ) %>%
  mutate(normal_ef = "LVEF >=50%")

rel_data_all_impaired_lvef_adjusted <-
  newdata %>%
  mutate(
    Estimate =
      predict(rel_mod_impaired_ef_adjusted, newdata = newdata)
  ) %>%
  mutate(normal_ef = "LVEF <50%")


##### Predict LOLE and MRST #####
var_type_ <- "ci" # Toggle CI on and off = c("n", "ci")

newdata <- expand_grid(
  c_age = c(50, 65, 80),
  rt_sex = factor(c("male", "female")),
  d_indexyear = 1991:2022
) %>%
  mutate(
    rt_agedays = c_age * 365.24,
    rt_year = parse_date_time(paste(d_indexyear, "0101"), "ymd"),
    c_sex = rt_sex
  )

newdata_lvef <- newdata %>% filter(d_indexyear > 1997)


lole_estimates_all_adjusted <- calc.LL(rel_mod_all_adjusted,
  newdata = newdata,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_all,
  var.type = var_type_ # CI´s takes a lot of time to compute
)

lole_estimates_normal_lvef_adjusted <- calc.LL(rel_mod_normal_ef_adjusted,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_normal_ef,
  var.type = var_type_ # CI´s takes a lot of time to compute
)
lole_estimates_impaired_lvef_adjusted <- calc.LL(rel_mod_impaired_ef_adjusted,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_impaired_ef,
  var.type = var_type_ # CI´s takes a lot of time to compute
)


MRL_estimates_all_adjusted <- calc.LL(rel_mod_all_adjusted,
  newdata = newdata,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_all,
  var.type = var_type_, # CI´s takes a lot of time to compute
  type = "mrl"
)

MRL_estimates_normal_lvef_adjusted <- calc.LL(rel_mod_normal_ef_adjusted,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_normal_ef,
  var.type = var_type_, # CI´s takes a lot of time to compute
  type = "mrl"
)

MRL_estimates_impaired_lvef_adjusted <- calc.LL(rel_mod_impaired_ef_adjusted,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_impaired_ef,
  var.type = var_type_, # CI´s takes a lot of time to compute
  type = "mrl"
)

# format all relative measures together
all_ <-
  newdata %>% bind_cols(
    lole_bind_rows(lole_estimates_all_adjusted),
    lole_bind_rows(MRL_estimates_all_adjusted)
  )

lvef_ <-
  newdata_lvef %>% bind_cols(
    lole_bind_rows(lole_estimates_normal_lvef_adjusted),
    lole_bind_rows(lole_estimates_impaired_lvef_adjusted),
    lole_bind_rows(MRL_estimates_normal_lvef_adjusted),
    lole_bind_rows(MRL_estimates_impaired_lvef_adjusted)
  )

lole_measures_adjusted <-
  all_ %>%
  left_join(lvef_) %>%
  mutate(
    Expected_all_adjusted.Estimate = MRL_estimates_all_adjusted.Estimate + lole_estimates_all_adjusted.Estimate,
    Expected_normal_lvef_adjusted.Estimate = lole_estimates_normal_lvef_adjusted.Estimate + MRL_estimates_normal_lvef_adjusted.Estimate,
    Expected_impaired_lvef_adjusted.Estimate = lole_estimates_impaired_lvef_adjusted.Estimate + MRL_estimates_impaired_lvef_adjusted.Estimate
  ) %>%
  pivot_longer(cols = lole_estimates_all_adjusted.Estimate:Expected_impaired_lvef_adjusted.Estimate, names_to = "what", values_to = "Value") %>%
  mutate(what = str_remove(what, ".ci")) %>% # Quick fix to handle formatting for CI.
  separate_wider_delim(what, delim = ".", names = c("type", "valuetype")) %>%
  mutate(
    normal_ef = case_when(
      str_detect(type, "_all") ~ "All",
      str_detect(type, "_normal") ~ "LVEF >=50%",
      str_detect(type, "_impaired") ~ "LVEF <50%"
    ),
    label = case_when(
      str_detect(type, "Expected") ~ "Expected residual lifetime",
      str_detect(type, "MRL") ~ "Mean residual lifetime",
      str_detect(type, "lole") ~ "Loss of lifetime",
    )
  )

# relative survival tables export
lole_measures_adjusted %>%
  saveRDS(., file = "output/LOLE_MEASURES_adjusted.rds")
