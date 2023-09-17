#### Setup ####
# Packages
library(tidyverse)
library(broom)
library(lubridate)

# Modeling pkgs
library(survival)
library(rstpm2)
library(cuRe) # LOLE functions


# Read data
source("R/import_data.R")
source("R/helper_functions.R")

#### Unadjusted LOLE/MRL ####

# Subset cases
study_population <- study_population %>%
  filter(case_or_control == "Case")

##### Setup for relative analyses #####
# Read custom ratetable objects
custom_ratetable <- readRDS(file = "data/customratetable_ALL.RDS")

# Get background hazard - using general.haz() function from CuRe
# yearly scale
study_population$background_hazard_yearly <- general.haz(
  time = study_population$c_timetodeath_days,
  rmap = list(
    age = study_population$rt_agedays,
    sex = study_population$rt_sex,
    year = study_population$rt_year
  ),
  data = study_population,
  ratetable = custom_ratetable,
  scale = 365.24
)


##### Relative models #####
# Stratifying by age, sex, index date
rel_mod_all <-
  study_population %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + nsx(d_indexyear, df = 2) + rt_sex,
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 7,
    bhazard = .$background_hazard_yearly
  )

# stratifying also by LVEF
rel_mod_normal_ef <-
  study_population %>%
  filter(normal_ef == ">=50%") %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + rt_sex + nsx(d_indexyear, df = 2),
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 7,
    bhazard = .$background_hazard_yearly
  )

rel_mod_impaired_ef <-
  study_population %>%
  filter(normal_ef == "<50%") %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + rt_sex + nsx(d_indexyear, df = 2),
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 6,
    bhazard = .$background_hazard_yearly
  )

# Prepare newdata grid for predictions
newdata <- expand_grid(
  c_timetodeath_years = 0:300 / 10,
  c_age = c(50, 60, 70, 80),
  rt_sex = factor(c("male", "female")),
  d_indexyear = c(2000, 2004, 2008, 2016)
)
# Predict [relative] survival
rel_data_all <-
  newdata %>%
  mutate(
    Estimate =
      predict(rel_mod_all, newdata = newdata)
  ) %>%
  mutate(normal_ef = "All")

rel_data_all_normal_lvef <-
  newdata %>%
  mutate(
    Estimate =
      predict(rel_mod_normal_ef, newdata = newdata)
  ) %>%
  mutate(normal_ef = "LVEF >=50%")

rel_data_all_impaired_lvef <-
  newdata %>%
  mutate(
    Estimate =
      predict(rel_mod_impaired_ef, newdata = newdata)
  ) %>%
  mutate(normal_ef = "LVEF <50%")

##### Predict LOLE and MRST #####
var_type_ <- "ci" # Toggle CI on and off = c("n", "ci") (CIs takes a lot of time)

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

newdata_lvef <- newdata %>%
  filter(d_indexyear > 1997) # Don't predict before LVEF measurements

lole_estimates_all <- calc.LL(rel_mod_all,
  newdata = newdata,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_ # CI´s takes a lot of time to compute
)

lole_estimates_normal_lvef <- calc.LL(rel_mod_normal_ef,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_ # CI´s takes a lot of time to compute
)
lole_estimates_impaired_lvef <- calc.LL(rel_mod_impaired_ef,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_ # CI´s takes a lot of time to compute
)


MRL_estimates_all <- calc.LL(rel_mod_all,
  newdata = newdata,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_, # CI´s takes a lot of time to compute
  type = "mrl"
)

MRL_estimates_normal_lvef <- calc.LL(rel_mod_normal_ef,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_, # CI´s takes a lot of time to compute
  type = "mrl"
)

MRL_estimates_impaired_lvef <- calc.LL(rel_mod_impaired_ef,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_, # CI´s takes a lot of time to compute
  type = "mrl"
)

# format all relative measures together
all_ <-
  newdata %>% bind_cols(
    lole_bind_rows(lole_estimates_all),
    lole_bind_rows(MRL_estimates_all)
  )

lvef_ <-
  newdata_lvef %>% bind_cols(
    lole_bind_rows(lole_estimates_normal_lvef),
    lole_bind_rows(lole_estimates_impaired_lvef),
    lole_bind_rows(MRL_estimates_normal_lvef),
    lole_bind_rows(MRL_estimates_impaired_lvef)
  )

lole_measures <-
  all_ %>%
  left_join(lvef_) %>%
  mutate(
    Expected_all.Estimate = MRL_estimates_all.Estimate + lole_estimates_all.Estimate,
    Expected_normal_lvef.Estimate = lole_estimates_normal_lvef.Estimate + MRL_estimates_normal_lvef.Estimate,
    Expected_impaired_lvef.Estimate = lole_estimates_impaired_lvef.Estimate + MRL_estimates_impaired_lvef.Estimate
  ) %>%
  pivot_longer(cols = lole_estimates_all.Estimate:Expected_impaired_lvef.Estimate, names_to = "what", values_to = "Value") %>%
  mutate(what = str_remove(what, ".ci")) %>% # Easy fix to handle formatting for CI.
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
lole_measures %>%
  saveRDS(., file = "output/LOLE_MEASURES.rds")
