#### Granular LVEF groups LOLE/MRL ####
##### Setup #####
# Packages
library(tidyverse)
library(broom)

# Modeling pkgs
library(rstpm2)
library(cuRe) # LOLE functions

# Read data
source("R/import_data.R")
source("R/helper_functions.R")

#### Setup for relative analyses ####
# Select cases
study_population <- study_population %>%
  filter(c_case == 1)

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

#### Models ####
rel_mod_gran_normal_ef_unadjusted <-
  study_population %>%
  filter(csub_lvef == "Normalt (>=50%)") %>%
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

rel_mod_gran_40_49_ef_unadjusted <-
  study_population %>%
  filter(csub_lvef == "Lätt nedsatt (40-49%)") %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + nsx(d_indexyear, df = 2) + rt_sex,
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 4,
    bhazard = .$background_hazard_yearly
  )

rel_mod_gran_30_39_ef_unadjusted <-
  study_population %>%
  filter(csub_lvef == "Måttligt nedsatt (30-39%)") %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + nsx(d_indexyear, df = 2) + rt_sex,
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 5,
    bhazard = .$background_hazard_yearly
  )

rel_mod_gran_below30ef_unadjusted <-
  study_population %>%
  filter(csub_lvef == "Kraftigt nedsatt (<30%)") %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + nsx(d_indexyear, df = 2) + rt_sex,
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 4,
    bhazard = .$background_hazard_yearly
  )

print("Unadjusted models OK")


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

# calc.LL
lole_estimates_gran_normalef_unadjusted <- calc.LL(rel_mod_gran_normal_ef_unadjusted,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_ # CI´s takes a lot of time to compute
)

lole_estimates_gran_40_49ef_unadjusted <- calc.LL(rel_mod_gran_40_49_ef_unadjusted,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_ # CI´s takes a lot of time to compute
)

lole_estimates_gran_30_39ef_unadjusted <- calc.LL(rel_mod_gran_30_39_ef_unadjusted,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_ # CI´s takes a lot of time to compute
)
lole_estimates_gran_below30ef_unadjusted <- calc.LL(rel_mod_gran_below30ef_unadjusted,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_ # CI´s takes a lot of time to compute
)



MRL_estimates_gran_normalef_unadjusted <- calc.LL(rel_mod_gran_normal_ef_unadjusted,
  type = "mrl",
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_ # CI´s takes a lot of time to compute
)


MRL_estimates_gran_40_49ef_unadjusted <- calc.LL(rel_mod_gran_40_49_ef_unadjusted,
  type = "mrl",
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_ # CI´s takes a lot of time to compute
)

MRL_estimates_gran_30_39ef_unadjusted <- calc.LL(rel_mod_gran_30_39_ef_unadjusted,
  type = "mrl",
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_ # CI´s takes a lot of time to compute
)
MRL_estimates_gran_below30ef_unadjusted <- calc.LL(rel_mod_gran_below30ef_unadjusted,
  type = "mrl",
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable,
  var.type = var_type_ # CI´s takes a lot of time to compute
)


# format all relative measures together
lole_measures_granular_unadjusted <-
  newdata_lvef %>%
  bind_cols(
    lole_bind_rows(lole_estimates_gran_normalef_unadjusted),
    lole_bind_rows(lole_estimates_gran_40_49ef_unadjusted),
    lole_bind_rows(lole_estimates_gran_30_39ef_unadjusted),
    lole_bind_rows(lole_estimates_gran_below30ef_unadjusted),
    lole_bind_rows(MRL_estimates_gran_normalef_unadjusted),
    lole_bind_rows(MRL_estimates_gran_40_49ef_unadjusted),
    lole_bind_rows(MRL_estimates_gran_30_39ef_unadjusted),
    lole_bind_rows(MRL_estimates_gran_below30ef_unadjusted)
  ) %>%
  mutate(
    Expected_normalef_unadjusted.Estimate = lole_estimates_gran_normalef_unadjusted.Estimate + MRL_estimates_gran_normalef_unadjusted.Estimate,
    Expected_40_49ef_unadjusted.Estimate = lole_estimates_gran_40_49ef_unadjusted.Estimate + MRL_estimates_gran_40_49ef_unadjusted.Estimate,
    Expected_30_39ef_unadjusted.Estimate = lole_estimates_gran_30_39ef_unadjusted.Estimate + MRL_estimates_gran_30_39ef_unadjusted.Estimate,
    Expected_below30ef_unadjusted.Estimate = lole_estimates_gran_below30ef_unadjusted.Estimate + MRL_estimates_gran_below30ef_unadjusted.Estimate,
  ) %>%
  pivot_longer(cols = lole_estimates_gran_normalef_unadjusted.Estimate:Expected_below30ef_unadjusted.Estimate, names_to = "what", values_to = "Value") %>%
  mutate(what = str_remove(what, ".ci")) %>% # Quick fix to handle formatting for CI.
  separate_wider_delim(what, delim = ".", names = c("type", "valuetype")) %>%
  mutate(
    normal_ef = case_when(
      str_detect(type, "_normal") ~ "Normal EF",
      str_detect(type, "_30_39") ~ "30-39%",
      str_detect(type, "_40_49") ~ "40-49%",
      str_detect(type, "below") ~ "Below 30%"
    ) %>% factor(),
    normal_ef = fct_relevel(normal_ef, "Below 30%"),
    label = case_when(
      str_detect(type, "Expected") ~ "Expected residual lifetime",
      str_detect(type, "MRL") ~ "Mean residual lifetime",
      str_detect(type, "lole") ~ "Loss of lifetime",
    )
  )

# relative survival tables export
lole_measures_granular_unadjusted %>%
  saveRDS(., file = "output/LOLE_MEASURES_unadjusted_granular.rds")


#### Adjusted models  ####
# Read custom ratetable objects
custom_ratetable_adjusted <- readRDS(file = "data/adjusted_ratetable_list.RDS")

# Get hazards from ratetable object
# yearly scale
study_population$background_hazard_yearly_adjusted_gran_normalEF <- general.haz(
  time = study_population$c_timetodeath_days,
  rmap = list(
    age = study_population$rt_agedays,
    sex = study_population$rt_sex,
    year = study_population$rt_year
  ),
  data = study_population,
  ratetable = custom_ratetable_adjusted$time_split_means_granular_normalEF,
  scale = 365.24
)
study_population$background_hazard_yearly_adjusted_gran_40_49EF <- general.haz(
  time = study_population$c_timetodeath_days,
  rmap = list(
    age = study_population$rt_agedays,
    sex = study_population$rt_sex,
    year = study_population$rt_year
  ),
  data = study_population,
  ratetable = custom_ratetable_adjusted$time_split_means_granular_40_49EF,
  scale = 365.24
)
study_population$background_hazard_yearly_adjusted_gran_30_39EF <- general.haz(
  time = study_population$c_timetodeath_days,
  rmap = list(
    age = study_population$rt_agedays,
    sex = study_population$rt_sex,
    year = study_population$rt_year
  ),
  data = study_population,
  ratetable = custom_ratetable_adjusted$time_split_means_granular_30_39EF,
  scale = 365.24
)
study_population$background_hazard_yearly_adjusted_gran_below30EF <- general.haz(
  time = study_population$c_timetodeath_days,
  rmap = list(
    age = study_population$rt_agedays,
    sex = study_population$rt_sex,
    year = study_population$rt_year
  ),
  data = study_population,
  ratetable = custom_ratetable_adjusted$time_split_means_granular_below30EF,
  scale = 365.24
)


#### Models####
rel_mod_gran_normal_ef_adjusted <-
  study_population %>%
  filter(csub_lvef == "Normalt (>=50%)") %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + nsx(d_indexyear, df = 2) + rt_sex,
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 4,
    bhazard = .$background_hazard_yearly_adjusted_gran_normalEF
  )

rel_mod_gran_40_49_ef_adjusted <-
  study_population %>%
  filter(csub_lvef == "Lätt nedsatt (40-49%)") %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + nsx(d_indexyear, df = 2) + rt_sex,
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 4,
    bhazard = .$background_hazard_yearly_adjusted_gran_40_49EF
  )


rel_mod_gran_30_39_ef_adjusted <-
  study_population %>%
  filter(csub_lvef == "Måttligt nedsatt (30-39%)") %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + nsx(d_indexyear, df = 2) + rt_sex,
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 4,
    bhazard = .$background_hazard_yearly_adjusted_gran_30_39EF
  )


rel_mod_gran_below30ef_adjusted <-
  study_population %>%
  filter(csub_lvef == "Kraftigt nedsatt (<30%)") %>%
  stpm2(Surv(c_timetodeath_years, c_dead) ~ nsx(c_age, df = 2) + nsx(d_indexyear, df = 2) + rt_sex,
    data = .,
    tvc = list(
      c_age = 2,
      d_indexyear = 2,
      rt_sex = 2
    ),
    df = 4,
    bhazard = .$background_hazard_yearly_adjusted_gran_below30EF
  )


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

lole_estimates_gran_normalef_adjusted <- calc.LL(rel_mod_gran_normal_ef_adjusted,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_granular_normalEF,
  var.type = var_type_ # CI´s takes a lot of time to compute
)


lole_estimates_gran_40_49ef_adjusted <- calc.LL(rel_mod_gran_40_49_ef_adjusted,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_granular_40_49EF,
  var.type = var_type_ # CI´s takes a lot of time to compute
)

lole_estimates_gran_30_39ef_adjusted <- calc.LL(rel_mod_gran_30_39_ef_adjusted,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_granular_30_39EF,
  var.type = var_type_ # CI´s takes a lot of time to compute
)
lole_estimates_gran_below30ef_adjusted <- calc.LL(rel_mod_gran_below30ef_adjusted,
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_granular_below30EF,
  var.type = var_type_ # CI´s takes a lot of time to compute
)


tictoc::tic()
MRL_estimates_gran_normalef_adjusted <- calc.LL(rel_mod_gran_normal_ef_adjusted,
  type = "mrl",
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_granular_normalEF,
  var.type = var_type_ # CI´s takes a lot of time to compute
)
tictoc::toc()

MRL_estimates_gran_40_49ef_adjusted <- calc.LL(rel_mod_gran_40_49_ef_adjusted,
  type = "mrl",
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_granular_40_49EF,
  var.type = var_type_ # CI´s takes a lot of time to compute
)

MRL_estimates_gran_30_39ef_adjusted <- calc.LL(rel_mod_gran_30_39_ef_adjusted,
  type = "mrl",
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_granular_30_39EF,
  var.type = var_type_ # CI´s takes a lot of time to compute
)
MRL_estimates_gran_below30ef_adjusted <- calc.LL(rel_mod_gran_below30ef_adjusted,
  type = "mrl",
  newdata = newdata_lvef,
  time = 0,
  rmap = list(age = rt_agedays, sex = rt_sex, year = rt_year),
  ratetable = custom_ratetable_adjusted$time_split_means_granular_below30EF,
  var.type = var_type_ # CI´s takes a lot of time to compute
)

# format all relative measures together
lole_measures_granular_adjusted <-
  newdata_lvef %>%
  bind_cols(
    lole_bind_rows(lole_estimates_gran_normalef_adjusted),
    lole_bind_rows(lole_estimates_gran_40_49ef_adjusted),
    lole_bind_rows(lole_estimates_gran_30_39ef_adjusted),
    lole_bind_rows(lole_estimates_gran_below30ef_adjusted),
    lole_bind_rows(MRL_estimates_gran_normalef_adjusted),
    lole_bind_rows(MRL_estimates_gran_40_49ef_adjusted),
    lole_bind_rows(MRL_estimates_gran_30_39ef_adjusted),
    lole_bind_rows(MRL_estimates_gran_below30ef_adjusted)
  ) %>%
  mutate(
    Expected_normalef_adjusted.Estimate = lole_estimates_gran_normalef_adjusted.Estimate + MRL_estimates_gran_normalef_adjusted.Estimate,
    Expected_40_49ef_adjusted.Estimate = lole_estimates_gran_40_49ef_adjusted.Estimate + MRL_estimates_gran_40_49ef_adjusted.Estimate,
    Expected_30_39ef_adjusted.Estimate = lole_estimates_gran_30_39ef_adjusted.Estimate + MRL_estimates_gran_30_39ef_adjusted.Estimate,
    Expected_below30ef_adjusted.Estimate = lole_estimates_gran_below30ef_adjusted.Estimate + MRL_estimates_gran_below30ef_adjusted.Estimate,
  ) %>%
  pivot_longer(
    cols = lole_estimates_gran_normalef_adjusted.Estimate:Expected_below30ef_adjusted.Estimate,
    names_to = "what", values_to = "Value"
  ) %>%
  mutate(what = str_remove(what, ".ci")) %>% # Quick fix to handle formatting for CI.
  separate_wider_delim(what, delim = ".", names = c("type", "valuetype")) %>%
  mutate(
    normal_ef = case_when(
      str_detect(type, "_normal") ~ "Normal EF",
      str_detect(type, "_30_39") ~ "30-39%",
      str_detect(type, "_40_49") ~ "40-49%",
      str_detect(type, "below") ~ "Below 30%"
    ) %>% factor(),
    normal_ef = fct_relevel(normal_ef, "Below 30%"),
    label = case_when(
      str_detect(type, "Expected") ~ "Expected residual lifetime",
      str_detect(type, "MRL") ~ "Mean residual lifetime",
      str_detect(type, "lole") ~ "Loss of lifetime",
    )
  )

# relative survival tables export
lole_measures_granular_adjusted %>%
  saveRDS(., file = "output/LOLE_MEASURES_adjusted_granular.rds")
