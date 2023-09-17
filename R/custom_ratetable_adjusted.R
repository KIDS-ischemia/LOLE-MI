#### create a custom ratetable where we adjust for comorbidity #####
# Code adapted from https://enochytchen.com/tutorials/ratetable/content/

# Packages
library(tidyverse)
library(broom)
library(patchwork)
library(lubridate)

# Modeling pkgs
library(survival)
library(rstpm2)
library(relsurv)

# Read data
source("R/import_data.R")
print("Data import OK")

# Get population data from HMD (popmort equivalent)
ratetable_sweden <- relsurv::transrate.hmd(
  "data/mltper_1x1.txt",
  "data/fltper_1x1.txt"
)

# Create dummy variables for sun2000 and disp_ink
study_population_dummy <-
  study_population %>%
  mutate(
    sun2000niva_n = case_when(is.na(sun2000niva_n) ~ "Saknas", TRUE ~ sun2000niva_n),
    c_disp_inkomst_kat = case_when(is.na(c_disp_inkomst_kat) ~ "Saknas", TRUE ~ c_disp_inkomst_kat)
  ) %>%
  # Create dummies for factor vars
  fastDummies::dummy_columns(select_columns = c("sun2000niva_n", "c_disp_inkomst_kat")) %>%
  # format time vars
  #' Approximate the date of birth,
  mutate(
    c_sex = as_factor(c_sex), ## as_factor preserves labels
    c_dead = as.numeric(c_dead),
    c_dob = d_indexdatum - c_age * 365.241,
    c_year = d_indexyear,
    c_exit = d_indexdatum + c_timetodeath_days
  ) %>%
  # Select variables for modeling
  select(
    c_case,
    csub_lvef,
    c_married,
    c_age,
    c_sex,
    c_utrikesfodd,
    kom_stroke,
    kom_njursvikt,
    kom_kol,
    kom_demens,
    kom_hsvikt,
    kom_diabetes,
    kom_perifer,
    kom_cancer,
    kom_dialys,
    TJ_hypertension,
    c_p2y12inh_6_months,
    c_asa_6_months,
    c_betablocker_6_months,
    c_acei_or_at2_6_months,
    c_arni_6_months,
    c_diuretics_6_months,
    c_kalciumantagonist_6_months,
    c_lipid_lowering_6_months,
    c_timetodeath_years,
    c_dead,
    c_missing_lakemedel,
    c_dob,
    c_year,
    c_exit,
    d_indexdatum,
    idnr,
    `sun2000niva_n_Eftergymnasial utbildning 3 år eller längre(exkl. forskarutbildning)`,
    `sun2000niva_n_Eftergymnasial utbildning kortare än 3 år`,
    `sun2000niva_n_Förgymnasial utbildning 9 år`,
    `sun2000niva_n_Förgymnasial utbildning kortare än 9 år`,
    sun2000niva_n_Forskarutbildning,
    `sun2000niva_n_Gymnasial utbildning 3 år`,
    `sun2000niva_n_Gymnasial utbildning högst 2-årig`,
    sun2000niva_n_Saknas,
    c_disp_inkomst_kat_Hög,
    c_disp_inkomst_kat_Låg,
    c_disp_inkomst_kat_Mellan,
    c_disp_inkomst_kat_Saknas,
  )

# select Comparators only
comparators_dummy <- study_population_dummy %>%
  filter(c_case == 0)

#' Split calendar time into 1-year intervals;
comp_dummy_split <- survSplit(Surv(decimal_date(d_indexdatum), decimal_date(c_exit), c_dead, type = "mstate") ~ .,
  data = comparators_dummy, cut = seq(1998, 2022, by = 1),
  event = "c_dead", episode = "period"
) %>%
  # changed word "censor" to 0, so to keep it consistent with original definition
  mutate(c_dead = as.numeric(ifelse(as.character(c_dead) == "censor",
    "0", as.character(c_dead)
  )))

#' For downstream analysis, we want age as primary time scale;
#' Calculate age at entry and age at exit
#' Also, it would be nice to have the period expressed as actual starting year
#' of the time interval; see ?survSplit for its definition
#' (i.e. 1 = before first interval)
comp_dummy_split <-
  mutate(
    comp_dummy_split,
    age_start = tstart - decimal_date(c_dob),
    age_stop = tstop - decimal_date(c_dob),
    period = 1998 + period - 2
  )


# Save the data
arrow::write_parquet(comp_dummy_split, "tmp/split_comp_adj.parquet")

#### Flexible parametric model ####
# Needs _lots_ of RAM.

print(paste0("Starting FPM at ", Sys.time()))
fpm_adjusted <-
  comp_dummy_split %>%
  stpm2(
    Surv(time = age_start, time2 = age_stop, event = c_dead == 1) ~ c_sex + nsx(period, df = 2) +
      c_married + c_utrikesfodd + kom_stroke + kom_njursvikt + kom_kol + kom_demens + kom_hsvikt + kom_diabetes + kom_perifer +
      kom_cancer + kom_dialys + TJ_hypertension + c_p2y12inh_6_months + c_asa_6_months + c_betablocker_6_months +
      c_acei_or_at2_6_months + c_arni_6_months + c_diuretics_6_months + c_kalciumantagonist_6_months + c_lipid_lowering_6_months +
      c_missing_lakemedel + `sun2000niva_n_Eftergymnasial utbildning kortare än 3 år` + `sun2000niva_n_Förgymnasial utbildning 9 år` +
      `sun2000niva_n_Förgymnasial utbildning kortare än 9 år` + sun2000niva_n_Forskarutbildning +
      `sun2000niva_n_Gymnasial utbildning 3 år` + `sun2000niva_n_Gymnasial utbildning högst 2-årig` +
      `sun2000niva_n_Saknas` + c_disp_inkomst_kat_Hög + c_disp_inkomst_kat_Låg + c_disp_inkomst_kat_Mellan,
    data = .,
    tvc = list(c_sex = 2, period = 2),
    df = 3
  )
saveRDS(fpm_adjusted, "tmp/mod_stpm2_adjusted_all.rds")

print("FPM finished and saved")

# read model object
fpm_adjusted  <-  readRDS("tmp/mod_stpm2_adjusted_all.rds")

#### Create list to save everything in ####
newdata_versions <- list()

### Make newdata based on case group
cases_dummy <- study_population_dummy %>%
  filter(c_case == 1) %>%
  select(
    c_case,
    csub_lvef,
    c_married,
    c_age,
    c_sex,
    c_utrikesfodd,
    kom_stroke,
    kom_njursvikt,
    kom_kol,
    kom_demens,
    kom_hsvikt,
    kom_diabetes,
    kom_perifer,
    kom_cancer,
    kom_dialys,
    TJ_hypertension,
    c_p2y12inh_6_months,
    c_asa_6_months,
    c_betablocker_6_months,
    c_acei_or_at2_6_months,
    c_arni_6_months,
    c_diuretics_6_months,
    c_kalciumantagonist_6_months,
    c_lipid_lowering_6_months,
    c_timetodeath_years,
    c_dead,
    c_missing_lakemedel,
    c_dob,
    c_year,
    c_exit,
    d_indexdatum,
    idnr,
    `sun2000niva_n_Eftergymnasial utbildning 3 år eller längre(exkl. forskarutbildning)`,
    `sun2000niva_n_Eftergymnasial utbildning kortare än 3 år`,
    `sun2000niva_n_Förgymnasial utbildning 9 år`,
    `sun2000niva_n_Förgymnasial utbildning kortare än 9 år`,
    sun2000niva_n_Forskarutbildning,
    `sun2000niva_n_Gymnasial utbildning 3 år`,
    `sun2000niva_n_Gymnasial utbildning högst 2-årig`,
    sun2000niva_n_Saknas,
    c_disp_inkomst_kat_Hög,
    c_disp_inkomst_kat_Låg,
    c_disp_inkomst_kat_Mellan,
    c_disp_inkomst_kat_Saknas,
  )

#' Split calendar time into 1-year intervals;
case_dummy_split <- survSplit(Surv(decimal_date(d_indexdatum), decimal_date(c_exit), c_dead, type = "mstate") ~ .,
  data = cases_dummy, cut = seq(1998, 2022, by = 1),
  event = "c_dead", episode = "period"
) %>%
  # changed word "censor" to 0, so to keep it consistent with original definition
  mutate(c_dead = as.numeric(ifelse(as.character(c_dead) == "censor",
    "0", as.character(c_dead)
  )))

# Format to match comparators
case_dummy_split <-
  mutate(
    case_dummy_split,
    age_start = tstart - decimal_date(c_dob),
    age_stop = tstop - decimal_date(c_dob),
    period = 1998 + period - 2
  ) %>%
  mutate(
    # floor ages for grouping
    age_stop = floor(age_stop),
    # format EF for grouping
    normal_ef = case_when(
      csub_lvef == "Normalt (>=50%)" ~ "normal",
      csub_lvef %in% c("Måttligt nedsatt (30-39%)", "Lätt nedsatt (40-49%)", "Kraftigt nedsatt (<30%)") ~ "impaired"
    )
  )


# Predict on case population
case_dummy_split$predicted_hazards <- predict(fpm_adjusted, newdata = case_dummy_split, type = "hazard")

# average Estimates
adjusted_hazards <- list()

adjusted_hazards$time_split_means_all <-
  case_dummy_split %>%
  select(c_sex, age_stop, period, predicted_hazards) %>%
  group_by(c_sex, age_stop, period) %>%
  summarise(mean_predicted_hazard = mean(predicted_hazards)) %>%
  ungroup() %>%
  mutate(
    normal_ef = "All",
    prob = exp(-mean_predicted_hazard)
  )

adjusted_hazard_rate_lvef <-
  case_dummy_split %>%
  select(c_sex, age_stop, period, normal_ef, predicted_hazards) %>%
  filter(!is.na(normal_ef)) %>%
  group_by(c_sex, age_stop, period, normal_ef) %>%
  summarise(mean_predicted_hazard = mean(predicted_hazards)) %>%
  ungroup() %>%
  mutate(prob = exp(-mean_predicted_hazard))

adjusted_hazards$time_split_means_normal_ef <-
  adjusted_hazard_rate_lvef %>%
  filter(normal_ef == "normal")

adjusted_hazards$time_split_means_impaired_ef <-
  adjusted_hazard_rate_lvef %>%
  filter(normal_ef == "impaired")

# granular EF
adjusted_hazard_rate_lvef_granular <-
  case_dummy_split %>%
  select(c_sex, age_stop, period, csub_lvef, predicted_hazards) %>%
  filter(!is.na(csub_lvef)) %>%
  group_by(c_sex, age_stop, period, csub_lvef) %>%
  summarise(mean_predicted_hazard = mean(predicted_hazards)) %>%
  ungroup() %>%
  mutate(prob = exp(-mean_predicted_hazard))

adjusted_hazards$time_split_means_granular_normalEF <-
  adjusted_hazard_rate_lvef_granular %>%
  filter(csub_lvef == "Normalt (>=50%)")
adjusted_hazards$time_split_means_granular_40_49EF <-
  adjusted_hazard_rate_lvef_granular %>%
  filter(csub_lvef == "Lätt nedsatt (40-49%)")
adjusted_hazards$time_split_means_granular_30_39EF <-
  adjusted_hazard_rate_lvef_granular %>%
  filter(csub_lvef == "Måttligt nedsatt (30-39%)")
adjusted_hazards$time_split_means_granular_below30EF <-
  adjusted_hazard_rate_lvef_granular %>%
  filter(csub_lvef == "Kraftigt nedsatt (<30%)")

# Label data frames
adjusted_hazards <- lapply(names(adjusted_hazards), function(name) {
  df <- adjusted_hazards[[name]]
  df$method <- name
  df
})

# join for convenience
adjusted_hazards <- bind_rows(adjusted_hazards) %>%
  # and rename
  rename(
    sex = c_sex,
    age = age_stop,
    year = period
  )


# format official swedish population values
general_sweden <- popEpi::ratetable_to_long_df(ratetable_sweden) %>%
  tibble() %>%
  rename(value_general_pop = value) %>%
  mutate(
    age = as.numeric(age),
    year = as.numeric(year),
    sex = ifelse(sex == "male", "Male", "Female"),
    value_general_pop = value_general_pop * 365.241
  ) # Convert from daily to yearly


# join together
adjusted_hazards <- adjusted_hazards %>%
  left_join(general_sweden, by = c("age" = "age", "sex" = "sex", "year" = "year")) %>%
  mutate(ratio = mean_predicted_hazard / value_general_pop) %>%
  # we remove vewry high ages as they are so high we cannot fit a model with them.
  filter(age <= 100)

# format values for modeling
adjusted_hazards <- adjusted_hazards %>% mutate(male = ifelse(sex == "Male", 1, 0))
# Back to list for modeling separately
adjusted_hazards_list <- split(adjusted_hazards, adjusted_hazards$method)

#### Fit a linear model on each subset ####
lin_reg_list <- lapply(adjusted_hazards_list, function(subset) {
  lm(
    ratio ~ nsx(age, df = 2, derivs = c(1, 1)) +
      nsx(year, df = 2, derivs = c(1, 1)) + nsx(value_general_pop, df = 2, derivs = c(1, 1)) +
      nsx(male, df = 1),
    data = subset
  )
})

# use HMD data as a grid to predict on (the size of the resulting table must be consistent to be able to format as ratetable later.)
prediction_grid <- setNames(
  replicate(length(names(lin_reg_list)),
    general_sweden %>% filter(year > 1990) %>% filter(age < 101) %>% # special high values for 110-year olds
      mutate(male = ifelse(sex == "Male", 1, 0)),
    simplify = FALSE
  ),
  names(lin_reg_list)
)

# And predict
# Predict using the models on each dataset and add predictions as a new column
for (i in seq_along(prediction_grid)) {
  dataset <- prediction_grid[[i]] # Get the dataset
  model <- lin_reg_list[[i]] # Get the corresponding linear model
  dataset$smoothed_ratio <- predict(model, newdata = dataset) # Predict and store the estimates in a new column
  prediction_grid[[i]] <- dataset # Update the dataset in the list
}

# back to long form
prediction_grid <- Map(function(df, method) {
  df %>% mutate(method = method)
}, prediction_grid, names(prediction_grid))
prediction_grid <- bind_rows(prediction_grid)

# add FPM data and ratios.
prediction_grid <- prediction_grid %>% left_join(adjusted_hazards)

# calc smoothed hazards
adjusted_hazards <- prediction_grid %>%
  mutate(smoothed_hazards = smoothed_ratio * value_general_pop)

####  create ratetable and final hazard dataset  ####
# make final hazard dataset
mortalityrates_custom_adjusted <-
  adjusted_hazards %>%
  select(age, year, sex, smoothed_hazards, method) %>%
  rename(yearly_hazards = smoothed_hazards) %>%
  mutate(
    daily_hazards = yearly_hazards / 365.241,
    prob = exp(-yearly_hazards)
  )

# to list again for looping over
mortalityrates_custom_adjusted_list <- split(mortalityrates_custom_adjusted, mortalityrates_custom_adjusted$method)

# and format as ratetable
adjusted_ratetable_list <-
  lapply(mortalityrates_custom_adjusted_list, function(subset) {
    to_rt <-
      general_sweden %>%
      left_join(subset) %>%
      mutate(
        value = ifelse(is.na(yearly_hazards), value_general_pop, yearly_hazards), # gives us 110-year olds hazards.
        prob = exp(-value)
      ) %>%
      filter(year > 1990)
    male <- to_rt %>%
      filter(sex == "Male") %>%
      select(age, year, prob) %>%
      pivot_wider(names_from = year, values_from = prob) %>%
      column_to_rownames(var = "age") %>%
      as.matrix()
    female <- to_rt %>%
      filter(sex == "Female") %>%
      select(age, year, prob) %>%
      pivot_wider(names_from = year, values_from = prob) %>%
      column_to_rownames(var = "age") %>%
      as.matrix()
    relsurv::transrate(men = male, women = female, yearlim = c(1991, 2021))
  })

# Rename list objects
names(adjusted_ratetable_list) <- names(mortalityrates_custom_adjusted_list)

# Check if they formatted correctly
lapply(adjusted_ratetable_list, is.ratetable)

# and save
saveRDS(adjusted_ratetable_list, file = "data/adjusted_ratetable_list.RDS")

print("Ratetables saved and OK")
