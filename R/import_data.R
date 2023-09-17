# Read data
study_population <- arrow::read_parquet("data/lole_mi_data.parquet")

# tidyverse
library(tidyverse)

# Create additional variables and format old ones
study_population <-
  study_population %>%
  # add some time to time var
  # mutate(c_timetodeath_days = ifelse(c_timetodeath_days == 0, 0.0001, c_timetodeath_days)) %>%
  # And add some variables
  mutate(
    c_pci_treatment = case_when(
      d_repertreatment == 2 ~ T,
      is.na(d_repertreatment) ~ NA,
      T ~ F
    ),
    c_coronary_angio = case_when(
      d_repertreatment %in% c(2, 4) ~ T,
      is.na(d_repertreatment) ~ NA,
      T ~ F
    ),
    c_LVEF = case_when(
      d_left_ventricular_function == 1 ~ "Normalt (>=50%)",
      d_left_ventricular_function == 2 ~ "Lätt nedsatt (40-49%)",
      d_left_ventricular_function == 3 ~ "Måttligt nedsatt (30-39%)",
      d_left_ventricular_function == 4 ~ "Kraftigt nedsatt (<30%)",
      T ~ NA_character_
    ) %>% as.factor(),
    d_cabg = as.logical(d_cabg),
    cpr_before_hospital = case_when(
      cpr_before_hospital == 1 ~ T,
      cpr_before_hospital == 0 ~ F,
      T ~ NA
    ),
    c_car_shock_at_arrival = case_when(
      cardiac_shock == 1 ~ T,
      cardiac_shock == 0 ~ F,
      T ~ NA
    ),
    c_length_of_stay = as.numeric(discharge_date - d_indexdatum),
    c_time_period = factor(case_when(
      d_indexyear < 2001 ~ "≤2000",
      d_indexyear < 2007 ~ "2001-2006",
      d_indexyear < 2011 ~ "2007-2010",
      T ~ "2011-2021"
    )),
    c_married = case_when(
      d_civil == "Gift/partner" ~ 1,
      TRUE ~ 0
    ),
    c_case = ifelse(case_or_control == "Case", 1, 0),
    normal_ef = case_when(
      csub_lvef == "Normalt (>=50%)" ~ ">=50%",
      is.na(csub_lvef) ~ NA_character_,
      T ~ "<50%"
    ) %>% factor(),
    age_group = case_when(
      c_age < 60 ~ "<60",
      c_age < 76 ~ "60-75",
      T ~ "76+"
    ),
    subgroup_all = paste0(c_sex, "s, age ", age_group),
    subgroup_with_ef = paste0(c_sex, "s, age ", age_group, ", LVEF ", normal_ef),
    subgroup_with_ef = ifelse(is.na(normal_ef), NA_character_, subgroup_with_ef),
  # Hypertension definition
    TJ_hypertension = case_when(
      kom_hypertoni == 1 ~ 1,
      c_kalciumantagonist_6_months == 1 ~ 1,
      (c_acei_or_at2_6_months + c_diuretics_6_months) > 1 & kom_hsvikt == 0 ~ 1,
      # c_betablocker_6_months == 1 && kom_hsvikt == 0 ~ 1, # No tachy diagnosis available rn
      TRUE ~ 0
    ),
  TJ_disp_ink_fam = 100 * case_when(
    !is.na(dispinkfam) ~ dispink, # to get in SEK * 100
    TRUE ~ dispinkfam04
  ),
    single_antiplatelet_discharge = ifelse(other_antiplatelet_discharge + aspirin_discharge == 1, 1, 0),
    dual_antiplatelet_discharge = ifelse(ifelse(other_antiplatelet_discharge > 0, 1, 0) + aspirin_discharge == 2, 1, 0),
    major_bleeding = ifelse(d_resuscitated_cardiac_arrest > 0, 1, 0),
    beta_blockers_discharge = ifelse(beta_blockers_discharge == 9, NA, beta_blockers_discharge),
    ace_inhibitors_discharge = ifelse(ace_inhibitors_discharge == 9, NA, ace_inhibitors_discharge),
    angiotensin_ii_block_discharge = case_when(
      angiotensin_ii_block_discharge == 9 ~ NA,
      angiotensin_ii_block_discharge == 2 ~ 1,
      T ~ angiotensin_ii_block_discharge
    ),
    ACEi_or_ARB = case_when(
      ace_inhibitors_discharge == 1 | angiotensin_ii_block_discharge %in% 1:2 ~ 1,
      is.na(ace_inhibitors_discharge) | is.na(angiotensin_ii_block_discharge) ~ NA_real_,
      T ~ 0
    ),
    statins_discharge = ifelse(statins_discharge == 9, NA, statins_discharge),
    d_myocardial_reinfarct_hospital = ifelse(d_myocardial_reinfarct_hospital == 9, NA, d_myocardial_reinfarct_hospital),
    d_cardiogenic_shock = ifelse(d_cardiogenic_shock == 9, NA, d_cardiogenic_shock),
    d_resuscitated_cardiac_arrest = case_when(
      d_resuscitated_cardiac_arrest == 9 ~ NA,
      d_resuscitated_cardiac_arrest == 8 ~ 1,
      T ~ d_resuscitated_cardiac_arrest
    ),
    # Vars to work with ratetable
    rt_sex = as.factor(ifelse(c_sex == "Male", "male", "female")),
    rt_agedays = c_age * 365.24,
    rt_year = d_indexdatum
  )

set.seed(999)