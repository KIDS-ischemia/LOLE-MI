#### Setup ####

# Packages
library(tidyverse)
library(naniar)
library(gtsummary)
library(gt)

# Read data
source("R/import_data.R")

#### 1. Baseline table ####
### Table one
study_population <-
  study_population %>%
  # format variables for table
  mutate(
    Age = c_age,
    Sex = c_sex,
    `Time period` = c_time_period,
    `Civil status` = case_when(
      d_civil == "Gift/partner" ~ "Married/reg. Partner",
      d_civil == "Ogift/ej partner" ~ "Unmarried",
      d_civil == "Frånskild" ~ "Divorced",
      d_civil == "Änka/änkling/efterlevande partner" ~ "Married/reg. Partner",
      T ~ NA
    ) %>%
      factor(levels = c("Married/reg. Partner", "Unmarried", "Divorced")),
    `Born abroad` = c_utrikesfodd,
    `Disposable income (per family, SEK/year)` = TJ_disp_ink_fam,
    `Disposable income strata` = case_when(c_disp_inkomst_kat == "Låg" ~ "Low",
                                            c_disp_inkomst_kat == "Mellan" ~ "Middle",
                                            c_disp_inkomst_kat == "Hög" ~ "High",
                                            T ~ NA_character_) %>% factor(levels = c("Low", "Middle", "High")),
    `Education level` = case_when(
      sun2000niva_n == "Förgymnasial utbildning kortare än 9 år" ~ "<9 years school",
      sun2000niva_n == "Förgymnasial utbildning 9 år" ~ "9 years school",
      sun2000niva_n == "Gymnasial utbildning högst 2-årig" ~ "2-year secondary school",
      sun2000niva_n == "Gymnasial utbildning 3 år" ~ "3 years secondary school",
      sun2000niva_n == "Eftergymnasial utbildning kortare än 3 år" ~ "<3 years university",
      sun2000niva_n == "Eftergymnasial utbildning 3 år eller längre(exkl. forskarutbildning)" ~ "≥3 years university",
      sun2000niva_n == "Forskarutbildning" ~ "Doctoral studies"
    ) %>%
      factor(levels = c("<9 years school", "9 years school", "2-year secondary school", "3 years secondary school", "<3 years university", "≥3 years university", "Doctoral studies")),
    `Previous stroke` = kom_stroke,
    `Kidney disease` = kom_njursvikt,
    COPD = kom_kol,
    Dementia = kom_demens,
    `Heart failure` = kom_hsvikt,
    Diabetes = kom_diabetes,
    `Peripheral vascular disease` = kom_perifer,
    Cancer = kom_cancer,
    Dialysis = kom_dialys,
    Hypertension = kom_hypertoni,
    `p2y12-inhibitor` = ifelse(c_missing_lakemedel == 1, NA, c_p2y12inh_6_months),
    c_asa_6_months = ifelse(c_missing_lakemedel == 1, NA, c_asa_6_months),
    `Betablocker ` = ifelse(c_missing_lakemedel == 1, NA, c_betablocker_6_months),
    `ACEi/ARB ` = ifelse(c_missing_lakemedel == 1, NA, c_acei_or_at2_6_months),
    `ARNi` = ifelse(c_missing_lakemedel == 1, NA, c_arni_6_months),
    `Lipid lowering drug ` = ifelse(c_missing_lakemedel == 1, NA, c_lipid_lowering_6_months),
    `Calcium antagonist ` = ifelse(c_missing_lakemedel == 1, NA, c_kalciumantagonist_6_months),
    `Diuretics ` = ifelse(c_missing_lakemedel == 1, NA, c_diuretics_6_months),
    tab1_strata = ifelse(c_case == 1, "Case", "Comparator")
  )
# Build table
table1 <-
  study_population %>%
  tbl_summary(
    include = c(
      "Age", "Sex", "Time period",
      "Civil status", "Born abroad", "Disposable income (per family, SEK/year)","Disposable income strata", "Education level",
      "Previous stroke", "Kidney disease", "COPD", "Dementia", "Heart failure", "Diabetes",
      "Peripheral vascular disease", "Cancer", "Dialysis", "Hypertension",
      "p2y12-inhibitor", "Betablocker ", "ACEi/ARB ",
      "ARNi", "Lipid lowering drug ", "Calcium antagonist ",
      "Diuretics ", "tab1_strata"
    ),
    by = tab1_strata,
    missing = "no" # don't list missing data separately
  ) %>%
  modify_header(label = " ") %>%
  bold_labels()

# Polish and save
table1 %>%
  as_gt() %>%
  tab_row_group(label = "", 1:9) %>%
  tab_row_group(label = "Socioeconomy", rows = 10:27) %>%
  tab_row_group(
    label = "Comorbidity",
    rows = 28:37
  ) %>%
  tab_row_group(
    label = "Medication dispensed within 6 months before admission",
    rows = 38:44
  ) %>%
  row_group_order(c("", "Socioeconomy", "Comorbidity", "Medication dispensed within 6 months before admission")) %>%
  gtsave("output/table_1.docx")

#### 2. additional tables and figures #####


# in-hospital characteristics over time
study_population <-
  study_population %>%
  # Group on 2-year intervals to make a little bit more readable
  mutate(
    year_group = paste0(floor(d_indexyear / 2) * 2, "-", floor(d_indexyear / 2) * 2 + 1),
    year_group = ifelse(d_indexyear == 2022, "2022", year_group),
    # Format variables for table
    STEMI = ifelse(csub_STEMI == "STEMI", 1, 0),
    `Coronary angiography` = d_angio,
    `PCI treatment` = d_pci,
    CABG = d_cabg,
    `Reinfarction during hospital stay` = d_myocardial_reinfarct_hospital,
    `Cardiogenic shock` = d_cardiogenic_shock,
    `Resuscitated cardiac arrest` = d_resuscitated_cardiac_arrest,
    `CPR before hospital` = cpr_before_hospital,
    `Cardiogenic shock on arrival` = c_car_shock_at_arrival,
    `Major Bleeding` = major_bleeding,
    `LVEF after MI` = case_when(
      csub_lvef == "Kraftigt nedsatt (<30%)" ~ "<30%",
      csub_lvef == "Måttligt nedsatt (30-39%)" ~ "30-39%",
      csub_lvef == "Lätt nedsatt (40-49%)" ~ "40-49%",
      csub_lvef == "Normalt (>=50%)" ~ "≥50%",
    ),
    `Single antiplatelet therapy` = single_antiplatelet_discharge,
    `Dual antiplatelet therapy` = dual_antiplatelet_discharge,
    `Beta blocker` = beta_blockers_discharge,
    `ACEi or ARB` = ACEi_or_ARB,
    `Statins` = statins_discharge
  )

tableS1 <-
  study_population %>%
  filter(case_or_control == "Case") %>%
  tbl_summary(
    include = c(
      "STEMI", "Coronary angiography",
      "PCI treatment", "CABG", "Reinfarction during hospital stay", "Cardiogenic shock", "Resuscitated cardiac arrest",
      "CPR before hospital", "Cardiogenic shock on arrival", "Major Bleeding", "LVEF after MI", "Single antiplatelet therapy",
      "Dual antiplatelet therapy", "Beta blocker", "ACEi or ARB", "Statins"
    ),
    by = year_group,
    missing = "no"
  ) 

# Format to tibble and polish (to easier fit on pdf page later)
table_S1_tibble <- 
tableS1 %>% 
  as_tibble() %>%
  add_row(`**Characteristic**` = "Medications at discharge", .after = 15)
table_S1_tibble[12:15, 2:5] <- "NA"
table_S1_tibble[c(11, 16), 2:8] <- ""
names(table_S1_tibble) <- str_remove_all(names(table_S1_tibble), coll("*"))
S1_varnames <- str_split_fixed(names(table_S1_tibble), ", N = ", n = 2)
names(table_S1_tibble) <- S1_varnames[,1]
N <- S1_varnames %>% t() %>% as.data.frame() %>% filter(row_number() == 2)
names(N) <- S1_varnames[,1]
N[1,1] <- "N ="
table_S1_tibble <- bind_rows(N, table_S1_tibble)
names(table_S1_tibble)[1]  <- " "
saveRDS(table_S1_tibble, "output/tableS1.RDS")



# Missing data
# Variables for both cases and comparators
vars_all <- c(
  "case_or_control",
  "Age", "Sex", 
  "Civil status", "Born abroad", "Disposable income (SEK/year)", "Education level",
  "Previous stroke", "Kidney disease", "COPD", "Dementia", "Heart failure", "Diabetes",
  "Peripheral vascular disease", "Cancer", "Dialysis", "Hypertension",
  "p2y12-inhibitor", "Betablocker ", "ACEi/ARB ",
  "ARNi", "Lipid lowering drug ", "Calcium antagonist ",
  "Diuretics ")
missingness_table_all <-
  study_population %>%
  select(all_of(vars_all)) %>%
  group_by(case_or_control) %>%
  miss_var_summary() %>%
  mutate(pct_miss = ifelse(pct_miss == 100, NA_real_, pct_miss))
saveRDS(missingness_table_all, file = "output/missingness_table_all.RDS")

# In-hospital variables (cases only)
vars_inhospital <-  c("STEMI", "Coronary angiography",
  "PCI treatment", "CABG", "Reinfarction during hospital stay", "Cardiogenic shock", "Resuscitated cardiac arrest",
  "CPR before hospital", "Cardiogenic shock on arrival", "Major Bleeding", "LVEF after MI", "Single antiplatelet therapy",
  "Dual antiplatelet therapy", "Beta blocker", "ACEi or ARB", "Statins"
)
missingness_table_inhospital <-
  study_population %>%
  filter(c_case == 1) %>%
  select(all_of(vars_inhospital)) %>%
  miss_var_summary() %>%
  mutate(pct_miss = ifelse(pct_miss == 100, NA_real_, pct_miss))
saveRDS(missingness_table_inhospital, file = "output/missingness_table_inhospital.RDS")

# year distribution
year_dist <-
  study_population %>%
  group_by(d_indexyear, case_or_control) %>%
  tally()
saveRDS(year_dist, file = "output/year_dist.RDS")

# Follow-up times
study_population %>% 
  group_by(c_case) %>% 
  summarise(sum_personyears = sum(c_timetodeath_years), median = median(c_timetodeath_years), 
            IQR_low = quantile(c_timetodeath_years, prob = .25), IQR_high = quantile(c_timetodeath_years, prob = .75)) %>%
  saveRDS(file = "output/followup_times.RDS")

# Number of events
study_population %>% 
  group_by(c_case) %>% 
  summarise(n_deaths = sum(c_dead), pct_dead = sum(c_dead) / n() * 100) %>%
  saveRDS(file = "output/n_deaths.RDS")
