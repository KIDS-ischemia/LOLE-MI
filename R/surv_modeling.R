#### Setup ####

# Packages
library(tidyverse)
library(broom)
library(ggsurvfit)
library(lubridate)

# Modeling pkgs
library(survival)
library(rstpm2)
library(cuRe) # LOLE functions

# Read data
source("R/import_data.R")

#### Survival modeling ####

# Flexible models - subgroups
list_subgroups_all <- list()
for (subgroup in unique(study_population$subgroup_all)) {
  print(paste0("processing ", subgroup))
  list_subgroups_all[[subgroup]] <- study_population %>%
    filter(subgroup_all == subgroup) %>%
    stpm2(Surv(c_timetodeath_years, c_dead) ~ c_case,
          data = .,
          tvc.formula = ~ c_case:nsx(log(c_timetodeath_years),
                                     knots = log(c(.001, 0.7, 15))
          ),
          df = 5
    )
}

list_subgroups_with_lvef <- list()
for (subgroup in unique(study_population$subgroup_with_ef)) {
  if (is.na(subgroup)) {
    next
  }
  print(paste0("processing ", subgroup))
  list_subgroups_with_lvef[[subgroup]] <- study_population %>%
    filter(subgroup_with_ef == subgroup) %>%
    stpm2(Surv(c_timetodeath_years, c_dead) ~ c_case,
          data = .,
          tvc.formula = ~ c_case:nsx(log(c_timetodeath_years),
                                     knots = log(c(.001, .7, 15))
          ),
          df = 5
    )
}

# Predict survival for the subgroup models (per month, 15 years)
newdata <- expand_grid(
  c_timetodeath_years = 0:(15 * 12) / 12, #
  c_case = 0:1
)

predicted_survival_all_without_ef <- list()
for (subgroup in unique(study_population$subgroup_all)) {
  predicted_survival_all_without_ef[[subgroup]] <-
    newdata %>%
    bind_cols(predict(list_subgroups_all[[subgroup]], newdata = newdata, se.fit = T)) %>%
    mutate(Subgroup = subgroup) %>%
    separate_wider_delim(Subgroup, delim = ", ", names = c("Sex", "Age group"), cols_remove = F) %>%
    mutate(LVEF = "All")
}

predicted_survival_all_with_lvef <- list()
for (subgroup in unique(study_population$subgroup_with_ef)) {
  if (is.na(subgroup)) {
    next
  }
  predicted_survival_all_with_lvef[[subgroup]] <-
    newdata %>%
    bind_cols(predict(list_subgroups_with_lvef[[subgroup]], newdata = newdata, se.fit = T)) %>%
    mutate(Subgroup = subgroup) %>%
    separate_wider_delim(Subgroup, delim = ", ", names = c("Sex", "Age group", "LVEF"), cols_remove = F)
}

# bind_row for later plotting
survival_curves <-
  bind_rows(
    bind_rows(predicted_survival_all_with_lvef),
    bind_rows(predicted_survival_all_without_ef)
  ) %>%
  mutate(case = ifelse(c_case == 1, "Cases", "Comparators") %>% factor())
saveRDS(survival_curves, file = "output/survival_curves.RDS")

# Predict mortality rates over time
newdata <- expand_grid(
  c_timetodeath_years = 0:(15 * 12) / 12, # 15 years at 1 month intervals
  c_case = 0:1
) %>%
  filter(c_timetodeath_years != 0) # Rstpm2 cannot handle rates at time = 0.

predicted_mortality_all_without_ef <- list()
for (subgroup in unique(study_population$subgroup_all)) {
  predicted_mortality_all_without_ef[[subgroup]] <-
    newdata %>%
    bind_cols(predict(list_subgroups_all[[subgroup]], newdata = newdata, se.fit = T, type = "hazard")) %>%
    mutate(Subgroup = subgroup) %>%
    separate_wider_delim(Subgroup, delim = ", ", names = c("Sex", "Age group"), cols_remove = F) %>%
    mutate(LVEF = "All")
}

predicted_mortality_all_with_lvef <- list()
for (subgroup in unique(study_population$subgroup_with_ef)) {
  if (is.na(subgroup)) {
    next
  }
  predicted_mortality_all_with_lvef[[subgroup]] <-
    newdata %>%
    bind_cols(predict(list_subgroups_with_lvef[[subgroup]], newdata = newdata, se.fit = T, type = "hazard")) %>%
    mutate(Subgroup = subgroup) %>%
    separate_wider_delim(Subgroup, delim = ", ", names = c("Sex", "Age group", "LVEF"), cols_remove = F)
}

# Save for mortality rate curves and table
bind_rows(
  bind_rows(predicted_mortality_all_without_ef),
  bind_rows(predicted_mortality_all_with_lvef)
) %>%
  saveRDS(., file = "output/mortalityrate_table.rds")


# Predict hazard ratios  over time
newdata <- expand_grid(
  c_timetodeath_years = 0:(15 * 12) / 12, #
  c_case = 0
) %>%
  filter(c_timetodeath_years != 0)

predicted_HR_all_without_ef <- list()
for (subgroup in unique(study_population$subgroup_all)) {
  predicted_HR_all_without_ef[[subgroup]] <-
    newdata %>%
    bind_cols(predict(list_subgroups_all[[subgroup]],
                      newdata = newdata, se.fit = T, type = "hr",
                      exposed = function(data) transform(data, c_case = 1)
    )) %>%
    mutate(Subgroup = subgroup) %>%
    separate_wider_delim(Subgroup, delim = ", ", names = c("Sex", "Age group"), cols_remove = F) %>%
    mutate(LVEF = "All")
}

predicted_HR_all_with_lvef <- list()
for (subgroup in unique(study_population$subgroup_with_ef)) {
  if (is.na(subgroup)) {
    next
  }
  predicted_HR_all_with_lvef[[subgroup]] <-
    newdata %>%
    bind_cols(predict(list_subgroups_with_lvef[[subgroup]],
                      newdata = newdata, se.fit = T, type = "hr",
                      exposed = function(data) transform(data, c_case = 1)
    )) %>%
    mutate(Subgroup = subgroup) %>%
    separate_wider_delim(Subgroup, delim = ", ", names = c("Sex", "Age group", "LVEF"), cols_remove = F)
}

# Save HR table
bind_rows(
  bind_rows(predicted_HR_all_without_ef),
  bind_rows(predicted_HR_all_with_lvef)
) %>%
  saveRDS(., file = "output/HR_table.rds")


# predict Hazard difference (Excess mortality)
# Predict hazard ratios  over time
newdata <- expand_grid(
  c_timetodeath_years = 0:(15 * 12) / 12, #
  c_case = 0
) %>%
  filter(c_timetodeath_years != 0)

predicted_MD_all_without_ef <- list()
for (subgroup in unique(study_population$subgroup_all)) {
  predicted_MD_all_without_ef[[subgroup]] <-
    newdata %>%
    bind_cols(predict(list_subgroups_all[[subgroup]],
                      newdata = newdata, se.fit = T, type = "hdiff",
                      exposed = function(data) transform(data, c_case = 1)
    )) %>%
    mutate(Subgroup = subgroup) %>%
    separate_wider_delim(Subgroup, delim = ", ", names = c("Sex", "Age group"), cols_remove = F) %>%
    mutate(LVEF = "All")
}


predicted_MD_all_with_lvef <- list()
for (subgroup in unique(study_population$subgroup_with_ef)) {
  if (is.na(subgroup)) {
    next
  }
  predicted_MD_all_with_lvef[[subgroup]] <-
    newdata %>%
    bind_cols(predict(list_subgroups_with_lvef[[subgroup]],
                      newdata = newdata, se.fit = T, type = "hdiff",
                      exposed = function(data) transform(data, c_case = 1)
    )) %>%
    mutate(Subgroup = subgroup) %>%
    separate_wider_delim(Subgroup, delim = ", ", names = c("Sex", "Age group", "LVEF"), cols_remove = F)
}

# bind and save table
bind_rows(
  bind_rows(predicted_MD_all_without_ef),
  bind_rows(predicted_MD_all_with_lvef)
) %>%
  saveRDS(., file = "output/MD_table.rds")
