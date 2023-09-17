#### Create a custom ratetable ####
# First part inspired by https://enochytchen.com/tutorials/ratetable/content/

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
source("data_import.R")

# Get population data from HMD
ratetable_sweden <- relsurv::transrate.hmd("data/mltper_1x1.txt", "data/fltper_1x1.txt")

# select Comparators only
comparators <- study_population %>%
  filter(c_case == 0)

# Approximate the date of birth,
# Will take care of derived variables (entry/exit years) later
comparators <- comparators %>%
  mutate(
    c_sex = as_factor(c_sex), ## as_factor preserves labels
    c_dead = as.numeric(c_dead),
    c_dob = d_indexdatum - c_age * 365.241,
    c_year = d_indexyear,
    c_exit = d_indexdatum + c_timetodeath_days
  ) %>%
  select(idnr, c_sex, c_dead, c_dob, c_year, c_age, c_exit, d_indexdatum)
str(comparators)
summary(comparators)

# Split calendar time into 1-year intervals;
comp_split <- survSplit(Surv(decimal_date(d_indexdatum), decimal_date(c_exit), c_dead, type = "mstate") ~ .,
  data = comparators, cut = seq(1991, 2022, by = 1),
  event = "c_dead", episode = "period"
) %>%
  # changed word "censor" to 0, so to keep it consistent with original definition
  mutate(c_dead = as.numeric(ifelse(as.character(c_dead) == "censor",
    "0", as.character(c_dead)
  )))

# Inspect: select the first 20 to take a look
head(comp_split, 20)

# For downstream analysis, we want age as primary time scale;
# Calculate age at entry and age at exit
# Also, it would be nice to have the period expressed as actual starting year
# of the time interval; see ?survSplit for its definition
#' (i.e. 1 = before first interval)
comp_split <-
  mutate(
    comp_split,
    age_start = tstart - decimal_date(c_dob),
    age_stop = tstop - decimal_date(c_dob),
    period_ = period,
    period = min(c_year) + period - 2
  )

####  Fit flexible model  ####

# Flexible parametric model
fpm <-
  comp_split %>%
  stpm2(Surv(time = age_start, time2 = age_stop, event = c_dead == 1) ~ c_sex + nsx(period, df = 2),
    data = .,
    tvc = list(c_sex = 2, period = 2),
    df = 3
  )
summary(fpm)

# # Save model as it takes so long to compute
saveRDS(fpm, "data/comp_split_model_ALL.rds", compress = F)

# # And read it if necessary
# fpm <- readRDS("data/comp_split_model.rds")


####  Predict from flexible model  ####
# Create a grid to predict everything on
comp_new <- expand_grid(
  c_sex = levels(comp_split$c_sex) %>% factor(),
  age_stop = 18:100,
  period = min(comp_split$period):max(comp_split$period)
)

# Populate the empty data frame with predicted hazards
comp_new$hazard <- predict(fpm, newdata = comp_new, type = "hazard")

comp_new <- comp_new %>%
  mutate(
    prob = exp(-hazard),
    age = age_stop,
    sex = c_sex
  ) %>%
  select(sex, age, period, prob, hazard)
str(comp_new)
comp_new

# add official swedish population values
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
comp_new <- comp_new %>%
  left_join(general_sweden, by = c("age" = "age", "sex" = "sex", "period" = "year")) %>%
  mutate(ratio = hazard / value_general_pop) %>%
  filter(age <= 100)

# We create a set with all rows from general_sweden, to make it easier to predict on later
replaced_hazards <-
  general_sweden %>%
  left_join(comp_new, by = c("age" = "age", "sex" = "sex", "year" = "period")) %>%
  filter(year >= 1980 & age <= 100) %>% # We do not need the full lenght
  select(-prob, -value_general_pop.y) %>% # drop columns
  rename(
    hazard_general_pop = value_general_pop.x,
    hazard_comparator_pop = hazard
  ) %>%
  mutate(
    male = ifelse(sex == "Male", 1, 0), # binary for the regression model
  )


####  Fit glm  ####
# fit a regression model
lin_reg <-
  replaced_hazards %>%
  lm(
    ratio ~ nsx(age, df = 2, derivs = c(1, 1)) +
      nsx(year, df = 2, derivs = c(1, 1)) + nsx(hazard_general_pop, df = 2, derivs = c(1, 1)) +
      nsx(male, df = 1),
    data = .
  )
summary(lin_reg)

# Add predictions from glm to data
replaced_hazards <-
  replaced_hazards %>%
  mutate(
    smoothed_ratios = predict(lin_reg, replaced_hazards),
    smoothed_hazards = smoothed_ratios * hazard_general_pop
  )


####  create ratetable and final hazard dataframe  ####
# make final hazard dataframe
mortalityrates_custom <-
  replaced_hazards %>%
  select(age, year, sex, smoothed_hazards) %>%
  rename(yearly_hazards = smoothed_hazards) %>%
  mutate(
    daily_hazards = yearly_hazards / 365.241,
    prob = exp(-yearly_hazards)
  )

# we'll use the general populations because it is (maybe) easier to format as a ratetable
to_rt <-
  general_sweden %>%
  left_join(mortalityrates_custom) %>%
  mutate(
    value = ifelse(is.na(yearly_hazards), value_general_pop, yearly_hazards),
    prob = exp(-value)
  )

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

# A ratetable containing the DAILY HAZARDS - as the HMD also does
custom_ratetable <- transrate(men = male, women = female, yearlim = c(1751, 2021))
is.ratetable(custom_ratetable) # TRUE

# Save for later
saveRDS(custom_ratetable, file = "data/customratetable_ALL.RDS")
