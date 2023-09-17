#### Plots ####

# Libraries
library(tidyverse)
library(patchwork)
library(extrafont)

# Set ggplot theme
theme_set(crgg::theme_ki_standard())
# load fonts
loadfonts()

# Mortality Rate curves
MortalityRate_curves_lvef <-
  readRDS(file = "output/mortalityrate_table.rds") %>%
  filter(LVEF != "All") %>%
  mutate(case = ifelse(c_case == 1, "Cases", "Comparators") %>% factor()) %>%
  ggplot(aes(x = c_timetodeath_years, y = Estimate * 1000, linetype = case, color = Sex, fill = Sex, ymin = lower * 1000, ymax = upper * 1000)) +
  geom_line() +
  geom_ribbon(alpha = .2, linewidth = .1) +
  facet_grid(cols = vars(LVEF), rows = vars(`Age group`)) +
  labs(
    caption = "Mortality rate per 1000 person-years",
    fill = "", linetype = "", color = ""
  ) +
  ylab("Mortality rate per 1000 person-years") +
  xlab("Time since MI or index date in years") +
  coord_cartesian(ylim = c(0, 250)) +
  crgg::scale_colour_ki_d() +
  crgg::scale_fill_ki_d()

# HR curves
HR_curves_lvef <-
  readRDS(file = "output/HR_table.rds") %>%
  filter(LVEF != "All") %>%
  ggplot(aes(x = c_timetodeath_years, y = Estimate, color = Sex, fill = Sex, ymin = lower, ymax = upper)) +
  geom_line() +
  geom_ribbon(alpha = .2, linewidth = .1) +
  facet_grid(cols = vars(LVEF), rows = vars(`Age group`)) +
  labs(
    caption = "Hazard Ratio over time",
    fill = "", linetype = "", color = ""
  ) +
  ylab("Hazard Ratio cases vs comparators") +
  xlab("Time since MI or index date in years") +
  coord_cartesian(ylim = c(0, 7), xlim = c(0, 15)) +
  crgg::scale_colour_ki_d() +
  crgg::scale_fill_ki_d()

# Combine and save
MR_HR_plot <- MortalityRate_curves_lvef + HR_curves_lvef + plot_annotation(tag_levels = 'A')
ggsave("output/figure_1.pdf", plot = MR_HR_plot, width = 12, height = 8)
# Embed font
embed_fonts("output/figure_1.pdf", outfile = "output/figure_1.pdf")


# Survival curves
survival_curves <-
  readRDS(file = "output/survival_curves.RDS") %>%
  mutate(LVEF = ifelse(LVEF == "All", "all", LVEF)) %>% # Capital "A" issue in pdf render...
  mutate(case = ifelse(c_case == 1, "Cases", "Comparators") %>% factor()) %>%
  ggplot(aes(x = c_timetodeath_years, y = Estimate, linetype = case, color = Sex, 
             fill = Sex, ymin = lower, ymax = upper)) +
  geom_line() +
  geom_ribbon(alpha = .2, linewidth = .1) +
  facet_grid(cols = vars(LVEF), rows = vars(`Age group`)) +
  labs(
    caption = "Survival predicted from models fit on subgroups (stratified on age, sex, all/LVEF)",
    fill = "", linetype = "", color = ""
  ) +
  ylab("Survival") +
  xlab("Time since MI or index date in years") +
  crgg::scale_colour_ki_d() +
  crgg::scale_fill_ki_d()
ggsave("output/figure_2.pdf", plot = survival_curves, width = 12, height = 8)
# Embed font
embed_fonts("output/figure_2.pdf", outfile = "output/figure_2.pdf")

# LOLE
LOLE_over_time <-
  readRDS("output/LOLE_MEASURES.rds") %>%
  mutate(normal_ef = case_when(normal_ef == "All" ~ "all", 
                               normal_ef == "Impaired EF" ~ "LVEF <50%",
                               normal_ef == "Normal EF" ~ "LVEF >=50%",
                               T ~ normal_ef
                               )) %>% 
  pivot_wider(values_from = Value, names_from = valuetype) %>%
  filter(str_detect(type, "lole")) %>%
  ggplot(aes(x = d_indexyear, y = Estimate, color = rt_sex, fill = rt_sex, ymin = lower, ymax = upper)) +
  geom_line() +
  geom_ribbon(alpha = .1, linewidth = 0) +
  facet_grid(rows = vars(normal_ef), cols = vars(c_age)) +
  labs(linetype = "", color = "") +
  crgg::theme_ki_standard() +
  crgg::scale_colour_ki_d() +
  crgg::scale_fill_ki_d() +
  labs(fill = "") +
  ylab("Loss of life expectancy (years)") +
  xlab("Calendar year") +
  coord_cartesian(ylim = c(0, 15))

LOLE_over_time_adjusted <-
  readRDS("output/LOLE_MEASURES_adjusted.rds")  %>%
  mutate(normal_ef = case_when(normal_ef == "All" ~ "all", 
                               normal_ef == "Impaired EF" ~ "LVEF <50%",
                               normal_ef == "Normal EF" ~ "LVEF >=50%",
                               T ~ normal_ef
  )) %>% 
  filter(str_detect(type, "lole")) %>%
  pivot_wider(values_from = Value, names_from = valuetype) %>%
  ggplot(aes(x = d_indexyear, y = Estimate, ymin = lower, ymax = upper, color = rt_sex, fill = rt_sex)) +
  geom_line() +
  geom_ribbon(alpha = .1, linewidth = 0) +
  facet_grid(rows = vars(normal_ef), cols = vars(c_age)) +
  labs(linetype = "", color = "") +
  crgg::theme_ki_standard() +
  crgg::scale_colour_ki_d() +
  crgg::scale_fill_ki_d() +
  labs(fill = "") +
  ylab("Loss of life expectancy (years)") +
  xlab("Calendar year") +
  coord_cartesian(ylim = c(0, 15))

LOLE_plots <- LOLE_over_time + LOLE_over_time_adjusted + plot_annotation(tag_levels = 'A')
ggsave("output/figure_3.pdf", plot = LOLE_plots, width = 12, height = 8)
# Embed font
embed_fonts("output/figure_3.pdf", outfile = "output/figure_3.pdf")


# MRL
MRL_over_time <-
  readRDS("output/LOLE_MEASURES.rds") %>%
  mutate(normal_ef = case_when(normal_ef == "All" ~ "all", 
                               normal_ef == "Impaired EF" ~ "LVEF <50%",
                               normal_ef == "Normal EF" ~ "LVEF >=50%",
                               T ~ normal_ef
  )) %>% 
  pivot_wider(values_from = Value, names_from = valuetype) %>%
  filter(str_detect(type, "Expected") | str_detect(type, "MRL")) %>%
  mutate(label = fct_relevel(label, "Mean residual lifetime")) %>%
  ggplot(aes(x = d_indexyear, y = Estimate, color = rt_sex, fill = rt_sex, ymin = lower, ymax = upper, linetype = label)) +
  geom_line() +
  geom_ribbon(alpha = .1, linewidth = 0) +
  facet_grid(rows = vars(normal_ef), cols = vars(c_age)) +
  labs(linetype = "", color = "") +
  crgg::theme_ki_standard() +
  crgg::scale_colour_ki_d() +
  crgg::scale_fill_ki_d() +
  labs(fill = "") +
  ylab("Residual lifetime (years)") +
  xlab("Calendar year") +
  coord_cartesian(ylim = c(0, 38))

MRL_over_time_adjusted <-
  readRDS("output/LOLE_MEASURES_adjusted.rds")  %>%
  mutate(normal_ef = case_when(normal_ef == "All" ~ "all", 
                               normal_ef == "Impaired EF" ~ "LVEF <50%",
                               normal_ef == "Normal EF" ~ "LVEF >=50%",
                               T ~ normal_ef
  )) %>% 
  filter(str_detect(type, "Expected") | str_detect(type, "MRL")) %>%
  pivot_wider(values_from = Value, names_from = valuetype) %>%
  mutate(label = fct_relevel(label, "Mean residual lifetime")) %>%
  ggplot(aes(x = d_indexyear, y = Estimate, ymin = lower, ymax = upper, color = rt_sex, fill = rt_sex, linetype = label)) +
  geom_line() +
  geom_ribbon(alpha = .2, linewidth = .2) +
  facet_grid(rows = vars(normal_ef), cols = vars(c_age)) +
  labs(linetype = "", color = "") +
  crgg::theme_ki_standard() +
  crgg::scale_colour_ki_d() +
  crgg::scale_fill_ki_d() +
  labs(fill = "") +
  ylab("Residual lifetime (years)") +
  xlab("Calendar year") +
  coord_cartesian(ylim = c(0, 38))

MRL_plots <- MRL_over_time + MRL_over_time_adjusted + plot_annotation(tag_levels = 'A')
ggsave("output/figure_4.pdf", plot = MRL_plots, width = 12, height = 8)
# Embed font
embed_fonts("output/figure_4.pdf", outfile = "output/figure_4.pdf")
