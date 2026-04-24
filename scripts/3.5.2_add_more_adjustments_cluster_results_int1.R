# Bootstrap results (from cluster)
# Created by: Yingyan Wu
# use hotdeck imputed datasets

#---- Package loading + options ----
rm(list = ls())
if (!require("pacman")) {
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}
p_load("here", "tidyverse", "broom", "magrittr")

# No scientific notation
options(scipen = 999)

# Paths
source(here::here("scripts", "0.paths.R"))

#---- Load the results ----
results_folder <- paste0(path_to_box, 
                         "Asian_Americans_dementia/Manuscripts/APOE4_ADRD/aa_apoe_adrd/", 
                         "output/R1_model_results/")
files <- list.files(path = results_folder, 
                    pattern = "^bootstraps_scenario_\\d+\\.csv$", full.names = TRUE)
results <- map_dfr(files, read_csv)

#---- PE, 95% CI ----
race_labels <- c("Overall", "Chinese", "Japanese", "Filipino", "Non-Latino White", "Asian American")
fu_yrs <- 0:18
RD_RR_time_wide <- results  %>%
  # pivot longer for calculation at each follow up years across model predicted risks
  pivot_longer(cols = !c(i_boot, scenario_num, race, apoe_prop),
               names_to = c(".value", "model", "year"),
               names_pattern = "(.*)_model_(\\d)_y(\\d+)") %>%
  mutate(model = case_when(model == 1 ~ 1,
                           model == 2 ~ 4, # following the order of main analysis models
                           model == 3 ~ 5)) %>%
  mutate(
    # Calculating excess risk
    # [P(Y1 = 1) - P(Y0 = 1)]/P(Y1 = 1)
    excsfrac = (mean_cif1 - mean_cif0)/mean_cif1,
    excsfrac_perc = excsfrac * 100,
    # Calculating PAR
    # [P(E=1)*(RR-1)/1 + P(E=1)* (RR-1)
    par = (apoe_prop * (RR - 1)) / (1 + apoe_prop * (RR - 1)),
    par_perc = par * 100) %>%
  pivot_wider(id_cols = c(i_boot, scenario_num, race, apoe_prop),
              names_from = c(model, year),
              values_from = !c(i_boot, scenario_num, race, apoe_prop, model, year),
              names_glue = "{.value}_model_{model}_y{year}")

# # Sanity check white
# RD_RR_time_wide %>%
#   filter(race == "wht", i_boot == 0) %>%
#   select(race, contains("y10")) %>%
#   select(contains("RR")) %>%
#   t()

temp1 <- RD_RR_time_wide %>% 
  filter(i_boot != 0) %>%
  select(-contains("int_RD"), -c(i_boot, scenario_num, apoe_prop)) %>%
  group_by(race) %>%
  reframe(across(everything(), ~quantile(.x, c(0.025, 0.975)))) %>%
  mutate(estimate = rep(c("p2.5th", "p97.5th"), 5)) %>%
  ungroup() %>%
  select(estimate, race, everything()) %>% 
  pivot_longer(cols = -c(estimate, race), 
               names_to = c(".value", "fu_yr"),
               names_pattern = "(.*)_y(.*)") %>%
  mutate(fu_yr = as.numeric(fu_yr)) %>%
  arrange(race, fu_yr, estimate)

temp2 <- RD_RR_time_wide %>% 
  filter(i_boot == 0) %>%
  select(-contains("int_RD"), -c(i_boot, scenario_num, apoe_prop)) %>%
  mutate(estimate = "pe") %>%
  select(estimate, race, everything()) %>%
  pivot_longer(cols = -c(estimate, race), 
               names_to = c(".value", "fu_yr"),
               names_pattern = "(.*)_y(.*)") %>%
  mutate(fu_yr = as.numeric(fu_yr)) %>%
  arrange(fu_yr, estimate)

race_labels <- c("Chinese", "Japanese", "Filipino", "Asian American", "Non-Latino White")
RD_RR_time_CI <- bind_rows(temp2, temp1) %>%
  mutate(race = case_when(race == "chn" ~ "Chinese",
                          race == "jpn" ~ "Japanese",
                          race == "phl" ~ "Filipino",
                          race == "wht" ~ "Non-Latino White",
                          race == "AA" ~ "Asian American"))

# # Sanity check wht
# RD_RR_time_CI %>%
#   filter(race == "Non-Latino White", fu_yr == 10) %>%
#   select(estimate, contains(c("RR", "RD")))

save(RD_RR_time_CI, file = paste0(results_folder, "/RD_RR_time_CI.RData"))

#---- Estimates ----
race_fac_labels <- c("Chinese", "Japanese", "Filipino", "Asian American", "Non-Latino White")
race_grp_fac <- c("Asian American ethnic groups", "All")
fy = 10
RD_RR_10yr_tib <- RD_RR_time_CI %>%
  mutate(race = factor(race, levels = race_labels)) %>%
  filter(fu_yr == fy) %>%
  pivot_longer(cols = -c(estimate, race, fu_yr),
               names_to = c(".value", "model"),
               names_pattern = "(.*)_model_(.*)") %>%
  select(race, estimate, model, RD_perc, RR, excsfrac_perc, par_perc) %>%
  mutate_at(vars(RD_perc, excsfrac_perc, par_perc), sprintf, fmt = '%#.1f') %>%
  mutate(RR = sprintf(fmt = "%#.2f", RR)) %>%
  pivot_wider(names_from = estimate, values_from = `RD_perc`:`par_perc`) %>%
  mutate(
    !!sym(paste0("RR (95% CI) at year ", fy)) :=
      paste0(`RR_pe`, "\n (", `RR_p2.5th`, ", ", 
             `RR_p97.5th`, ")"),
    !!sym(paste0("RD % (95% CI) at year ", fy, "(percentage)")) := 
      paste0(`RD_perc_pe`, "\n (", `RD_perc_p2.5th`, ", ", 
             `RD_perc_p97.5th`, ")"),
    !!sym(paste0("Excess fraction % (95% CI) at year ", fy)) :=
      paste0(`excsfrac_perc_pe`, "\n (", `excsfrac_perc_p2.5th`, ", ", 
             `excsfrac_perc_p97.5th`, ")"),
    !!sym(paste0("PAR % (95% CI) at year ", fy)) :=
      paste0(`par_perc_pe`, "\n (", `par_perc_p2.5th`, ", ", 
             `par_perc_p97.5th`, ")"),) %>%
  select(race, model, contains("(95% CI)")) %>%
  pivot_wider(names_from = model, values_from = contains("(95% CI)")) %>%
  mutate(race = factor(race, levels = race_fac_labels),
         race_grp = ifelse(race %in% c("Non-Latino White", "Asian American"), "All", 
                           "Asian American ethnic groups") %>%
           factor(., levels = race_grp_fac)) %>%
  select(race_grp, race, everything()) %>% 
  arrange(race_grp, race)

writexl::write_xlsx(RD_RR_10yr_tib,
                    here::here("output", "tables", "R1_RD_RR_int1_10yr_add_covar.xlsx"))
