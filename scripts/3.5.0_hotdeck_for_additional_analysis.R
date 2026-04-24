# Hotdeck imputation for R1

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "tidyverse", "magrittr", "VIM")
# Do not load "summarytools" as it contains conflicted functions as dplyr

#No scientific notation
options(scipen = 999)

# Paths
source(here::here("scripts", "0.paths.R"))

#---- Load the data ----
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/aa_apoe_tte_selected_e4all.RData"))

#---- check the missingness ----
aa_apoe_tte_selected %>%
  select(edu_3, usaborn_rev, income_pp, prev_htn_flag, sr_depress, prev_diab_flag) %>%
  DataExplorer::plot_missing()

#---- check age dist ----
aa_apoe_tte_selected %>% 
  ggplot(aes(x = survey_age_r)) + 
  geom_density()

aa_apoe_tte_selected %>% reframe(median(survey_age_r)) 

aa_apoe_tte_selected <- aa_apoe_tte_selected %>% 
  mutate(age_le69 = ifelse(survey_age_r <= 69, 1, 0))

#---- Hot deck imputation for covariates ----
# for missingness in education, nativity, household income, hypertension, depression, and diabetes
# Hotdeck on race/ethnicity, sex, age group(above and below overall median)
set.seed(62283)
impute_data <- 
  hotdeck(aa_apoe_tte_selected, 
          variable = c("edu_3", "usaborn_rev", "income_pp"), 
          ord_var = c("ethnicity_rev", "female", "age_le69"))

#---- join in long data ----
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                   "analysis_data/aa_apoe_long_tte_selected_e4all.RData"))

aa_apoe_long_tte_selected <- aa_apoe_long_tte_selected %>% 
  select(-c(edu_3, usaborn_rev, income_pp)) %>% 
  left_join(impute_data %>% select(subjid, edu_3, usaborn_rev, income_pp,
                                   prev_htn_flag, prev_diab_flag), by = "subjid") 
save(aa_apoe_long_tte_selected, 
     file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                   "analysis_data/aa_apoe_long_tte_selected_e4all_hotdeck.RData"))

# Save imputed data
aa_apoe_tte_selected <- impute_data
save(aa_apoe_tte_selected, 
     file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                   "analysis_data/aa_apoe_tte_selected_e4all_hotdeck.RData"))

