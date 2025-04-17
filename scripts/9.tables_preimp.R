# Tables
# Created by: Yingyan Wu
# Nov.16.2023
# Using non-imputed datasets

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "dplyr", "haven", "labelled", 
       "gtsummary")
#No scientific notation
options(scipen = 999)

# Paths
source(here::here("scripts", "0.paths.R"))

#---- Load the data ----
# load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
#             "analysis_data/aa_apoe_tte_ehr.RData"))
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/aa_apoe_tte_selected.RData"))

#---- FU time ----
aa_apoe_tte_selected %>%
  mutate(fu_time = round(main_dem_v1_fu_time)) %>%
  dplyr::count(fu_time) %>%
  mutate(prop = prop.table(n)*100)

#---- Table 1 ----
mean(aa_apoe_tte_selected$survey_age) # 70.72584
sd(aa_apoe_tte_selected$survey_age) # 7.364147
mean(aa_apoe_tte_selected$main_dem_v1_fu_time) # 10.7
sd(aa_apoe_tte_selected$main_dem_v1_fu_time) # 3.85

data_for_table_1 <- aa_apoe_tte_selected %>%
  mutate(income_pp = income_pp/10000) %>%
  select(survey_age, female, ethnicity_rev, edu_3, usaborn_rev, 
         marital_2, income_pp, generalhealth_3, ehr_ht_median, 
         sr_bmi, sr_depress, prev_htn_flag, prev_diab_flag, 
         prev_combined_stroke_flag,
         main_dem_v1_end_type, main_dem_v1_fu_time, 
         apoe_y) %>%
  set_variable_labels(apoe_y = "APOE e4 carrier",
                      survey_age = "Survey age, years",
                      female = "Female",
                      edu_3 = "Education",
                      usaborn_rev = "US born",
                      marital_2 = "Married or living as married",
                      income_pp = "Household-adjusted income, 10,000 dollars",
                      generalhealth_3 = "General health (self-reported)",
                      ehr_ht_median = "EHR height, in",
                      sr_bmi = "Self-reported BMI",
                      sr_depress = "Self-reported depression",
                      prev_htn_flag = "Hypertension diagnosis pre survey",
                      prev_diab_flag = "Diabetes diagnosis pre survey",
                      prev_combined_stroke_flag = "Stroke diagnosis pre survey",
                      main_dem_v1_end_type = "End of follow up event",
                      # sr_dem_kin = "Family member conditions: Dementia/Alzheimer's Disease",
                      main_dem_v1_fu_time = "Follow up time",
                      ethnicity_rev = "Race/ethnicity") %>%
  set_value_labels(apoe_y = c("APOE-e4+" = 1, "APOE-e4-" = 0),
                   edu_3 = c("Less than high school" = 1, 
                             "High school and some college" = 2,
                             "College and above" = 3),
                   marital_2 = c("Yes" = 1, "No" = 0),
                   generalhealth_3 = c("Excellent/Very Good" = 3,
                                       "Good" = 2,
                                       "Fair/Poor" = 1),
                   main_dem_v1_end_type = 
                     c("Dementia" = "DEMENTIA",
                       "Death" = "DEATH",
                       "Administratively Censored" = "ADMIN CENSORED",
                       "End of Membership" = "END OF MEMBERSHIP"),
                   ethnicity_rev = c(
                     "South Asian" = 1,
                     "Chinese" = 2,
                     "Japanese" = 3,
                     "Korean" = 4,
                     "Filipino" = 5,
                     "Vietnamese" = 6,
                     "Other Southeast Asian" = 7,
                     "Any Pacific Islander" = 8,
                     "White" = 9,
                     "Multiple Asian ethn" = 10))
                   # sr_dem_kin = c("Yes" = 1, "No" = 0))

table_1_res <- data_for_table_1 %>%
  filter(ethnicity_rev %in% c(2, 3, 5, 9)) %>%
  modify_if(is.labelled, to_factor) %>%
  tbl_strata(
    strata = ethnicity_rev,
    .tbl_fun =
      ~ .x %>%
      tbl_summary(type = all_continuous() ~ "continuous",
                  digits = list(all_continuous() ~ 1,
                                all_categorical() ~ c(0, 1)),
                  statistic = list(all_categorical() ~ "{n} ({p})",
                                   all_continuous() ~ "{mean} ({sd})"),
                  by = apoe_y,
                  missing = "ifany",
                  missing_text = "Missing") %>%
      # add_overall %>%
      modify_header(label = "") %>%
      modify_spanning_header(starts_with("stat_") ~ "**APOE**")) %>%
  modify_header(all_stat_cols() ~ 
                  "**{level}**, \nN = {n} ({style_percent(p, digits = 1)}%)") %>%
  bold_labels()

table_1_res %>% as_flex_table()

table_1_res %>% 
  as_hux_xlsx(here::here("output", "tables", "table1_e2e4exclu.xlsx"))

#---- Table 2. APOE * race ----
apoe_fortable <- aa_apoe_tte_selected %>%
  filter(ethnicity_rev != 1) %>%
  mutate(`Race/Ethnicity` = case_when(ethnicity_rev == 2 ~ "Chinese",
                                      ethnicity_rev == 3 ~ "Japanese",
                                      ethnicity_rev == 5 ~ "Filipino",
                                      ethnicity_rev == 9 ~ "White"),
         apoe = case_when(apoe == "e2e2" ~ "e2/e2",
                          apoe == "e3e2" ~ "e2/e3",
                          apoe == "e4e2" ~ "e2/e4",
                          apoe == "e3e3" ~ "e3/e3",
                          apoe == "e3e4" ~ "e3/e4",
                          apoe == "e4e4" ~ "e4/e4"))

apoe_fortable %<>%
  dplyr::select(`Race/Ethnicity`, apoe) %>%
  rbind(apoe_fortable %>%
          dplyr::mutate(`Race/Ethnicity` = case_when(
            ethnicity_rev %in% c(2, 3, 5) ~ "Asian")) %>%
          filter(!is.na(`Race/Ethnicity`)) %>%
          select(`Race/Ethnicity`, apoe)) %>%
  mutate(`Race/Ethnicity` = factor(`Race/Ethnicity`, 
                                   levels = c("White", "Asian", 
                                              "Chinese", "Japanese", "Filipino")))

apoe_table <- apoe_fortable %>%
  tbl_summary(by = `Race/Ethnicity`,
              statistic = list(all_categorical() ~ "{n} ({p})"),
              digits = list(all_categorical() ~ c(0, 1)))

apoe_table

apoe_table %>%
  as_hux_xlsx(here::here("output", "tables", "table2_apoe.xlsx"))


