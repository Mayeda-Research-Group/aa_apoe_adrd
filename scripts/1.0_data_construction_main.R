# APOE data construction
# Created by: Yingyan Wu
# Jun.20.2023

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "dplyr", "haven", "labelled", 
       "gtsummary")
# Do not load "summarytools" as it contains conflicted functions as dplyr
#No scientific notation
options(scipen = 999)

# Paths
source(here::here("scripts", "0.paths.R"))

#---- Paths ----
path_to_cleaned_data<- paste0(path_to_box,
                              "Asian_Americans_dementia_data/analysis_data_tables/")
path_to_raw_data <- paste0(path_to_box,
                           "Asian_Americans_dementia_data/raw_data_tables/")

#---- Load data ----
aa_adrd_tte <- read_sas(paste0(path_to_cleaned_data, 
                               "aa_adrd_cardiometabolic_tte.sas7bdat"))
apoe_data <- read_csv(paste0(path_to_raw_data,
                             "apoe_data/apoe_10049.csv"))
# Merge the data
aa_apoe_tte <- apoe_data %>%
  select(-starts_with("chr")) %>% # Move the local ancestry attempted estimates
  left_join(aa_adrd_tte, by = c("iid" = "SUBJID", "female" = "FEMALE")) %>%
  dplyr::rename(SUBJID = iid) %>%
  filter(!is.na(STUDY)) # Restrict to those who filled out CMHS or RPGEH
# 52190-51679 = 511
colnames(aa_apoe_tte) <- tolower(colnames(aa_apoe_tte))

#---- Check the data ----
with(aa_apoe_tte, table(sex, female, useNA = "ifany"))
with(apoe_data, table(sex, female, useNA = "ifany"))
with(aa_apoe_tte, table(apoe, useNA = "ifany"))
with(aa_apoe_tte, summarytools::freq(ethnicity_rev, useNA = "ifany"))

#---- APOE ----
aa_apoe_tte %<>%
  mutate(apoe_y = case_when(str_detect(apoe, "e3e4|e4e4") ~ 1,
                            TRUE ~ 0))
with(aa_apoe_tte, table(apoe, apoe_y, useNA = "ifany"))
summarytools::freq(aa_apoe_tte$apoe_y)

#---- Education ----
aa_apoe_tte %<>%
  mutate(edu_3 = case_when(education_rev %in% c(1, 2) ~ 1,
                           education_rev %in% c(3, 4) ~ 2,
                           education_rev %in% c(5, 6) ~ 3))

#---- Marital status ----
# "Never Married" = 1, "Married or living as married" = 2,
# "Separted/Divorced" = 3, "Widowed" = 4
aa_apoe_tte %<>%
  mutate(marital_2 = case_when(maritalstatus == 2 ~ 1,
                               maritalstatus %in% c(1, 3, 4) ~ 0))

#---- General Health ----
# "Excellent" = 5, "Very Good" = 4, "Good" = 3, "Fair" = 2, "Poor" = 1
aa_apoe_tte %<>%
  mutate(generalhealth_3 = case_when(generalhealth %in% c(1, 2) ~ 1,
                                     generalhealth == 3 ~ 2,
                                     generalhealth %in% c(4, 5) ~ 3))

#---- Smoking status ----
# Ever/never smoker
aa_apoe_tte %<>%
  mutate(smoking_2 = case_when(smoking_status %in% c(2, 3) ~ 1,
                               smoking_status == 1 ~ 0))

#---- Alcohol ----
# 2010 Dietary Guidelines for Americans:
# Moderate alcohol consumption: <= 1 per day for women, <= 2 per day for men
# Heavy/high-risk drinking: >= 3 any day/ >= 7 per week for women
# >= 4 any day / >= 14 per week for men
aa_apoe_tte %<>%
  mutate(alcohol_3 = case_when(
    alcohol_drinksperweek_rev == 0 & alcohol_binge == 1 ~ 0,
    ((alcohol_drinksperweek_rev >= 7) & female == 1)|
      ((alcohol_drinksperweek_rev >= 14) & female == 0)|alcohol_binge >= 2 ~ 2,
    !is.na(alcohol_drinksperweek_rev) | alcohol_binge == 1 | 
      !is.na(alcohol_daysperweek) ~ 1,
    TRUE ~ NA))
# # Sanity check
# summarytools::freq(aa_apoe_tte$alcohol_3)
# with(aa_apoe_tte, table(alcohol_3, alcohol_binge, useNA ="ifany"))
# with(aa_apoe_tte %>% filter(is.na(alcohol_3)), 
#      table(alcohol_drinksperweek_rev, female, useNA = "ifany"))
# with(aa_apoe_tte %>% filter(is.na(alcohol_3)), 
#      table(alcohol_daysperweek, useNA = "ifany"))
# with(aa_apoe_tte, table(alcohol_3, alcohol_drinksperweek_rev, useNA ="ifany"))

#---- Physical activity ----
# Can't really use AHA guideline thresholds as the second smallest value is 74.25
female_median <- with(aa_apoe_tte %>% filter(female == 1 & totalmet != 0), 
                      quantile(totalmet, 0.5, na.rm = T))
male_median <- with(aa_apoe_tte %>% filter(female == 0 & totalmet != 0), 
                    quantile(totalmet, 0.5, na.rm = T))

aa_apoe_tte %<>%
  mutate(pa_3 = case_when(totalmet == 0 ~ 0,
                          (female == 1 & totalmet <= female_median) |
                            (female == 0 & totalmet <= male_median) ~ 1,
                          is.na(totalmet) ~ NA_real_,
                          TRUE ~ 2))
# # Sanity check
# summarytools::freq(aa_apoe_tte$pa_3)
# with(aa_apoe_tte, table(pa_3, totalmet, useNA = "ifany"))

#---- Asian ----
aa_apoe_tte %<>%
  mutate(ethn_asian = case_when(ethnicity_rev == 9 ~ 0,
                                !is.na(ethnicity_rev) ~ 1))
# Sanity check
with(aa_apoe_tte, table(ethn_asian, ethnicity_rev, useNA = "ifany"))

#---- Censor 90 + ----
# We are using the imputed 90 + age and change the censored 90+ type to admin censored
aa_apoe_tte %<>%
  mutate(main_dem_v1_end_type =
           case_when(main_dem_v1_end_type == "CENSORED 90+" ~ "ADMIN CENSORED",
                     TRUE ~ main_dem_v1_end_type))

# Sanity check
with(aa_apoe_tte, table(main_dem_v1_end_type, main_dem_v1_sample, useNA = "ifany"))
with(aa_apoe_tte, summary(main_dem_v1_end_age))

#---- Filter participants ----
# n = 485 not in Main dementia definition sample and non dementia at survey
aa_apoe_tte %>% filter(!(main_dem_v1_sample == 1)) %>% nrow()
with(aa_apoe_tte %>% filter(!(main_dem_v1_sample == 1)), summary(main_dem_v1_dem_flag))
  
aa_apoe_tte %<>% filter(main_dem_v1_sample == 1) # n = 51194

# n = 0 subjects not aged 60 and above
aa_apoe_tte %<>% mutate(survey_age_r = round(survey_age))
aa_apoe_tte %>% filter(!survey_age_r >= 60) %>% nrow()

# n = 4636 not in main Asian/white group
with(aa_apoe_tte, table(ethnicity_rev, main_dem_v1_end_dem_flag, useNA = "ifany"))
with(aa_apoe_tte, table(main_dem_v1_end_dem_flag, apoe_y, ethnicity_rev, useNA = "ifany"))
table(aa_apoe_tte$ethnicity_rev, useNA = "ifany")
aa_apoe_tte %>% filter(!ethnicity_rev %in% c(2, 3, 5, 9)) %>% nrow()
aa_apoe_tte %>% filter(!ethnicity_rev %in% c(2, 3, 5, 9)) %>%
  with(table(ethnicity_rev, useNA = "ifany"))
aa_apoe_tte %>% filter(ethnicity_rev %in% c(2, 3, 5, 9)) %>% nrow()
# Since South Asian does not have dementia, drop south asian (ethnicity_rev == 1)
aa_apoe_tte %<>%
  filter(ethnicity_rev %in% c(2, 3, 5, 9)) # n = 46558
# filter(ethnicity_rev %in% c(1, 2, 3, 5, 9))
with(aa_apoe_tte, table(ethnicity_rev, useNA = "ifany"))

# n = 1040 with APOE e2e4 alleles
aa_apoe_tte %<>% 
  mutate(apoe_e2 = case_when(str_detect(apoe, "e2") ~ 1,
                             TRUE ~ 0),
         apoe_e2e4 = case_when(str_detect(apoe, "e4e2") ~ 1,
                               TRUE ~ 0))
with(aa_apoe_tte, table(apoe_e2e4, useNA = "ifany"))
with(aa_apoe_tte, table(apoe_e2e4, ethnicity_rev, useNA = "ifany"))
aa_apoe_tte %<>% filter(apoe_e2e4 == 0) # n = 45518

# #---- Save the data ----
save(aa_apoe_tte,
     file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                   "analysis_data/aa_apoe_tte.RData"))

#---- Descriptives ----
#---- Genetic ancestry ----
aa_apoe_tte %>%
  select(subjid, ethnicity_rev, contains("global")) %>%
  pivot_longer(cols = starts_with("global"),
               names_to = "global_ancestry",
               names_prefix = "global_",
               values_to = "est") %>%
  ggplot() +
  geom_histogram(aes(x = est, color = global_ancestry, fill = global_ancestry), alpha = 0.3) +
  theme_bw() +
  facet_wrap(~ethnicity_rev, scales = "free",
             labeller =
               labeller(ethnicity_rev = c("1" = "South Asian",
                                          "2" = "Chinese", "3" = "Japanese",
                                          "5" = "Filipino", "9" = "White")))
fig_list <- list()
for (ethn in c(1, 2, 3, 5, 9)){
  ethn_label <- case_when(ethn == 1 ~ "South Asian",
                          ethn == 2 ~ "Chinese",
                          ethn == 3 ~ "Japanese",
                          ethn == 5 ~ "Filipino",
                          ethn == 9 ~ "White")
  fig_list[[ethn]] <- aa_apoe_tte %>%
    filter(ethnicity_rev == ethn) %>%
    select(subjid, ethnicity_rev, contains("global")) %>%
    pivot_longer(cols = starts_with("global"),
                 names_to = "global_ancestry",
                 names_prefix = "global_",
                 values_to = "est") %>%
    ggplot() +
    geom_histogram(aes(x = est, color = global_ancestry, fill = global_ancestry), alpha = 0.3) +
    facet_wrap(~global_ancestry, scales = "free") +
    theme_bw() +
    labs(title = paste0("Within ", ethn_label))
}
fig_list[[2]]
fig_list[[3]]
fig_list[[5]]
fig_list[[9]]

with(aa_apoe_tte %>% filter(ethnicity_rev == 2), summary(global_eu))
with(aa_apoe_tte %>% filter(ethnicity_rev == 3), summary(global_eu))

#---- Race/ethnicity * pcAnalysisgroup ----
with(aa_apoe_tte %>% mutate(
  ethn = case_when(ethnicity_rev == 2 ~ "Chinese",
                   ethnicity_rev  == 3 ~ "Japanese",
                   ethnicity_rev  == 5 ~ "Filipino",
                   ethnicity_rev  == 9 ~ "White")), 
  table(ethn, pcanalysisgroup, useNA = "ifany"))

with(aa_apoe_tte %>% mutate(
  ethn = case_when(ethnicity_rev == 2 ~ "Chinese",
                   ethnicity_rev  == 3 ~ "Japanese",
                   ethnicity_rev  == 5 ~ "Filipino",
                   ethnicity_rev  == 9 ~ "White")), 
  table(ethn, kit, useNA = "ifany"))

#---- Check missingness ----
aa_apoe_tte %>%
  select(survey_age, female, ethnicity_rev, edu_3, usaborn_rev, 
         marital_2, income_pp, generalhealth_3, ehr_ht_median, 
         sr_bmi, sr_depress, smoking_2, alcohol_3, pa_3,
         main_dem_v1_end_type, main_dem_v1_fu_time, 
         apoe_y, contains("global_")) %>%
  DataExplorer::plot_missing()

missing_summary <- aa_apoe_tte %>%
  select(survey_age, female, ethnicity_rev, edu_3, usaborn_rev, 
         marital_2, income_pp, generalhealth_3, ehr_ht_median, 
         sr_bmi, sr_depress, smoking_2, alcohol_3, pa_3,
         main_dem_v1_end_type, main_dem_v1_fu_time, 
         apoe_y, contains("global_")) %>%
  DataExplorer::profile_missing() %>%
  arrange(pct_missing) %>%
  mutate(pct_missing = paste0(round(pct_missing*100, 2), "%")) %>%
  print(n = Inf)

#---- APOE descriptives before filtering ----
# Run scripts till line 135 and then this chunk for descriptives
with(aa_apoe_tte, table(apoe_y, useNA = "ifany"))

apoe_fortable <- aa_apoe_tte %>%
  mutate(`Race/Ethnicity` = case_when(
    ethnicity_rev == 1 ~ "South Asian",
    ethnicity_rev == 2 ~ "Chinese",
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
  dplyr::select(`Race/Ethnicity`, apoe, apoe_y)

apoe_table <- apoe_fortable %>%
  filter(`Race/Ethnicity` == "South Asian") %>%
  tbl_summary(by = `Race/Ethnicity`)

apoe_table

#---- OLD -----
# # Raw table 1 without any restriction
# table_1 <- data_for_table_1 %>%
#   modify_if(is.labelled, to_factor) %>%
#   # tbl_strata(
#   #   strata = ethnicity_rev,
#   #   .tbl_fun =
#   #     ~ .x %>%
#       tbl_summary(type = all_continuous() ~ "continuous",
#                   digits = list(all_continuous() ~ 1,
#                                 all_categorical() ~ c(0, 1)),
#                   statistic = list(all_categorical() ~ "{n}({p}%)",
#                                    all_continuous() ~ "{mean} ({sd})"),
#                   by = apoe_y,
#                   missing = "ifany",
#                   missing_text = "Missing") %>%
#       add_overall %>%
#       modify_header(label = "") %>%
#       modify_spanning_header(starts_with("stat_") ~ "**APOE**") %>%
#   bold_labels()
# 
# table_1 %>% as_flex_table() %>%
#   flextable::save_as_docx(
#     path = here::here("output", "tables", "table1_prelim.docx"))
# #---- Death age 90+ cleanning ----
# # # Imputed death, not imputd dementia
# # view(aa_apoe_tte %>% filter(main_dem_v1_90flag == 0 & death_90flag == 1) %>% select(
# #   death_age, main_dem_v1_dem_age, last_membership_end_age, main_dem_v1_end_type,
# # ) %>% mutate(death_dem = death_age - main_dem_v1_dem_age))
# # # Imputed dementia, not imputed death
# # view(aa_apoe_tte %>% filter(main_dem_v1_90flag == 1 & death_90flag == 0) %>% select(
# #   death_age, main_dem_v1_dem_age, last_membership_end_age, main_dem_v1_end_type,
# # ) %>% mutate(death_dem = death_age - main_dem_v1_dem_age))
# # # Imputed dementia and imputed death
# # view(aa_apoe_tte %>% filter(main_dem_v1_90flag == 1 & death_90flag == 1) %>% select(
# #   subjid, death_age, main_dem_v1_dem_age, last_membership_end_age, main_dem_v1_end_type,
# # ) %>% mutate(death_dem = death_age - main_dem_v1_dem_age))
# 
# aa_apoe_tte %<>%
#   mutate(death_age_v2 = case_when(
#     main_dem_v1_90flag == 1 & death_90flag == 1 &
#       death_age - main_dem_v1_dem_age < 0 ~ main_dem_v1_dem_age,
#     TRUE ~ death_age))
# # # Sanity check
# # view(aa_apoe_tte %>% filter(main_dem_v1_90flag == 1 & death_90flag == 1) %>% select(
# #   subjid, death_age, death_age_v2, main_dem_v1_dem_age, last_membership_end_age, 
# #   main_dem_v1_end_type) %>% mutate(death_dem = death_age - main_dem_v1_dem_age))
# 
# #---- Censored 90+ cleaning ----
# # Do not need censored 90+ in this study
# aa_apoe_tte %>% 
#   filter(main_dem_v1_end_type == "CENSORED 90+") %>%
#   select(subjid, death_age, main_dem_v1_dem_age, membership_end_age, 
#          last_membership_end_age, main_dem_v1_end_type, main_dem_v1_end_age) %>% 
#   summary()
# 
# aa_apoe_tte %<>%
#   mutate(main_dem_v1_end_type_v2 = 
#            case_when(main_dem_v1_end_type == "CENSORED 90+" ~ "ADMIN CENSORED",
#                      TRUE ~ main_dem_v1_end_type),
#          adm_censor_age_v2 = 
#            case_when(main_dem_v1_end_type == "CENSORED 90+" ~ main_dem_v1_end_age,
#                      TRUE ~ adm_censor_age))
# 
# #---- End of membership & admin censoring cleaning ----
# # aa_apoe_tte %>% 
# #   filter(main_dem_v1_end_type_v2 == "END OF MEMBERSHIP") %>%
# #   select(subjid, death_age_v2, main_dem_v1_dem_age, membership_end_age, last_membership_end_age, 
# #          main_dem_v1_end_type_v2, adm_censor_age, adm_censor, main_dem_v1_end_age) %>% view()
# # aa_apoe_tte %>% 
# #   filter(main_dem_v1_end_type_v2 == "ADMIN CENSORED") %>%
# #   select(subjid, death_age_v2, main_dem_v1_dem_age, membership_end_age, last_membership_end_age, 
# #          main_dem_v1_end_type_v2, adm_censor_age, adm_censor, main_dem_v1_end_age) %>% view()
# 
# aa_apoe_tte %<>%
#   mutate(main_dem_v1_end_type_v2 = case_when(
#     !is.na(main_dem_v1_dem_age) ~ "DEMENTIA",
#     !is.na(death_age_v2) ~ "DEATH",
#     TRUE ~ main_dem_v1_end_type_v2),
#     main_dem_v1_end_age_v2 = case_when(
#       main_dem_v1_end_type_v2 == "DEMENTIA" ~ main_dem_v1_dem_age,
#       main_dem_v1_end_type_v2 == "ADMIN CENSORED" ~ adm_censor_age_v2,
#       main_dem_v1_end_type_v2 == "DEATH" ~ death_age,
#       main_dem_v1_end_type_v2 == "END OF MEMBERSHIP" & 
#         !is.na(last_membership_end_age) ~ last_membership_end_age,
#       main_dem_v1_end_type_v2 == "END OF MEMBERSHIP" & 
#         is.na(last_membership_end_age) ~ membership_end_age),
#     main_dem_v1_fu_time_v2 = main_dem_v1_end_age_v2 - survey_age)
# 
# # # Sanity check
# # table(aa_apoe_tte$main_dem_v1_end_type_v2, useNA = "ifany")
# # summary(aa_apoe_tte$main_dem_v1_end_age_v2)
# # with(aa_apoe_tte, table(main_dem_v1_end_type_v2, main_dem_v1_dem_flag, useNA = "ifany"))
# # with(aa_apoe_tte, table(main_dem_v1_end_type_v2, death_flag, useNA = "ifany"))
# # with(aa_apoe_tte %>% filter(main_dem_v1_end_age_v2 == main_dem_v1_end_age),
# #      table(main_dem_v1_end_type, main_dem_v1_end_type_v2))
# with(aa_apoe_tte %>% filter(main_dem_v1_end_age_v2 != main_dem_v1_end_age),
#      table(main_dem_v1_end_type, main_dem_v1_end_type_v2))
# aa_apoe_tte %>% filter(main_dem_v1_end_age_v2 != main_dem_v1_end_age) %>%
#   select(subjid, contains("main_dem_v1_end")) %>% view()