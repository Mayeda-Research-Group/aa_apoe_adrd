# Descriptives
# Created by: Yingyan Wu
# Nov.16.2023
# Using non-imputed datasets

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}
p_load("here", "readr", "tidyverse", "magrittr", "dplyr", "haven", "labelled", 
       "gtsummary", "cobalt")
#No scientific notation
options(scipen = 999)

# Paths
source(here::here("scripts", "0.paths.R"))

#---- Load the data ----
# load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
#             "analysis_data/aa_apoe_tte_ehr.RData"))
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/aa_apoe_tte_selected_e4all.RData"))

#---- FU time ----
aa_apoe_tte_selected %>%
  mutate(fu_time = round(main_dem_v1_fu_time)) %>%
  dplyr::count(fu_time) %>%
  mutate(prop = prop.table(n)*100)

#---- Table 1 ----
mean(aa_apoe_tte_selected$survey_age) # 70.71
sd(aa_apoe_tte_selected$survey_age) # 7.36
mean(aa_apoe_tte_selected$main_dem_v1_fu_time) # 10.7
sd(aa_apoe_tte_selected$main_dem_v1_fu_time) # 3.85

data_for_table_1 <- aa_apoe_tte_selected %>%
  mutate(income_pp = income_pp/10000) %>%
  select(survey_age, female, ethnicity_rev, edu_3, usaborn_rev, 
         marital_2, income_pp, generalhealth_3, # ehr_ht_median, sr_bmi, 
         sr_depress, prev_htn_flag, prev_diab_flag, 
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
                      # ehr_ht_median = "EHR height, in",
                      # sr_bmi = "Self-reported BMI",
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
                       "End of Membership" = "END OF MEMBERSHIP",
                       "Administratively Censored" = "ADMIN CENSORED"),
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
                  statistic = list(all_categorical() ~ "{p}",
                                   # all_categorical() ~ "{n} ({p})",
                                   all_continuous() ~ "{mean} ({sd})"),
                  by = apoe_y,
                  missing = "no") %>%
      # missing = "ifany",
      # missing_text = "Missing") %>%
      # add_overall %>%
      modify_header(label = "") %>%
      modify_spanning_header(starts_with("stat_") ~ "**APOE**")) %>%
  modify_header(all_stat_cols() ~ 
                  "**{level}**, \nN = {n} ({style_percent(p, digits = 1)}%)") %>%
  bold_labels()

table_1_res %>% as_flex_table()

table_1_res %>% 
  as_hux_xlsx(here::here("output", "tables", "table1_e4all.xlsx"))

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
  as_hux_xlsx(here::here("output", "tables", "table2_apoe_e4all.xlsx"))

#---- Table X. genetic ancestry dist ----
color_palette <- c("#7B5AA3", "#ED9DB2", "#54B663", "#4B4B4B")
aa_apoe_tte_selected  %>%
  select(ethnicity_rev, global_eu, global_ea) %>%
  pivot_longer(cols = c(global_eu, global_ea), 
               names_to = "type", values_to = "global_anc_prop",
               names_prefix = "global_") %>%
  ggplot() +
  geom_density(aes(x = global_anc_prop, color = as.factor(ethnicity_rev),
                   linetype = type, ..scaled..)) +
  scale_color_manual(values = color_palette) +
  scale_linetype_manual(name = "Global ancestry",
                        labels = c("East Asian", "European"),
                        values = c("solid", "dashed")) +
  facet_wrap(~ethnicity_rev, scales = "free",
             labeller =
               labeller(ethnicity_rev = c("2" = "Chinese", "3" = "Japanese",
                                          "5" = "Filipino", "9" = "Non-Latino White"))) +
  theme_bw() +
  labs(x = "Global ancestry proportion", y = "Density",
       shape = "Models")+
  guides(color = "none") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 11))
# save the figure
ggsave(here::here("output", "figures", "dist_genetic_ancestry_prop_e4all.png"),
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

#---- Table 3. Death and no Died ----
data_for_table_3 <- aa_apoe_tte_selected %>%
  mutate(income_pp = income_pp/10000,
         death_fu = case_when(main_dem_v1_end_type == "DEATH" ~ 1,
                              main_dem_v1_end_type != "DEATH" ~ 0)) %>%
  select(apoe_y, survey_age, female, ethnicity_rev, edu_3, usaborn_rev, 
         marital_2, income_pp, generalhealth_3, # ehr_ht_median, sr_bmi, 
         sr_depress, prev_htn_flag, prev_diab_flag, 
         prev_combined_stroke_flag,
         death_fu, main_dem_v1_fu_time) %>%
  set_variable_labels(apoe_y = "APOE e4 carrier",
                      survey_age = "Survey age, years",
                      female = "Female",
                      edu_3 = "Education",
                      usaborn_rev = "US born",
                      marital_2 = "Married or living as married",
                      income_pp = "Household-adjusted income, 10,000 dollars",
                      generalhealth_3 = "General health (self-reported)",
                      # ehr_ht_median = "EHR height, in",
                      # sr_bmi = "Self-reported BMI",
                      sr_depress = "Self-reported depression",
                      prev_htn_flag = "Hypertension diagnosis pre survey",
                      prev_diab_flag = "Diabetes diagnosis pre survey",
                      prev_combined_stroke_flag = "Stroke diagnosis pre survey",
                      death_fu = "Death during follow-up",
                      # sr_dem_kin = "Family member conditions: Dementia/Alzheimer's Disease",
                      main_dem_v1_fu_time = "Follow up time",
                      ethnicity_rev = "Race/ethnicity") %>%
  set_value_labels(apoe_y = c("Yes" = 1, "No" = 0),
                   death_fu = c("Yes" = 1, "No" = 0),
                   edu_3 = c("Less than high school" = 1, 
                             "High school and some college" = 2,
                             "College and above" = 3),
                   marital_2 = c("Yes" = 1, "No" = 0),
                   generalhealth_3 = c("Excellent/Very Good" = 3,
                                       "Good" = 2,
                                       "Fair/Poor" = 1),
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

table_3_res <- data_for_table_3 %>%
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
                  by = death_fu,
                  missing = "no") %>%
      # missing = "ifany",
      # missing_text = "Missing") %>%
      # add_overall %>%
      modify_header(label = "") %>%
      modify_spanning_header(starts_with("stat_") ~ "**APOE**")) %>%
  modify_header(all_stat_cols() ~ 
                  "**{level}**, \nN = {n} ({style_percent(p, digits = 1)}%)") %>%
  bold_labels()

table_3_res %>% as_flex_table()

table_3_res %>% 
  as_hux_xlsx(here::here("output", "tables", "table3_death_table1_e4all.xlsx"))

# ##---- Covariate balance plot ----
# data_for_death_covbal <- aa_apoe_tte_selected %>%
#   mutate(income_pp = income_pp/10000,
#          death_fu = case_when(main_dem_v1_end_type == "DEATH" ~ 1,
#                               main_dem_v1_end_type != "DEATH" ~ 0),
#          across(c(edu_3, generalhealth_3), as.factor)) %>%
#   select(survey_age, female, ethnicity_rev, edu_3, usaborn_rev, 
#          marital_2, income_pp, generalhealth_3, # ehr_ht_median, sr_bmi, 
#          sr_depress, prev_htn_flag, prev_diab_flag, 
#          prev_combined_stroke_flag,
#          death_fu, main_dem_v1_fu_time, 
#          apoe_y)
# 
# cov_bal_vars <- c("survey_age", "female", "edu_3", "usaborn_rev", 
#                   "marital_2", "income_pp", "generalhealth_3", # ehr_ht_median, sr_bmi, 
#                   "sr_depress", "prev_htn_flag", "prev_diab_flag", 
#                   "prev_combined_stroke_flag",
#                   "main_dem_v1_fu_time", 
#                   "apoe_y")
# 
# smd_cal <- function(data, ethn, treat){
#   # data = data_for_death_covbal
#   # treat = "death_fu"
#   data <- data %>% filter(ethnicity_rev == ethn) %>%
#     as.data.frame()
#   smd <- cobalt::col_w_smd(
#     mat = data %>% select(all_of(cov_bal_vars)),
#     treat = data[, treat],
#     std = TRUE,
#     abs = F,
#     s.d.denom = "control")# (mean(treated)-mean(control))/sd(control)
#   # # (mean(died)-mean(not died))/sd(not died)
#   return(smd)
# }
# 
# race_val <- c(2, 3, 5, 9)
# smd_race_death <- 
#   lapply(race_val, function(x) smd_cal(data_for_death_covbal, x, "death_fu")) %>%
#   bind_rows() %>%
#   t() %>% as.data.frame() %>%
#   mutate(var = rownames(.)) %>%
#   select(var, everything()) %>%
#   set_colnames(c("var", "Chinese", "Japanese", "Filipino", "White")) %>%
#   pivot_longer(cols = c(Chinese, Japanese, Filipino, White), 
#                names_to = "ethn", values_to = "mean") %>%
#   mutate(var_labels = 
#            case_when(var == "survey_age" ~ "Survey age",
#                      var == "female" ~ "Female",
#                      var == "edu_3_1" ~ "Less than high school",
#                      var == "edu_3_2" ~ "High school or some college",
#                      var == "edu_3_3" ~ "College degree",
#                      var == "usaborn_rev" ~ "US-born",
#                      var == "marital_2" ~ "Married or living as married",
#                      var == "income_pp" ~ "Household-adjusted income",
#                      var == "generalhealth_3_3" ~ "General health: Excellent/Very Good",
#                      var == "generalhealth_3_2" ~ "General health: Good",
#                      var == "generalhealth_3_1" ~ "General health: Fair/Poor",
#                      var == "sr_depress" ~ "Self-reported history of depression",
#                      var == "prev_htn_flag" ~ "Hypertension diagnosis",
#                      var == "prev_diab_flag" ~ "Diabetes diagnosis",
#                      var == "prev_combined_stroke_flag" ~ "Stroke diagnosis",
#                      var == "main_dem_v1_fu_time" ~ "Follow up time",
#                      var == "apoe_y" ~ "APOE e4 carriers"))
# 
# color_palette <- c("#7B5AA3", "#ED9DB2", "#54B663", "#4B4B4B")
# smd_race_death %>%
#   ggplot(aes(x = mean, y = var_labels, group = ethn, color = ethn)) +
#   geom_point(size = 1) +
#   geom_vline(xintercept = 0) +
#   geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", 
#              color = "dark grey") +
#   labs(x = "SMD\n(noncensor-censor)/SD[noncensor]", 
#        y = NULL) +
#   theme_bw()
#   # theme(
#   #   axis.title.y = element_text(size = 10),
#   #   legend.position = "none", 
#   #   aspect.ratio = 5 / 3)
# # xlim(-2.1, 1.1)

