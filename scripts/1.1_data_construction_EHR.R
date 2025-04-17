# Diagnosis & EHR data construction
# Created by: Yingyan Wu
# Adpated from Juliet Zhou's code for stroke dementia project
# Sep.12.2023

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
path_to_cleaned_data <- paste0(path_to_box,
                               "Asian_Americans_dementia_data/analysis_data_tables/")
path_to_raw_data <- paste0(path_to_box,
                           "Asian_Americans_dementia_data/raw_data_tables/")

#---- Load data ----
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/aa_apoe_tte.RData"))

# SUBJID and SURVEY_AGE
apoe_tte_data_survey_age <- aa_apoe_tte %>% 
  select(subjid, survey_age_r)

#--- Diabetes, HBP ----
# Prevalent diabetes at survey
# First recorded hypertension Dx date if available
ehr_tib <- aa_apoe_tte %>% 
  mutate(prev_diab_flag = case_when(diab_reg_age < survey_age ~ 1,
                                    TRUE ~ 0), 
         # whether entry into diabetes registry is 90+ censored
         diab_reg_90plus_flag = case_when(diab_reg_age == 90 ~ 1,
                                          TRUE ~ 0),
         prev_htn_flag = case_when(first_htn_dx_age < survey_age ~ 1,
                                   TRUE ~ 0), 
         htn_dx_90plus_flag = case_when(first_htn_dx_age == 90 ~ 1,
                                        TRUE ~ 0)) %>% 
  select(subjid, survey_age, survey_age_r,
         diab_reg_age, prev_diab_flag, diab_reg_90plus_flag,
         first_htn_dx_age, prev_htn_flag, htn_dx_90plus_flag)

#---- Stroke ----
stroke_data <- read_sas(paste0(path_to_box, "Asian_Americans_dementia_data/",
                               "raw_data_tables/Membership_EHR_mortality_data/",
                               "c10049_all_obs_dx.sas7bdat"))
colnames(stroke_data) <- tolower(colnames(stroke_data))
# basic cleaning
stroke_sub <- stroke_data %>% 
  mutate(subjid = as.numeric(subjid),
         # calculate age at event
         stroke_age = ifelse(
           as.Date(dx_date) == "1600-01-01", 90,
           interval(as.Date("1900-01-01"), as.Date(dx_date)) / years(1))) %>% 
  # exclude subjects not in APOE subset
  filter(subjid %in% aa_apoe_tte$subjid) %>% 
  select(subjid, stroke_age, everything(), -dx_date) %>% 
  arrange(subjid, stroke_age)

# convert to long format all stroke subtype variables
stroke_long <- stroke_sub %>% 
  pivot_longer(cols = acvd:tia,
               names_to = "stroke_subtype",
               values_to = "indicator")

# filter rows that do have a stroke subtype
stroke_long <- stroke_long %>% filter(!is.na(indicator))

# summarize the min age of stroke, for each id and each subtype
stroke_long_minage <- stroke_long %>% 
  group_by(subjid, stroke_subtype) %>% 
  dplyr::summarize(min_age = min(stroke_age)) %>% 
  ungroup()

## convert to a wide format
stroke_wide_minage <- stroke_long_minage %>% 
  pivot_wider(names_from = stroke_subtype,
              values_from = min_age)

# 1. Create a combined category for ACVD and ischemic stroke,
# but prioritize ischemic stroke time, not the min between these 2
stroke_wide_minage <- stroke_wide_minage %>% 
  mutate(isd_acvd = case_when(!is.na(isd) ~ isd,
                              is.na(isd) & !is.na(acvd) ~ acvd))

# 2. For hstroke subtype, if their min age of hstroke is equal to the min age of 
# isd, then prioritize isd and not hstroke, this should not affect any variables
# that combine hstroke and ischemic
stroke_wide_minage %>% count(!is.na(hstroke))
stroke_wide_minage %>% count(!is.na(isd), !is.na(hstroke)) 
stroke_wide_minage %>% count(isd == hstroke)  

stroke_wide_minage <- stroke_wide_minage %>% 
  mutate(hstroke_noties = case_when(
    is.na(isd) & is.na(hstroke) ~ NA_real_,
    (!is.na(isd) & !is.na(hstroke)) & isd == hstroke ~ NA_real_,
    TRUE ~ hstroke))

stroke_wide_minage %>% count(!is.na(hstroke), !is.na(hstroke_noties))
stroke_wide_minage %>% count(!is.na(hstroke_noties))

# 3. Create a combined stroke, no TIA and prioritize hierarchy of ISD > HStroke > ACVD
# this is the main stroke variable
stroke_wide_minage <- stroke_wide_minage %>% 
  mutate(combined_stroke = case_when(!is.na(isd) ~ isd,
                                     is.na(isd) & !is.na(hstroke) ~ hstroke,
                                     is.na(isd) & is.na(hstroke) ~ acvd))


## 4. Create a combined stroke that excludes ACVD
stroke_wide_minage <- stroke_wide_minage %>% 
  mutate(combined_isd_hstroke = case_when(!is.na(isd) ~ isd,
                                          is.na(isd) & !is.na(hstroke) ~ hstroke))

## separately, derive earliest age at any stroke
# icih: Iatrogenic cerebrovascular infarction or hemorrhage
stroke_wide_minage <- stroke_wide_minage %>% 
  mutate(any_stroke = pmin(acvd, isd, tia, hstroke, icih, na.rm = TRUE)) 

ehr_tib <- stroke_wide_minage %>% 
  select(subjid, combined_stroke, any_stroke) %>% 
  right_join(ehr_tib, by = "subjid") 

ehr_tib %<>%
  mutate(prev_combined_stroke_flag = case_when(combined_stroke < survey_age ~ 1, 
                                               TRUE ~ 0), # no dx ever or no dx by survey
         prev_any_stroke_flag = case_when(any_stroke < survey_age ~ 1, 
                                          TRUE ~ 0)) %>% # no dx ever or no dx by survey
  rename(combined_stroke_age = combined_stroke,
         any_stroke_age = any_stroke)

# Sanity check
with(ehr_tib, table(prev_combined_stroke_flag, useNA = "ifany"))

#---- Actute MI ----
## this is derived the same way as the stroke variables that Taylor did 
## i.e. flag and age for first prevalent and first incident events
## this is the case for the rest of the EHR variables in this script
ami_hypogly <- read_sas(paste0(path_to_box, "Asian_Americans_dementia_data/",
                               "raw_data_tables/Membership_EHR_mortality_data/",
                               "c10049_all_obs_dx20220531.sas7bdat"))
colnames(ami_hypogly) <- tolower(colnames(ami_hypogly))
# there could be multiple ami records for each subject 
ami_data <- ami_hypogly %>% 
  # use AMI indicators only
  filter(ami == 1) %>% select(-hypogly) %>% 
  mutate(subjid = as.numeric(subjid),
         # calculate age at event 
         age = ifelse(
           as.Date(dx_date) == "1600-01-01", 90, 
           interval(as.Date("1900-01-01"), as.Date(dx_date)) / years(1))) %>% 
  arrange(subjid, age) %>% 
  left_join(apoe_tte_data_survey_age, by = "subjid") %>% 
  # subset to subjects in the analytical data
  filter(subjid %in% aa_apoe_tte$subjid) %>% 
  mutate(prev_ami_flag = ifelse(age < survey_age_r, 1, 0)) %>% 
  group_by(subjid, prev_ami_flag) %>% 
  # slices the earliest events, one from prevalent dx's, one from incident dx's,
  slice_head() %>% # since we have arranged by age already
  ungroup()

ami_data_wide <- ami_data %>% 
  select(subjid, age, prev_ami_flag) %>% 
  pivot_wider(
    names_from = prev_ami_flag, 
    values_from = age, 
    names_prefix = "first_ami_age_"
  ) %>% 
  dplyr::rename(
    first_prev_ami_age = first_ami_age_1,
    first_inc_ami_age = first_ami_age_0
  ) 

ehr_tib %<>% left_join(ami_data_wide, by = "subjid") %>% 
  mutate(
    first_prev_ami_flag = ifelse(is.na(first_prev_ami_age), 0, 1), 
    # Since dataset is restricted to < 90 yo at survey, prev ami will not be at 90+
    first_inc_ami_flag = ifelse(is.na(first_inc_ami_age), 0, 1), 
    first_inc_ami_90plus_flag = case_when(first_inc_ami_age == 90 ~ 1,
                                          TRUE ~ 0))
#---- CHF and PVD ----
# first recorded Dx for congestive heart failure (CHF)
# first recorded Dx for peripheral vascular disease (PVD)
chf_pvd <- read_sas(paste0(path_to_raw_data, 
                           "Membership_EHR_mortality_data/",
                           "c10049_cvd20220531.sas7bdat"))
colnames(chf_pvd) <- tolower(colnames(chf_pvd))
cvd_tib <- chf_pvd %>%
  select(subjid, contains("chf"), contains("pvd")) %>% 
  mutate(subjid = as.numeric(subjid),
         # calculate age at event 
         chf_dx_age = ifelse(
           as.Date(first_chf_dx_date) == "1600-01-01", 90,
           interval(as.Date("1900-01-01"), as.Date(first_chf_dx_date))/years(1)),
         pvd_dx_age = ifelse(
           as.Date(first_pvd_dx_date) == "1600-01-01", 90,
           interval(as.Date("1900-01-01"), as.Date(first_pvd_dx_date))/years(1))
  ) %>% 
  select(subjid, chf_dx_age, pvd_dx_age)

ehr_tib %<>%
  left_join(cvd_tib, by = "subjid") %>%
  mutate(prev_chf_flag = case_when(chf_dx_age < survey_age ~ 1,
                                   TRUE ~ 0), 
         chf_dx_90plus_flag = case_when(chf_dx_age == 90 ~ 1,
                                        TRUE ~ 0),
         prev_pvd_flag = case_when(pvd_dx_age < survey_age ~ 1,
                                   TRUE ~ 0),
         pvd_dx_90plus_flag = case_when(pvd_dx_age == 90 ~ 1,
                                        TRUE ~ 0))
#---- IHD and cancer ----
# first recorded Dx for cancer, and
# first recorded Dx for ischemic heart disease (IHD) 
ihd_cancer <- read_sas(paste0(path_to_raw_data, 
                              "Membership_EHR_mortality_data/",
                              "c10049_cvd_20230630.sas7bdat")) 
colnames(ihd_cancer) <- tolower(colnames(ihd_cancer))

ihd_cancer_tib <- ihd_cancer %>% 
  mutate(subjid = as.numeric(subjid),
         # calculate age at event 
         cancer_dx_age = ifelse(
           as.Date(first_cancer_dx_date) == "1600-01-01", 90, 
           interval(as.Date("1900-01-01"), as.Date(first_cancer_dx_date)) / years(1)),
         ihd_dx_age = ifelse(
           as.Date(first_ihd_dx_date) == "1600-01-01", 90, 
           interval(as.Date("1900-01-01"), as.Date(first_ihd_dx_date)) / years(1))) %>%
  select(subjid, cancer_dx_age, ihd_dx_age) 

ehr_tib %<>%
  left_join(ihd_cancer_tib, by = "subjid") %>%
  mutate(prev_cancer_flag = case_when(cancer_dx_age < survey_age ~ 1,
                                      TRUE ~ 0), 
         cancer_dx_90plus_flag = case_when(cancer_dx_age == 90 ~ 1,
                                           TRUE ~ 0),
         prev_ihd_flag = case_when(ihd_dx_age < survey_age ~ 1,
                                   TRUE ~ 0), 
         ihd_dx_90plus_flag = case_when(ihd_dx_age == 90 ~ 1,
                                        TRUE ~ 0))
# # Sanity check
# with(ehr_tib, summary(cancer_dx_age))
# with(ehr_tib, summary(ihd_dx_age)

#---- dyslipidemia ----
# first recorded Dx for dyslipidemia
dyslip <- read_sas(paste0(path_to_raw_data, 
                          "Membership_EHR_mortality_data/",
                          "c10049_dyslip_20230630.sas7bdat")) 
colnames(dyslip) <- tolower(colnames(dyslip))
# diffdf::diffdf(chf_pvd %>% select(subjid) %>% mutate(subjid = as.numeric(subjid)) %>%
#                  arrange(subjid),
#                dyslip %>% select(subjid) %>% mutate(subjid = as.numeric(subjid)) %>%
#                  arrange(subjid))

hcl_tib <- dyslip %>% 
  mutate(subjid = as.numeric(subjid),
         # calculate age at event 
         dyslip_dx_age = ifelse(
           as.Date(first_dyslip_dx_date) == "1600-01-01", 90, 
           interval(as.Date("1900-01-01"), as.Date(first_dyslip_dx_date)) / years(1))) %>% 
  select(subjid, dyslip_dx_age)

ehr_tib %<>%
  left_join(hcl_tib, by = "subjid") %>%
  mutate(prev_dyslip_flag = case_when(dyslip_dx_age < survey_age ~ 1,
                                      TRUE ~ 0), 
         dyslip_dx_90plus_flag = case_when(
           dyslip_dx_age == 90 ~ 1,
           TRUE ~ 0))

#---- Time-varying continuous measure ----
#---- height, weight, BMI ----
# ht and wt are medians in each year of age (empirical scale)
ht_wt_data <- read_sas(paste0(path_to_raw_data, 
                              "Membership_EHR_mortality_data/",
                              "c10049_ht_wt.sas7bdat"))
colnames(ht_wt_data) <- tolower(colnames(ht_wt_data))
# filter for subjects that are in the tte dataset 
bmi_tib <- ht_wt_data %>% filter(subjid %in% aa_apoe_tte$subjid)
length(unique(bmi_tib$subjid)) # n = 45229

# there are 289 subjects with no ht or wt data from EHR
# we will impute these values, with the help of self-reported BMI info
length(unique(aa_apoe_tte$subjid))  - length(unique(bmi_tib$subjid)) 

bmi_tib %<>% 
  mutate(subjid = as.numeric(subjid)) %>% 
  arrange(subjid, age) %>% 
  mutate(bmi = wt_median / ht_median^2 * 703) %>% 
  filter(!is.na(bmi))

# ht and wt data are also 90+ censored
# i.e. measurements taken at or after age 90 are not distinguishable
with(bmi_tib %>% filter(age >= 90) %>% group_by(subjid, age) %>%
       dplyr::summarize(n = n()), table(n, useNA = "ifany"))
# bmi_tib %>% group_by(subjid, age) %>% mutate(n = n()) %>% filter(n > 1) %>% view()
# n = 3964
id_bmi_90plus <- bmi_tib %>% filter(age >= 90) %>% group_by(subjid, age) %>%
  dplyr::summarize(n = n()) %>% filter(n > 1) %>% pull(subjid)
set.seed("62283")
id_sub <- sample(id_bmi_90plus, size = 50, replace = F)
bmi_tib %>% 
  filter(subjid %in% id_sub) %>%
  filter(age >= 90) %>%
  group_by(subjid, age) %>% mutate(n = 1:n()) %>%
  ggplot(aes(x = n, y = bmi, group = subjid, color = factor(subjid))) +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "none")

nrow(bmi_tib) # n = 484153
bmi_tib %>% distinct(subjid, age) %>% nrow() # n = 472716

with(bmi_tib %>% filter(age >= 90) %>% group_by(subjid, age) %>%
       dplyr::summarize(n = n()) %>% filter(n > 1), summary(n))

# merge in survey age rounded by yr
bmi_tib %<>% left_join(apoe_tte_data_survey_age, by = "subjid")

##---- baseline ----
bmi_baseline <- bmi_tib %>% filter(age == survey_age_r) %>%
  arrange(subjid, age, bmi) %>%
  group_by(subjid) %>%
  slice(1) %>%
  select(-age) %>%
  dplyr::rename_at(vars(-subjid, -survey_age_r), function(x) paste0("bl_", x))

ehr_tib %<>% left_join(bmi_baseline %>% select(-survey_age_r), by = "subjid")

##---- followup ----
bmi_followup <- bmi_tib %>% filter(age >= survey_age_r)

#---- SBP ----
# derivation of baseline and time-varying BP from EHR is similar to BMI
bp_data <- read_sas(paste0(path_to_raw_data, 
                           "Membership_EHR_mortality_data/",
                           "c10049_bp.sas7bdat"))
colnames(bp_data) <- tolower(colnames(bp_data))
# filter for subjects that are in the tte dataset 
bp_tib <- bp_data %>% filter(subjid %in% aa_apoe_tte$subjid)
length(unique(bp_tib$subjid)) # n = 45335

# there are 183 subjects with no bp data from EHR
length(unique(aa_apoe_tte$subjid)) - length(unique(bp_tib$subjid))

# # there are "categorical" bp's in the dataset
# table(bp_tib$categorical_dbp, useNA = "ifany")
# table(bp_tib$categorical_sbp, useNA = "ifany")

bp_tib %<>%
  mutate(
    subjid = as.numeric(subjid),
    # merge categorical bp's with continuous
    continuous_dbp = ifelse(is.na(continuous_dbp), categorical_dbp, continuous_dbp),
    continuous_sbp = ifelse(is.na(continuous_sbp), categorical_sbp, continuous_sbp),
    # calculate age at event 
    age = ifelse(as.Date(measure_date) == "1600-01-01", 90, 
                 interval(as.Date("1900-01-01"), as.Date(measure_date)) / years(1)),
    age_yr = round(age),
    bp_flag_90plus = case_when(age == 90 ~ 1, TRUE ~ 0)) %>% 
  dplyr::rename(dbp = continuous_dbp, sbp = continuous_sbp) %>% 
  select(-categorical_dbp, -categorical_sbp, -measure_date, -enctype) %>% 
  arrange(subjid, age)

with(bp_tib %>% filter(age_yr == 90), table(bp_flag_90plus, useNA = "ifany"))

bp_median_before90 <- bp_tib %>% 
  filter(bp_flag_90plus == 0) %>%
  group_by(subjid, age_yr) %>% # summarise median bp's 
  dplyr::summarise(
    dbp_median = median(dbp), 
    sbp_median = median(sbp)) %>% 
  ungroup() %>% 
  dplyr::rename(dbp = dbp_median, sbp = sbp_median)

# merge in survey age rounded by yr
bp_median_before90 <- bp_median_before90 %>% 
  left_join(apoe_tte_data_survey_age, by = "subjid")

id_bp_90plus <- bp_tib %>% filter(age >= 90) %>% group_by(subjid, age) %>%
  dplyr::summarize(n = n()) %>% filter(n > 1) %>% pull(subjid)
set.seed("62283")
id_sub <- sample(id_bp_90plus, size = 50, replace = F)
bp_tib %>% 
  filter(subjid %in% id_sub) %>%
  filter(age >= 90) %>%
  group_by(subjid, age) %>% mutate(n = 1:n()) %>%
  ggplot(aes(x = n, y = sbp, group = subjid, color = factor(subjid))) +
  geom_line() +
  theme_minimal() +
  theme(legend.position = "none")

with(bp_tib %>% filter(age >= 90) %>% group_by(subjid, age) %>%
       dplyr::summarize(n = n()) %>% filter(n > 1), table(n))
with(bp_tib %>% filter(age >= 90) %>% group_by(subjid, age) %>%
       dplyr::summarize(n = n()) %>% filter(n > 1), summary(n))
id_weird <- bp_tib %>% filter(age >= 90) %>% group_by(subjid, age) %>%
  dplyr::summarize(n = n()) %>% filter(n > 100) %>% pull(subjid)
bp_tib %>% filter(age >= 90 & subjid %in% id_weird) %>% view()

##---- baseline ----
bp_baseline <- bp_median_before90 %>% filter(age_yr == survey_age_r) %>%
  select(-age_yr) %>%
  dplyr::rename_at(vars(-subjid, -survey_age_r), function(x) paste0("bl_", x))

length(unique(bp_baseline$subjid))
# add into tte_add 
ehr_tib %<>% left_join(bp_baseline %>% select(-survey_age_r), by = "subjid")

##---- followup ----
bp_followup_before90 <- bp_median_before90 %>% filter(age_yr >= survey_age_r)

# #---- Save the data ----
save(ehr_tib, file =
       paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
              "analysis_data/ehr_bl.RData"))
save(bmi_followup, file =
       paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
              "analysis_data/bmi_followup.RData"))
save(bp_followup_before90, file=
       paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
              "analysis_data/bp_followup_before90.RData"))

#---- OLD ----
# # Because round has a weird behavior
# # > round(60.500000000000001)
# # [1] 60
# # > round(60.50000000000001)
# # [1] 61
# bp_tib %>% mutate(age_char = as.character(age)) %>% 
#   filter(str_detect(age_char, fixed(".500000"))) %>% view()
##---- BMI ----
###---- 90+ rules ----
# bmi_test <- bmi_tib %>%
#   left_join(aa_apoe_tte %>% 
#               select(subjid, death_age_v2, main_dem_v1_dem_age, 
#                      membership_end_age, last_membership_end_age, 
#                      main_dem_v1_end_age_v2,
#                      main_dem_v1_end_type_v2), by = "subjid")
# bmi_test %>% filter(age >= 90) %>% view()
# bmi_test %>% filter(age >= 90 & main_dem_v1_end_type_v2 == "END OF MEMBERSHIP")
# # For those who got the event before 90, exclude their bmi records after 90
# bmi_test %<>% filter(!(main_dem_v1_end_age_v2 < 90 & age >= 90))
# length(unique(bmi_test$subjid)) # 46091
# 
# bmi_test_90plus <- bmi_test %>% filter(age >= 90) %>%
#   group_by(subjid, age) %>%
#   mutate(n = 1:n(), n_total = n()) %>% ungroup()
# with(bmi_test_90plus %>% select(subjid, n_total) %>% distinct(), 
#      table(n_total, useNA = "ifany"))
# 
# bmi_test_90plus %>% filter(n_total == 1) %>% view()
# 
# bmi_test %>% filter(age >= 90) %>% view()
# bmi_test %>% filter(age >= 90 & main_dem_v1_end_type_v2 == "DEMENTIA") %>% view()
# bmi_test %>% filter(age >= 90 & main_dem_v1_end_type_v2 == "DEATH") %>% view()
# 
# table(bmi_test$n_total, useNA = "ifany")
##---- SBP ----
###---- 90+ rules ----
# bp_test <- bp_tib %>%
#   left_join(aa_apoe_tte %>% 
#               select(subjid, death_age_v2, main_dem_v1_dem_age, 
#                      membership_end_age, last_membership_end_age, 
#                      main_dem_v1_end_age_v2,
#                      main_dem_v1_end_type_v2), by = "subjid")
# # For those who got the event before 90, exclude their bp records after 90
# bp_test %<>% filter(!(main_dem_v1_end_age_v2 < 90 & age >= 90))
# length(unique(bp_test$subjid)) # 46470
# 
# bp_test_90plus <- bp_test %>% filter(age >= 90) %>%
#   group_by(subjid, age) %>%
#   mutate(n = 1:n(), n_total = n()) %>% ungroup()
# with(bp_test_90plus %>% select(subjid, n_total) %>% distinct(), 
#      table(n_total, useNA = "ifany"))
# 
# bp_test_90plus %>% filter(n_total == 1) %>% view()
# 
# bp_test %>% filter(age >= 90) %>% view()
# bp_test_90plus %>% filter(main_dem_v1_end_type_v2 == "DEMENTIA") %>% view()
# bp_test_90plus %>% filter(main_dem_v1_end_type_v2 == "DEATH") %>% view()
