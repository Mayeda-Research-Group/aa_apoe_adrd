# Long data construction
# Created by: Yingyan Wu
# Jun.23.2023
# Adapted from Yixuan (Juliet) Zhou's long dataset constrcution code from stroke dementia project

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")){
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("here", "readr", "tidyverse", "magrittr", "dplyr", "haven", "labelled", 
       "gtsummary", "survival")
#No scientific notation
options(scipen = 999)

# Paths
source(here::here("scripts", "0.paths.R"))

#---- Load the data ----
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/aa_apoe_tte_e4all.RData"))
# EHR baseline
load(paste0(path_to_box, "asian_americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/ehr_bl_e4all.RData"))
# EHR long
load(paste0(path_to_box, "asian_americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/bp_followup_before90_e4all.RData"))
load(paste0(path_to_box, "asian_americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/bmi_followup_e4all.RData"))

#---- pre-processing the data ----
aa_apoe_tte_ehr <- aa_apoe_tte %>%
  left_join(., ehr_tib, by = c("subjid", "survey_age", "diab_reg_age",
                               "first_htn_dx_age", "survey_age_r"))

#---- 90 + rules ----
# Carry covariates at 89 over to 90 and 90 +
# Keep the end age (90+ imputation)
aa_apoe_tte_ehr %<>%
  dplyr::mutate(across(c(main_dem_v1_end_age, combined_stroke_age), round, 
                       .names = "{.col}_r"),
                max_fu_yr = main_dem_v1_end_age_r - survey_age_r,
                event_endmem = ifelse(main_dem_v1_end_type == "END OF MEMBERSHIP", 1, 0),
                event_death = ifelse(main_dem_v1_end_type == "DEATH", 1, 0), 
                event_dem = ifelse(main_dem_v1_end_type == "DEMENTIA", 1, 0),
                start = -0.1)

colSums(is.na(aa_apoe_tte_ehr %>% select(contains("event"))))
# there are 233 subjects whose max_fu_yr is 0 because of rounding
aa_apoe_tte_ehr %>% filter(max_fu_yr == 0) %>% nrow()
table(aa_apoe_tte_ehr$max_fu_yr == 0, aa_apoe_tte_ehr$ethnicity_rev, useNA = "ifany")
table(aa_apoe_tte_ehr$max_fu_yr == 0, aa_apoe_tte_ehr$ethnicity_rev, useNA = "ifany") %>% 
  prop.table(2)
aa_apoe_tte_ehr %>% filter(max_fu_yr == 0) %>%
  mutate(fu_yr = main_dem_v1_end_age - survey_age) %>%
  ggplot() +
  geom_histogram(aes(x = fu_yr*12)) +
  scale_x_continuous(breaks = seq(0, 12, 1), name = "Follow up time (month)")
with(aa_apoe_tte_ehr %>% filter(max_fu_yr == 0) %>%
       mutate(fu_yr = main_dem_v1_end_age - survey_age), 
     table(fu_yr > 0.5, useNA = "ifany"))


# remove these subjects from analysis
# (or not, they shouldn't contribute to tte analysis anyways)
aa_apoe_tte_ehr %<>% filter(max_fu_yr != 0)
nrow(aa_apoe_tte_ehr) # n = 46325 (removed 233 subjects)

with(aa_apoe_tte_ehr, table(ethnicity_rev, useNA = "ifany"))
#---- Long dataset ----
tte_vars_tbl <- aa_apoe_tte_ehr %>%
  select(subjid, survey_age_r, main_dem_v1_end_age_r, 
         main_dem_v1_end_type, max_fu_yr, 
         starts_with("event"), start) %>%
  arrange(subjid)
max(tte_vars_tbl$max_fu_yr) # subjects are followed for up to 18 yrs

# this is the main data we will merge additional info onto 
# For Adniministrative censoring
long_tte_data <- survSplit(
  data = tte_vars_tbl, 
  cut = seq(0, max(tte_vars_tbl$max_fu_yr), 1), 
  start = "start", 
  end = "max_fu_yr", 
  event = "event_endmem"
) %>% 
  dplyr::rename(fu_yr = max_fu_yr)

# For Death
long_Death_data <- survSplit(
  data = tte_vars_tbl, 
  cut = seq(0, max(tte_vars_tbl$max_fu_yr), 1), 
  start = "start", 
  end = "max_fu_yr", 
  event = "event_death"
) %>% 
  dplyr::rename(fu_yr = max_fu_yr)

# For Event (dementia)
long_Dem_data <- survSplit(
  data = tte_vars_tbl, 
  cut = seq(0, max(tte_vars_tbl$max_fu_yr), 1), 
  start = "start", 
  end = "max_fu_yr", 
  event = "event_dem"
) %>% 
  dplyr::rename(fu_yr = max_fu_yr)

# add in correct tv indicators 
long_tte_data$event_death <- long_Death_data$event_death
long_tte_data$event_dem <- long_Dem_data$event_dem

# Check end of membership, death and dementia
with(long_tte_data, table(event_death, event_dem, useNA = "ifany"))
with(long_tte_data, table(event_death, event_endmem, useNA = "ifany"))
with(long_tte_data, table(event_dem, event_endmem, useNA = "ifany"))
endmem_id <- long_tte_data %>% filter(event_endmem == 1) %>% pull(subjid) %>% 
  unique() 
with(long_tte_data %>% filter(subjid %in% endmem_id), table(event_dem))
with(long_tte_data %>% filter(subjid %in% endmem_id), table(event_death))
death_id <- long_tte_data %>% filter(event_death == 1) %>% pull(subjid) %>% 
  unique() 
with(long_tte_data %>% filter(subjid %in% death_id), table(event_dem))

long_tte_data %<>%
  mutate(
    # if censored due to end of membership, then death and dem status are missing
    event_death = ifelse(event_endmem == 1, NA, event_death),
    event_dem = ifelse(event_endmem == 1, NA, event_dem),
    # if died, then dem status is missing
    event_dem = ifelse(event_death == 1, NA, event_dem)
  )

# Sanity check
# summary(long_tte_data$fu_yr) #maximum 18 years
long_tte_data %>% filter(subjid == 928480746) %>% View()
# we see that the last interval starts at yr 17 and ends at yr 18

#---- Time-invarying Covariates and conditions ----
cov_tbl <- aa_apoe_tte_ehr %>% 
  select(subjid, 
         # Time-invariant covariates
         female, ethnicity_rev, 
         edu_3, usaborn_rev, marital_2, income_pp, smoking_2,
         pa_3, alcohol_3, sr_depress,
         # Exposure
         contains("apoe"), starts_with("global"), 
         contains("pc"), contains("analysisgroup"),
         # Health conditions
         contains(c("diab", "htn", "combined_stroke")),
         htn_dx_flag, diab_reg_age:dyslip_dx_90plus_flag, -survey_age_r, 
         -main_dem_v1_end_type) %>%
  mutate(across(ends_with("age"), round))

aa_apoe_long_tte <- cov_tbl %>% right_join(long_tte_data, by = "subjid")
colnames(aa_apoe_long_tte)
colnames(aa_apoe_long_tte %>% select(ends_with("_age")))

#---- EHR conditions ----
aa_apoe_long_tte %<>%
  mutate(diab = case_when(diab_reg_flag == 0 ~ 0,
                          diab_reg_age == 90 ~ 0,
                          diab_reg_age > survey_age_r + fu_yr ~ 0,
                          diab_reg_age <= survey_age_r + fu_yr ~ 1),
         hbp = case_when(is.na(htn_dx_flag) ~ 0,
                         first_htn_dx_age == 90 ~ 0,
                         first_htn_dx_age > survey_age_r + fu_yr ~ 0,
                         first_htn_dx_age <= survey_age_r + fu_yr ~ 1),
         combined_stroke = case_when(is.na(combined_stroke_age) ~ 0,
                                     combined_stroke_age == 90 ~ 0,
                                     combined_stroke_age > survey_age_r + fu_yr ~ 0,
                                     combined_stroke_age <= survey_age_r + fu_yr ~ 1),
         # baseline ami
         prev_ami = first_prev_ami_flag,
         # time varying ami
         inc_ami = case_when(first_inc_ami_flag == 0 ~ 0,
                             first_inc_ami_age == 90 ~ 0,
                             first_inc_ami_age > survey_age_r + fu_yr ~ 0,
                             first_inc_ami_age <= survey_age_r + fu_yr ~ 1),
         chf = case_when(is.na(chf_dx_age) ~ 0,
                         chf_dx_age == 90 ~ 0,
                         chf_dx_age > survey_age_r + fu_yr ~ 0,
                         chf_dx_age <= survey_age_r + fu_yr ~ 1),
         pvd = case_when(is.na(pvd_dx_age) ~ 0,
                         pvd_dx_age == 90 ~ 0,
                         pvd_dx_age > survey_age_r + fu_yr ~ 0,
                         pvd_dx_age <= survey_age_r + fu_yr ~ 1),
         cancer = case_when(is.na(cancer_dx_age) ~ 0,
                            cancer_dx_age == 90 ~ 0,
                            cancer_dx_age > survey_age_r + fu_yr ~ 0,
                            cancer_dx_age <= survey_age_r + fu_yr ~ 1),
         ihd = case_when(is.na(ihd_dx_age) ~ 0,
                         ihd_dx_age == 90 ~ 0,
                         ihd_dx_age > survey_age_r + fu_yr ~ 0,
                         ihd_dx_age <= survey_age_r + fu_yr ~ 1),
         dyslip = case_when(is.na(dyslip_dx_age) ~ 0,
                            dyslip_dx_age == 90 ~ 0,
                            dyslip_dx_age > survey_age_r + fu_yr ~ 0, 
                            dyslip_dx_age <= survey_age_r + fu_yr ~ 1)) 

# Check NAs
colSums(is.na(aa_apoe_long_tte %>% 
                select(subjid, hbp, diab, prev_ami, inc_ami, chf, pvd, cancer, combined_stroke,
                       ihd, dyslip)))
with(aa_apoe_long_tte %>% filter(prev_diab_flag == 1), table(diab, useNA = "ifany"))
with(aa_apoe_long_tte %>% filter(prev_chf_flag == 1), table(chf, useNA = "ifany"))
with(aa_apoe_long_tte %>% filter(prev_combined_stroke_flag == 1), table(combined_stroke, useNA = "ifany"))
# Drop variables
aa_apoe_long_tte %<>%
  select(-contains(c("diab_", "htn_", "prev_ami_", "inc_ami_", "_dx_", "prev_")))
colnames(aa_apoe_long_tte)

#---- Time-varying covariates ----
##---- BMI ----
bmi_followup_before90 <- bmi_followup %>% filter(age != 90) %>%
  mutate(fu_yr = age - survey_age_r) %>%
  select(subjid, bmi, fu_yr)

bp_followup_before90 <- bp_followup_before90 %>%
  mutate(fu_yr = age_yr - survey_age_r) %>%
  select(subjid, sbp, fu_yr)

aa_apoe_long_tte %<>%
  left_join(bmi_followup_before90, by = c("subjid", "fu_yr"))

aa_apoe_long_tte %<>%
  left_join(bp_followup_before90, by = c("subjid", "fu_yr"))

# # Check NAs
# aa_apoe_long_tte %>% filter(is.na(bmi) & survey_age_r + fu_yr >= 90) %>% 
#   select(subjid, bmi, survey_age_r, fu_yr) %>% view()
# 
# aa_apoe_long_tte %>% filter(is.na(bmi) & survey_age_r + fu_yr < 90) %>% 
#   select(subjid, bmi, survey_age_r, fu_yr) %>% view()
# 
# aa_apoe_long_tte %>% select(subjid, bmi, survey_age_r, fu_yr) %>%
#   view()

#---- Time interval: 2 years for the long dataset ----
# Check missingness
colSums(is.na(aa_apoe_long_tte))
# aa_apoe_long_tte %>% filter(is.na(event_death)) %>% 
#   select(subjid, death_age_v2:event_endmem) %>% view()


##---- bmi ----
bmi_wide <- bmi_followup_before90 %>%
  pivot_wider(names_from = fu_yr, values_from = bmi,
              names_glue = "bmi_{fu_yr}") 

for (i in 1:18){
  bmi_wide %<>%
    mutate(!!paste0("bmi_y", i) := 
             case_when(is.na(!!sym(paste0("bmi_", i))) ~ !!sym(paste0("bmi_", i-1)),
                       !is.na(!!sym(paste0("bmi_", i))) ~ !!sym(paste0("bmi_", i))))
}

bmi_long <- bmi_wide %>%
  dplyr::rename(bmi_y0 = bmi_0) %>%
  select(subjid, contains("bmi_y")) %>%
  pivot_longer(!subjid, names_to = "fu_yr",
               values_to = "bmi", names_prefix = "bmi_y") %>%
  mutate(fu_yr = as.numeric(fu_yr))

##---- sbp ----
sbp_wide <- bp_followup_before90 %>%
  arrange(subjid, fu_yr) %>%
  pivot_wider(names_from = fu_yr, values_from = sbp, names_glue = "sbp_{fu_yr}")

for (i in 1:18){
  sbp_wide %<>%
    mutate(!!paste0("sbp_y", i) := 
             case_when(is.na(!!sym(paste0("sbp_", i))) ~ !!sym(paste0("sbp_", i-1)),
                       !is.na(!!sym(paste0("sbp_", i))) ~ !!sym(paste0("sbp_", i))))
}

sbp_long <- sbp_wide %>%
  dplyr::rename(sbp_y0 = sbp_0) %>%
  select(subjid, contains("sbp_y")) %>%
  pivot_longer(!subjid, names_to = "fu_yr",
               values_to = "sbp", names_prefix = "sbp_y") %>%
  mutate(fu_yr = as.numeric(fu_yr))

##---- Merge ----
aa_apoe_long_tte_temp <- aa_apoe_long_tte %>%
  select(-bmi, -sbp) %>%
  left_join(bmi_long, by = c("subjid", "fu_yr")) %>%
  left_join(sbp_long, by = c("subjid", "fu_yr"))

# Sanity check
aa_apoe_long_tte_temp %>% select(contains(".")) %>% colnames()

# aa_apoe_long_tte %>% select(subjid, fu_yr, bmi, sbp) %>% view()
with(aa_apoe_long_tte_temp %>% filter(fu_yr == 0), table(event_dem, useNA = "ifany"))
table(aa_apoe_tte_ehr$max_fu_yr, useNA = "ifany")

# Two year interval cut
# id_odd_fu <- aa_apoe_tte_ehr %>% filter(max_fu_yr %% 2 == 1) %>% pull(subjid)
# view(aa_apoe_long_tte %>% filter(subjid %in% id_odd_fu) %>% 
#        select(subjid, fu_yr, contains("event"), main_dem_v1_end_age_v2_r, 
#               main_dem_v1_end_type_v2))
# 
# aa_apoe_long_tte_int2 <- aa_apoe_long_tte_temp %>%
#   mutate(fu_yr = case_when(fu_yr == main_dem_v1_end_age_r - survey_age_r &
#                              fu_yr %% 2 == 1 ~ fu_yr + 1,
#                            TRUE ~ fu_yr)) %>%
#   filter(fu_yr %% 2 == 0)

# Sanity check
# The event numbers should not change!
# For wide dataset, treat the event variables the same way as the long dataset
with(aa_apoe_tte_ehr %>% mutate(
  event_death = ifelse(event_endmem == 1, NA, event_death),
  event_dem = ifelse(event_endmem == 1, NA, event_dem),
  event_dem = ifelse(event_death == 1, NA, event_dem)), 
  table(event_dem, max_fu_yr, useNA = "ifany"))
with(aa_apoe_long_tte %>% select(subjid, event_dem, fu_yr) %>% 
       group_by(subjid) %>% 
       arrange(subjid, fu_yr) %>%
       slice_tail() %>% ungroup(), table(event_dem, fu_yr, useNA = "ifany"))

# aa_apoe_long_tte %<>%
#   filter(fu_yr %% 2 == 0 | fu_yr == main_dem_v1_end_age_v2_r - survey_age_r)

#---- Select variables ----
aa_apoe_tte_selected <- aa_apoe_tte_ehr %>% 
  select(-contains("sapc"), -contains("lapc"), -contains("afpc")) %>%
  # remove South Asian, latino, African PCs
  mutate(apoe_e2 = case_when(str_detect(apoe, "e2") ~ 1,
                             TRUE ~ 0),
         apoe_e2e4 = case_when(str_detect(apoe, "e4e2") ~ 1,
                               TRUE ~ 0))

aa_apoe_long_tte_selected <- aa_apoe_long_tte %>% 
  select(-contains("sapc"), -contains("lapc"), -contains("afpc")) %>%
  # remove South Asian, latino, African PCs
  mutate(apoe_e2 = case_when(str_detect(apoe, "e2") ~ 1,
                             TRUE ~ 0),
         apoe_e2e4 = case_when(str_detect(apoe, "e4e2") ~ 1,
                               TRUE ~ 0))

#---- Sensitivity analysis dataset: excluding e2e4 ----
aa_apoe_tte_selected_exclue2e4 <- aa_apoe_tte_selected %>%
  filter(apoe_e2e4 == 0)
aa_apoe_long_tte_selected_exclue2e4 <- aa_apoe_long_tte_selected %>%
  filter(apoe_e2e4 == 0)

#---- Save the data ----
# save(aa_apoe_tte_ehr,
#      file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
#                    "analysis_data/aa_apoe_tte_ehr_e4all.RData"))
save(aa_apoe_tte_selected,
     file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                   "analysis_data/aa_apoe_tte_selected_e4all.RData"))
save(aa_apoe_long_tte,
     file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                   "analysis_data/aa_apoe_long_tte_e4all.RData"))
# save(aa_apoe_long_tte_int2,
#      file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
#                    "analysis_data/aa_apoe_long_tte_int2.RData"))
save(aa_apoe_long_tte_selected,
     file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                   "analysis_data/aa_apoe_long_tte_selected_e4all.RData"))

# Exclude e2e4
save(aa_apoe_tte_selected_exclue2e4,
     file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                   "analysis_data/aa_apoe_tte_selected_exclue2e4.RData"))
save(aa_apoe_long_tte_selected_exclue2e4,
     file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                   "analysis_data/aa_apoe_long_tte_selected_exclue2e4.RData"))
#---- OLD ----
# #---- Dementia * age 90 ----
# table(aa_apoe_tte$main_dem_v1_90flag, useNA = "ifany")
# with(aa_apoe_tte %>% filter(main_dem_v1_90flag == 1), 
#      summary(main_dem_v1_dem_age, useNA = "ifany"))
# with(aa_apoe_tte %>% filter(main_dem_v1_90flag == 1), 
#      summary(main_dem_v1_end_age, useNA = "ifany"))
# 
# with(aa_apoe_tte %>% filter(main_dem_v1_90flag == 1), 
#      table(main_dem_v1_end_type, useNA = "ifany"))
# with(aa_apoe_tte %>% filter(main_dem_v1_90flag == 1 & 
#                               main_dem_v1_end_type == "DEMENTIA"), 
#      summary(dem90_surv_age_agetrunc, useNA = "ifany"))
# with(aa_apoe_tte %>% filter(main_dem_v1_90flag == 1 & 
#                               main_dem_v1_end_type == "DEMENTIA"), 
#      summary(main_dem_v1_end_age, useNA = "ifany"))
# with(aa_apoe_tte %>% filter(main_dem_v1_90flag == 1 & 
#                               main_dem_v1_end_type == "DEMENTIA"), 
#      table(dem90_surv_age_agetrunc == main_dem_v1_end_age, useNA = "ifany"))