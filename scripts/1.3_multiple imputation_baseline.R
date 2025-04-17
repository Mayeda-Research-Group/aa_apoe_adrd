# Multiple imputation for baseline characteristics
# Created by: Yingyan Wu
# Mar.19.2025

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")) {
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("haven", "tidyverse", "magrittr", "mice", "DataExplorer")

options(scipen = 999, digits = 8)

# Paths
source(here::here("scripts", "0.paths.R"))

#---- Load the data ----
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/aa_apoe_tte_selected.RData"))

# DataExplorer::profile_missing(aa_apoe_tte_selected %>% 
#                                 select(-subjid)) %>%
#   mutate(pct_missing = pct_missing*100) %>%
#   arrange(pct_missing) %>%
#   print(n = Inf)

# ---- prepare pre-imputation dataset ----
aa_apoe_tte_pre_mi <- aa_apoe_tte_selected %>% 
  # select variables
  select(subjid, 
         apoe, main_dem_v1_end_type, main_dem_v1_fu_time, 
         # demographics
         survey_age, female, ethnicity_rev, education_rev, usaborn_rev, 
         maritalstatus, sizeofhh, income, income_pp,
         # Health related
         generalhealth, ehr_ht_median, 
         sr_bmi, sr_depress, prev_htn_flag, prev_diab_flag, prev_combined_stroke_flag,
         sr_dem_kin,
         # Behavioral
         smoking_status, alcohol_3, pa_3)

# ---- check missingness ----
DataExplorer::plot_missing(aa_apoe_tte_pre_mi %>% select(-subjid))

missing_summary <- DataExplorer::profile_missing(aa_apoe_tte_pre_mi %>% 
                                                   select(-subjid)) %>%
  mutate(pct_missing = pct_missing*100) %>%
  arrange(pct_missing) %>%
  print(n = Inf)

ordered_var_list <- paste0(missing_summary$feature)

# # Assess missingness by race/ethnicity (Supplement)
# # note that this is before subsetting on pre-survey membership criteria
# # therefore total sample size is larger compared to Table 1
# 
# missing_x_ethn <- data.frame(
#   varname = impute.var.list, 
#   miss_w = NA, # White
#   miss_a = NA, # Asian
#   miss_c = NA, # Chinese
#   miss_f = NA, # Filipino
#   miss_j = NA, # Japanese
#   miss_sa = NA # South Asian
# )
# 
# col_ethn <- names(missing_x_ethn)[-1]
# row.names(missing_x_ethn) <- impute.var.list
# 
# conds <- c(
#   "ASIAN == 0", # White
#   "ASIAN == 1", # Asian
#   "ETHNICITY_REV ==  2", # Chinese
#   "ETHNICITY_REV ==  5", # Filipino
#   "ETHNICITY_REV ==  3", # Japanese
#   "ETHNICITY_REV ==  1" # South Asian
# )
# 
# for (i in impute.var.list){
#   subset_data <- tte_data_pre_mi %>% 
#     filter(PRESURV5YR_SAMPLE == 1) %>% 
#     select(all_of(i), ETHNICITY_REV, ASIAN)
#   for (j in 1:6) {
#     missing_x_ethn[i, col_ethn[j]] <- subset_data %>% 
#       filter(!!parse_expr(conds[j])) %>% 
#       summarise(pctmiss = sum(is.na(get(i))) / n() * 100) %>% as.numeric()
#   }
# }
# 
# table(tte_data_pre_mi$ASIAN)
# table(tte_data_pre_mi$ETHNICITY_REV)
# 
# missing_x_ethn
# missing_x_ethn[order(missing_x_ethn$miss_a), ]
# # output table 
# # write.xlsx(missing_x_ethn, 
# #            file = paste0(path_to_box, 
# #                          "Asian_Americans_dementia/Manuscripts/AA_ADRD_diabetes/",
# #                          "Code/Cleaned_Scripts/output/supplement/",
# #                          "imp_summary_missing_x_ethn.xlsx"))

#---- Data prep ----
# prep data by ordering by missingness
impute_data <- aa_apoe_tte_pre_mi[, ordered_var_list] %>%
  as.data.frame() %>% print()
ncol(impute_data)
colnames(impute_data)
str(impute_data)

# Variable types
variable_length_mat <- impute_data %>%
  map_dfr(function(x) length(unique(x))) %>%
  pivot_longer(cols = everything(),
               names_to = "Var",
               values_to = "length")
table(variable_length_mat$length, useNA = "ifany")
# view(variable_length_mat %>% filter(length <= 30 & !Var %in% binary_vars))

(binary_vars <- 
    variable_length_mat[variable_length_mat$length %in% c(2, 3), "Var"][["Var"]])

ordinal_vars <- c("apoe", "education_rev", "income", "sizeofhh", "smoking_status", 
                  "pa_3", "alcohol_3")

cat_vars <- c("ethnicity_rev", "maritalstatus", "generalhealth", "main_dem_v1_end_type")

(cont_vars <- variable_length_mat %>%
    filter(!Var %in% c(binary_vars, ordinal_vars, cat_vars)) %>% pull(Var))

# Check
length(c(cont_vars, binary_vars, cat_vars, ordinal_vars))
ncol(impute_data)
sum(duplicated(c(cont_vars, binary_vars, cat_vars, ordinal_vars)))

# Check types
impute_data %>%
  map_dfc(function(x) class(x)[1]) %>%
  t() %>% print()

# Set classes by type
impute_data %<>%
  mutate(across(all_of(binary_vars), \(x) as.factor(x)),
         across(all_of(ordinal_vars), \(x) factor(x, ordered = T)),
         across(all_of(cat_vars), \(x) as.factor(x)),
         across(all_of(cont_vars), \(x) as.numeric(x)))

#---- Initiate imputation (PMM) ----
ini_imp <- mice(impute_data, maxit = 0,
                defaultMethod = c("pmm", "pmm", "pmm", "pmm"),
                seed = 62283)
# ini_imp$loggedEvents
(meth <- ini_imp$method)
(pred <- ini_imp$predictorMatrix)

#change predictor matrix so income_pp doesn't predict income or sizeofhh
pred[c("income", "sizeofhh"), "income_pp"] <- 0
pred

#---- Run imputation ----
tictoc::tic()
imp_pmm_20 <- mice(impute_data, m = 20, maxit = 10, 
                   pred = pred, meth = meth,
                   n.burn = 10,
                   defaultMethod = c("pmm", "pmm", "pmm", "pmm"),
                   seed = 62283)
tictoc::toc()
# 10 iterations 20 imp: 10 min

imp_pmm_20$loggedEvents# null
plot(imp_pmm_20)
densityplot(imp_pmm_20)

save(imp_pmm_20, file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                               "analysis_data/imp_pmm_bl_20.RData"))

#---- Stack imputed dataset ----
load(file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                   "analysis_data/imp_pmm_bl_20.RData"))

imp_num = 20
imp_temp <- list()
for (i in 1:imp_num){
  imp_temp[[i]] <- complete(get(paste0("imp_pmm_", imp_num)),
                            action = i)
  imp_temp[[i]] <- cbind(aa_apoe_tte_pre_mi %>% select(subjid), 
                         imp_temp[[i]][, colnames(impute_data)])
  imp_temp[[i]][, "imp"] <- i
}
imp_pmm_stacked <- do.call(rbind, imp_temp)

colnames(imp_pmm_stacked)

# Reconstruct the cleaned variables
imp_pmm_stacked %<>%
  mutate(edu_3 = as.factor(case_when(education_rev %in% c(1, 2) ~ 1,
                           education_rev %in% c(3, 4) ~ 2,
                           education_rev %in% c(5, 6) ~ 3)),
         marital_2 = case_when(maritalstatus == 2 ~ 1,
                               maritalstatus %in% c(1, 3, 4) ~ 0),
         generalhealth_3 = case_when(generalhealth %in% c(1, 2) ~ 1,
                                     generalhealth == 3 ~ 2,
                                     generalhealth %in% c(4, 5) ~ 3),
         sizeofhh = as.numeric(as.character(sizeofhh)),
         smoking_2 = case_when(smoking_status %in% c(2, 3) ~ 1,
                               smoking_status == 1 ~ 0),
         apoe_y = case_when(str_detect(apoe, "e3e4|e4e4") ~ 1,
                            TRUE ~ 0)) %>%
  select(-education_rev, -maritalstatus, -generalhealth, -smoking_status,
         -sizeofhh, -income, -apoe)

#---- Merge datasets with wide and long pre imputed datasets ----
# Wide
aa_apoe_tte_selected_imp <- imp_pmm_stacked %>%
  left_join(aa_apoe_tte_selected %>% 
              select(-colnames(imp_pmm_stacked %>% select(-subjid, -imp))),
            by = "subjid")

# Long
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/aa_apoe_long_tte_selected.RData"))
aa_apoe_long_tte_selected_imp <- imp_pmm_stacked %>%
  left_join(aa_apoe_long_tte_selected %>% 
              select(-colnames(
                imp_pmm_stacked %>% 
                  select(-subjid, -imp, -main_dem_v1_fu_time, -survey_age, 
                         -prev_htn_flag, -prev_diab_flag, -prev_combined_stroke_flag,
                         -ehr_ht_median, -sr_bmi, -sr_dem_kin))),
            by = "subjid",
            relationship = "many-to-many")
# Sanity check
with(aa_apoe_tte_selected_imp, table(imp, useNA = "ifany"))
with(aa_apoe_long_tte_selected_imp, table(imp, useNA = "ifany"))
# diffdf::diffdf(aa_apoe_long_tte_selected_imp %>% filter(imp == 1) %>% arrange(subjid),
#                aa_apoe_long_tte_selected %>% arrange(subjid))

#---- Save imputed data ----
saveRDS(aa_apoe_tte_selected_imp, 
        file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                      "analysis_data/imp_pmm_20_stacked_wide.RDS"))
saveRDS(aa_apoe_long_tte_selected_imp, 
        file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                      "analysis_data/imp_pmm_20_stacked_long.RDS"))
