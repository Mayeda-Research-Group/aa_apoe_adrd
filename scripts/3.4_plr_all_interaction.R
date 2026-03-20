# Pooled logistic regression -- local scripts
# Created by: Yingyan Wu
# use nonimputed datasets

#---- Package loading ----
library(tidyverse)
library(magrittr)
library(purrr)
library(openxlsx)
# No scientific notation
options(scipen = 999)

#---- Load the data ----
# Paths (local)
source(here::here("scripts", "0.paths.R"))
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/aa_apoe_long_tte_selected_e4all.RData"))
# 
# # Paths cluster
# project_folder <- paste0("/u/home/y/yingyanw/aa_apoe_adrd/")
# load(paste0(project_folder, "data/aa_apoe_long_tte_selected_e4all.RData"))

##---- Formulas ----
formula_1 <- "event_dem ~ apoe_y*poly(fu_yr, 2, raw = T) + poly(survey_age_r, 2, raw = T) + female +
apoe_y*as.factor(ethnicity_rev)"
pcs <- paste0(# paste0("eupc", 1:10, collapse = " + "), " + ",
              paste0("eapc", 1:6, collapse = " + "))

model_formulas <- list(
  # Model 1: age, sex
  formula_1,
  # Model 2: + EU Ancestry + EA Ancestry
  formula_2 = paste0(formula_1, " + ", "global_eu + global_ea"),
  # Model 3: + PCs (all relevant PCs)
  formula_3 = paste0(formula_1, " + ", pcs)
)

#---- 2. Pooled-logistic regression models ----
plr_models <- map(model_formulas, ~ glm(as.formula(.x),
                        data = aa_apoe_long_tte_selected %>% 
                          filter(ethnicity_rev != 9),
                        family = binomial("logit")))



names(plr_models) <- paste0("plr_model_", seq_along(plr_models))

#---- Get p values from the interaction term between APOE and ethnicity ----

lapply(plr_models, \(x) broom::tidy(x)) %>% 
  write.xlsx(here("output", "tables", "plr_interaction.xlsx"))
