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

#---- Data preparation ----
aa_apoe_long_tte_selected %<>%
  mutate(ethnicity_rev = factor(ethnicity_rev, levels = c(9, 2, 3, 5)))

#---- Pooled-logistic regression models ----
##---- Asian aggregated ----
formula_1 <- "event_dem ~ apoe_y*poly(fu_yr, 2, raw = T) + poly(survey_age_r, 2, raw = T) + female +
apoe_y*ethn_asian"
pcs <- # paste0(
  paste0("eupc", 1:10, collapse = " + ")# , " + ")
# paste0("eapc", 1:6, collapse = " + "))

model_formulas <- list(
  # Model 1: age, sex
  formula_1,
  # Model 2: + EU Ancestry + EA Ancestry
  formula_2 = paste0(formula_1, " + ", "global_eu + global_ea"),
  # Model 3: + PCs (all relevant PCs)
  formula_3 = paste0(formula_1, " + ", pcs)
)

plr_models_aggr <- map(model_formulas, ~ glm(as.formula(.x),
                                        data = aa_apoe_long_tte_selected,
                                        family = binomial("logit")))

names(plr_models_aggr) <- paste0("plr_model_", seq_along(plr_models_aggr))

plr_res_aggr <- lapply(plr_models_aggr, \(x) broom::tidy(x)) %>%
  bind_rows(.id = "model")

##---- Asian ethnicity ----
formula_1 <- "event_dem ~ apoe_y*poly(fu_yr, 2, raw = T) + poly(survey_age_r, 2, raw = T) + female +
apoe_y*ethnicity_rev"
pcs <- # paste0(
  paste0("eupc", 1:10, collapse = " + ")# , " + ")
# paste0("eapc", 1:6, collapse = " + "))

model_formulas <- list(
  # Model 1: age, sex
  formula_1,
  # Model 2: + EU Ancestry + EA Ancestry
  formula_2 = paste0(formula_1, " + ", "global_eu + global_ea"),
  # Model 3: + PCs (all relevant PCs)
  formula_3 = paste0(formula_1, " + ", pcs)
)

plr_models_ethn <- map(model_formulas, ~ glm(as.formula(.x),
                                        data = aa_apoe_long_tte_selected,
                                        family = binomial("logit")))

names(plr_models_ethn) <- paste0("plr_model_", seq_along(plr_models_ethn))

plr_res_ethn <- lapply(plr_models_ethn, \(x) broom::tidy(x)) %>%
  bind_rows(.id = "model")

#---- Get p values from the interaction term between APOE and ethnicity ----
plr_res <- plr_res_aggr %>% rbind(plr_res_ethn)
# plr_res %>% 
#   write.xlsx(here::here("output", "tables", "plr_interaction_results.xlsx"))
plr_res <- read.xlsx(here::here("output", "tables", "plr_interaction_results.xlsx"))

plr_res_int_p_value <- plr_res %>%
  filter(str_detect(term, "apoe_y:ethnicity_rev") |
           str_detect(term, "apoe_y:ethn_asian")) %>%
  select(model, term, p.value) %>%
  mutate(p.value_fmt = as.character(round(p.value, 3)))

plr_res_int_p_value_fmt <- plr_res_int_p_value %>%
  select(model, term, p.value_fmt) %>%
  mutate(model = case_when(model == "plr_model_1" ~ "Model 1",
                           model == "plr_model_2" ~ "Model 2",
                           model == "plr_model_3" ~ "Model 3")) %>%
  pivot_wider(names_from = model, values_from = p.value_fmt) %>%
  mutate(ethn = case_when(term == "apoe_y:ethn_asian" ~ "Asian",
                          str_detect(term, "rev2") ~ "Chinese",
                          str_detect(term, "rev3") ~ "Japanese",
                          str_detect(term, "rev5") ~ "Filipino"),
         group = case_when(ethn == "Asian" ~ "All",
                           TRUE ~ "Asian American ethnic groups")) %>%
  arrange(desc(group)) %>%
  select(group, ethn, starts_with("Model"))

write.xlsx(plr_res_int_p_value_fmt, 
           here::here("output", "tables", "plr_interaction_p_values.xlsx"))