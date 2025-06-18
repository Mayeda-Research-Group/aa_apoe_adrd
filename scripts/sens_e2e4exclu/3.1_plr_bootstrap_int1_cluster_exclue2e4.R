# Pooled logistic regression CI -- Cluster scripts
# Created by: Yingyan Wu
# Feb.05.2024
# use nonimputed datasets

#---- Argument ----
# Read in the arguments listed in the:
# R CMD BATCH --no-save --no-restore "--args scenario_num=$SGE_TASK_ID"  
## expression:

args = (commandArgs(TRUE))

# Check to see if arguments are passed and set default values if not.
# If so, parse the arguments. (I only have one argument here.)
if (length(args) == 0) {
  print("No arguments supplied.")
  ##supply default values
  scenario_num <- 1
} else {
  for (i in 1:length(args)) {
    eval(parse(text = args[[i]])) # parse each argument (only have 1 here)
  }
}

#---- Package loading ----
library(tidyverse)
library(magrittr)
# No scientific notation
options(scipen = 999)

#---- Load the data ----
# Paths (local)
# source(here::here("scripts", "0.paths.R"))
# load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
#             "analysis_data/aa_apoe_long_tte_selected_exclue2e4.RData"))

# Paths cluster
project_folder <- paste0("/u/home/y/yingyanw/aa_apoe_adrd/")
load(paste0(project_folder, "data/aa_apoe_long_tte_selected_exclue2e4.RData"))

#---- plr function ----
plr_func <- function(data, b, formula, bootstrap = TRUE){
  # data = aa_apoe_long_tte_selected_exclue2e4
  # i_boot = 1
  # b = 1
  # formula = model_formulas
  # bootstrap = i_boot != 0
  #---- 0. Data preparation ----
  # Filter according to race/ethnicity groups
  data_filtered <- data %>% 
    filter(ethnicity_rev %in% as.integer(unlist(strsplit(subsample, ","))))
  # Obtain sampling with replacement ids within the specific race/ethnicity groups
  if (bootstrap == TRUE) {
    set.seed(b)
    boot_ids <- sample(unique(data_filtered$subjid), 
                       length(unique(data_filtered$subjid)), 
                       replace = TRUE)
  } else {
    boot_ids <- unique(data_filtered$subjid)
  }
  
  boot_data <- tibble(subjid = boot_ids, 
                      # a new ID is created to distinguish between these repeats. 
                      new_id = 1:length(boot_ids)) %>% 
    left_join(data, by = "subjid", relationship = "many-to-many") %>% 
    select(-subjid)
  
  #---- 1. Pooled-logistic regression models ----
  tryCatch({
    plr_model_1 <<- glm(as.formula(formula[[1]]),
                        data = boot_data, family = binomial("logit"))
    w_glm_1 <<- c(w_glm_1, 0)
  }, warning = function(w) {
    w_glm_1 <<- c(w_glm_1, 1)
    plr_model_1 <<- glm(as.formula(formula[[1]]),
                        data = boot_data, family = binomial("logit"))
  })
  
  tryCatch({
    plr_model_2 <<- glm(as.formula(formula[[2]]),
                        data = boot_data, family = binomial("logit"))
    w_glm_2 <<- c(w_glm_2, 0)
  }, warning = function(w) {
    w_glm_2 <<- c(w_glm_2, 1)
    plr_model_2 <<- glm(as.formula(formula[[2]]),
                        data = boot_data, family = binomial("logit"))
  })
  
  tryCatch({
    plr_model_3 <<- glm(as.formula(formula[[3]]),
                        data = boot_data, family = binomial("logit"))
    w_glm_3 <<- c(w_glm_3, 0)
  }, warning = function(w) {
    w_glm_3 <<- c(w_glm_3, 1)
    plr_model_3 <<- glm(as.formula(formula[[3]]),
                        data = boot_data, family = binomial("logit"))
  })
  
  plr_models <- list(plr_model_1, plr_model_2, plr_model_3)
  #---- 2. Predict survival ----
  cloned_data <- boot_data %>%
    select(new_id, survey_age_r, female, ethnicity_rev, apoe_y,
           global_eu, global_ea, paste0("eupc", 1:10),
           paste0("eapc", 1:6)) %>%
    unique() %>%
    expand_grid(fu_yr = 0:18)
  
  for (m in 1:3){
    cloned_data %<>%
      dplyr::mutate(!!paste0("pdem1_", "model_", m) := 
                      predict(plr_models[[m]], 
                              newdata = cloned_data %>% mutate(apoe_y = 1),
                              type = "response"),
                    !!paste0("pdem0_", "model_", m) := 
                      predict(plr_models[[m]], 
                              newdata = cloned_data %>% mutate(apoe_y = 0),
                              type = "response"))
  }
  
  
  #---- 3. Implement the cross-product Kaplan Meier method ----
  RD_RR_time <- tibble(.rows = length(0:18))
  for (m in 1:3){
    RD_RR_time %<>% cbind(
      cloned_data %>%
        dplyr::arrange(new_id, fu_yr) %>%
        dplyr::group_by(new_id) %>%
        # Cumulative survival probs
        dplyr::mutate(
          p1 = !!sym(paste0("pdem1_", "model_", m)),
          p0 = !!sym(paste0("pdem0_", "model_", m)),
          s1 = 1 - p1,
          s0 = 1 - p0,
          cum_s1 = cumprod(s1),
          cum_s0 = cumprod(s0),
          cum_s1_lag = lag(cum_s1),
          cum_s0_lag = lag(cum_s0),
          cif1 = ifelse(fu_yr == 0, p1, p1 * cum_s1_lag),
          cif0 = ifelse(fu_yr == 0, p0, p0 * cum_s0_lag),
          cum_cif1 = cumsum(cif1),
          cum_cif0 = cumsum(cif0)) %>%
        dplyr::group_by(fu_yr) %>%
        dplyr::summarise(
          # mean survival probs at each follow up time
          mean_cums1 = mean(cum_s1, na.rm = T),
          mean_cums0 = mean(cum_s0, na.rm = T),
          mean_cif1 = mean(cum_cif1, na.rm = T),
          mean_cif0 = mean(cum_cif0, na.rm = T),
          # mean risk at each follow up time
          mean_cif1_perc = mean_cif1*100,
          mean_cif0_perc = mean_cif0*100,
          # RD at each follow up time
          RD = mean_cif1 - mean_cif0,
          RD_perc = mean_cif1_perc - mean_cif0_perc,
          # RR at each follow up time
          RR = mean_cif1/mean_cif0) %>%
        dplyr::ungroup() %>%
        select(-fu_yr) %>%
        dplyr::rename_all(function(x) paste0(x, "_model_", m)))
  }
  
  RD_RR_time %<>%
    mutate(fu_yr = 0:18,
           race = race_varlabel) %>%
    select(race, fu_yr, everything())
  
  RD_RR_time_wide <- RD_RR_time %>%
    pivot_wider(names_from = fu_yr,
                names_sep = "_y",
                values_from = !c(race, fu_yr))
  
  return(RD_RR_time_wide)
}

#----- scenario tibble ----
n_bootstrap_per_job = 200
n_bootstrap = 1000
scenario_tib <- expand_grid(
  race_varlabels = c("overall", "chn", "jpn", "phl", "wht", "AA"),
  boot_i = c(paste0(seq(1, n_bootstrap-n_bootstrap_per_job + 1, 
                        by = n_bootstrap_per_job), "-",
                    seq(n_bootstrap_per_job, n_bootstrap, 
                        by = n_bootstrap_per_job)))) %>%
  mutate(subsamples = case_when(race_varlabels == "overall" ~ "2, 3, 5, 9",
                                race_varlabels == "chn" ~ "2",
                                race_varlabels == "jpn" ~ "3",
                                race_varlabels == "phl" ~ "5",
                                race_varlabels == "wht" ~ "9",
                                race_varlabels == "AA" ~ "2, 3, 5")) %>%
  separate(boot_i, into = c("boot_i_start", "boot_i_end"), sep = "-") %>%
  mutate_at(vars(boot_i_start), ~ifelse(boot_i_start == 1, 0, .))

# Make sure that scenario_num is not greater than all the scenario index
if (scenario_num > nrow(scenario_tib)) {
  print("out of bound senario number")
  # supply default value
  scenario_num <- 1
}

# print values just to make sure:
print(scenario_num)
scenario <- scenario_tib[scenario_num, ]
race_varlabel <- scenario$race_varlabels
subsample <- scenario$subsamples
boot_i_start <- scenario$boot_i_start
boot_i_end <- scenario$boot_i_end

print(race_varlabel)
print(subsample)
print(boot_i_start)
print(boot_i_end)

w_glm_1 <- w_glm_2 <- w_glm_3 <- numeric()

##---- Formulas ----
formula_1 <- "event_dem ~ apoe_y*poly(fu_yr, 2, raw = T) + poly(survey_age_r, 2, raw = T) + female"
pcs <- case_when(race_varlabel %in% c("overall", "wht") ~ 
                   paste0("eupc", 1:10, collapse = " + "), 
                 race_varlabel %in% c("chn", "jpn", "phl", "AA") ~ 
                   paste0("eapc", 1:6, collapse = " + "))

model_formulas <- list(
  # Model 1: age, sex
  formula_1,
  # Model 2: + EU Ancestry + EA Ancestry
  formula_2 = paste0(formula_1, " + ", "global_eu + global_ea"),
  # Model 3: + PCs (all relevant PCs)
  formula_3 = paste0(formula_1, " + ", pcs))

#---- Bootstrap ----
for (i_boot in boot_i_start:boot_i_end) {
  print(i_boot)
  # i_boot <- 0 # point estimate, using full dataset, no resampling
  # i_boot <- 1 # and so on: bootstrap
  start <- Sys.time() # record starting time
  result <- plr_func(aa_apoe_long_tte_selected_exclue2e4, i_boot,  
                     formula = model_formulas,
                     bootstrap = i_boot != 0)
  out <- tibble(i_boot = i_boot, 
                scenario_num = scenario_num) %>%
    cbind(., result)
  end <- Sys.time() # record ending time
  saveRDS(out, paste0(project_folder, "/model_results/plr_bootstrap/",
                      paste0(race_varlabel,
                             "_RD_RR_time_wide_", i_boot, ".RDS")))
  warning_summary <- tibble(
    i_boot = i_boot,
    race = race_varlabel,
    model = model_formulas, 
    scenario = scenario_num,
    duration = difftime(end, start, units = 'mins'),
    # 3 rows but it actually was within one run
    w_glm_1 = sum(w_glm_1), 
    w_glm_2 = sum(w_glm_2),
    w_glm_3 = sum(w_glm_3))
  
  saveRDS(warning_summary,
          paste0(project_folder, "/model_results/warnings/",
                 paste0(race_varlabel,"_warning_", i_boot, ".RDS")))
}
