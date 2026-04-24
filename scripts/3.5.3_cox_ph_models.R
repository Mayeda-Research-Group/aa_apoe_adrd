# Cox proportional models
# Created by: Yingyan Wu
# use hotdeck imputed datasets

#---- Package loading + options ----
rm(list = ls())
if (!require("pacman")) {
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}
p_load("here", "tidyverse", "broom", "magrittr", "survival")

# No scientific notation
options(scipen = 999)

# Paths
source(here::here("scripts", "0.paths.R"))

#---- Load the data ----
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/aa_apoe_tte_selected_e4all_hotdeck.RData"))

#---- preprocess the data ----
aa_apoe_tte_selected %<>%
  mutate(logT = log(main_dem_v1_fu_time))

#---- COX PH models ----
race_vals <- c(2, 3, 5, 9)
all_asian_vals <- c(2, 3, 5)
subsample <- list(race_vals, 2, 3, 5, 9, all_asian_vals)
race_labels <- c("Overall", "Chinese", "Japanese", "Filipino", "Non-Latino White", "Asian American")
race_varlabels <- c("overall", "chn", "jpn", "phl", "wht", "AA")

##---- Test the PH assumption ----
p_load(survminer)
diag <- survfit(Surv(main_dem_v1_fu_time, event_dem) ~ apoe_y, 
                data = aa_apoe_tte_selected)
plot(diag, col = c("red", "blue"), 
     fun = "cloglog",
     main = "Testing the proportional hazard assumption",
     xlab = "time",
     ylab = "log(-log(survival))")

##---- Cox PH formulas----
# age as time scale
formula_1 <- "Surv(survey_age, main_dem_v1_end_age_r, event_dem) ~ 
apoe_y + female"
model_formulas <- list(
  ###---- Model 1: age, sex ----
  formula_1 = rep(formula_1, length(race_labels)),
  ###---- Model 2: + EU Ancestry ----
  formula_2 = rep(paste0(formula_1, " + ", "global_eu + global_ea"), length(race_labels)),
  ###---- Model 3: + PCs (all relevant PCs) ----
  formula_3s = c(
    # Overall
    paste0(formula_1, " + ", paste0("eupc", 1:10, collapse = " + ")),
    # CHN, JPN, PHL
    rep(paste0(formula_1, " + ", paste0("eapc", 1:6, collapse = " + ")), 3),
    # WHT
    paste0(formula_1, " + ", paste0("eupc", 1:10, collapse = " + ")),
    # AA
    paste0(formula_1, " + ", paste0("eapc", 1:6, collapse = " + "))),
  ###---- Model 4: + edu, usborn, income ----
  formula_4 = rep(paste0(formula_1, " + ", "edu_3 + usaborn_rev + income_pp"), length(race_labels)),
  ###---- Model 5: model 4 + hbp + diab + depression ----
  formula_5 = rep(paste0(formula_1, " + ", "edu_3 + usaborn_rev + income_pp + prev_htn_flag + prev_diab_flag + sr_depress"), length(race_labels)))

##---- COX PH models ----
coxph_models <- list()
coxph_results <- list()
for (m in 1:length(model_formulas)){
  temp <- vector(mode = "list", length = length(race_labels))
  temp_results <- vector(mode = "list", length = length(race_labels))
  for (r in 1:length(race_labels)){
    temp[[r]] <- coxph(
      as.formula(model_formulas[[m]][[r]]),
      data = aa_apoe_tte_selected %>% 
        filter(ethnicity_rev %in% subsample[[r]]))
    temp_results[[r]] <- tidy(temp[[r]], conf.int = T) %>%
      mutate(race = race_labels[[r]],
             HR = exp(estimate),
             p2.5th = exp(conf.low),
             p97.5th = exp(conf.high)) %>%
      select(term, race, HR, everything())
  }
  names(temp) <- race_varlabels
  names(temp_results) <- race_varlabels
  coxph_models[[m]] <- temp
  coxph_results[[m]] <- temp_results
}

save(coxph_results,
     file = paste0(path_to_box, 
                   "Asian_Americans_dementia/Manuscripts/APOE4_ADRD/aa_apoe_adrd/", 
                   "output/R1_model_results/", "coxph_results_e4all_R1.RData"))

##---- HR table ----
load(paste0(path_to_box, 
            "Asian_Americans_dementia/Manuscripts/APOE4_ADRD/aa_apoe_adrd/", 
            "output/R1_model_results/", "coxph_results_e4all_R1.RData"))
HR_forplot <- tibble()

apoe_prev_ethn <- aa_apoe_tte_selected %>%
  mutate(race = "Overall") %>% # overall
  rbind(., aa_apoe_tte_selected %>%
          mutate(race = case_when(ethnicity_rev == 2 ~ "Chinese",
                                  ethnicity_rev == 3 ~ "Japanese",
                                  ethnicity_rev == 5 ~ "Filipino",
                                  ethnicity_rev == 9 ~ "Non-Latino White")),
        aa_apoe_tte_selected %>%
          filter(ethnicity_rev %in% c(2, 3, 5)) %>%
          mutate(race = "Asian American")) %>% # Asian American 
  group_by(race) %>%
  summarize(prev_apoe_y = mean(apoe_y))

HR_table <- tibble(race = race_labels)

for (m in 1:length(model_formulas)){
  HR_forplot %<>% rbind(bind_rows(coxph_results[[m]]) %>% 
                          filter(term %in% "apoe_y") %>%
                          mutate(model = m) %>%
                          select(model, race, HR, p2.5th, p97.5th))
  HR_table %<>% cbind(
    bind_rows(coxph_results[[m]]) %>% 
      filter(term %in% "apoe_y") %>%
      left_join(apoe_prev_ethn, by = "race") %>%
      mutate(!!paste0("PAR Model ", m) := 
               (prev_apoe_y * (HR - 1))/(1 + prev_apoe_y * (HR - 1))) %>%
      mutate(across(c(HR, p2.5th, p97.5th), 
                    ~ sprintf(fmt = '%#.2f', .x))) %>%
      mutate(!!paste0("Model ",m) := 
               paste0(HR, " (", p2.5th, ", ", p97.5th, ")")) %>%
      select(contains("Model")))
}

HR_table_fmt <- HR_table %>%
  filter(race != "Overall") %>%
  mutate(race = factor(race, levels = c("Non-Latino White", "Asian American",
                                        "Chinese", "Japanese", "Filipino")),
         race_grp = ifelse(race %in% c("Non-Latino White", "Asian American"), "All", 
                           "Asian American ethnic groups")) %>%
  relocate(race_grp, race, starts_with("Model"), starts_with("PAR")) %>%
  arrange(race_grp, race)

writexl::write_xlsx(HR_table_fmt, here::here("output", "tables",
                                             "coxph_e4all_R1.xlsx"))

# ##---- HR plot ----
# color_palette <- c("#4B4B4B", "#B69833", "#7B5AA3", "#ED9DB2", "#54B663")
# HR_forplot %>%
#   filter(race != "Overall") %>%
#   mutate(race = factor(race, levels = c("Non-Latino White", 
#                                         "Asian American",
#                                         "Chinese", "Japanese", "Filipino")),
#          race_grp = ifelse(race %in% c("Non-Latino White", "Asian American"), "All", 
#                            "Asian American ethnic groups")) %>%
#   ggplot(aes(x = race, y = HR, color = race, shape = as.factor(model))) +
#   geom_errorbar(aes(ymin = p2.5th, ymax = p97.5th), width = 0.2,
#                 position = position_dodge((width = 0.3))) +
#   geom_point(position = position_dodge(width = 0.3), fill = "white") +
#   scale_shape_manual(values = c(21, 24, 22)) +
#   scale_y_continuous(breaks = seq(0, 5.5, 0.5)) +
#   # , position = "right") +
#   facet_grid(~ race_grp, scales = "free", space = "free",
#              switch = "y") +
#   scale_color_manual(values = color_palette) +
#   theme_bw() +
#   geom_hline(aes(yintercept = 1), linetype = "dashed") +
#   labs(x = element_blank(), y = "Hazard Ratio (95% CI)",
#        shape = "Models")+
#   guides(color = "none") +
#   theme(legend.position = "bottom",
#         axis.text = element_text(size = 10),
#         axis.title = element_text(size = 11),
#         strip.text = element_text(size = 11))
# 
# ggsave(file = here::here("output", "figures", "efigure3_coxph_e4all_model123.png"),
#        device = "png", width = 7, height = 5, units = "in", dpi = 300)

#---- Secondary: F/U time as timescale ----
##---- Cox PH formulas----
# age as time scale
formula_1 <- "Surv(main_dem_v1_fu_time, event_dem) ~ 
apoe_y + female + survey_age"
model_formulas <- list(
  # Model 1: age, sex
  formula_1 = rep(formula_1, length(race_labels)),
  # Model 2: + EU Ancestry
  formula_2 = rep(paste0(formula_1, " + ", "global_eu + global_ea"), length(race_labels)),
  # Model 3: + PCs (all relevant PCs)
  formula_3s = c(
    # Overall
    paste0(formula_1, " + ", paste0("eupc", 1:10, collapse = " + ")),
    # CHN, JPN, PHL
    rep(paste0(formula_1, " + ", paste0("eapc", 1:6, collapse = " + ")), 3),
    # WHT
    paste0(formula_1, " + ", paste0("eupc", 1:10, collapse = " + ")),
    # AA
    paste0(formula_1, " + ", paste0("eapc", 1:6, collapse = " + "))))

##---- COX PH models ----
coxph_models_futime <- list()
coxph_results_futime <- list()
for (m in 1:length(model_formulas)){
  temp <- vector(mode = "list", length = length(race_labels))
  temp_results <- vector(mode = "list", length = length(race_labels))
  for (r in 1:length(race_labels)){
    temp[[r]] <- coxph(
      as.formula(model_formulas[[m]][[r]]),
      data = aa_apoe_tte_selected %>% 
        filter(ethnicity_rev %in% subsample[[r]]))
    temp_results[[r]] <- tidy(temp[[r]], conf.int = T) %>%
      mutate(race = race_labels[[r]],
             HR = exp(estimate),
             p2.5th = exp(conf.low),
             p97.5th = exp(conf.high)) %>%
      select(term, race, HR, everything())
  }
  names(temp) <- race_varlabels
  names(temp_results) <- race_varlabels
  coxph_models_futime[[m]] <- temp
  coxph_results_futime[[m]] <- temp_results
}

save(coxph_results_futime,
     file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                   "model_results/", "coxph_results_futime_e4all.RData"))

#---- HR table ----
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
            "model_results/", "coxph_results_futime_e4all.RData"))
HR_forplot <- tibble()

apoe_prev_ethn <- aa_apoe_tte_selected %>%
  mutate(race = "Overall") %>% # overall
  rbind(., aa_apoe_tte_selected %>%
          mutate(race = case_when(ethnicity_rev == 2 ~ "Chinese",
                                  ethnicity_rev == 3 ~ "Japanese",
                                  ethnicity_rev == 5 ~ "Filipino",
                                  ethnicity_rev == 9 ~ "Non-Latino White")),
        aa_apoe_tte_selected %>%
          filter(ethnicity_rev %in% c(2, 3, 5)) %>%
          mutate(race = "Asian American")) %>% # Asian American 
  group_by(race) %>%
  summarize(prev_apoe_y = mean(apoe_y))

HR_table <- tibble(race = race_labels)

for (m in 1:length(model_formulas)){
  HR_forplot %<>% rbind(bind_rows(coxph_results_futime[[m]]) %>% 
                          filter(term %in% "apoe_y") %>%
                          mutate(model = m) %>%
                          select(model, race, HR, p2.5th, p97.5th))
  HR_table %<>% cbind(
    bind_rows(coxph_results_futime[[m]]) %>% 
      filter(term %in% "apoe_y") %>%
      left_join(apoe_prev_ethn, by = "race") %>%
      mutate(!!paste0("PAR Model ", m) := 
               (prev_apoe_y * (HR - 1))/(1 + prev_apoe_y * (HR - 1))) %>%
      mutate(across(c(HR, p2.5th, p97.5th), 
                    ~ sprintf(fmt = '%#.2f', .x))) %>%
      mutate(!!paste0("Model ",m) := 
               paste0(HR, " (", p2.5th, ", ", p97.5th, ")")) %>%
      select(contains("Model")))
}

HR_table_fmt <- HR_table %>%
  filter(race != "Overall") %>%
  mutate(race = factor(race, levels = c("Non-Latino White", "Asian American",
                                        "Chinese", "Japanese", "Filipino")),
         race_grp = ifelse(race %in% c("Non-Latino White", "Asian American"), "All", 
                           "Asian American ethnic groups")) %>%
  relocate(race_grp, race, starts_with("Model"), starts_with("PAR")) %>%
  arrange(race_grp, race)

writexl::write_xlsx(HR_table_fmt, here::here("output", "tables", "exploratory_analyses",
                                             "coxph_futime_e4all.xlsx"))
