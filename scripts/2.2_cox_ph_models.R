# Cox proportional models
# Created by: Yingyan Wu
# Nov.20.2023
# use nonimputed datasets

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
            "analysis_data/aa_apoe_tte_selected.RData"))

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
    paste0(formula_1, " + ", paste0("eapc", 1:6, collapse = " + "))))

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
     file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
                   "model_results/", "coxph_results_09162024.RData"))

#---- HR table ----
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
            "model_results/", "coxph_results_09162024.RData"))
HR_forplot <- tibble()
HR_table <- tibble(race = race_labels)
for (m in 1:length(model_formulas)){
  HR_forplot %<>% rbind(bind_rows(coxph_results[[m]]) %>% 
                          filter(term %in% "apoe_y") %>%
                          mutate(model = m) %>%
                          select(model, race, HR, p2.5th, p97.5th))
  HR_table %<>% cbind(
    bind_rows(coxph_results[[m]]) %>% 
      filter(term %in% "apoe_y") %>%
      mutate_if(is.numeric, sprintf, fmt = '%#.2f') %>%
      mutate(!!paste0("Model ",m) := 
               paste0(HR, " (", p2.5th, ", ", p97.5th, ")")) %>%
      select(contains("Model")))
}

HR_table_fmt <- HR_table %>%
  filter(race != "Overall") %>%
  mutate(race = factor(race, levels = c("Non-Latino White", 
                                        "Asian American",
                                        "Chinese", "Japanese", "Filipino")),
         race_grp = ifelse(race %in% c("Non-Latino White", "Asian American"), "All", 
                           "Asian American ethnic groups")) %>%
  relocate(race_grp, race) %>%
  arrange(race_grp, race)

writexl::write_xlsx(HR_table_fmt, here::here("output", "tables",
                                             "coxph_exclue2e4_09162024.xlsx"))

#---- HR plot ----
color_palette <- c("#4B4B4B", "#B69833", "#7B5AA3", "#ED9DB2", "#54B663")
HR_forplot %>%
  filter(race != "Overall") %>%
  mutate(race = factor(race, levels = c("Non-Latino White", 
                                        "Asian American",
                                        "Chinese", "Japanese", "Filipino")),
         race_grp = ifelse(race %in% c("Non-Latino White", "Asian American"), "All", 
                           "Asian American ethnic groups")) %>%
  ggplot(aes(x = race, y = HR, color = race, shape = as.factor(model))) +
  geom_errorbar(aes(ymin = p2.5th, ymax = p97.5th), width = 0.2,
                position = position_dodge((width = 0.3))) +
  geom_point(position = position_dodge(width = 0.3), fill = "white") +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_y_continuous(breaks = seq(0, 5.5, 0.5)) +
  # , position = "right") +
  facet_grid(~ race_grp, scales = "free", space = "free",
             switch = "y") +
  scale_color_manual(values = color_palette) +
  theme_bw() +
  geom_hline(aes(yintercept = 1), linetype = "dashed") +
  labs(x = element_blank(), y = "Hazard Ratio (95% CI)",
       shape = "Models")+
  guides(color = "none") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 11))
# theme(legend.position = "none")

ggsave(file = here::here("output", "figures", "efigure_coxph_exclue2e4_model123.png"),
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

# #---- OLD ----
# ###---- restricted at 13 years ----
# coxph_models_y13 <- vector(mode = "list", length = length(model_formulas))
# coxph_results_y13 <- vector(mode = "list", length = length(model_formulas))
# for (m in 1:3){
#   temp <- vector(mode = "list", length = length(race_labels))
#   temp_results <- vector(mode = "list", length = length(race_labels))
#   for (r in 1:length(race_labels)){
#     temp[[r]] <- coxph(
#       as.formula(model_formulas[[m]][[r]]),
#       data = aa_apoe_tte_selected_y13 %>% 
#         filter(ethnicity_rev %in% subsample[[r]]))
#     temp_results[[r]] <- tidy(temp[[r]], conf.int = T) %>%
#       mutate(race = race_labels[[r]],
#              HR = exp(estimate),
#              p2.5th = exp(conf.low),
#              p97.5th = exp(conf.high)) %>%
#       select(term, race, HR, everything())
#   }
#   names(temp) <- race_varlabels
#   names(temp_results) <- race_varlabels
#   coxph_models_y13[[m]] <- temp
#   coxph_results_y13[[m]] <- temp_results
# }
# 
# save(coxph_results_y13, 
#      file = paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
#                    "model_results/", "coxph_results_y13.RData"))
# 
# HR_forplot_y13 <- tibble()
# HR_table_y13 <- tibble(race = race_labels)
# for (m in 1:3){
#   HR_forplot_y13 %<>% rbind(bind_rows(coxph_results_y13[[m]]) %>% 
#                               filter(term %in% "apoe_y") %>%
#                               mutate(model = m) %>%
#                               select(model, race, HR, p2.5th, p97.5th))
#   HR_table_y13 %<>% cbind(
#     bind_rows(coxph_results_y13[[m]]) %>% filter(term %in% "apoe_y") %>%
#       mutate_if(is.numeric, sprintf, fmt = '%#.2f') %>%
#       mutate(!!paste0("Model ",m, " HR (95% CI)") := 
#                paste0(HR, " (", p2.5th, ", ", p97.5th, ")")) %>%
#       select(contains("Model")))
# }
# 
# #---- HR table ----
# load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
#             "model_results/", "coxph_results.RData"))
# HR_forplot <- tibble()
# HR_table <- tibble(race = race_labels)
# for (m in 1:3){
#   HR_forplot %<>% rbind(bind_rows(coxph_results[[m]]) %>% 
#                           filter(term %in% "apoe_y") %>%
#                           mutate(model = m) %>%
#                           select(model, race, HR, p2.5th, p97.5th))
#   HR_table %<>% cbind(
#     bind_rows(coxph_results[[m]]) %>% filter(term %in% "apoe_y") %>%
#       mutate_if(is.numeric, sprintf, fmt = '%#.2f') %>%
#       mutate(!!paste0("Model ",m, " HR (95% CI)") := 
#                paste0(HR, " (", p2.5th, ", ", p97.5th, ")")) %>%
#       select(contains("Model")))
# }
# 
# 
# 
# HR_table
# HR_table_y13
# HR_table_y13_y18 <- HR_table_y13 %>% 
#   cbind(HR_table %>% select(-race))
# HR_table_y13_y18 <-
#   rbind(setNames(as.data.frame(c("", rep("Follow up 13 years", ncol(HR_table) - 1), 
#                                  rep("Follow up 18 years", ncol(HR_table) - 1)) %>% t()),
#                  names(HR_table_y13_y18)),
#         HR_table_y13_y18)
# 
# writexl::write_xlsx(HR_table_y13_y18,
#                     here::here("output", "tables", "coxph_exclue2e4.xlsx"))
# 
# #---- HR plot ----
# HR_forplot_y13_y18 <- HR_forplot_y13 %>% 
#   mutate(fu_yrs = "Follow up 13 years") %>%
#   bind_rows(HR_forplot %>% mutate(fu_yrs = "Follow up 18 years"))
# 
# HR_forplot_y13_y18 %>%
#   filter(race != "Overall") %>%
#   mutate(race = factor(race, levels = race_labels),
#          # model = factor(model),
#          group_labels = case_when(
#            race %in% c("White", "Asian") ~ "All",
#            TRUE ~ "Asian ethnic groups")) %>%
#   ggplot(aes(x = race, y = HR, color = race, shape = as.factor(model))) +
#   geom_errorbar(aes(ymin = p2.5th, ymax = p97.5th), width = 0.2,
#                 position = position_dodge((width = 0.3))) +
#   geom_point(position = position_dodge(width = 0.3), fill = "white") +
#   scale_shape_manual(values = c(21, 24, 22)) +
#   scale_y_continuous(breaks = seq(0, 5.5, 0.5), position = "right") +
#   facet_grid(fu_yrs ~ group_labels, scales = "free", space = "free",
#              switch = "y") +
#   theme_bw() +
#   geom_hline(aes(yintercept = 1), linetype = "dashed") +
#   labs(x = element_blank(), y = "HR",
#        shape = "Models")+
#   guides(color = "none")
# # theme(legend.position = "none")
# 
# ggsave(file = here::here("output", "figures", "efigure3_coxph_exclue2e4_y1318.png"), 
#        device = "png", width = 7, height = 5, units = "in", dpi = 300)
# 
# 
# ##---- year specific HR plot ----
# HR_forplot %>%
#   filter(race != "Overall") %>%
#   mutate(race = factor(race, levels = race_labels),
#          # model = factor(model),
#          group_labels = case_when(
#            race %in% c("White", "Asian") ~ "All",
#            TRUE ~ "Asian ethnic groups")) %>%
#   ggplot(aes(x = race, y = HR, color = race, shape = as.factor(model))) +
#   geom_errorbar(aes(ymin = p2.5th, ymax = p97.5th), width = 0.2,
#                 position = position_dodge((width = 0.3))) +
#   geom_point(position = position_dodge((width = 0.3))) +
#   scale_y_continuous(breaks = seq(0, 4.2, 0.5)) +
#   facet_grid(. ~group_labels, scales = "free", space = "free") +
#   theme_bw() +
#   theme(legend.position = "none")
# 
# ggsave(file = here::here("output", "figures", "coxph_exclue2e4.png"), 
#        device = "png", width = 7, height = 5, units = "in", dpi = 300)
# 
# HR_forplot_y13 %>%
#   filter(race != "Overall") %>%
#   mutate(race = factor(race, levels = race_labels),
#          # model = factor(model),
#          group_labels = case_when(
#            race %in% c("White", "Asian") ~ "All",
#            TRUE ~ "Asian ethnic groups")) %>%
#   ggplot(aes(x = race, y = HR, color = race, shape = as.factor(model))) +
#   geom_errorbar(aes(ymin = p2.5th, ymax = p97.5th), width = 0.2,
#                 position = position_dodge((width = 0.3))) +
#   geom_point(position = position_dodge((width = 0.3))) +
#   scale_y_continuous(breaks = seq(0, 4, 0.5)) +
#   facet_grid(. ~group_labels, scales = "free", space = "free") +
#   theme_bw() +
#   theme(legend.position = "none")
# ggsave(file = here::here("output", "figures", "coxph_exclue2e4_y13.png"), 
#        device = "png", width = 7, height = 5, units = "in", dpi = 300)
