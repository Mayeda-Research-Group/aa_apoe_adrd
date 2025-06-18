# Descriptive statistics for baseline and cumulative incidence
# Created by: Yingyan Wu
# use nonimputed datasets

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")) {
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}
p_load("here", "tidyverse", "dplyr","formattable", "writexl", "readxl", "gridExtra",
       "magrittr", "survival")

#No scientific notation
options(scipen = 999)

# Paths
source(here::here("scripts", "0.paths.R"))

#----load the dataset ----
load(paste0(path_to_box, "Asian_Americans_dementia_data/aa_apoe_dementia/",
            "analysis_data/aa_apoe_tte_selected_e4all.RData"))
aa_apoe_tte_selected %<>% 
  mutate(apoe_y = factor(apoe_y, levels = c(1, 0)),
         ethnicity_rev = factor(ethnicity_rev,
                                levels = c(9, 2, 3, 5),
                                labels = c("Non-Latino White", 
                                           "Chinese", "Japanese", "Filipino")))

#----Density plot of baseline age by race/ethnicity APOE----
aa_apoe_tte_selected %>%
  ggplot(aes(x = survey_age, group = apoe_y, color = apoe_y)) +
  geom_density() +
  facet_grid(~ethnicity_rev, scales = "free_y") +
  scale_color_manual(
    name = NULL,
    values = c("#3271AE", "#E5A84B"),
    labels = c(expression(paste(italic("APOE-"), italic(epsilon), italic("4"), " carriers")),
               "Non-carriers")) +
  labs(x = "Age", y = "Frequency") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave(here::here("output", "figures", "fig_baseline_age_ethn_e4all.png"),
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

#---- count of deaths after dementia ----
aa_apoe_tte_selected %<>%
  mutate(death_after_dem_flag = 
           case_when(main_dem_v1_end_dem_flag == 1 & death_flag == 1 ~ 1,
                     main_dem_v1_end_dem_flag == 1 & death_flag == 0 ~ 0)) 


death_after_dem_tib <- aa_apoe_tte_selected %>%
  group_by(ethnicity_rev, apoe_y) %>%
  dplyr::summarize(n = n(),
                   n_dem = sum(main_dem_v1_end_dem_flag),
                   prop_dem = n_dem/n*100,
                   n_death_after_dem = sum(death_after_dem_flag, na.rm = T),
                   prop_death_after_dem = n_death_after_dem/n*100) %>%
  mutate(n_dem_prop = paste0(n_dem, " (", round(prop_dem, 1), "%)"),
         n_death_after_dem_prop = paste0(n_death_after_dem, " (", 
                                         round(prop_death_after_dem, 1), "%)")) %>%
  select(ethnicity_rev, apoe_y, n, n_dem_prop, n_death_after_dem_prop) %>%
  set_colnames(c("Race/ethnicity", "APOE-e4 carriers", "N", "Dementia cases (%)", 
                 "Death cases after dementia (%)"))

# write_xlsx(death_after_dem_tib, here::here("output", "tables",, 
#                                            "table_death_counts_after_dem.xlsx"))

#---- Cumulative incidence using Aalen-Johansen estimator (left truncation, competing risk) ----
aa_apoe_tte_selected %<>%
  mutate(event_multi =  factor(case_when(event_dem == 1 & event_death == 0 ~ 1,
                                         event_death == 1 ~ 2,
                                         TRUE ~ 0),
                               labels = c("s(0)", 
                                          "dementia", "death")),
         starttime = 0,
         main_dem_v1_end_time = main_dem_v1_end_age - survey_age)

#---- Follow up time as time scale ----
##---- Age adjusted multi state model ----
###---- Age weights ----
# p(apoe|age) by different ethnicity groups
pS_denom_model_ethn <- glm(apoe_y ~ survey_age * ethnicity_rev, data = aa_apoe_tte_selected, 
                           family = "binomial")
pS_num_model_ethn <- glm(apoe_y ~ ethnicity_rev, data = aa_apoe_tte_selected, 
                         family = "binomial")
aa_apoe_tte_selected %<>%
  mutate(
    pred_pS_denom_ethn = predict(pS_denom_model_ethn, type = "response", newdata = aa_apoe_tte_selected),
    pred_pS_num_ethn = predict(pS_num_model_ethn, type = "response", newdata = aa_apoe_tte_selected),
    pS_denom_ethn = ifelse(apoe_y == 1, pred_pS_denom_ethn, 1 - pred_pS_denom_ethn),
    pS_num_ethn = ifelse(apoe_y == 1, pred_pS_num_ethn, 1 - pred_pS_num_ethn),
    age_wt_ethn = pS_num_ethn/pS_denom_ethn) %>%
  group_by(ethnicity_rev) %>%
  mutate(age_wt_ethn_trunc = ifelse(age_wt_ethn > quantile(age_wt_ethn, 0.99),
                                    quantile(age_wt_ethn, 0.99), age_wt_ethn)) %>%
  ungroup()

# Sanity check
with(aa_apoe_tte_selected, tapply(age_wt_ethn, ethnicity_rev, summary))
with(aa_apoe_tte_selected, tapply(age_wt_ethn_trunc, ethnicity_rev, summary))

# p(apoe|age) pooled ethnicity groups
pS_denom_model <- glm(apoe_y ~ survey_age, data = aa_apoe_tte_selected, 
                      family = "binomial")
pS_num_model <- glm(apoe_y ~ 1, data = aa_apoe_tte_selected, 
                    family = "binomial")
aa_apoe_tte_selected %<>%
  mutate(
    pred_pS_denom = predict(pS_denom_model, type = "response", newdata = aa_apoe_tte_selected),
    pred_pS_num = predict(pS_num_model, type = "response", newdata = aa_apoe_tte_selected),
    pS_denom = ifelse(apoe_y == 1, pred_pS_denom, 1 - pred_pS_denom),
    pS_num = ifelse(apoe_y == 1, pred_pS_num, 1 - pred_pS_num),
    age_wt = pS_num/pS_denom,
    age_wt_trunc = ifelse(age_wt > quantile(age_wt, 0.99),
                          quantile(age_wt, 0.99), age_wt)) %>%
  ungroup()

# Sanity check
summary(aa_apoe_tte_selected$age_wt_trunc)
with(aa_apoe_tte_selected, tapply(age_wt, ethnicity_rev, summary))
with(aa_apoe_tte_selected, tapply(age_wt_trunc, ethnicity_rev, summary))

##---- Aalen-Johansen estimator ----
aa_apoe_tte_selected_aj <- aa_apoe_tte_selected %>%
  rbind(., aa_apoe_tte_selected %>% 
          filter(ethnicity_rev %in% c("Chinese", "Filipino", "Japanese")) %>%
          mutate(ethnicity_rev = "Asian American",
                 subjid = subjid*1000))

AJ_fit_age_adj <- 
  survfit(Surv(time = starttime, time2 = main_dem_v1_end_time, 
               event_multi, type = "mstate") ~ 
            strata(ethnicity_rev, apoe_y), data = aa_apoe_tte_selected_aj,
          weights = age_wt_trunc,
          cluster = subjid, id = subjid)

AJ_tib_age_adj <- broom::tidy(AJ_fit_age_adj) %>%
  mutate(strata = stringr::str_remove(strata, "strata\\(ethnicity_rev, apoe_y\\)="),
         ethnicity = factor(str_split_i(strata, ",", 1), 
                            levels = levels(aa_apoe_tte_selected_aj$ethnicity_rev)),
         apoe_y = factor(str_split_i(strata, ", ", 2),
                         levels = levels(aa_apoe_tte_selected$apoe_y)),
         state = case_when(state == "(s0)" ~ "Alive without dementia",
                           state == "death" ~ "Death without dementia",
                           state == "dementia" ~ "Dementia")) %>%
  mutate(ethnicity = factor(ethnicity, levels = c("Non-Latino White", "Asian American", 
                                                  "Chinese", "Japanese", "Filipino")),
         apoe_y = factor(apoe_y, levels = c(1, 0)),
         ethn_grp = ifelse(ethnicity %in% c("Non-Latino White", "Asian American"), "All", 
                           "Asian American ethnic groups"))

##---- 10 year dementia Cumulative incidence ----
color_palette <- c("#4B4B4B", "#B69833", "#7B5AA3", "#ED9DB2", "#54B663")
AJ_tib_age_adj %>% 
  mutate(across(c("estimate", "conf.low", "conf.high"), function(x) x*100)) %>%
  filter(state == "Dementia", time >= 10) %>%
  arrange(ethnicity, apoe_y, time) %>%
  group_by(ethnicity, apoe_y) %>%
  slice(1) %>%
  ggplot(aes(x = ethnicity, y = estimate, 
             fill = ethnicity, alpha = apoe_y)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(
    stat = "identity",
    aes(ymin = conf.low, ymax = conf.high), width = .2,
    position = position_dodge(.9)) + 
  # scale_y_continuous(breaks = seq(0, 30, by = 5)) +
  scale_alpha_discrete(range = c(1, 0.5)) +
  scale_fill_manual(values = color_palette) +
  #theme(axis.text.x = element_text(angle = 45,  vjust = 1, hjust = 1)) +
  theme_bw() + 
  theme(legend.position = "none",
        axis.text = element_text(size = 11),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 11)) +
  labs(
    # title = "Age-adjusted incidence rates by race/ethnicity and diabetes", 
    x = element_blank(), y = "Age-adjusted cumulative incidence at 10 years of follow-up") +
  # scale_fill_brewer(palette = "Paired") +
  facet_grid(. ~ ethn_grp, scales = "free", space = "free")

ggsave(here::here("output", "figures", "figure_age_adj_ci_dem_ethn_10yrs_e4all.png"),
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

##---- 10 year dementia and dementia-free death cumulative incidence table ----
AJ_tib_age_adj %>%
  filter(state != "Alive without dementia", ethnicity == "Filipino") %>%
  view()
# Filipino APOE carriers does not have time 13, jumped from 12.8 to 14. Use 12.8 estimate

AJ_age_dem_death_table <- AJ_tib_age_adj %>%
  mutate(across(c("estimate", "conf.low", "conf.high"), function(x) x*100),
         apoe_y = factor(apoe_y, levels = c(1, 0), labels = c("Carriers", "Non-carriers")),
         state = factor(state, levels = c("Alive without dementia", "Dementia", 
                                          "Death without dementia")),
         across(c(estimate, conf.high, conf.low), ~sprintf(., fmt = '%#.1f')),
         pe_ci = paste0(estimate, " \n(", conf.low, ", ", conf.high, ")"),
         yr_grp = case_when(
           ethnicity == "Filipino" & apoe_y == "Carriers" & time >= 12.79 ~ 13,
           TRUE ~ floor(time))) %>%
  filter(state != "Alive without dementia", yr_grp %in% c(5, 8, 10, 13)) %>%
  arrange(state, ethnicity, apoe_y, time) %>%
  group_by(state, ethnicity, apoe_y, yr_grp) %>%
  slice(1) %>%
  ungroup() %>%
  select(yr_grp, state, ethnicity, apoe_y, pe_ci) %>%
  pivot_wider(names_from = c(ethnicity, apoe_y),
              values_from = c(pe_ci),
              names_sep = "_") 

write_xlsx(AJ_age_dem_death_table, 
           here::here("output", "tables",
                      "table_age_adj_cum_inc_dem_death_581013_e4all.xlsx"))

##---- Cumulative incidence plot Dementia only ----
AJ_tib_age_adj %>% 
  mutate(across(c("estimate", "conf.low", "conf.high"), function(x) x*100)) %>%
  # cut the y axis at 50, top coded conf.high at 50
  mutate(conf.high_tc = case_when(conf.high >= 50 ~ 50,
                                  TRUE ~ conf.high)) %>%
  filter(state == "Dementia") %>%
  ggplot(aes(x = time, y = estimate, group = apoe_y)) +
  geom_line(aes(color = apoe_y)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high_tc, fill = apoe_y), alpha = 0.2) +
  facet_grid(cols = vars(ethnicity), scales = "free", space = "free") +
  scale_color_manual(
    name = NULL,
    values = c("#3271AE", "#E5A84B"),
    labels = c(expression(paste(italic("APOE-"), italic(epsilon), italic("4"), " carriers")),
               "Non-carriers")) +
  scale_fill_manual(
    name = NULL,
    values = c("#3271AE", "#E5A84B"),
    labels = c(expression(paste(italic("APOE-"), italic(epsilon), italic("4"), " carriers")),
               "Non-carriers")) +
  scale_y_continuous(limits = c(0, 50), breaks = seq(0, 50, 10)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Follow-up time (years)",
       y = "Age-adjusted cumulative incidence (%)") 

# for footnotes
AJ_tib_age_adj %>% 
  mutate(across(c("estimate", "conf.low", "conf.high"), function(x) x*100)) %>%
  # cut the y axis at 50, top coded conf.high at 50
  mutate(conf.high_tc = case_when(conf.high >= 50 ~ 50,
                                  TRUE ~ conf.high)) %>%
  filter(state == "Dementia", ethnicity == "Filipino", conf.high_tc < conf.high) %>%
  select(time, conf.high) %>%
  apply(., 2, summary)

ggsave(here::here("output", "figures", "fig_cuminc_AJ_dem_ethns_withCI_e4all.png"),
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

##---- Cumulative incidence plot dementia-free death only ----
AJ_tib_age_adj %>% 
  mutate(across(c("estimate", "conf.low", "conf.high"), function(x) x*100)) %>%
  # # cut the y axis at 50, top coded conf.high at 50
  # mutate(conf.high_tc = case_when(conf.high >= 50 ~ 50,
  #                                 TRUE ~ conf.high)) %>%
  filter(state == "Death without dementia") %>%
  ggplot(aes(x = time, y = estimate, group = apoe_y)) +
  geom_line(aes(color = apoe_y), linetype = "dashed") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = apoe_y), alpha = 0.2) +
  facet_grid(cols = vars(ethnicity), scales = "free", space = "free") +
  scale_color_manual(
    name = NULL,
    values = c("#3271AE", "#E5A84B"),
    labels = c(expression(paste(italic("APOE-"), italic(epsilon), italic("4"), " carriers")),
               "Non-carriers")) +
  scale_fill_manual(
    name = NULL,
    values = c("#3271AE", "#E5A84B"),
    labels = c(expression(paste(italic("APOE-"), italic(epsilon), italic("4"), " carriers")),
               "Non-carriers")) +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 10)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Follow-up time (years)",
       y = "Age-adjusted cumulative incidence (%)") 

ggsave(here::here("output", "figures", "fig_cuminc_AJ_death_ethns_withCI_e4all.png"),
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

#---- OLD ----
##---- Cumulative incidence plot dem and death ----
# crude_AJ_tib %>% 
#   filter(state %in% c("Dementia", "Death without dementia"), age <= 95) %>%
#   ggplot(aes(x = age, y = estimate, group = interaction(state, apoe_y))) + 
#   geom_line(aes(color = apoe_y, linetype = state)) +
#   facet_grid(cols = vars(ethnicity)) +
#   scale_color_manual(name = NULL,
#                      values = c("#3271AE", "#E5A84B"),
#                      labels = c("Women", "Men")) +
#   scale_linetype_manual(name = NULL,
#                         values = c(2, 1),
#                         labels = c("Death without dementia", "Dementia")) +
#   theme_bw() +
#   theme(legend.position = "bottom") +
#   labs(x = "Age (years)",
#        y = "Cumulative incidence (%)") 
# 
# ggsave(here::here("02_figs", "fig_cuminc_AJ_dem_death_ethns.png"),
#        device = "png", width = 7, height = 5, units = "in", dpi = 300)