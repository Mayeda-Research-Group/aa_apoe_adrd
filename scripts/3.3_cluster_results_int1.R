# Bootstrap results (from cluster)
# Created by: Yingyan Wu
# Nov.7.2023
# use nonimputed datasets

#---- Package loading + options ----
rm(list = ls())
if (!require("pacman")) {
  install.packages("pacman", repos = 'http://cran.us.r-project.org')
}
p_load("here", "tidyverse", "broom", "magrittr")

# No scientific notation
options(scipen = 999)

# Paths
source(here::here("scripts", "0.paths.R"))

#---- Load the results ----
race_varlabel <- c("overall", "chn", "jpn", "phl", "wht", "AA")
for (r in race_varlabel){
  temp <- tibble()
  for (b in 0:1000){
    temp <- bind_rows(
      temp, readRDS(paste0(path_to_box, 
                           "Asian_Americans_dementia_data/aa_apoe_dementia/", 
                           "model_results/hoffman2/plr_bootstrap/",
                           paste0(r, "_RD_RR_time_wide_", b, ".RDS"))))
  }
  assign(paste0(r, "_RD_RR_time_wide"), temp)
}
rm(temp)

# # warnings
# race_varlabel <- c("overall", "chn", "jpn", "phl", "wht", "AA")
# warnings <- tibble()
# for (r in race_varlabel){
#   for (b in 0:1000){
#     warnings <- bind_rows(
#       warnings, readRDS(paste0(path_to_box, 
#                                "Asian_Americans_dementia_data/aa_apoe_dementia/", 
#                                "model_results/hoffman2/warnings/", 
#                                paste0(r, "_warning_", b, ".RDS"))))
#   }
# }

# # Check the distribution of estimates with and without warnings
# warnings_filtered <- warnings %>%
#   filter(w_glm_1 == 0, w_glm_2 == 0, w_glm_3 == 0)
# 
# warnings_filtered %>%
#   group_by(race, model) %>%
#   summarise(n = n())
# 
# jpn_nowarnings <- warnings_filtered %>%
#   filter(race == "jpn") %>%
#   mutate(new_id = paste0(scenario, "_", i_boot)) %>%
#   pull(unique(new_id))
# 
# chn_nowarnings <- warnings_filtered %>%
#   filter(race == "chn") %>%
#   mutate(new_id = paste0(scenario, "_", i_boot)) %>%
#   pull(unique(new_id))
# 
# phl_nowarnings <- warnings_filtered %>%
#   filter(race == "phl") %>%
#   mutate(new_id = paste0(scenario, "_", i_boot)) %>%
#   pull(unique(new_id))
# 
# jpn_RD_RR_time_wide %>%
#   mutate(new_id = paste0(scenario_num, "_", i_boot),
#          flag = case_when(new_id %in% jpn_nowarnings ~ "no warnings",
#                           TRUE ~ "warnings")) %>%
#   ggplot(aes(x = RR_model_3_y10, color = flag)) +
#   geom_density()
# 
# chn_RD_RR_time_wide %>%
#   mutate(new_id = paste0(scenario_num, "_", i_boot),
#          flag = case_when(new_id %in% chn_nowarnings ~ "no warnings",
#                           TRUE ~ "warnings")) %>%
#   ggplot(aes(x = RR_model_3_y10, color = flag)) +
#   geom_density()
# 
# phl_RD_RR_time_wide %>%
#   mutate(new_id = paste0(scenario_num, "_", i_boot),
#          flag = case_when(new_id %in% phl_nowarnings ~ "no warnings",
#                           TRUE ~ "warnings")) %>%
#   ggplot(aes(x = RR_model_3_y10, color = flag)) +
#   geom_density()

#---- PE, 95% CI ----
race_labels <- c("Overall", "Chinese", "Japanese", "Filipino", "Non-Latino White", "Asian American")
RD_RR_time_CI <- tibble()
for (i in 1:length(race_varlabel)){
  temp <- get(paste0(race_varlabel[[i]], "_RD_RR_time_wide")) %>% 
    filter(i_boot != 0) %>%
    select(-race, -contains("int_RD"), -i_boot, -scenario_num) %>%
    reframe(across(everything(), ~quantile(.x, c(0.025, 0.975)))) %>%
    mutate(estimate = c("p2.5th", "p97.5th"),
           race = race_labels[[i]]) %>%
    select(estimate, race, everything()) %>% 
    pivot_longer(cols = -c(estimate, race), 
                 names_to = c(".value", "fu_yr"),
                 names_pattern = "(.*)_y(.*)") %>%
    mutate(fu_yr = as.numeric(fu_yr)) %>%
    arrange(fu_yr, estimate)
  temp2 <- get(paste0(race_varlabel[[i]], "_RD_RR_time_wide")) %>% 
    filter(i_boot == 0) %>%
    select(-contains("int_RD"), -i_boot, -scenario_num) %>%
    mutate(estimate = "pe",
           race = race_labels[[i]]) %>%
    select(estimate, race, everything()) %>%
    pivot_longer(cols = -c(estimate, race), 
                 names_to = c(".value", "fu_yr"),
                 names_pattern = "(.*)_y(.*)") %>%
    mutate(fu_yr = as.numeric(fu_yr)) %>%
    arrange(fu_yr, estimate)
  RD_RR_time_CI <- bind_rows(RD_RR_time_CI, temp2, temp)
}

remove(temp, temp2)
save(RD_RR_time_CI, file = paste0(path_to_box, 
                                  "Asian_Americans_dementia_data/aa_apoe_dementia/", 
                                  "model_results/hoffman2/RD_RR_time_CI.RData"))

#---- Estimates ----
load(paste0(path_to_box, 
            "Asian_Americans_dementia_data/aa_apoe_dementia/", 
            "model_results/hoffman2/RD_RR_time_CI.RData"))
for (fy in c(10, 13, 18)){
  temp <- RD_RR_time_CI %>%
    mutate(race = factor(race, levels = race_labels)) %>%
    filter(fu_yr == fy) %>%
    mutate_if(is.numeric, sprintf, fmt = '%#.2f') %>%
    pivot_longer(cols = -c(estimate, race, fu_yr),
                 names_to = c(".value", "model"),
                 names_pattern = "(.*)_model_(.*)") %>%
    select(race, estimate, model, RD_perc, RR) %>%
    pivot_wider(names_from = estimate, values_from = `RD_perc`:`RR`) %>%
    mutate(
      !!sym(paste0("RD (95% CI) at year ", fy, "(percentage)")) := 
        paste0(`RD_perc_pe`, " (", `RD_perc_p2.5th`, ", ", 
               `RD_perc_p97.5th`, ")"),
      !!sym(paste0("RR (95% CI) at year ", fy)) :=
        paste0(`RR_pe`, " (", `RR_p2.5th`, ", ", 
               `RR_p97.5th`, ")"),) %>%
    select(race, model, contains("(95% CI)"))
  assign(paste0("RD_RR_", fy, "_tib"), temp)
}
RD_RR_101318_tib <- plyr::join_all(list(RD_RR_10_tib, RD_RR_13_tib, RD_RR_18_tib), 
                                   type = "left",
                                   by = c("race", "model"))

##---- RD and RR 13, 18 years of follow up  table ----
race_fac_labels <- c("Non-Latino White", "Asian American", "Chinese", "Japanese", "Filipino")
RD_RR_101318_wide <- RD_RR_101318_tib %>%
  filter(race != "Overall") %>%
  pivot_wider(names_from = model, values_from = contains("(95% CI)")) %>%
  mutate(race = factor(race, levels = race_fac_labels),
         race_grp = ifelse(race %in% c("Non-Latino White", "Asian American"), "All", 
                           "Asian American ethnic groups")) %>%
  select(race_grp, race, everything()) %>% 
  arrange(race_grp, race)

RD_RR_101318_wide %<>%
  rbind(setNames(
    as.data.frame(
      c("", "", rep("Follow up 10 years", (ncol(RD_RR_101318_wide) - 1)/3),
        rep("Follow up 13 years", (ncol(RD_RR_101318_wide) - 1)/3), 
        rep("Follow up 18 years", (ncol(RD_RR_101318_wide) - 1)/3)) %>% t()),
    names(RD_RR_101318_wide)),
    .)

writexl::write_xlsx(RD_RR_101318_wide,
                    here::here("output", "tables", "RD_RR_int1_exclue2e4_101318.xlsx"))

##---- RR and RD whisker plots ----
RD_RR_time_CI_model_long <- RD_RR_time_CI %>%
  filter(race != "Overall") %>%
  pivot_longer(cols = -c(estimate, race, fu_yr),
               names_to = c(".value", "model"),
               names_pattern = "(.*)_model_(.*)")

RD_RR_time_forerrorbarplot <- RD_RR_time_CI_model_long %>%
  select(estimate, race, model, fu_yr, RD_perc, RR) %>%
  mutate(race = factor(race, levels = race_fac_labels),
         race_facet = ifelse(race %in% c("Non-Latino White", "Asian American"), "All", 
                             "Asian American ethnic groups")) %>%
  pivot_wider(names_from = estimate, values_from = `RD_perc`:`RR`)

color_palette <- c("#4B4B4B", "#B69833", "#7B5AA3", "#ED9DB2", "#54B663")
###---- Model 1 at 10 years of follow up----
hline_tib <- tibble(scale = c("RR", "RD_perc"), ref = c(1, 0))
RD_RR_time_forerrorbarplot %>%
  filter(model == 1 & !race %in% c("Overall") & fu_yr == 10) %>%
  pivot_longer(RD_perc_pe:RR_p97.5th, 
               names_to = c("scale", ".value"),
               names_pattern = "(.*)_(.*)") %>%
  ggplot(aes(x = race, y = pe, color = race, shape = as.factor(fu_yr))) +
  geom_point(position = position_dodge((width = 0.3))) +
  geom_errorbar(aes(ymin = p2.5th, ymax = p97.5th), width = 0.2,
                position = position_dodge((width = 0.3))) +
  facet_grid(rows = vars(scale), # vars(fct_rev(scale)),
             cols = vars(race_facet),
             scales = "free", space = "free_x", switch = "y",
             labeller = as_labeller(
               c(`RD_perc` = "Risk Difference % (95% CI)", 
                 `RR` = "Risk Ratio (95% CI)",
                 `All` = "All",
                 `Asian American ethnic groups` = "Asian American ethnic groups"))) +
  scale_y_continuous(position = "right") +
  scale_color_manual(values = color_palette) +
  geom_hline(data = hline_tib, aes(yintercept = ref), linetype = "dashed") +
  theme_bw() +
  labs(x = element_blank(), y = element_blank(),
       shape = "Follow up years") +
  # title = expression("Risk ratio and risk difference(%) at 10 years of follow up for the association between "*italic("APOE-"*epsilon*"4")*" and dementia")) +
  theme(legend.position = "none",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 11))

ggsave(file = here::here("output", "figures", "RD_RR_model1_exclue2e4_10yrs.png"), 
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

###---- model 1 at 10, 13 and 18 years of follow up ----
# hline_tib <- tibble(scale = c("Risk Ratio", "Risk Difference (%)"), ref = c(1, 0))
# RD_RR_time_forerrorbarplot %>%
#   filter(model == 1 & !race %in% c("Overall") & fu_yr %in% c(10, 13, 18)) %>%
#   pivot_longer(RD_perc_pe:RR_p97.5th, 
#                names_to = c("scale", ".value"),
#                names_pattern = "(.*)_(.*)") %>%
#   mutate(scale = case_when(scale == "RD_perc" ~ "Risk Difference (%)",
#                            scale == "RR" ~ "Risk Ratio")) %>%
#   ggplot(aes(x = race, y = pe, color = race, shape = as.factor(fu_yr))) +
#   geom_point(position = position_dodge((width = 0.3))) +
#   geom_errorbar(aes(ymin = p2.5th, ymax = p97.5th), width = 0.2,
#                 position = position_dodge((width = 0.3))) +
#   facet_grid(rows = vars(fct_rev(scale)),
#              cols = vars(race_facet),
#              scales = "free", space = "free_x", switch = "y") +
#   scale_y_continuous(position = "right") +
#   geom_hline(data = hline_tib, aes(yintercept = ref), linetype = "dashed") +
#   theme_bw() +
#   labs(x = element_blank(), y = element_blank(),
#        shape = "Follow up years") +
#   guides(color = "none")
# 
# ggsave(file = here::here("output", "figures", "RD_RR_model1_exclue2e4_all_fuyrs_left.png"), 
#        device = "png", width = 7, height = 5, units = "in", dpi = 300)

###---- Model 1-3 at 10 years of follow up ----
hline_tib <- tibble(scale = c("RR", "RD_perc"), ref = c(1, 0))
RD_RR_time_forerrorbarplot %>%
  filter(fu_yr == 10 & !race %in% c("Overall")) %>%
  pivot_longer(RD_perc_pe:RR_p97.5th, 
               names_to = c("scale", ".value"),
               names_pattern = "(.*)_(.*)") %>%
  ggplot(aes(x = race, y = pe, color = race, shape = as.factor(model))) +
  geom_errorbar(aes(ymin = p2.5th, ymax = p97.5th), width = 0.2,
                position = position_dodge((width = 0.3))) +
  geom_point(position = position_dodge((width = 0.3)), fill = "white") +
  facet_grid(rows = vars(scale), # vars(fct_rev(scale)),
             cols = vars(race_facet),
             scales = "free", space = "free_x", switch = "y",
             labeller = as_labeller(
               c(`RD_perc` = "Risk Difference % (95% CI)", 
                 `RR` = "Risk Ratio (95% CI)",
                 `All` = "All",
                 `Asian American ethnic groups` = "Asian American ethnic groups"))) +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_color_manual(values = color_palette) +
  scale_y_continuous(position = "right") +
  geom_hline(data = hline_tib, aes(yintercept = ref), linetype = "dashed") +
  theme_bw() +
  labs(x = element_blank(), y = element_blank(),
       shape = "Models")+
  # title = expression("Risk Difference (%) and Risk Ratio at 13 years of follow up for the association between "*italic("APOE-"*epsilon*"4")*" and dementia")) +
  guides(color = "none") +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 11),
        strip.text = element_text(size = 11))

ggsave(file = here::here("output", "figures", "RD_RR_10fuyrs_exclue2e4_all_models.png"), 
       device = "png", width = 7, height = 5, units = "in", dpi = 300)

#---- Plot curves over time ----
##---- risk over time ----
# Model 1
RD_RR_time_CI_model_long %>% 
  filter(!race %in% c("Overall")) %>%
  mutate(race = factor(race, levels = race_fac_labels),
         model = paste0("Model ", model)) %>%
  select(model, race, fu_yr, estimate, mean_cif1_perc, mean_cif0_perc) %>%
  pivot_longer(mean_cif1_perc:mean_cif0_perc, names_to = "apoe_y", 
               values_to = "mean_cif_perc", names_pattern = "mean_cif(.)_perc") %>%
  pivot_wider(names_from = estimate, values_from = mean_cif_perc, names_prefix = "mean_cif_") %>%
  mutate(apoe_y = factor(apoe_y, levels = c(1, 0))) %>%
  ggplot(aes(x = fu_yr, y = mean_cif_pe, group = apoe_y)) + 
  geom_line(aes(color = race, linetype = apoe_y)) + 
  geom_ribbon(aes(ymin = mean_cif_p2.5th, ymax = mean_cif_p97.5th,
                  fill = race, alpha = apoe_y)) +
  facet_grid(cols = vars(race), rows = vars(model), 
             scales = "free", space = "free_x") +
  scale_alpha_discrete(range = c(0.4, 0.2)) +
  guides(color = "none", fill = "none", alpha = "none") +
  scale_color_manual(values = color_palette) +
  scale_fill_manual(values = color_palette) +
  scale_linetype_manual(
    name = NULL,
    values = c("solid", "dashed"),
    labels = c(expression(paste(italic("APOE-"), italic(epsilon), italic("4"), " carriers")),
               "Non-carriers")) +
  scale_y_continuous(limits = c(0, 65), breaks = seq(0, 60, 10)) +
  scale_x_continuous(limits = c(0, 18),
                     breaks = c(0, 5, 10, 15, 18)) +
  theme_bw() +
  theme(legend.position = "bottom",
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15)) +
  labs(x = "Follow-up time (years)",
       y = "Cumulative incidence (%)") 

ggsave(file = here::here("output", "figures", "adjusted_cumulative_incidence_curves_exclue2e4.png"), 
       device = "png", width = 12, height = 9, units = "in", dpi = 300)

#---- OLD ----
##---- Adjusted survival curves ----
# # Model 1
# RD_RR_time_CI_model_long %>% 
#   filter(model == 1 & !race %in% c("Overall")) %>%
#   mutate(race = factor(race, levels = race_fac_labels)) %>%
#   select(race, fu_yr, model, estimate, mean_cums1, mean_cums0) %>%
#   pivot_longer(mean_cums1:mean_cums0, names_to = "APOE4 e4", values_to = "mean_cums", 
#                names_prefix = "mean_cums") %>%
#   pivot_wider(names_from = estimate, values_from = mean_cums, names_prefix = "mean_cums_") %>%
#   mutate(`APOE4 e4` = ifelse(`APOE4 e4` == 1, "Yes", "No"),
#          `APOE4 e4` = factor(`APOE4 e4`, levels = c("Yes", "No"))) %>%
#   ggplot(aes(x = fu_yr)) + 
#   geom_line(aes(y = mean_cums_pe, colour = race, linetype = `APOE4 e4`)) + 
#   geom_ribbon(aes(ymin = mean_cums_p2.5th, ymax = mean_cums_p97.5th,
#                   group = `APOE4 e4`, fill = race), alpha = 0.2) +
#   labs(
#     title = "Adjusted survival curves",
#     x = "Years",
#     y = "Survival probability",
#     colour = "Race/Ethnicity",
#     fill = "Race/Ethnicity") + 
#   scale_x_continuous(limits = c(0, 18), 
#                      breaks = seq(0, 18, 2)) +
#   facet_wrap(~race) +
#   # scale_y_continuous(limits = c(0, 1),  breaks=seq(0, 1, 0.2)) +
#   theme_bw()
# # theme(legend.position="bottom")
# 
# ggsave(file = here::here("output", "figures", "adjusted_survival_curves_exclue2e4.png"), 
#        device = "png", width = 7, height = 5, units = "in", dpi = 300)

##---- RD over time ----
# # Model 1
# RD_RR_time_CI_model_long %>% 
#   filter(model == 1 & !race %in% c("Overall")) %>%
#   mutate(race = factor(race, levels = race_fac_labels)) %>%
#   select(race, fu_yr, model, estimate, RD_perc) %>%
#   pivot_wider(names_from = estimate, values_from = RD_perc, 
#               names_prefix = "RD_perc_") %>%
#   ggplot(aes(x = fu_yr)) + 
#   geom_line(aes(y = RD_perc_pe, color = race, group = model)) +
#   geom_ribbon(aes(ymin = RD_perc_p2.5th, ymax = RD_perc_p97.5th,
#                   fill = race), alpha = 0.2) +
#   facet_wrap(~race) +
#   labs(title = "Risk difference (in percentage) over time",
#        x = "Years",
#        y = "RD (%)",
#        color = "Race/Ethnicity",
#        fill = "Race/Ethnicity") + 
#   scale_x_continuous(limits = c(0, 18), 
#                      breaks = seq(0, 18, 2)) +
#   theme_bw()
# 
# ggsave(file = here::here("output", "figures", "RD_curves_exclue2e4.png"), 
#        device = "png", width = 7, height = 5, units = "in", dpi = 300)

##---- RR over time ----
# RD_RR_time_CI_model_long %>% 
#   filter(model == 1 & !race %in% c("Overall")) %>%
#   mutate(race = factor(race, levels = race_fac_labels)) %>%
#   select(race, fu_yr, model, estimate, RR) %>%
#   pivot_wider(names_from = estimate, values_from = RR, 
#               names_prefix = "RR_") %>%
#   ggplot(aes(x = fu_yr)) + 
#   geom_line(aes(y = RR_pe, color = race, group = model)) +
#   geom_ribbon(aes(ymin = RR_p2.5th, ymax = RR_p97.5th,
#                   fill = race), alpha = 0.2) +
#   facet_wrap(~race, scale = "free") +
#   labs(title = "Risk Ratio over time",
#        x = "Years",
#        y = "RR",
#        color = "Race/Ethnicity",
#        fill = "Race/Ethnicity") + 
#   scale_x_continuous(limits = c(0, 18), 
#                      breaks = seq(0, 18, 2)) +
#   theme_bw()
# ggsave(file = here::here("output", "figures", "RR_curves_exclue2e4.png"), 
#        device = "png", width = 7, height = 5, units = "in", dpi = 300)


##---- RR and RD whisker plots ----
# RD_RR_time_model1_forerrorbarplot <- RD_RR_time_model1 %>%
#   select(estimate, race, fu_yr, RD_perc, RR) %>%
#   mutate(race_facet = case_when(
#     race %in% c("Overall", "Asian", "White") ~ "All",
#     TRUE ~ "Asian ethnic groups"),
#     race = factor(race, levels = race_fac_labels) %>%
#   pivot_wider(names_from = estimate, values_from = `RD_perc`:`RR`)
# 
# RD_RR_time_model1_forerrorbarplot %>%
#   filter(fu_yr == 13 & race != "Overall") %>%
#   ggplot(aes(x = race, y = RR_pe, color = race)) +
#   geom_point(position = position_dodge((width = 0.3))) +
#   geom_errorbar(aes(ymin = RR_p2.5th, ymax = RR_p97.5th), width = 0.2,
#                 position = position_dodge()) +
#   facet_grid(cols = vars(race_facet), scales = "free_x", space="free") +
#   geom_hline(yintercept = 1) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(x = "Race/Ethnicity", y = "Risk Ratio at 13 years of follow up")
# 
# ggsave(file = here::here("output", "figures", "RR_13fuyrs_exclue2e4.png"), 
#        device = "png", width = 7, height = 5, units = "in", dpi = 300)
# 
# RD_RR_time_model1_forerrorbarplot %>%
#   filter(fu_yr == 13 & race != "Overall") %>%
#   ggplot(aes(x = race, y = RD_perc_pe, color = race)) +
#   geom_point(position = position_dodge((width = 0.3))) +
#   geom_errorbar(aes(ymin = RD_perc_p2.5th, ymax = RD_perc_p97.5th), width = 0.2,
#                 position = position_dodge()) +
#   facet_grid(cols = vars(race_facet), scales = "free_x", space = "free") +
#   geom_hline(yintercept = 0) +
#   theme_bw() +
#   theme(legend.position = "none") +
#   labs(x = "Race/Ethnicity", y = "Risk Difference (%) at 13 years of follow up")
# 
# ggsave(file = here::here("output", "figures", "RD_perc_13fuyrs_exclue2e4.png"), 
#        device = "png", width = 7, height = 5, units = "in", dpi = 300)
# 
# hline_tib <- tibble(scale = c("Risk Ratio", "Risk Difference (%)"), ref = c(1, 0))
# RD_RR_time_model1_forerrorbarplot %>%
#   filter(fu_yr == 13 & !race %in% c("Overall")) %>%
#   pivot_longer(RD_perc_pe:RR_p97.5th, 
#                names_to = c("scale", ".value"),
#                names_pattern = "(.*)_(.*)") %>%
#   mutate(scale = case_when(scale == "RD_perc" ~ "Risk Difference (%)",
#                            scale == "RR" ~ "Risk Ratio")) %>%
#   ggplot(aes(x = race, y = pe, color = race)) +
#   geom_point(position = position_dodge((width = 0.3))) +
#   geom_errorbar(aes(ymin = p2.5th, ymax = p97.5th), width = 0.2,
#                 position = position_dodge()) +
#   facet_grid(rows = vars(fct_rev(scale)),
#              cols = vars(race_facet),
#              scales = "free", space = "free_x", switch = "y") +
#   scale_y_continuous(position = "right") +
#   geom_hline(data = hline_tib, aes(yintercept = ref), linetype = "dashed") +
#   theme_bw() +
#   labs(x = element_blank(), y = element_blank(), 
#        title = expression("Risk Difference (%) and Risk Ratio at 13 years of follow up for the association between "*italic("APOE-"*epsilon*"4")*" and dementia")) +
#   theme(legend.position = "none",
#         plot.title = element_text(hjust = 0.5, size = 10, face = "bold"))
# 
# ggsave(file = here::here("output", "figures", "RD_RR_13fuyrs_exclue2e4_title.jpg"), 
#        device = "jpg", width = 7, height = 5, units = "in", dpi = 300)
