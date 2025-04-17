# IPCW for PLR (use baseline characteristics)
# Yingyan Wu
# Mar 19 2025
# use imputed dataset

#---- Package loading + options ----
rm(list = ls())

if (!require("pacman")) {
  install.packages("pacman", repos='http://cran.us.r-project.org')
}

p_load("haven", "tidyverse", "magrittr", "dplyr", "cobalt", "survival", "splines")

options(scipen = 999, digits = 8)

# Paths
source(here::here("scripts", "0.paths.R"))

#---- Load the data ----
aa_apoe_long_tte_selected_imp <- 
  readRDS(paste0(path_to_box, 
                 "Asian_Americans_dementia_data/aa_apoe_dementia/",
                 "analysis_data/imp_pmm_20_stacked_long.RDS")) %>%
  mutate(event_cens = case_when(event_endmem == 1|event_death == 1 ~ 1,
                                event_endmem == 0 & event_death == 0 ~ 0)) %>%
  # for stratified CIF
  mutate(
    survey_age_cat = cut(survey_age_r, c(60, 70, 80, 90)), 
    survey_age_cat = str_sub(survey_age_cat, 2, 6) %>% str_replace(",", "-"),
    sr_bmi_cat_25 = ifelse(sr_bmi < 25, "<25", ">25"), 
    sr_bmi_cat_30 = ifelse(sr_bmi < 30, "<30", ">30"),
    income_pp_cat = cut(income_pp, c(0, 20000, 40000, 60000, 80000, 100000, Inf),
                        dig.lab = 6),
    income_pp_cat = str_replace(income_pp_cat, "\\(", "") %>% 
      str_replace("\\]", "") %>% str_replace(",", "-"),
    ehr_ht_median_cat = cut(ehr_ht_median, c(57, 64, 70, 76)),
    ehr_ht_median_cat = str_sub(ehr_ht_median_cat, 2, 6) %>% str_replace(",", "-"))

num_imp = 20
cov_bal_vars <- c("apoe_y", "female", "survey_age", "usaborn_rev", "edu_3", 
                  "income_pp", "generalhealth_3", "smoking_2", "prev_htn_flag",
                  "prev_combined_stroke_flag", "prev_diab_flag", "sr_bmi", 
                  "sr_depress", "alcohol_3", "pa_3", "ehr_ht_median")

# # Check the class of the variables
aa_apoe_long_tte_selected_imp %>%
  select(all_of(cov_bal_vars)) %>%
  map_dfr(class) %>%
  t()

with(aa_apoe_long_tte_selected_imp1, table(fu_yr, ethnicity_rev, event_death, useNA = "ifany"))

#---- CIF for events stratified ----
plot_km <- function(data = data_for_cif, event = "death", stratify_by, 
                    weight = F, weightvar = NULL) {
  # # Test
  # data <- data_for_cif
  # stratify_by <- "education_4"
  
  formula <- paste0("Surv(time = start, time2 = fu_yr, event = event_", event, ") ~ ", 
                    stratify_by, " + ethnicity_rev")
  if (weight == FALSE){
    km_object <- survfit(
      as.formula(formula), data = data, cluster = subjid
    )
  } else if (weight == TRUE){
    km_object <- survfit(
      as.formula(formula), data = data, cluster = subjid, weights = get(weightvar),
    )
  }
  
  
  p <- km_object %>% 
    broom::tidy() %>% 
    mutate(
      stratum = str_split_i(strata, ",", i = 1) %>% str_split_i("[^<]=", i = 2),
      # stroke = ifelse(str_detect(strata, "1"), "Yes", "No"), 
      ethn = str_split_i(strata, ",", i = 2) %>% str_split_i("[^<]=", i = 2) %>% str_trim(),
      cif = 1 - estimate, 
      cif.conf.low = 1 - conf.low, 
      cif.conf.high = 1 - conf.high, 
    ) %>% 
    ggplot(aes(time, cif)) +
    geom_line(aes(color = stratum), linewidth = 0.7) +
    geom_ribbon(aes(ymin = cif.conf.low, ymax = cif.conf.high, group = stratum), alpha = 0.1) +
    scale_x_continuous(
      # limits = c(0, 17),
      breaks = seq(0, 18, 3), minor_breaks = 0:18) +
    labs(y = paste0("Cumulative Incidence for ", event), color = stratify_by) +
    # labs(title = title_text) + 
    theme_bw() + 
    facet_wrap(~ethn, labeller = as_labeller(c("2" = "Chinese", 
                                               "3" = "Japanese", 
                                               "5" = "Filipino", 
                                               "9" = "White")))
  
  # ggsave(here::here("output", "figures", "crude_km_stratified", 
  #                   paste0(stratify_by, "_", event, ".png")),
  #        device = "png", width = 7, height = 5, units = "in", dpi = 300,
  #        plot = p)
  
  return(p)
}

##---- Stratified by baseline covariates ----
data_for_cif <- aa_apoe_long_tte_selected_imp %>% filter(imp == 1) 
# # Sanity check
# table(data_for_cif$income_pp_cat, useNA = "ifany")
# table(data_for_cif$ehr_ht_median_cat, useNA = "ifany")

# Death
plot_km(stratify_by = "survey_age_cat")
plot_km(stratify_by = "female")
plot_km(stratify_by = "smoking_2")
plot_km(stratify_by = "generalhealth_3")
plot_km(stratify_by = "income_pp_cat")
plot_km(stratify_by = "prev_htn_flag")
plot_km(stratify_by = "prev_combined_stroke_flag")
plot_km(stratify_by = "prev_diab_flag")
plot_km(stratify_by = "usaborn_rev")
plot_km(stratify_by = "edu_3")
plot_km(stratify_by = "pa_3")

# # Not different for all race/ethnicity
# plot_km(stratify_by = "sr_bmi_cat_25")
# plot_km(stratify_by = "sr_bmi_cat_30")
# plot_km(stratify_by = "ehr_ht_median_cat")
# plot_km(stratify_by = "sr_depress")
# plot_km(stratify_by = "alcohol_3")

# End membership
plot_km(event = "endmem", stratify_by = "survey_age_cat")
plot_km(event = "endmem", stratify_by = "female")
plot_km(event = "endmem", stratify_by = "generalhealth_3")
plot_km(event = "endmem", stratify_by = "prev_htn_flag")
plot_km(event = "endmem", stratify_by = "prev_combined_stroke_flag")
plot_km(event = "endmem", stratify_by = "prev_diab_flag")
plot_km(event = "endmem", stratify_by = "usaborn_rev")
plot_km(event = "endmem", stratify_by = "edu_3")
plot_km(event = "endmem", stratify_by = "sr_depress")
plot_km(event = "endmem", stratify_by = "alcohol_3")
plot_km(event = "endmem", stratify_by = "sr_bmi_cat_25")
plot_km(event = "endmem", stratify_by = "sr_bmi_cat_30")

# # Not different for all race/ethnicity
# plot_km(event = "endmem", stratify_by = "smoking_2")
# plot_km(event = "endmem", stratify_by = "pa_3")
# plot_km(event = "endmem", stratify_by = "income_pp_cat")

#---- Distribution of covariates by death ----
aa_apoe_long_tte_selected_imp1 %>% filter(subjid == 594947172) %>% 
  select(subjid, fu_yr, contains("event_")) %>% View()

aa_apoe_long_tte_selected_imp1 <- aa_apoe_long_tte_selected_imp %>%
  filter(imp == 1)

aa_apoe_long_tte_selected_imp1 %>%
  filter(fu_yr == 3, !is.na(event_death)) %>% 
  ggplot() +
  geom_density(aes(x = survey_age, color = as.factor(event_death))) +
  facet_wrap(~ethnicity_rev, 
             labeller = as_labeller(c("2" = "Chinese", 
                                      "3" = "Japanese", 
                                      "5" = "Filipino", 
                                      "9" = "White")))

aa_apoe_long_tte_selected_imp1 %>%
  filter(fu_yr == 10, !is.na(event_death)) %>% 
  ggplot() +
  geom_density(aes(x = survey_age, color = as.factor(event_death))) +
  facet_wrap(~ethnicity_rev, 
             labeller = as_labeller(c("2" = "Chinese", 
                                      "3" = "Japanese", 
                                      "5" = "Filipino", 
                                      "9" = "White")))


with(aa_apoe_long_tte_selected_imp1, table(fu_yr, ethnicity_rev, event_death,
                                           useNA = "ifany"))
with(aa_apoe_long_tte_selected_imp1 %>%
       filter(fu_yr == 3), table(ethnicity_rev, event_death, useNA = "ifany"))

with(aa_apoe_long_tte_selected_imp1 %>%
       filter(fu_yr == 10), table(ethnicity_rev, event_death, useNA = "ifany"))

with(aa_apoe_long_tte_selected_imp1 %>%
       filter(fu_yr == 3), table(event_death, event_endmem, useNA = "ifany"))

#---- Functions ----
wt_cov_bal <- function(longdata, ethn_vals,  wt = FALSE, cens_event = "death",
                       fml = NULL, glmnet = FALSE) {
  # # testing arguments
  # longdata <- aa_apoe_long_tte_selected_imp
  # ethn_vals <- 2 # Chinese
  # # for death as censoring event
  # wt = T
  # cens_event = "death"
  # fml <- c(# time term and treatment
  #   "event_death ~ fu_yr + I(fu_yr ^ 2) + apoe_y",
  #   # sociodemo
  #   " + female + survey_age + usaborn_rev + edu_3 + income_pp ")
  # # wt = F
  # # # for death as censoring event
  # # cens_event = "death"
  # # fml <- NULL
  
  # subset to a specific ethnicity
  longdata_ethn <- longdata %>% filter(ethnicity_rev %in% ethn_vals)
  
  original_vars <- colnames(longdata_ethn)
  
  #---- Model PS & weights ----
  if (wt == TRUE){
    for (i in 1:num_imp) {
      longdata_ethn_imp <- longdata_ethn %>%
        filter(imp %in% i) %>%
        as.data.frame()
      
      # if (glmnet == F){
      # denominator
      plrFit <- glm(as.formula(paste(fml, collapse = "")), 
                    data = longdata_ethn_imp, family = binomial())
      
      # numerator for stabilized weights
      plrFit_num <- glm(as.formula(fml[1]), 
                        data = longdata_ethn_imp, family = binomial())
      
      longdata_ethn_imp <- longdata_ethn_imp %>% 
        mutate(pS = 1 - predict(plrFit, newdata = longdata_ethn_imp, 
                                type = "response"), 
               pS_num = 1 - predict(plrFit_num, newdata = longdata_ethn_imp, 
                                    type = "response")) %>%
        group_by(subjid) %>% 
        mutate(denom_cum = cumprod(pS), 
               num_cum = cumprod(pS_num)) %>% 
        ungroup()
      
      pred_wt_vars <- c("denom_cum", "num_cum")
      
      longdata_ethn[which(longdata_ethn$imp == i), pred_wt_vars] <- 
        longdata_ethn_imp %>% select(all_of(pred_wt_vars))
      
      ##---- Attempt for glmnet ----
      #   # PENDING
      # } else if (glmnet == T){
      #   response <- longdata_ethn_imp %>% select(paste0("event_", cens_event)) %>%
      #     as.matrix()
      #   create_intx <- as.formula(response ~ .^2)
      #   pred_mat <- longdata_ethn_imp %>% select(all_of(cov_bal_vars))
      #   predictors_intx <- model.matrix(create_intx, pred_mat)[, -1]
      #   
      #   elasticnet_second <- glmnet::glmnet(x = predictors_intx,
      #                                       y = response,
      #                                       family = "binomial"(link = "logit"),
      #                                       nlambda = 50)
      #   
      #   cv_second <- glmnet::cv.glmnet(x = predictors_intx,
      #                                  y = response,
      #                                  family = "binomial"(link = "logit"),
      #                                  type.measure = "mse",
      #                                  nfolds = 20)
      #   
      #   cv_second$lambda.min
      #   cv_second$lambda.1se
      #   
      #   coef(elasticnet_second, s = cv_second$lambda.min)
      #   
      #   # Generate IOWs
      #   # Main effects Only
      #   
      #   # harmonized$predprob_glmnet <- predict(elasticnet_second, newx = predictors_intx[1:1350,], type = "response", s = cv_second$lambda.min)
      # }
      
    }
    
    longdata_ethn <- longdata_ethn %>% 
      mutate(wt = 1 / denom_cum, # unstabilized
             wt_stab = num_cum / denom_cum, # stabilized
             # Truncate stabilized censoring weights
             wt_stab_trunc = ifelse(wt_stab > quantile(wt_stab, 0.99), 
                                    quantile(wt_stab, 0.99), wt_stab))
    
    # exclude intermediate variables
    longdata_ethn <- longdata_ethn %>% select(all_of(original_vars),
                                              starts_with("wt"))
  }
  
  #---- Calculate covariate balance -----
  for (yr in fu_yrs_cov){
    for (i in 1:num_imp){
      temp_dat1 <- longdata_ethn %>% 
        filter(imp %in% i, fu_yr == yr,
               !is.na(!!sym(paste0("event_", cens_event))))
      
      if (wt == FALSE){
        temp_dat1 <- temp_dat1 %>%
          mutate(wt_stab_trunc = 1)
      }
      
      temp <- cobalt::col_w_smd(
        mat = temp_dat1 %>% select(all_of(cov_bal_vars)),
        treat = temp_dat1[, paste0("event_", cens_event)],
        weights = temp_dat1[, "wt_stab_trunc"], 
        std = TRUE,
        abs = F,
        s.d.denom = "control") # (mean(treated)-mean(control))/sd(control)
      # # (mean(censor)-mean(noncensor))/sd(noncensor)
      
      # store std.eff.sz for each imputation and for each year
      if (i == 1 & yr == fu_yrs_cov[1]){
        std_eff_sz <- expand_grid(fu_yr = fu_yrs_cov, var = names(temp))
      }
      # reverse the sign
      # (mean(noncensor)-mean(censor))/sd(noncensor)
      std_eff_sz[std_eff_sz$fu_yr == yr, paste0("imp_", i)] <- -temp
    }
  }
  # calculate mean std effect size across imputations
  std_eff_sz[, "mean"] <- rowMeans(std_eff_sz %>% select(starts_with("imp")))
  
  return(
    if(wt == T){list(longdata_ethn = longdata_ethn, 
                     std_eff_sz = std_eff_sz,
                     plr_model = plrFit)
    } else if (wt == F){std_eff_sz}
  )
}

##---- Covbal plot ----
cov_bal_yr_plot <- function(data){
  p <- data %>%
    ggplot(aes(x = mean, y = var)) +
    geom_point(size = 1) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", 
               color = "dark grey") +
    labs(x = "SMD\n(noncensor-censor)/SD[noncensor]", 
         y = NULL) +
    facet_wrap(~fu_yr) +
    theme_bw() +
    theme(
      axis.title.y = element_text(size = 10),
      legend.position = "none", 
      aspect.ratio = 5 / 3) +
    xlim(-2.1, 1.1)
  return(p)
}

cov_bal_yr_compare_plot <- function(data_unw, data_wtd, censor_type){
  data <- data_unw %>% mutate(type = "Unweighted") %>%
    rbind(data_wtd %>% mutate(type = "Weighted"))
  
  p <- data %>%
    ggplot(aes(x = mean, y = var)) +
    geom_point(aes(shape = type), size = 1, color = "#0072B2") +
    scale_shape_manual(values = c(21, 16)) +
    geom_vline(xintercept = 0) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", 
               color = "dark grey") +
    labs(x = "SMD\n(noncensor-censor)/SD[noncensor]", 
         y = paste0("IPCW for ", censor_type)) +
    facet_wrap(~fu_yr) +
    theme_bw() +
    theme(aspect.ratio = 5/3,
          panel.grid.minor.x = element_blank(),
          text = element_text(size = 10, color = "black"),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"), #, face = "bold"),
          legend.position = "none") +
    xlim(-2, 1.1)
  
  return(p)
}

#----- Weight dev wide dataset ----
# Test on one imputation
aa_apoe_tte_selected_imp1 <- 
  aa_apoe_long_tte_selected_imp1 %>% select(subjid:main_dem_v1_end_age_r) %>%
  distinct() %>%
  mutate(event_death = case_when(main_dem_v1_end_type == "DEATH" ~ 1,
                                 main_dem_v1_end_type == "ADMIN CENSORED" ~ 0),
         event_endmem = case_when(main_dem_v1_end_type == "END OF MEMBERSHIP" ~ 1,
                                  main_dem_v1_end_type == "ADMIN CENSORED" ~ 0))


ipcw_coxph <- function(widedata, ethn_vals,  cens_event = "death",
                       fml = NULL){
  # # testing arguments
  # widedata = aa_apoe_long_tte_selected_imp1
  # ethn_vals = 2
  # fml = ch1_death
  
  data_forcox <- widedata %>% filter(ethnicity_rev %in% ethn_vals)
  coxph_model_ethn <- coxph(as.formula(fml), 
                            data = data_forcox)
  pred.death_coxph <- data_forcox %>%
    filter(!is.na(get(paste0("event_", cens_event)))) %>%
    mutate(p_surv_denom = predict(coxph_model_ethn, newdata = ., 
                                  type = "survival"),
           ipcw = 1/p_surv_denom,
           ipcw_trunc = ifelse(ipcw > quantile(ipcw, 0.99),
                               quantile(ipcw, 0.99), ipcw)) %>%
    select(subjid, ethnicity_rev, p_surv_denom, ipcw, ipcw_trunc)
  
  data_forcox <- data_forcox %>%
    left_join(pred.death_coxph, by = c("subjid", "ethnicity_rev"))
  
  return(list(coxph_model_ethn = coxph_model_ethn,
              data_forcox = data_forcox,
              pred.death_coxph = pred.death_coxph))
}

aa_apoe_tte_selected_imp1_coxph <- aa_apoe_tte_selected_imp1 %>%
  mutate(main_dem_v1_fu_time = case_when(main_dem_v1_fu_time >= 12 ~ 12,
                                         TRUE ~ main_dem_v1_fu_time))

fml0_coxph <-  "Surv(time = main_dem_v1_fu_time, event = event_death) "
base_age_term <- "ns(survey_age, df = 4)"
ch1_death <- paste0(fml0_coxph, "~",
                    "+ female + ", base_age_term, " + smoking_2 + generalhealth_3 +
               prev_htn_flag + prev_combined_stroke_flag + prev_diab_flag")
# ch1_death <- paste0(fml0_coxph, "~",
#                     "+ female + ", base_age_term, " + generalhealth_3 + 
#                prev_htn_flag")

test_ch1_death <- ipcw_coxph(aa_apoe_tte_selected_imp1_coxph, 2, "death", ch1_death)
summary(test_ch1_death$coxph_model_ethn)


plot_km(data = aa_apoe_long_tte_selected_imp1 %>% filter(ethnicity_rev == 2, fu_yr <= 12) %>%
          left_join(test_ch1_death$pred.death_coxph),
        weight = T, weightvar = "ipcw_trunc",
        stratify_by = "survey_age_cat") +
  labs(title = "weighted (IPCW from cox ph)")
plot_km(data = data_for_cif %>% filter(ethnicity_rev == 2) %>%
          filter(fu_yr <= 12),
        stratify_by = "survey_age_cat") +
  labs(title = "observed")

pred.death_coxph_long <- aa_apoe_long_tte_selected_imp1 %>% 
  left_join(pred.death_coxph, by = c("subjid", "ethnicity_rev"))

plot_km(data = pred.death_coxph_long,
        weight = T, weightvar = "ipcw_trunc",
        stratify_by = "survey_age_cat")
plot_km(data = data_for_cif %>% filter(ethnicity_rev == 2),
        stratify_by = "female")
plot_km(data = pred.death_coxph_long,
        weight = T, weightvar = "ipcw",
        stratify_by = "female")

# coxph_fit <- survfit(coxph_model, newdata = aa_apoe_tte_selected_imp1_chn)
# plot(coxph_fit)
km_pred.death_coxph <- survfit(Surv(time = main_dem_v1_fu_time, 
                                    event = pred_death_event) ~ 1,
                               data = pred.death_coxph, 
                               cluster = subjid)
plot(km_pred.death_coxph)
plot(survfit(Surv(time = start,  time2 = fu_yr, event = event_death) ~ 1, 
             data = aa_apoe_long_tte_selected_imp1 %>% filter(ethnicity_rev == 2), 
             cluster = subjid))



#---- Check covariate balance by ethnicity ----
# time period
# Based on the distribution of the censoring event (must be binary)
fu_yrs_cov <- 3:12
# time term
fml0 <-  "event_death ~ ns(fu_yr, df = 4)"
base_age_term <- "ns(survey_age, df = 4)"

##---- Chinese ----
###---- Unweighted ----
unwtd_ch_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                             wt = F, cens_event = "death",
                             ethn_vals = 2)
unwtd_ch_endmem <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                              wt = F, cens_event = "endmem",
                              ethn_vals = 2)
cov_bal_yr_plot(unwtd_ch_death)
cov_bal_yr_plot(unwtd_ch_endmem)

###---- Weight 1----
ch1_death <- c(fml0, 
               "+ female + ", base_age_term, " + smoking_2 + generalhealth_3 + 
               prev_htn_flag + prev_combined_stroke_flag + prev_diab_flag")

wtd_ch1_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp %>% filter(fu_yr <= 12), 
                            ethn_vals = 2, 
                            wt = T, cens_event = "death", ch1_death)
cov_bal_yr_compare_plot(unwtd_ch_death, 
                        wtd_ch1_death$std_eff_sz, "death")
broom::tidy(wtd_ch1_death$plr_model) %>%
  mutate_if(is.numeric, round, 4)

pred.death <- aa_apoe_long_tte_selected_imp %>%
  # Use the last imputation for prediction
  filter(imp == 20, ethnicity_rev == 2, fu_yr <= 12) %>%
  # predicted hazard
  mutate(p1 = predict(wtd_ch1_death$plr_model, 
                      newdata = ., 
                      type = "response"),
         p0 = predict(plrFit_num, newdata = ., type = "response")) %>%
  group_by(subjid) %>%
  arrange(fu_yr) %>%
  mutate(
    s1 = 1 - p1,
    s0 = 1 - p0,
    cum_s1 = cumprod(s1),
    cum_s0 = cumprod(s0),
    cif_1 = 1 - cum_s1,
    cif_0 = 1 - cum_s0) %>%
  ungroup() %>%
  mutate(pred_event_death1 = rbinom(n(), 
                                    size = 1, 
                                    prob = cif_1),
         pred_event_death0 = rbinom(n(),
                                    size = 1, 
                                    prob = cif_0)) 

pred.death1 <- pred.death %>%
  group_by(subjid) %>%
  mutate(cum_pred_event_death1 = cumsum(pred_event_death1),
         cum_pred_event_death0 = cumsum(pred_event_death0)) %>%
  ungroup() %>%
  select(subjid, fu_yr, 
         p1, p0, 
         cif_1, cif_0,
         pred_event_death1, pred_event_death0,
         cum_pred_event_death1, cum_pred_event_death0)
# filter(pred_event_death1 == 1)

pred.death_subj <- pred.death1 %>%
  select(-contains("0")) %>%
  filter(cum_pred_event_death1 == 1) %>%
  group_by(subjid) %>%
  arrange(fu_yr) %>%
  slice(1)

pred.death_subj <- pred.death_subj %>%
  rbind(pred.death1 %>%
          select(-contains("0")) %>%
          filter(cum_pred_event_death1 == 0, !subjid %in% pred.death_subj$subjid) %>%
          group_by(subjid) %>%
          arrange(desc(fu_yr)) %>%
          slice(1))


model1 <- survfit(Surv(time = fu_yr, event =  cum_pred_event_death1) ~ 1, 
                  data = pred.death_subj, cluster = subjid)
plot(model1)
title(main = "Survival curve (predicted death event, Chinese)")

model_km <- survfit(Surv(time = start,  time2 = fu_yr, event = event_death) ~ 1, 
                    data = aa_apoe_long_tte_selected_imp %>% filter(imp == 20, ethnicity_rev == 2, fu_yr <= 12), 
                    cluster = subjid)
plot(model_km)
title(main = "Survival curve (observed, Chinese)")

summary(pred.death$p1)
summary(pred.death$p0)
summary(pred.death$cif_1)
summary(pred.death$cif_0)

pred.death %>%
  group_by(fu_yr) %>%
  summarise(
    mean_cif_1 = mean(cif_1, na.rm = TRUE),
    mean_cif_0 = mean(cif_0, na.rm = TRUE)) %>%
  ggplot() +
  geom_line(aes(x = fu_yr, y = mean_cif_1 * 100), 
            color = "#0072B2", size = 1) +
  geom_line(aes(x = fu_yr, y = mean_cif_0 * 100,
                color = "red")) +
  scale_x_continuous(
    limits = c(0, 18),
    breaks = seq(0, 18, 3), 
    minor_breaks = seq(0, 18, 1))

pred.death %>%
  filter(subjid == 8050719474) %>%
  ggplot() +
  geom_line(aes(x = fu_yr, 
                y = cif_1 * 100), 
            color = "#0072B2", size = 1) +
  geom_line(aes(x = fu_yr, 
                y = cif_0 * 100,
                color = "red"))



with(aa_apoe_long_tte_selected_imp1, table(fu_yr, ethnicity_rev, event_death, useNA = "ifany"))

# Weighted CIF did not change.
plot_km(data = data_for_cif %>% filter(ethnicity_rev == 2),
        stratify_by = "survey_age_cat")
plot_km(data = wtd_ch1_death$longdata_ethn %>% filter(imp == 1),
        weight = T, weightvar = "wt_stab_trunc",
        stratify_by = "survey_age_cat")
plot_km(data = data_for_cif %>% filter(ethnicity_rev == 2),
        stratify_by = "female")
plot_km(data = wtd_ch1_death$longdata_ethn %>% filter(imp == 1),
        weight = T, weightvar = "wt_stab_trunc",
        stratify_by = "female")

# Check weights
wt_sum_ch1 <- wtd_ch1_death$longdata_ethn %>% group_by(imp) %>% 
  summarize(sum_wt = sum(wt))
wt_sum_ch1

with(wtd_ch1_death$longdata_ethn %>% filter(imp == 1),
     summary(wt))
with(wtd_ch1_death$longdata_ethn %>% filter(imp == 1),
     summary(wt_stab_trunc))

with(wtd_ch1_death$longdata_ethn, table(imp, useNA = "ifany"))


###---- Weight 2----
# ch2_death <- c(fml0, 
#                " + ", base_age_term, " + generalhealth_3")
# 
# wtd_ch2_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
#                             ethn_vals = 2, 
#                             wt = T, cens_event = "death", ch2_death)
# cov_bal_yr_compare_plot(unwtd_ch_death, 
#                         wtd_ch2_death$std_eff_sz, "death")
# broom::tidy(wtd_ch2_death$plr_model) %>%
#   mutate_if(is.numeric, round, 4)

###---- Weight 3 ----
# Full model with no interaction terms
cov_terms_ipcw <- c("apoe_y", "female", base_age_term, "usaborn_rev", "edu_3", 
                    "income_pp", "generalhealth_3", "smoking_2", "prev_htn_flag",
                    "prev_combined_stroke_flag", "prev_diab_flag", "sr_bmi", 
                    "sr_depress", "alcohol_3", "pa_3", "ehr_ht_median")
ch3_death <- c(fml0, " + ",
               paste(cov_terms_ipcw, collapse = " + "))

wtd_ch3_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                            ethn_vals = 2, 
                            wt = T, cens_event = "death", ch3_death)
cov_bal_yr_compare_plot(unwtd_ch_death, 
                        wtd_ch3_death$std_eff_sz, "death")
broom::tidy(wtd_ch3_death$plr_model) %>%
  mutate_if(is.numeric, round, 4) %>% print(n = Inf)

# Weighted CIF did not change.
plot_km(data = wtd_ch3_death$longdata_ethn %>% filter(imp == 1),
        weight = T, weightvar = "wt_stab_trunc",
        stratify_by = "survey_age_cat")
plot_km(data = data_for_cif %>% filter(ethnicity_rev == 2),
        stratify_by = "survey_age_cat")

plot_km(data = data_for_cif %>% filter(ethnicity_rev == 2),
        stratify_by = "female")
plot_km(data = wtd_ch3_death$longdata_ethn %>% filter(imp == 1),
        weight = T, weightvar = "wt_stab_trunc",
        stratify_by = "female")
plot_km(data = wtd_ch3_death$longdata_ethn %>% filter(imp == 1),
        weight = T, weightvar = "wt_stab_trunc",
        stratify_by = "generalhealth_3")

###---- Weight 4 ----
cov_terms_ipcw1 <- c(base_age_term, "edu_3", 
                     "generalhealth_3", "prev_htn_flag")
ch4_death <- c(fml0, " + ",
               paste(cov_terms_ipcw1, collapse = " + "))

wtd_ch4_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                            ethn_vals = 2, 
                            wt = T, cens_event = "death", ch4_death)
cov_bal_yr_compare_plot(unwtd_ch_death, 
                        wtd_ch4_death$std_eff_sz, "death")
broom::tidy(wtd_ch4_death$plr_model) %>%
  mutate_if(is.numeric, round, 4) %>% print(n = Inf)
# Weighted CIF did not change.
plot_km(data = wtd_ch4_death$longdata_ethn %>% filter(imp == 1),
        weight = T, weightvar = "wt_stab_trunc",
        stratify_by = "survey_age_cat")
plot_km(data = data_for_cif %>% filter(ethnicity_rev == 2),
        stratify_by = "female")
plot_km(data = wtd_ch4_death$longdata_ethn %>% filter(imp == 1),
        weight = T, weightvar = "wt_stab_trunc",
        stratify_by = "female")
plot_km(data = wtd_ch4_death$longdata_ethn %>% filter(imp == 1),
        weight = T, weightvar = "wt_stab_trunc",
        stratify_by = "generalhealth_3")

###---- Weight 5 ----
cov_terms_ipcw2 <- c(base_age_term, "edu_3", paste0(base_age_term, " * generalhealth_3"),
                     paste0(base_age_term, " * ehr_ht_median"), "usaborn_rev:prev_htn_flag",
                     "usaborn_rev:prev_diab_flag",  "edu_3:prev_combined_stroke_flag",
                     "income_pp:prev_htn_flag", "generalhealth_3:prev_htn_flag",
                     "prev_combined_stroke_flag:pa_3", "alcohol_3:pa_3")
ch5_death <- c(fml0, " + ",
               paste(cov_terms_ipcw2, collapse = " + "))

wtd_ch5_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                            ethn_vals = 2, 
                            wt = T, cens_event = "death", ch5_death)
cov_bal_yr_compare_plot(unwtd_ch_death, 
                        wtd_ch5_death$std_eff_sz, "death")
broom::tidy(wtd_ch5_death$plr_model) %>%
  mutate_if(is.numeric, round, 4) %>% print(n = Inf)
# Weighted CIF did not change.
plot_km(data = wtd_ch5_death$longdata_ethn %>% filter(imp == 1),
        weight = T, weightvar = "wt_stab_trunc",
        stratify_by = "survey_age_cat")
plot_km(data = data_for_cif %>% filter(ethnicity_rev == 2),
        stratify_by = "female")
plot_km(data = wtd_ch5_death$longdata_ethn %>% filter(imp == 1),
        weight = T, weightvar = "wt_stab_trunc",
        stratify_by = "female")


# End of follow up
ch1_endmem <- c(fml0,
                # sociodemo
                " + female + survey_age + edu_3 + income_pp",
                " + smoking_2 + prev_htn_flag + prev_combined_stroke_flag +
    prev_diab_flag + sr_bmi + sr_depress + alcohol_3 + pa_3 + ehr_ht_median")
cov_bal_yr_compare_plot(unwtd_ch0$std_eff_sz_endmem, 
                        wtdev_ch2$std_eff_sz_endmem, "end membership")

##---- Japanese ----
###---- Unweighted ----
unwtd_jp_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                             wt = F, cens_event = "death",
                             ethn_vals = 3)
unwtd_jp_endmem <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                              wt = F, cens_event = "endmem",
                              ethn_vals = 3)
cov_bal_yr_plot(unwtd_jp_death)
cov_bal_yr_plot(unwtd_jp_endmem)

###---- Weight 1----
jp1_death <- c(fml0, 
               "+ female + survey_age + smoking_2 + generalhealth_3 + 
               prev_htn_flag + prev_combined_stroke_flag + prev_diab_flag")

wtd_jp1_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                            ethn_vals = 3, 
                            wt = T, cens_event = "death", jp1_death)
cov_bal_yr_compare_plot(unwtd_ch_death, 
                        wtd_jp1_death$std_eff_sz, "death")
broom::tidy(wtd_jp1_death$plr_model)

##---- Filipino ----
###---- Unweighted ----
unwtd_fl_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                             wt = F, cens_event = "death",
                             ethn_vals = 5)
unwtd_fl_endmem <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                              wt = F, cens_event = "endmem",
                              ethn_vals = 5)
cov_bal_yr_plot(unwtd_fl_death)
cov_bal_yr_plot(unwtd_fl_endmem)

##---- AA ----
###---- Unweighted ----
unwtd_AA_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                             wt = F, cens_event = "death",
                             ethn_vals = c(2, 3, 5))
unwtd_AA_endmem <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                              wt = F, cens_event = "endmem",
                              ethn_vals = c(2, 3, 5))
cov_bal_yr_plot(unwtd_AA_death)
cov_bal_yr_plot(unwtd_AA_endmem)

##---- White ----
###---- Unweighted ----
unwtd_wh_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                             wt = F, cens_event = "death",
                             ethn_vals = 9)
unwtd_wh_endmem <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                              wt = F, cens_event = "endmem",
                              ethn_vals = 9)
cov_bal_yr_plot(unwtd_wh_death)
cov_bal_yr_plot(unwtd_wh_endmem)
# Endmem pretty balanced! test on adding survey age to the weight worsen the balance a lot
# no weights for White endmem

###---- Weight 1----
wh1_death <- c(fml0, "+ survey_age * female + generalhealth_3")
wtd_wh1_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                            ethn_vals = 9, 
                            wt = T, cens_event = "death", wh1_death)
cov_bal_yr_compare_plot(unwtd_wh_death, 
                        wtd_wh1_death$std_eff_sz, "death")

wh1b_death <- c(fml0, "+ survey_age * female + generalhealth_3 + prev_htn_flag + income_pp + pa_3")
wtd_wh1b_death <- wt_cov_bal(longdata = aa_apoe_long_tte_selected_imp, 
                             ethn_vals = 9, 
                             wt = T, cens_event = "death", wh1b_death)
cov_bal_yr_compare_plot(unwtd_wh_death, 
                        wtd_wh1b_death$std_eff_sz, "death")


#---- **OLD** ----
#---- Function model death and endmem in one function ----
# wt_cov_bal <- function(longdata, wt = FALSE, ethn_vals, fml_death = NULL, 
#                        fml_endmem = NULL) {
#   # # testing arguments
#   # longdata <- aa_apoe_long_tte_selected_imp
#   # ethn_vals <- 2 # Chinese
#   # wt = T
#   # # for death as censoring event
#   # fml_death <- c(# time term and treatment
#   #   "event_death ~ fu_yr + I(fu_yr ^ 2) + apoe_y",
#   #   # sociodemo
#   #   " + female + survey_age + usaborn_rev + edu_3 + income_pp ")
#   # # for end of membership as censoring event
#   # fml_endmem <- c(# time term and treatment
#   #   "event_endmem ~ fu_yr + I(fu_yr ^ 2) + apoe_y",
#   #   # sociodemo
#   #   " + female + survey_age + usaborn_rev + edu_3 + income_pp")
#   
#   # subset to a specific ethnicity
#   longdata_ethn <- longdata %>% filter(ethnicity_rev %in% ethn_vals)
#   
#   original_vars <- colnames(longdata_ethn)
#   
#   if (wt == TRUE){
#     for (i in 1:num_imp) {
#       longdata_ethn_imp <- longdata_ethn %>%
#         filter(imp %in% i) %>%
#         as.data.frame()
#       
#       #---- Death as censoring event ----
#       # denominator
#       plrFit_Death <- glm(as.formula(paste0(fml_death, collapse = "")), 
#                           data = longdata_ethn_imp, family = binomial())
#       
#       # numerator for stabilized weights
#       plrFit_Death_num <- glm(as.formula(fml_death[1]), 
#                               data = longdata_ethn_imp, family = binomial())
#       
#       longdata_ethn_imp <- longdata_ethn_imp %>% 
#         mutate(predDeath = 1 - predict(plrFit_Death, newdata = longdata_ethn_imp, 
#                                        type = "response"), 
#                predDeath_num = 1 - predict(plrFit_Death_num, newdata = longdata_ethn_imp, 
#                                            type = "response")) %>%
#         group_by(subjid) %>% 
#         mutate(denom_cum_death = cumprod(predDeath), 
#                num_cum_death = cumprod(predDeath_num)) %>% 
#         ungroup()
#       
#       #---- End of membership censoring ----
#       # denominator
#       plrFitC <- glm(as.formula(paste0(fml_endmem, collapse = "")), 
#                      data = longdata_ethn_imp, family = binomial())
#       
#       # numerator
#       plrFitCnum <- glm(as.formula(fml_endmem[1]), 
#                         data = longdata_ethn_imp, family = binomial())
#       
#       longdata_ethn_imp <- longdata_ethn_imp %>% 
#         mutate(predC = 1 - predict(plrFitC, newdata = longdata_ethn_imp, type = "response"),
#                predC_num = 1 - predict(plrFitCnum, newdata = longdata_ethn_imp, 
#                                        type = "response")) %>% 
#         group_by(subjid) %>% 
#         mutate(denom_cum_endmem = cumprod(predC), 
#                num_cum_endmem = cumprod(predC_num)) %>% 
#         ungroup()
#       
#       pred_wt_vars <- c("denom_cum_endmem", "num_cum_endmem",
#                         "denom_cum_death", "num_cum_death")
#       
#       longdata_ethn[which(longdata_ethn$imp == i), pred_wt_vars] <- 
#         longdata_ethn_imp %>% select(all_of(pred_wt_vars))
#     }
#     
#     #---- Calculate weights ----
#     longdata_ethn <- longdata_ethn %>% 
#       mutate(# censoring due to end of membership
#         wt_endmem = 1 / denom_cum_endmem, # unstabilized
#         wt_endmem_stab = num_cum_endmem / denom_cum_endmem, # stabilized
#         # Truncate stabilized censoring weights
#         wt_endmem_stab_trunc = ifelse(wt_endmem_stab > quantile(wt_endmem_stab, 0.99), 
#                                       quantile(wt_endmem_stab, 0.99), wt_endmem_stab),
#         
#         # censoring due to death
#         wt_death = 1 / denom_cum_death,
#         wt_death_stab = num_cum_death / denom_cum_death,
#         # Truncate stabilized censoring weights
#         wt_death_stab_trunc = ifelse(wt_death_stab > quantile(wt_death_stab, 0.99), 
#                                      quantile(wt_death_stab, 0.99), wt_death_stab),
#         
#         # combined censoring weights
#         wt_cens = wt_endmem * wt_death, 
#         wt_cens_stab = wt_endmem_stab * wt_death_stab,
#         # Truncate stabilized censoring weights
#         wt_cens_stab_trunc = ifelse(wt_cens_stab > quantile(wt_cens_stab, 0.99), 
#                                     quantile(wt_cens_stab, 0.99), wt_cens_stab))
#     
#     # exclude intermediate variables
#     longdata_ethn <- longdata_ethn %>% select(all_of(original_vars),
#                                               starts_with("wt_"))
#   }
#   
#   #---- Calculate covariate balance -----
#   for (cens in c("death", "endmem")){
#     for (yr in fu_yrs_cov){
#       for (i in 1:num_imp){
#         temp_dat1 <- longdata_ethn %>% 
#           filter(imp %in% i, fu_yr == yr,
#                  !is.na(!!sym(paste0("event_", cens))))
#         
#         if (wt == FALSE){
#           temp_dat1 <- temp_dat1 %>%
#             mutate(!!paste0("wt_", cens, "_stab_trunc") := 1)
#         }
#         
#         temp <- cobalt::col_w_smd(
#           mat = temp_dat1 %>% select(all_of(cov_bal_vars)),
#           treat = temp_dat1[, paste0("event_", cens)],
#           weights = temp_dat1[, paste0("wt_", cens, "_stab_trunc")], 
#           std = TRUE,
#           abs = F,
#           s.d.denom = "control") # (mean(treated)-mean(control))/sd(control)
#         # # (mean(censor)-mean(noncensor))/sd(noncensor)
#         
#         # store std.eff.sz for each imputation and for each year
#         if (i == 1 & yr == fu_yrs_cov[1]){
#           std_eff_sz <- expand_grid(fu_yr = fu_yrs_cov, var = names(temp))
#           # reverse the sign
#           # (mean(noncensor)-mean(censor))/sd(noncensor)
#           std_eff_sz[std_eff_sz$fu_yr == yr, paste0("imp_", i)] <- -temp
#         } else {
#           # reverse the sign
#           std_eff_sz[std_eff_sz$fu_yr == yr, paste0("imp_", i)] <- -temp
#         }
#       }
#     }
#     # calculate mean std effect size across imputations
#     std_eff_sz[, "mean"] <- rowMeans(std_eff_sz %>% select(starts_with("imp")))
#     assign(paste0("std_eff_sz_", cens), std_eff_sz)
#   }
#   
#   return(list(longdata_ethn_wt = 
#                 longdata_ethn %>% select(subjid, imp, starts_with("wt_")), 
#               std_eff_sz_death = std_eff_sz_death, 
#               std_eff_sz_endmem = std_eff_sz_endmem))
# }
