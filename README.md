# aa_apoe_adrd

This repo contains all the code necessary (folder scripts) to replicate data construction and analyses in the paper: 

Association of APOE-ε4 genotype and risk of dementia over 10 years in a cohort of Asian American and non-Latino White older adults in California

## Scripts description
*R scripts for data construction, statistical analysis, cluster submission, and generating tables and figures were in the `scripts` folder. The R scripts are named in the order they should be run. Functions created by the authors were sourced where needed.*
* **0.paths.R**: specifies paths to data folders. 
* **1.0._data_construction_main.R** and **1.0.1_data_construction_EHR.R**: cleans and checks the distributions for genetic, survey-reported and EHR measures.
* **1.2_long_data_construction.R**: constructs the time-to-event dataset in a long format for analysis. 
* **2.1_descriptives_cumulative_incidence.R**: calculates cumulative incidence at specific follow-up years (Aalen–Johansen estimators).
* **2.2_cox_ph_models.R**: fit the Cox PH model for Asian aggregated and disaggregated and output results
* **3.1_plr_bootstrap_int1_cluster.R**: conducts main analysis: G-computed pooled logistic regression model with bootstrap on the Hoffman cluster. The script contains a function for bootstrapping data, fitting a pooled logistic regression Q model, creating clone datasets, obtaining potential outcomes in the cloned dataset, and calculating probabilities, survival probabilities, cumulative incidence, risk, and RD and RR. 
* **3.2_submission.sh**: submission script for `3.1_plr_bootstrap_int1_cluster.R` to Hoffman cluster
* **3.3_cluster_results_int1.R**: cleans analysis results and prepares tables and figures.
* **9.tables_preimp.R**: generates descriptive tables
