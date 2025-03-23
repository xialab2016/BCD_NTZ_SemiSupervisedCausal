# BCD & NTZ Analysis 
## Dominic DiSanto
## Updated October.17.2023



# Packages
if("pacman" %in% installed.packages()){
  library(pacman)
}else{
  install.packages("pacman")
  library(pacman)
}

cat("Initial wd:", getwd())
if(Sys.getenv("HOME")!="/Users/jdomi") setwd("/n/data1/hsph/biostat/celehs/lab/SHARE/MS_Causal")
cat("Set wd:", getwd())


pacman::p_load(rms, ordinal, ordinalNet
               , foreach, doParallel, doRNG
               , lubridate
               , stringr
               , here) # necessary to load shared set-up
summarize = dplyr::summarize # loading over Hmisc's summarize function 

# Cluster Submission 
  # lazily check if on my machine, otherwise set to cluster wd 

bash_args = commandArgs(trailingOnly = T)
cat("Bash Arguments:", bash_args)
.n_perturb = as.numeric(bash_args[1])
cat(".n_perturb (first bash argument) arg is", .n_perturb)

if("3yr" %in% bash_args){
  .look_yrs = 3
}else{
  .look_yrs = 2
}

# .prefix = paste0("perturb_1300_", .look_yrs, "yr_Increase")  
.prefix = bash_args[3]
if(is.na(.prefix)) .prefix = paste0("perturb_ATE_", .look_yrs, "yrs_Increase_"
                                    , gsub("-|[ ]|:|\\.", "_", lubridate::now())
                                    )

.outcome_var = bash_args[4]
# if(is.na(.outcome_var) | .outcome_var %nin% c("Y_change", "Y_NegChange")) stop("Specify outcome, Y_change or Y_NegChange")

# Parameters #### 
  # for any given DMT, either the construct boolean must be T/TRUE/1 or a string filename supplied  
  # if constructing cohorts, PLEASE edit construction 
  # filenames presented without extensions

  .results_subdir <- "BCD_NTZ"
  
  .bcd_filename_3yr <- "bcd_allMS_3yr_2024_03_07.csv" #"bcd_all_3yr_2023_12_20.csv"
  .bcd_filename_2yr <- "bcd_allMS_2yr_2024_03_07.csv" #"bcd_all_2yr_2023_12_20.csv"
  
  .ntz_filename_3yr <- "ntz_allMS_3yr_2024_03_07.csv" #"ntz_all_3yr_2023_12_20.csv"  
  .ntz_filename_2yr <- "ntz_allMS_2yr_2024_03_07.csv" #"ntz_all_2yr_2023_12_20.csv"  

  .ehr_cor_threshold <- 0.1
  DMT_reference <- "NTZ" # set to 0 in treatment assignment vector
  
  
# Data Creation/Import #### 
  .skip_build <- T  # skips cohort_building.R in DMT script 
  source(here("Code", "driver_scripts", "bcd.R"))
  bcd_3yr <- read_csv(here("Data", "created_data", "BCD", .bcd_filename_3yr))
  bcd_2yr <- read_csv(here("Data", "created_data", "BCD", .bcd_filename_2yr))

  .skip_build <- T  # skips cohort_building.R in DMT script 
  source(here("Code", "driver_scripts", "ntz.R"))
  ntz_3yr <- read_csv(here("Data", "created_data", "NTZ", .ntz_filename_3yr))
  ntz_2yr <- read_csv(here("Data", "created_data", "NTZ", .ntz_filename_2yr))

  bcd_ntz_3yr <- rbind(bcd_3yr, ntz_3yr)
  bcd_ntz_2yr <- rbind(bcd_2yr, ntz_2yr)

  
# Creating Shared Cohorts #### 
  bcd_ntz_3yr_joint <- bcd_ntz_3yr %>% 
    inner_join(bcd_ntz_2yr[,"PATIENT_NUM"], "PATIENT_NUM")
  
  bcd_ntz_2yr_joint <- bcd_ntz_2yr %>% 
    inner_join(bcd_ntz_3yr[,"PATIENT_NUM"], "PATIENT_NUM")

  # sample summary 
  bcd_ntz_3yr_joint %>% mutate(`PDDS Outcome` = case_when(PDDS_Change==1 ~ "Increase/Progression"
                                                         , PDDS_NegChange==1 ~ "Decrease/Improvement"
                                                         , PDDS_Change==0 & PDDS_NegChange ==0 ~ "No Sustained Change"
                                                         , T ~ "Insufficient Data (Unlabelled)")
                              ) %>% 
    mutate(Labelled = ifelse(`PDDS Outcome`!="Insufficient Data (Unlabelled)", "Labelled", "Unlabelled")) %>% 
    count(Labelled, `PDDS Outcome`) %>% mutate(`%` = 100*n / sum(n)) %>% xtable::xtable()
  
  bcd_ntz_2yr_joint %>% mutate(`PDDS Outcome` = case_when(PDDS_Change==1 ~ "Increase/Progression"
                                                          , PDDS_NegChange==1 ~ "Decrease/Improvement"
                                                          , PDDS_Change==0 & PDDS_NegChange ==0 ~ "No Sustained Change"
                                                          , T ~ "Insufficient Data (Unlabelled)")
  ) %>% 
    mutate(Labelled = ifelse(`PDDS Outcome`!="Insufficient Data (Unlabelled)", "Labelled", "Unlabelled")) %>% 
    count(Labelled, `PDDS Outcome`) %>% mutate(`%` = 100*n / sum(n)) %>% xtable::xtable()
  
# Load general imputation models of interest for latent scores 
  # lasso_cauchit = readRDS(here("Data", "PDDSModels", "lasso_cauchit_tune_EHR_NLP.RDS"))
  # lasso_logit = readRDS(here("Data", "PDDSModels", "lasso_logit_tune_EHR_NLP.RDS"))
  # ordinal_KG_NLP = readRDS(here("Data", "PDDSModels", "ordinal_cauchit_KG_ONCE.RDS"))
  ridge_12mo_EHR_NLP = readRDS(here("Data", "PDDSModels", "lookback12mo_KGclean_ridge_logit_tune_EHR_NLP.RDS"))
  ridge_All_EHR_NLP = readRDS(here("Data", "PDDSModels", "lookbackAll_KGclean_ridge_logit_tune_EHR_NLP.RDS"))
  
# Analytic wrangling #### 
  source(here("Code", "analysis_scripts", "00_causal_core_Robust_Perturbed.R"))
  source(here("Code", "analysis_scripts", "01_analytic_wrangling.R"))
  source(here("Code", "analysis_scripts", "01_analytic_wrangling_postDMT.R"))
  
  if("3yr" %in% bash_args){
    .df = bcd_ntz_3yr_joint
  }else{
    .df = bcd_ntz_2yr_joint
  }

  
  analytic_object <- analytic_wrangling(.analytic_df = .df
                                        , pre_process = T
                                        , efficacy="high"
                                        , .EHR_trim = 0.95 # also NLP trim
                                        , .NLP = T
                                        , .NLP_all = F
                                        , .EHR_all = F
                                        , .RXNORM_Exclude_All = T
                                        , .DMT_reference = DMT_reference
                                        , .diff_window = dmonths(6)
                                        , .lookback_window = Inf
                                        )
  
  X_lookforward = analytic_wrangling_postDMT(object = .df
                                             , ehr_ft_list = analytic_object[["EHR_fts"]]
                                             , nlp_ft_list = analytic_object[["NLP_fts"]]
                                             , pre_process = T
                                             , .lag = lubridate::dyears(.look_yrs) 
                                             ) # my own identifier warning, taken care of in-house (in-function)
  
  # ATE Estimate 
  set.seed(239579235) 
  ATE <- AIPW_robust_analytic(analytic_object = analytic_object
                              , .analytic_object_1yr = X_lookforward # misnomer, can take any lookforward, I just didn't bother to rename and ensure correctness
                              , .outcome = .outcome_var
                              # , pweights = NULL
                              , pweights = rep(1, nrow(analytic_object$X))
                              , .lookforward_window = dyears(.look_yrs)
                              , cv.outcome = T
                              , crump.trim = F
                              )

  .results = function(x) {cat("ATE" = round(x$ATE, 4), "\nAccuracy = ", round(x$Impute_Acc, 4), "\nAUC = ", round(x$Impute_AUC, 4))   }
  .results(ATE)

  # Perturbed Estimates
  doParallel::registerDoParallel()
  
  start = Sys.time() 
  # .n_perturb = 10
  ATE_perturb = foreach(icount(.n_perturb), .errorhandling="pass") %dorng% {
         AIPW_robust_analytic(analytic_object = analytic_object
                                , .analytic_object_1yr = X_lookforward
                                , .outcome = .outcome_var
                                , pweights = 4*rbeta(nrow(analytic_object$X), 0.5, 1.5)
                                , .lookforward_window = dyears(.look_yrs)
                                , cv.outcome = F
                                )
    }
  end = Sys.time()

  
  saveRDS(object = list(ATE_full = ATE$ATE # ATE[names(ATE)!="ImputeModel"]
                        , Outcome_ImputeModel = ATE$ImputeModel
                        , Outcome_ImputeModel_Coefs = list(Trt = ATE$Trt.coef, OC = ATE$OC.coef)
                        , Outcome_ImputeModel_Perf = c(Accuracy = ATE$Impute_Acc, AUC = ATE$Impute_AUC)
                        , FullModel_IntOnly = c(Trt_KS_IntOnly = ATE$Trt_KS_IntOnly, Oucome_KS_IntOnly = ATE$Outcome_KS_IntOnly)
                        , ATE_Perturbed = unlist(lapply(ATE_perturb, `[[`, 1))
                        , ATE_Trt_KS_IntOnly = unlist(lapply(ATE_perturb, `[[`, 2))
                        , ATE_Outcome_KS_IntOnly = unlist(lapply(ATE_perturb, `[[`, 3))
                        # , ATE_PerturbWeights = lapply(ATE_perturb, `[`, 2)
                        , seeds = attr(ATE_perturb, "rng")
                        , Runtime_Real = end-start
                        , N_Cores_Used = getDoParWorkers()
                        , Cohort_Yrs = .look_yrs
                        , Outcome = .outcome_var
                        )
          , file = here("Results", "BCD_NTZ", "ATE", "Sensitivity", paste0(.prefix, "_noTrim_SensitivityAnalysis.RDS")) 
          )
  