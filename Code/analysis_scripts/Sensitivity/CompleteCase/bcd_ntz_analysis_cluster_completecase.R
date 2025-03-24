# BCD & NTZ Analysis 
## Dominic DiSanto

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
               , dplyr
               , stringr
               , here) # necessary to load shared set-up

summarize = dplyr::summarize # loading over Hmisc's summarize function 




which(is.na(unlist(lapply(ATE_perturb, `[[`, 1))))

.outcome = "Y_avg"
.outcome_yrs = 3
.n_perturb = 200


# bash_args = commandArgs(trailingOnly = T)
# .n_perturb = as.numeric(bash_args[1])
# .outcome_yrs = bash_args[2]
# .outcome = bash_args[3]

.outcome_paste = case_when(.outcome=="Y_change" ~ "Increase"
                           , .outcome=="Y_NegChange" ~ "Decrease"
                           , .outcome=="Y_avg" ~ "Average")
.outcome.df = case_when(.outcome=="Y_change" ~ "PDDS_Change"
                        , .outcome == "Y_NegChange" ~ "PDDS_NegChange"
                        , T ~ "Average_PDDS")

.prefix = paste0("perturb_ATE_", .outcome_yrs, "yrs_", .outcome_paste)

cat(.outcome, "//", .outcome.df, "\n", .prefix, "\nperturb=", .n_perturb, ", outcome.yrs=", .outcome_yrs)


# Parameters #### 
  # for any given DMT, either the construct boolean must be T/TRUE/1 or a string filename supplied  
  # if constructing cohorts, PLEASE edit construction 
  # filenames presented without extensions

  .results_subdir <- "BCD_NTZ"

  .bcd_filename_3yr <- "bcd_allMS_3yr_2024_10_07.csv" #"bcd_all_3yr_2023_12_20.csv"
  .bcd_filename_2yr <- "bcd_allMS_2yr_2024_10_07.csv" #"bcd_all_2yr_2023_12_20.csv"
  
  .ntz_filename_3yr <- "ntz_allMS_3yr_2024_10_07.csv" #"ntz_all_3yr_2023_12_20.csv"  
  .ntz_filename_2yr <- "ntz_allMS_2yr_2024_10_07.csv" #"ntz_all_2yr_2023_12_20.csv"  
  
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
  
  # tmp = read_excel(here("Data", "Gold Standard Labels", "PRT Dx 2024-02-09.xlsx"))
  
  
# Load general imputation models of interest for latent scores 
  # lasso_cauchit = readRDS(here("Data", "PDDSModels", "lasso_cauchit_tune_EHR_NLP.RDS"))
  # lasso_logit = readRDS(here("Data", "PDDSModels", "lasso_logit_tune_EHR_NLP.RDS"))
  # ordinal_KG_NLP = readRDS(here("Data", "PDDSModels", "ordinal_cauchit_KG_ONCE.RDS"))
  ridge_12mo_EHR_NLP = readRDS(here("Data", "PDDSModels", "lookback12mo_KGclean_ridge_logit_tune_EHR_NLP.RDS"))
  ridge_All_EHR_NLP = readRDS(here("Data", "PDDSModels", "lookbackAll_KGclean_ridge_logit_tune_EHR_NLP.RDS"))
  
# Analytic wrangling #### 
  source(here("Code", "analysis_scripts", "Sensitivity", "CompleteCase", "00_completecase_causalcore.R"))
  source(here("Code", "analysis_scripts", "01_analytic_wrangling.R"))
  source(here("Code", "analysis_scripts", "01_analytic_wrangling_postDMT.R"))
  
  
  if(.outcome_yrs==2){
    .df = bcd_ntz_2yr_joint %>% 
      filter(!is.na(.data[[.outcome.df]]))
  }else{
    .df = bcd_ntz_3yr_joint %>% 
      filter(!is.na(.data[[.outcome.df]]))
  }
  .look_yrs = 1

  
  analytic_object <- analytic_wrangling(.analytic_df = .df
                                        , pre_process = T
                                        , efficacy="high"
                                        , .EHR_trim = 0.90 # also NLP trim
                                        , .NLP = T
                                        , .NLP_all = F
                                        , .EHR_all = F
                                        , .RXNORM_Exclude_All = T
                                        , .DMT_reference = DMT_reference
                                        , .diff_window = dmonths(6)
                                        , .lookback_window = Inf
                                        )
  
  # ATE Estimate 
  set.seed(239579235) 
  ATE <- AIPW_robust_analytic(analytic_object = analytic_object
                              , .outcome = .outcome
                              , pweights = rep(1, nrow(analytic_object$X))
                              )

  .results = function(x) {
    .ate = ifelse("ATE" %in% names(x), x$ATE, x$ATE_full)
    cat("ATE" = round(.ate, 4), "\n Crump Labeled =", ATE$Crump.Labeled.Removed, "Unlabeled =", ATE$Crump.Unlabeled.Removed)
    }
  .results(ATE)
  
  # Perturbed Estimates
  doParallel::registerDoParallel()
  
  start = Sys.time() 
  # .n_perturb = 10
  ATE_perturb = foreach(icount(.n_perturb), .errorhandling="stop") %dorng% {
         AIPW_robust_analytic(analytic_object = analytic_object
                                , .outcome = .outcome
                                , pweights = 4*rbeta(nrow(analytic_object$X), 0.5, 1.5)
                                )
    }
  end = Sys.time()

  
  # RNGkind("L'Ecuyer-CMRG")
  # .Random.seed = attr(ATE_perturb, "rng")[which(is.na(unlist(lapply(ATE_perturb, `[[`, 1))))][[1]]
  # 
  # tmp = AIPW_robust_analytic(analytic_object = analytic_object
  #                      , .outcome = .outcome
  #                      , pweights = 4*rbeta(nrow(analytic_object$X), 0.5, 1.5)
  # )
  # 
  # .results(tmp)  

  saveRDS(object = list(ATE_full = ATE$ATE # ATE[names(ATE)!="ImputeModel"]
                        , Outcome_ImputeModel_Coefs = list(Trt = ATE$Trt.coef, OC = ATE$OC.coef)
                        , FullModel_IntOnly = c(Trt_KS_IntOnly = ATE$Trt_KS_IntOnly, Oucome_KS_IntOnly = ATE$Outcome_KS_IntOnly)
                        , DiPS = ATE$DiPS.notrim 
                        , ATE_Perturbed = unlist(lapply(ATE_perturb, `[[`, 1))
                        , ATE_Trt_KS_IntOnly = unlist(lapply(ATE_perturb, `[[`, 2))
                        , ATE_Outcome_KS_IntOnly = unlist(lapply(ATE_perturb, `[[`, 3))
                        # , ATE_PerturbWeights = lapply(ATE_perturb, `[`, 2)
                        , seeds = attr(ATE_perturb, "rng")
                        , Runtime_Real = end-start
                        , N_Cores_Used = getDoParWorkers()
                        , Cohort_Yrs = .outcome_yrs
                        )
          , file = here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", paste0(.prefix, "_ATE_BCD_NTZ.RDS")) 
          )
  
  message(paste0("Output stored at Sensitivity/CompleteCase/", .prefix, "_ATE_BCD_NTZ.RDS"))