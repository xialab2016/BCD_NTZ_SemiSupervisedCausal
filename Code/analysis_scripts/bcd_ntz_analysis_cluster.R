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

if("2yr" %in% bash_args){
  .outcome_yrs = 2
}else{
  .outcome_yrs = 3
}


# .prefix = paste0("perturb_1300_", .outcome_yrs, "yr_Increase")  
.prefix = bash_args[3]
if(is.na(.prefix)) .prefix = paste0("perturb_ATE_", .outcome_yrs, "yrs_Increase_"
                                    , gsub("-|[ ]|:|\\.", "_", lubridate::now())
                                    )

# Parameters #### 
  # for any given DMT, either the construct boolean must be T/TRUE/1 or a string filename supplied  
  # if constructing cohorts, PLEASE edit construction 
  # filenames presented without extensions

  .results_subdir <- "BCD_NTZ"

  # .bcd_filename_3yr <- "bcd_allMS_3yr_2024_10_07.csv" #"bcd_all_3yr_2023_12_20.csv"
  # .bcd_filename_2yr <- "bcd_allMS_2yr_2024_10_07.csv" #"bcd_all_2yr_2023_12_20.csv"
  # 
  # .ntz_filename_3yr <- "ntz_allMS_3yr_2024_10_07.csv" #"ntz_all_3yr_2023_12_20.csv"
  # .ntz_filename_2yr <- "ntz_allMS_2yr_2024_10_07.csv" #"ntz_all_2yr_2023_12_20.csv"

  .bcd_filename_3yr <- "bcd_allMS_3yr_2025_02_05.csv" #"bcd_all_3yr_2023_12_20.csv"
  .bcd_filename_2yr <- "bcd_allMS_2yr_2025_02_05.csv" #"bcd_all_2yr_2023_12_20.csv"

  .ntz_filename_3yr <- "ntz_allMS_3yr_2025_02_05.csv" #"ntz_all_3yr_2023_12_20.csv"
  .ntz_filename_2yr <- "ntz_allMS_2yr_2025_02_05.csv" #"ntz_all_2yr_2023_12_20.csv"


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
  
  # write.csv(bcd_ntz_3yr_joint, here("Data", "created_data", "BCD_NTZ", "bcd_ntz_3yr_joint.csv"))
  # write.csv(bcd_ntz_2yr_joint, here("Data", "created_data", "BCD_NTZ", "bcd_ntz_2yr_joint.csv"))
  
  # tmp = read_excel(here("Data", "Gold Standard Labels", "PRT Dx 2024-02-09.xlsx"))
  
  
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
  
  if("2yr" %in% bash_args){
    .df = bcd_ntz_2yr_joint
    .outcome_yrs = 2
  }else{
    .df = bcd_ntz_3yr_joint
    .outcome_yrs = 3
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
  
  
  X_lookforward = analytic_wrangling_postDMT(object = .df
                                             , ehr_ft_list = analytic_object[["EHR_fts"]]
                                             , nlp_ft_list = analytic_object[["NLP_fts"]]
                                             , pre_process = T
                                             , .lag = lubridate::dyears(.look_yrs) 
                                             ) # my own identifier warning, taken care of in-house (in-function)
  
  # ATE Estimate 
  RNGkind("L'Ecuyer")
  set.seed(239579235) 
  ATE <- AIPW_robust_analytic(analytic_object = analytic_object
                              , .analytic_object_1yr = X_lookforward # misnomer, can take any lookforward, I just didn't bother to rename and ensure correctness
                              , .outcome = "Y_change"
                              # , pweights = NULL
                              , pweights = rep(1, nrow(analytic_object$X))
                              , .lookforward_window = dyears(.look_yrs)
                              , cv.outcome = T
                              )

  .results = function(x) {cat("ATE" = round(x$ATE, 4), "\nAccuracy = ", round(x$Impute_Acc, 4), "\nAUC = ", round(x$Impute_AUC, 4)
                              , "\nCrump Removed: Labeled =", x$Crump.Labeled.Removed, "Unlabeled =", x$Crump.Unlabeled.Removed)   }
  .results(ATE)
  
  # Perturbed Estimates
  doParallel::registerDoParallel()
  
  start = Sys.time() 
  # .n_perturb = 5
  ATE_perturb = foreach(icount(.n_perturb), .errorhandling="pass") %dorng% {
         AIPW_robust_analytic(analytic_object = analytic_object
                                , .analytic_object_1yr = X_lookforward
                                , .outcome = "Y_change"
                                , pweights = 4*rbeta(nrow(analytic_object$X), 0.5, 1.5)
                                , .lookforward_window = dyears(.look_yrs)
                                , cv.outcome = F
                                )
    }
  end = Sys.time()

  # summary(unlist(lapply(ATE_perturb, `[[`, 1)))
  # quantile(unlist(lapply(ATE_perturb, `[[`, 1)), c(0.025, 0.975))
  
  saveRDS(object = list(ATE_full = ATE$ATE # ATE[names(ATE)!="ImputeModel"]
                        , Outcome_ImputeModel = ATE$ImputeModel
                        , Outcome_ImputeModel_Coefs = list(Trt = ATE$Trt.coef, OC = ATE$OC.coef)
                        , Outcome_ImputeModel_Perf = c(Accuracy = ATE$Impute_Acc, AUC = ATE$Impute_AUC)
                        , FullModel_IntOnly = c(Trt_KS_IntOnly = ATE$Trt_KS_IntOnly, Oucome_KS_IntOnly = ATE$Outcome_KS_IntOnly)
                        , Crump.Labeled.Removed = ATE$Crump.Labeled.Removed
                        , Crump.Unlabeled.Removed = ATE$Crump.Unlabeled.Removed
                        , DiPS = ATE$DiPS.notrim 
                        , Y_dag = ATE$Y_dag 
                        , ATE_Perturbed = unlist(lapply(ATE_perturb, `[[`, 1))
                        , ATE_Trt_KS_IntOnly = unlist(lapply(ATE_perturb, `[[`, 2))
                        , ATE_Outcome_KS_IntOnly = unlist(lapply(ATE_perturb, `[[`, 3))
                        # , ATE_PerturbWeights = lapply(ATE_perturb, `[`, 2)
                        , seeds = attr(ATE_perturb, "rng")
                        , Runtime_Real = end-start
                        , N_Cores_Used = getDoParWorkers()
                        , Cohort_Yrs = .outcome_yrs
                        )
          , file = here("Results", "BCD_NTZ", "ATE", paste0(.prefix, "_ATE_BCD_NTZ.RDS")) 
          )
  

  
# # Useful code to reference stored objects
  # # Example working with random seed from %dorng% to regenerate ATE and weights 
  # .Random.seed = attr(ATE_perturb, "rng")[[9]]
  # tmp = 4*rbeta(nrow(analytic_object$X), 0.5, 1.5)
  # identical(tmp, ATE_perturb[[9]][[1]])
  # x = AIPW_robust_analytic(analytic_object = analytic_object
  #                      , .analytic_object_1yr = X_lookforward
  #                      , .outcome = "Y_change"
  #                      , pweights = tmp 
  #                      , .lookforward_window = dyears(.outcome_yrs)
  #                      , cv.outcome = F
  #                      )
  # x
  # ATE_perturb[[9]][[2]]
  # 
  
  
# # Example workthrough with model saved in ATE object   
  # 
  # Y_test_class = predict(ATE$ImputeModel$Model
  #                        , newx = as.matrix(ATE$ImputeModel$data[ATE$ImputeModel$test_outcome_indx,])
  #                        , s="lambda.min", type = "class")
  # Y_test_prob = predict(ATE$ImputeModel$Model
  #                        , newx = as.matrix(ATE$ImputeModel$data[ATE$ImputeModel$test_outcome_indx,])
  #                        , s="lambda.min", type = "response")
  # 
  # test_acc = sum(Y_test_class == ATE$ImputeModel$outcome[ATE$ImputeModel$test_outcome_indx]) / length(ATE$ImputeModel$test_outcome_indx)
  # 
  # test_auc = survival::concordancefit(y = ATE$ImputeModel$outcome[ATE$ImputeModel$test_outcome_indx]
  #                                     , x = Y_test_prob)$concordance
