# PDDS Modelling
### Dominic DiSanto

# Set-Up #### 
library(lubridate) 

.prefix = "lookback12mo_KGclean_"
.lookback_months = dmonths(12)
setwd("/n/data1/hsph/biostat/celehs/lab/SHARE/PDDS_Imputation")

bash_args = commandArgs(trailingOnly = T)
cat("Bash Arguments (link, alpha[, .nlambda]):", bash_args)
.link_fn = bash_args[1]
.alpha = as.numeric(bash_args[2])
if(is.na(.alpha)){
  stop("No .alpha supplied")
}else if(.alpha==0){
  .reg = "ridge"
}else if(.alpha==1){
  .reg = "lasso"
}else if(!is.na(.alpha)){
  paste0("elNet", .alpha)
}

.nlambda = 20 # package default

if("pacman" %in% installed.packages()[,1]){
  library(pacman)
}else{
  install.packages("pacman")
  library(pacman)
}

p_load(here, caret
       , corrplot 
       , rms, ordinal, ordinalNet # linear methods
       , psych # wkappa 
       , DescTools # SomersDelta
       , readr
)



select <- dplyr::select # Resolve conflict with MASS ridge function
summarize <- dplyr::summarize # Resolve conflict with Hmisc
union = dplyr::union # silent conflict with lubridate
`%nin%` = Negate(`%in%`)

source(here("Code", "analytic_wrangling.R"))

cohort <- read_csv(here("Data", "created_data", "pdds_cohort_12molookback_2023_12_20.csv")
                   , col_types = cols(PATIENT_NUM = "character"))  %>%
  mutate(score = case_when(score==8 ~ 7
                           , T ~ score))


# possible removal of causal patients 
if("keepcausal" %nin% tolower(bash_args)){
  .txt = "Removing PDDS scores used in causal analysis from imputation model"
  print(.txt)
  message(.txt)
  
  .bcd_filename_3yr <- "bcd_allMS_3yr_2024_03_07.csv" 
  .bcd_filename_2yr <- "bcd_allMS_2yr_2024_03_07.csv" 
  
  .ntz_filename_3yr <- "ntz_allMS_3yr_2024_03_07.csv" 
  .ntz_filename_2yr <- "ntz_allMS_2yr_2024_03_07.csv" 
  
  bcd_3yr <- read_csv(here("Data", "created_data", "causal_data", .bcd_filename_3yr)) %>% 
    filter(!is.na(PDDS_Change))
  
  bcd_2yr <- read_csv(here("Data", "created_data", "causal_data", .bcd_filename_2yr)) %>% 
  filter(!is.na(PDDS_Change))
  
  ntz_3yr <- read_csv(here("Data", "created_data", "causal_data", .ntz_filename_3yr)) %>% 
  filter(!is.na(PDDS_Change))
  
  ntz_2yr <- read_csv(here("Data", "created_data", "causal_data", .ntz_filename_2yr)) %>% 
  filter(!is.na(PDDS_Change))
  
  .remove_df = bcd_3yr %>% union(bcd_2yr) %>% union(ntz_3yr) %>% union(ntz_2yr) %>% 
    filter(!is.na(PDDS_Change)) %>% 
    mutate(PATIENT_NUM = as.character(PATIENT_NUM)
           , RemoveTag = 1) %>% 
    distinct(PATIENT_NUM, DMT_Study_Start, RemoveTag)
  
  .causalIDs = unique(c(bcd_3yr$PATIENT_NUM, ntz_3yr$PATIENT_NUM
                        , bcd_2yr$PATIENT_NUM, ntz_2yr$PATIENT_NUM)
                      )
  
  cohortRemoved = .countdf %>% 
    filter(is.na(RemoveTag) | (RemoveTag == 1 & DMT_Study_Start >= date))
  
  cohort = cohortRemoved
}



analytic_output_KG <- analytic_wrangling(.analytic_df = cohort
                                         , pre_process = T
                                         , .EHR_trim = 0.9
                                         , .NLP = F
                                         , .EHR_all = F
                                         , .NLP_all = F
                                         , .import_NLP = F
                                           , .lookback_window = .lookback_months
                                           )

analytic_output_KG_ONCE <- analytic_wrangling(.analytic_df = cohort
                                              , pre_process = T
                                              , .EHR_trim = 0.9
                                              , .NLP = T
                                              , .EHR_all = F
                                              , .NLP_all = F
                                              , .import_NLP = F
                                              , .lookback_window = .lookback_months
                                              )


analytic_output_EHR <- analytic_wrangling(.analytic_df = cohort
                                              , pre_process = T
                                              , .EHR_trim = 0.9
                                              , .NLP = F
                                              , .EHR_all = T
                                              , .NLP_all = F
                                              , .import_NLP = F
                                              , .lookback_window = .lookback_months
                                              )


analytic_output_EHR_NLP <- analytic_wrangling(.analytic_df = cohort
                                              , pre_process = T
                                              , .EHR_trim = 0.9
                                              , .NLP = T
                                              , .EHR_all = T
                                              , .NLP_all = T
                                              , .import_NLP = F
                                              , .lookback_window = .lookback_months
                                              )


# Model Fitting ####

## Test/Train Split ####
set.seed(826178)

cv_ids = sample(unique(analytic_output_KG$analytic_df$PATIENT_NUM)
                , size = 0.8*length(unique(analytic_output_KG$analytic_df$PATIENT_NUM))
                , replace = F)

cv_indx = analytic_output_KG$analytic_df$PATIENT_NUM %in% cv_ids
test_indx = analytic_output_KG$analytic_df$PATIENT_NUM %nin% cv_ids

# ## KG ####
if("KG" %in% bash_args){
  print("KG"); message("KG")
  .tmp_df <- as.data.frame(cbind(analytic_output_KG$X[cv_indx,]
                                 , score = analytic_output_KG$Y[cv_indx])
  )

  .tmp_x <- subset(analytic_output_KG$X, select = -c(Disease_Subtype_1, Disease_Subtype_2)) %>%
    mutate_if(is.character, as.factor) %>%
    mutate_if(is.factor, as.numeric) %>%
    as.matrix()


  ordinalNet_KG <- ordinalNetTune(x = .tmp_x[cv_indx,]
                              , y = as.factor(analytic_output_KG$Y[cv_indx])
                              , alpha = .alpha
                              , nLambda = .nlambda
                              , family = "cumulative"
                              , link = .link_fn
                              , nFolds = 10
                              )


  saveRDS(ordinalNet_KG, here("Results", "Models", paste0(.prefix, .reg, "_", .link_fn, "_tune_KG.RDS")))
  saveRDS(analytic_output_KG, here("Results", "Models", paste0(.prefix, .reg, "_", .link_fn, "_analytic_output_KG.RDS")))
}

# ## KG+ONCE ####
if ("KG_ONCE" %in% bash_args){
  print("KG_ONCE"); message("KG_ONCE")
  .tmp_df <- as.data.frame(cbind(analytic_output_KG_ONCE$X[cv_indx,]
                                 , score = analytic_output_KG_ONCE$Y[cv_indx])
  )

  .tmp_x <- subset(analytic_output_KG_ONCE$X, select = -c(Disease_Subtype_1, Disease_Subtype_2)) %>%
    mutate_if(is.character, as.factor) %>%
    mutate_if(is.factor, as.numeric) %>%
    as.matrix()

  ordinalNet_KG_ONCE <- ordinalNetTune(x = .tmp_x[cv_indx,]
                              , y = as.factor(analytic_output_KG_ONCE$Y[cv_indx])
                              , alpha = .alpha
                              , nLambda = .nlambda
                              , family = "cumulative"
                              , link = .link_fn
                              , nFolds = 10
                              )

  saveRDS(ordinalNet_KG_ONCE, here("Results", "Models", paste0(.prefix, .reg, "_", .link_fn, "_tune_KG_ONCE.RDS")))
  saveRDS(analytic_output_KG_ONCE, here("Results", "Models", paste0(.prefix, .reg, "_", .link_fn, "_analytic_output_KG_ONCE.RDS")))
}


# EHR ####
if("EHR" %in% bash_args){
  print("EHR"); message("EHR")
  .tmp_df <- as.data.frame(cbind(analytic_output_EHR$X[cv_indx,]
                                 , score = analytic_output_EHR$Y[cv_indx])
  )

  .tmp_x <- subset(analytic_output_EHR$X, select = -c(Disease_Subtype_1, Disease_Subtype_2)) %>%
    mutate_if(is.character, as.factor) %>%
    mutate_if(is.factor, as.numeric) %>%
    as.matrix()

  ordinalNet_EHR <- ordinalNetTune(x = .tmp_x[cv_indx,]
                              , y = as.factor(analytic_output_EHR$Y[cv_indx])
                              , alpha = .alpha
                              , nLambda = .nlambda
                              , family = "cumulative"
                              , link = .link_fn
                              , nFolds = 10
                              )

  saveRDS(ordinalNet_EHR, here("Results", "Models", paste0(.prefix, .reg, "_", .link_fn, "_tune_EHR.RDS")))
  saveRDS(analytic_output_EHR, here("Results", "Models", paste0(.prefix, .reg, "_", .link_fn, "_analytic_output_EHR.RDS")))
}


## EHR+NLP ####
if("EHR_NLP" %in% bash_args){
  print("EHR_NLP"); message("EHR_NLP")
  .tmp_df <- as.data.frame(cbind(analytic_output_EHR_NLP$X[cv_indx,]
                                 , score = analytic_output_EHR_NLP$Y[cv_indx])
  )

  .tmp_x <- subset(analytic_output_EHR_NLP$X, select = -c(Disease_Subtype_1, Disease_Subtype_2)) %>%
    mutate_if(is.character, as.factor) %>%
    mutate_if(is.factor, as.numeric) %>%
    as.matrix()

  ordinalNet_EHR_NLP <- ordinalNetTune(x = .tmp_x[cv_indx,]
                                   , y = as.factor(analytic_output_EHR_NLP$Y[cv_indx])
                                   , alpha = .alpha
                                   , nLambda = .nlambda
                                   , family = "cumulative"
                                   , link = .link_fn
                                   , nFolds = 10
                                   )
  print("EHR_NLP model completed")
  saveRDS(ordinalNet_EHR_NLP, here("Results", "Models", paste0(.prefix, .reg, "_", .link_fn, "_tune_EHR_NLP.RDS")))
  saveRDS(analytic_output_EHR_NLP, here("Results", "Models", paste0(.prefix, .reg, "_", .link_fn, "_analytic_output_EHR_NLP.RDS")))
}

