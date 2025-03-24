
# Set-Up #### 
  library(pacman)
  p_load(here, ggplot2, Matrix, dplyr, RColorBrewer, cowplot, scales
         , forcats)
  old <- options(pillar.sigfig = 5)
  .save = F
  
  model = readRDS(here("Results", "BCD_NTZ", "ATE", "SI_3yr_FullModel_Coefs.RDS"))
  model$OC.coef
  
  coef.df = tibble(Coefficient = model$OC.coef
         , Variable = names(model$OC.coef)) %>% mutate(Model = "Outcome") %>% 
    union_all(
          tibble(Coefficient = model$Trt.coef
                 , Variable = names(model$Trt.coef)
                 ) %>% mutate(Model = "Treatment")
          ) %>% 
    filter(Variable!="(Intercept)") %>% 
    group_by(Model) %>% 
    mutate(Coefficient.Std = Coefficient / max(abs(Coefficient)))
  
  .wilcox.fn = function(x, data, group.var){
    .fm = as.formula(paste0(x, "~", group.var))
    data.frame(var=x
               , p.val=wilcox.test(.fm, data)$p.value
    )
  } 
  
  .chisq.fn = function(x, data, group.var){
    if(anyNA(data$x)){
      message("Removing missing observations from ")
      data = data %>% filter(!is.na() & !is.na())
    }
    data.frame(var = x
               , p.val =  chisq.test(data[[x]], data[[group.var]])[['p.value']]
    )
  }
  
  
## Data Import ####

# .bcd_filename_3yr <- "bcd_allMS_3yr_2024_10_07.csv" 
# .bcd_filename_2yr <- "bcd_allMS_2yr_2024_10_07.csv" 
# 
# .ntz_filename_3yr <- "ntz_allMS_3yr_2024_10_07.csv"   
# .ntz_filename_2yr <- "ntz_allMS_2yr_2024_10_07.csv"  

.bcd_filename_3yr <- "bcd_allMS_3yr_2025_02_05.csv" 
.bcd_filename_2yr <- "bcd_allMS_2yr_2025_02_05.csv" 

.ntz_filename_3yr <- "ntz_allMS_3yr_2025_02_05.csv"   
.ntz_filename_2yr <- "ntz_allMS_2yr_2025_02_05.csv"  
  
  
.ehr_cor_threshold <- 0.1
DMT_reference <- "NTZ" # set to 0 in treatment assignment vector


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


## Analytic Wrangling #####
bcd_ntz_3yr_joint <- bcd_ntz_3yr %>% 
  inner_join(bcd_ntz_2yr[,"PATIENT_NUM"], "PATIENT_NUM") %>% 
  mutate(`PDDS Outcome` = case_when(PDDS_Change==1 ~ "Increase/Progression"
                                    , PDDS_NegChange==1 ~ "Decrease/Improvement"
                                    , PDDS_Change==0 & PDDS_NegChange ==0 ~ "No Sustained Change"
                                    , T ~ "Insufficient Data (Unlabeled)")
  )

bcd_ntz_2yr_joint <- bcd_ntz_2yr %>% 
  inner_join(bcd_ntz_3yr[,"PATIENT_NUM"], "PATIENT_NUM") %>% 
  mutate(`PDDS Outcome` = case_when(PDDS_Change==1 ~ "Increase/Progression"
                                    , PDDS_NegChange==1 ~ "Decrease/Improvement"
                                    , PDDS_Change==0 & PDDS_NegChange ==0 ~ "No Sustained Change"
                                    , T ~ "Insufficient Data (Unlabeled)")
  )


# Adding specific BCD medications for Table 1 #### 
.dmt_IDs <- bcd_ntz_2yr_joint %>% distinct(PATIENT_NUM, id_participant) %>% mutate(PATIENT_NUM=as.character(PATIENT_NUM))
drugs <- data.frame(RxNormID = c("1876366", "121191", "712566")
                    , Type_Generic = c("Ocrelizumab", "Rituximab", "Ofatumumab")
                    , Type_BrandName = c("Ocrevus", "Rituxan", "Kesimpta"))
source(here("Code", "01_dmt_cleaning.R"))

.tmp = bcd_ntz_2yr_joint %>% # equivalent using 2yr or 3yr cohort 
  left_join(ehr_rxnorm_sub %>% select(patient_num, treatment, start_date) %>% distinct() %>% 
              mutate(treatment=tolower(gsub("BCD:", "", treatment)))
            , by=c("PATIENT_NUM"="patient_num"
                   , "DMT_Study_Start"="start_date")
            )

.tmp3.reg = .tmp %>% 
  left_join(registry_dmt_sub %>% select(id_participant, treatment=type_generic, start) %>% distinct()
            , by=c("id_participant"
                   , "DMT_Study_Start"="start"
                   )
            ) %>% 
  mutate(treatment=coalesce(treatment.x, treatment.y)) %>% 
  arrange(PATIENT_NUM) %>% group_by(PATIENT_NUM) %>% 
  filter(row_number()==1) %>% 
  ungroup() %>% 
  select(-c(treatment.x, treatment.y))

# lazily doing for both 2 and 3 year, allows use of either in subsequent code 
bcd_ntz_3yr_joint = bcd_ntz_3yr_joint %>% 
  left_join(.tmp3.reg %>% select(PATIENT_NUM, treatment), "PATIENT_NUM") %>% 
  mutate(treatment = ifelse(is.na(treatment), "NTZ", paste0("BCD - ", treatment) )) 

bcd_ntz_2yr_joint = bcd_ntz_2yr_joint %>%
  left_join(.tmp3.reg %>% select(PATIENT_NUM, treatment), "PATIENT_NUM") %>% 
  mutate(treatment = ifelse(is.na(treatment), "NTZ", paste0("BCD - ", treatment) )) 

# dim(.tmp3)
# dim(bcd_ntz_3yr_joint)
# all.equal(bcd_ntz_3yr_joint[,1:27], .tmp3[,1:27])


# .tmp3.reg %>% arrange(PATIENT_NUM) %>% group_by(PATIENT_NUM) %>% filter(n()>1) %>% select(PATIENT_NUM, DMT_Study_Start, treatment)
  # one patient has two records on same day (rituximab, ocrelizumab)

# .tmp3.reg %>% filter(DMT=="BCD") %>% count(treatment) %>% mutate(p=n/sum(n))
# .tmp3.reg %>% filter(DMT=="NTZ") %>% count(treatment) # confirmed, all NA as expected


# Load general imputation models of interest for latent scores 
# lasso_cauchit = readRDS(here("Data", "PDDSModels", "lasso_cauchit_tune_EHR_NLP.RDS"))
# lasso_logit = readRDS(here("Data", "PDDSModels", "lasso_logit_tune_EHR_NLP.RDS"))
# ordinal_KG_NLP = readRDS(here("Data", "PDDSModels", "ordinal_cauchit_KG_ONCE.RDS"))
ridge_12mo_EHR_NLP = readRDS(here("Data", "PDDSModels", "lookback12mo_KGclean_ridge_logit_tune_EHR_NLP.RDS"))
ridge_All_EHR_NLP = readRDS(here("Data", "PDDSModels", "lookbackAll_KGclean_ridge_logit_tune_EHR_NLP.RDS"))

# Analytic wrangling # 
source(here("Code", "analysis_scripts", "00_causal_core_Robust_Perturbed.R"))
source(here("Code", "analysis_scripts", "01_analytic_wrangling.R"))
source(here("Code", "analysis_scripts", "01_analytic_wrangling_postDMT.R"))


analytic_object_3yr <- analytic_wrangling(.analytic_df = bcd_ntz_3yr_joint
                                          , pre_process = F
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

X_lookforward_3yr = analytic_wrangling_postDMT(object = bcd_ntz_3yr_joint
                                               , ehr_ft_list = analytic_object_3yr[["EHR_fts"]]
                                               , nlp_ft_list = analytic_object_3yr[["NLP_fts"]]
                                               , pre_process = F
                                               , .lag = lubridate::dyears(3) 
) # my own identifier warning, taken care of in-house (in-function)

analytic_object_2yr <- analytic_wrangling(.analytic_df = bcd_ntz_2yr_joint
                                          , pre_process = F
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

X_lookforward_2yr = analytic_wrangling_postDMT(object = bcd_ntz_2yr_joint
                                               , ehr_ft_list = analytic_object_2yr[["EHR_fts"]]
                                               , nlp_ft_list = analytic_object_2yr[["NLP_fts"]]
                                               , pre_process = F
                                               , .lag = lubridate::dyears(2) 
) # my own identifier warning, taken care of in-house (in-function)


source(here("Code", "analysis_scripts", "00_causal_core_Robust_Perturbed.R"))
ridge_12mo_EHR_NLP = readRDS(here("Data", "PDDSModels", "lookback12mo_KGclean_ridge_logit_tune_EHR_NLP.RDS"))
ridge_All_EHR_NLP = readRDS(here("Data", "PDDSModels", "lookbackAll_KGclean_ridge_logit_tune_EHR_NLP.RDS"))


## Labelling/Presentation Edits ####

pdds.baseline.12mo = pdds_impute(ridge_12mo_EHR_NLP
                                 , .analytic_X = analytic_object_2yr$X %>%  
                                   mutate(across(-c(White_NonHispanic, subject_sex_ch)
                                                 , .fns = ~ (.x-mean(.x)) / sd(.x)))
                                 , .print = F, type = "response")

pdds.baseline.All = pdds_impute(ridge_All_EHR_NLP
                                , .analytic_X = analytic_object_2yr$X %>%  
                                  mutate(across(-c(White_NonHispanic, subject_sex_ch)
                                                , .fns = ~ (.x-mean(.x)) / sd(.x)))
                                , .print = F, type = "response")


pdds.LF.12mo.2yr = pdds_impute(ridge_12mo_EHR_NLP
                               , .analytic_X = X_lookforward_2yr %>%  
                                 mutate(across(-c(White_NonHispanic, subject_sex_ch)
                                               , .fns = ~ (.x-mean(.x)) / sd(.x)))
                               , .print = F, type = "response")

pdds.LF.All.2yr = pdds_impute(ridge_All_EHR_NLP
                              , .analytic_X = X_lookforward_2yr %>%  
                                mutate(across(-c(White_NonHispanic, subject_sex_ch)
                                              , .fns = ~ (.x-mean(.x)) / sd(.x)))
                              , .print = F, type = "response")

pdds.LF.12mo.3yr = pdds_impute(ridge_12mo_EHR_NLP
                               , .analytic_X = X_lookforward_3yr %>%  
                                 mutate(across(-c(White_NonHispanic, subject_sex_ch)
                                               , .fns = ~ (.x-mean(.x)) / sd(.x)))
                               , .print = F, type = "response")

pdds.LF.All.3yr = pdds_impute(ridge_All_EHR_NLP
                              , .analytic_X = X_lookforward_3yr %>%  
                                mutate(across(-c(White_NonHispanic, subject_sex_ch)
                                              , .fns = ~ (.x-mean(.x)) / sd(.x)))
                              , .print = F, type = "response")


# KG objects 
ehr_ft_list = rbind(once_disability_ehr, once_ms_ehr)

once_disability_ehr %>% filter(tolower(substr(term, 1, 3))=="rxn") %>% arrange(term)
once_ms_ehr %>% filter(tolower(substr(term, 1, 3))=="rxn") %>% arrange(term)
sort(paste0("RXNORM:", StdEff_Exclude_RxNormIDs))

colnames(analytic_object_3yr$X)[which(substr(colnames(analytic_object_3yr$X), 1, 3)=="RXN")]

cbind(StdEff_Exclude_RxNormIDs
      , paste0("RXNORM:", StdEff_Exclude_RxNormIDs) %in% colnames(analytic_object_3yr$X)
        )



data = data.frame(PATIENT_NUM = bcd_ntz_2yr_joint$PATIENT_NUM
                  , Registry = ifelse(!is.na(bcd_ntz_2yr_joint$id_participant), 1, 0)
                  # , DMT = ifelse(analytic_object_2yr$A, "BCD", "NTZ")
                  , DMT = bcd_ntz_2yr_joint$treatment 
                  , analytic_object_2yr$X
                  , Utilization.3yr = analytic_object_3yr$X$Utilization_Total
                  # , pdds.baseline.12mo = pdds.baseline.12mo
                  , pdds.baseline.All = pdds.baseline.All
                  , pdds.LF.12mo.2yr = pdds.LF.12mo.2yr
                  , pdds.LF.All.2yr = pdds.LF.All.2yr
                  , pdds.LF.12mo.3yr = pdds.LF.12mo.3yr
                  , pdds.LF.All.3yr = pdds.LF.All.3yr
                  , SI_2yr = analytic_object_2yr$Y_change
                  , SD_2yr = analytic_object_2yr$Y_NegChange
                  , SI_3yr = analytic_object_3yr$Y_change
                  , SD_3yr = analytic_object_3yr$Y_NegChange
                  ) %>% as_tibble()

data$Years_First_MS_PheCode_to_Study_DMT <- data$Disease_Duration/365
data$Years_FirstPheCode_to_Study_DMT <- pmax(data$Followup_Duration/365, 0)

data$White_NonHispanic_factor <- factor(data$White_NonHispanic
                                        , labels = c("Non-White/Hispanic", "White and Non-Hispanic"))
data$subject_sex_ch_factor = factor(data$subject_sex_ch, labels = c("Female", "Male"))

label(data$Age_at_DMT_Initiation) <- "Age (at DMT Start)"
label(data$Years_First_MS_PheCode_to_Study_DMT) <- "Disease Duration (Years)"
label(data$Years_FirstPheCode_to_Study_DMT) <- "Follow-up Duration (Years)"
label(data$subject_sex_ch) <- "Sex"  
label(data$White_NonHispanic_factor) <- "Race & Ethnicity"  
label(data$Utilization_Total) <- "Utilization at 2-years (No. of Codified+Narrative Codes)"
label(data$Utilization.3yr) <- "Utilization at 3-years (No. of Codified+Narrative Codes)"
# label(data$pdds.baseline.12mo) = "Expected PDDS Baseline (12 mo-Lookback model)"
label(data$pdds.baseline.All) = "Expected PDDS Baseline (All-Lookback model)"
label(data$pdds.LF.12mo.2yr) = "Expected PDDS 2yr-lookforward (12 mo-Lookback model)"
label(data$pdds.LF.All.2yr) = "Expected PDDS 2yr-lookforward (All-Lookback model)"
label(data$pdds.LF.12mo.3yr) = "Expected PDDS 3yr-lookforward (12 mo-Lookback model)"
label(data$pdds.LF.All.3yr) = "Expected PDDS 3yr-lookforward (All-Lookback model)"
# label(data$Days_Earliest_DMT_to_Study_DMT) = ""

data$PDDSOutcome_2yr = bcd_ntz_2yr_joint$`PDDS Outcome`
data$PDDSOutcome_3yr = bcd_ntz_3yr_joint$`PDDS Outcome`

# (Tbl 1) Summary Table #### 


meanSD = function(x){
  m = round(mean(x), 3)
  sd = round(sd(x), 3)
  
  m = formatC(mean(x), format = 'f', flag='0', digits = 2)
  sd = formatC(sd(x), format = 'f', flag='0', digits = 2)
  
  paste0(m, " (", sd, ")")
}

npct = function(x, cat=1){
  .p = mean(x==cat)
  n = sum(x==cat)
  paste0(format(n, big.mark=",")
         , " ("
         , scales::percent(.p, accuracy = 0.01)
         , ")")
}

.tbl = as_tibble(data) %>% 
  union_all(data %>% mutate(Registry=-1)) %>% 
  mutate(Registry=case_when(
    Registry==-1 ~ "Overall"
    , Registry==0 ~ "Non-Registry"
    , Registry==1 ~ "Registry"
  )
         ) %>% 
  group_by(Registry) %>% 
  summarize(n = paste0("(n=", n(), ")") 
            , `BCD - ocrelizumab: n (%)` = npct(DMT, "BCD - ocrelizumab")
            , `BCD - ofatumumab: n (%)` = npct(DMT, "BCD - ofatumumab")
            , `BCD - rituximab: n (%)` = npct(DMT, "BCD - rituximab")
            , `NTZ: n (%)` = npct(DMT, "NTZ")
            , `Age at Target Treatment Initiation: years, mean (SD)` = meanSD(Age_at_DMT_Initiation)
            , `Disease Duration: years, mean (SD)` = meanSD(Years_First_MS_PheCode_to_Study_DMT)
            , `Follow-up Duration` = meanSD(Years_FirstPheCode_to_Study_DMT)
            , `Self-Reported Gender (Male): n (%)` = npct(subject_sex_ch)
            , `Self-Reported Gender (Female): n (%)` = npct(1-subject_sex_ch)
            , `Race/Ethnicity - White, Non-Hispanic: n (%)` = npct(White_NonHispanic)
            , `Race/Ethnicity - Non-white or Hispanic: n (%)` = npct(1-White_NonHispanic)
            , `Utilization: number of codes, mean (SD) ` = meanSD(Utilization_Total)
            # , `Baseline PDDS Risk (12 mo.)` = meanSD(pdds.baseline.12mo)
            , `Baseline PDDS Risk (All): mean (SD)` = meanSD(pdds.baseline.All)
            , `Year 1 PDDS Risk: mean (SD)` = meanSD(pdds.LF.12mo.3yr)
            , `Sustained Worsening: n (%)` = npct(PDDSOutcome_3yr, "Increase/Progression")
            , `Sustained Improvement: n (%)` = npct(PDDSOutcome_3yr, "Decrease/Improvement")
            , `No Sustained Change: n (%)` = npct(PDDSOutcome_3yr, "No Sustained Change")
            , `Unlabeled: n (%)` = npct(PDDSOutcome_3yr, "Insufficient Data (Unlabeled)")
  ) %>% 
  # t() %>% as_tibble()
  # pivot_longer(values_to = "Value", cols=everything()) %>% 
  pivot_longer(values_to = "Value", cols=n:`Unlabeled: n (%)`) %>%
  # mutate() %>% 
  arrange(name!="n", name %in% c("Unlabeled: n (%)"
                                  , "Sustained Improvement: n (%)"
                                  , "Sustained Worsening: n (%)"
                                 , "No Sustained Change: n (%)")
          , name) %>% 
  pivot_wider(id_cols = name, values_from = Value, names_from=Registry) %>% 
  select(name, Overall, Registry, `Non-Registry`)

  if(.save) write.csv(.tbl, here("Results", "BCD_NTZ", "Table_CohortSummary.csv"))

  ### P-values for Main Table 1 #### 
    data.test = as_tibble(data) %>% 
      mutate(Registry=case_when(
        Registry==0 ~ "Non-Registry"
        , Registry==1 ~ "Registry"
        )
        )

    wilcox.vars = c('Age_at_DMT_Initiation', 'Years_First_MS_PheCode_to_Study_DMT'
                    , 'Years_FirstPheCode_to_Study_DMT', 'Utilization_Total'
                    , 'pdds.baseline.All', "pdds.LF.All.3yr"
                    )
    chisq.vars = c('subject_sex_ch', 'White_NonHispanic'
                   , 'PDDSOutcome_3yr'
                   , 'DMT')
    
    .wilcox.df = Reduce(`rbind`
                        , lapply(wilcox.vars, .wilcox.fn, data=data.test, group.var="Registry")
                        )
    .chisq.df = Reduce(`rbind`
                       , lapply(chisq.vars, .chisq.fn, data=data.test, group.var="Registry")
                       )
    
    .pval.df = as_tibble(union(.wilcox.df, .chisq.df)) %>% 
      mutate(var = case_when(
        var == 'Age_at_DMT_Initiation' ~ 'Age at Target Treatment Initiation: years, mean (SD)'
        , var == 'Years_First_MS_PheCode_to_Study_DMT' ~ 'Disease Duration: years, mean (SD)'
        , var == 'Years_FirstPheCode_to_Study_DMT' ~ 'Follow-up Duration'
        , var == 'Utilization_Total' ~ 'Utilization: number of codes, mean (SD) '
        , var == 'pdds.baseline.All' ~ 'Baseline PDDS Risk (All): mean (SD)'
        , var == 'pdds.LF.All.3yr' ~ 'Year 1 PDDS Risk: mean (SD)'
        , var == 'subject_sex_ch' ~ 'Self-Reported Gender (Male): n (%)'
        , var == 'White_NonHispanic' ~ 'Race/Ethnicity - White, Non-Hispanic: n (%)'
        , var == 'DMT' ~ 'BCD - ocrelizumab: n (%)'
        , var == 'PDDSOutcome_3yr' ~ 'No Sustained Change: n (%)'
      )) %>% 
      rename(name=var) %>% 
      mutate(p.val = formatC(p.val, format = 'f', flag='0', digits = 3)) %>% 
      mutate(p.val = case_when(p.val=="0.000" ~ "<0.001"
                               , T ~ p.val)
      )
    
    
    .tbl1.pval = .tbl %>% 
      left_join(.pval.df, by=c('name')) %>% 
      mutate(order = 
               case_when(str_detect(name, 'Age at') ~ 1
                         , str_detect(name, 'Male') ~ 2
                         , str_detect(name, 'Female') ~ 3
                         , str_detect(name, 'White,') ~ 4
                         , str_detect(name, 'Non-white') ~ 5
                         , str_detect(name, 'Disease Duration') ~ 6
                         , str_detect(name, 'Follow-up') ~ 7
                         , str_detect(name, "BCD") ~ 7.2 
                         , str_detect(name, "NTZ") ~ 7.4
                         , str_detect(name, "Utilization") ~ 8
                         , str_detect(name, 'Baseline') ~ 9
                         , str_detect(name, 'Year 1') ~ 10
                         , str_detect(name, ': n') ~ 11
                         , T ~ Inf) 
      ) %>% 
      arrange(order)
    
    
    if(.save) write.csv(.tbl1.pval
                        , here("Results", "BCD_NTZ", "Table_CohortSummary_PValues.csv"))
    

    
# Figure by DMT #### 

  plot.df = data %>% 
      mutate(DMT = ifelse(str_detect(DMT, "BCD"), "BCD", "NTZ")) %>% 
    select(-starts_with(c("SI", "SD"))) %>% 
    union_all(data %>%
                mutate(DMT = "Overall") %>% 
                select(-starts_with(c("SI", "SD")))
              ) %>%
    group_by(DMT) %>%
    # group_by(Registry) %>%
    summarize(n = n()
              , Age = mean(Age_at_DMT_Initiation)
              , Disease_Duration = mean(Years_First_MS_PheCode_to_Study_DMT)
              , Followup_Duration = mean(Years_FirstPheCode_to_Study_DMT)
              , Female = mean(1-subject_sex_ch)
              , White_NonHispanic = mean(White_NonHispanic)
              , Utilization_Total = mean(Utilization_Total)
              , pdds.baseline.12mo = mean(pdds.baseline.12mo)
              , pdds.baseline.All = mean(pdds.baseline.All)
              , pdds.LF.12mo.2y = mean(pdds.LF.12mo.2yr)
              , pdds.LF.All.2yr = mean(pdds.LF.12mo.2yr)
              , pdds.LF.12mo.3yr = mean(pdds.LF.12mo.3yr)
              , pdds.LF.All.3yr = mean(pdds.LF.All.3yr)
              )
    
  
  plot.df[1, 3:ncol(plot.df)] = plot.df[1, 3:ncol(plot.df)] / plot.df[3, 3:ncol(plot.df)] - 1
  plot.df[2, 3:ncol(plot.df)] = plot.df[2, 3:ncol(plot.df)] / plot.df[3, 3:ncol(plot.df)] - 1 
  
  coef.df.wide = coef.df %>% 
    pivot_wider(names_from = "Model", values_from=c(Coefficient, Coefficient.Std))
  

    plt.dmt = plot.df %>% 
      select(-c(pdds.LF.All.3yr, pdds.LF.12mo.2y, pdds.LF.All.2yr)) %>% 
      pivot_longer(names_to = "Variable", values_to = "Value"
                   , cols = -c(DMT, n)) %>% 
      filter(Variable!="pdds.baseline.12mo") %>% 
      mutate(Variable = 
               case_when(Variable == "Age" ~ "Age at Target Treatment Initiation"
                         , Variable == "Disease_Duration" ~ "Disease Duration" 
                         , Variable == "Followup_Duration" ~ "Follow-up Duration" 
                         , Variable == "pdds.LF.12mo.3yr" ~ "Year 1 PDDS Risk" 
                         , Variable == "Utilization_Total" ~ "Baseline, Total Healthcare Utilization" 
                         , Variable == "White_NonHispanic" ~ "Race and Ethnicity (White, Non-Hispanic)" 
                         , Variable == "Female" ~ "Self-Reported Gender (Women)"
                         , Variable == "pdds.baseline.All" ~ "Baseline PDDS Risk" 
                         , T ~ Variable
                         )
      ) %>% 
      # inner_join(coef.df.wide, by="Variable") %>% 
      filter(DMT!="Overall") %>% 
      mutate(DMT = paste0(DMT, "\n(n=", format(n, big.mark=',', trim=T), ")")) %>% 
      mutate(Variable = 
               fct_relevel(Variable, c("Age at Target Treatment Initiation"
                                       , "Self-Reported Gender (Women)"
                                       , "Race and Ethnicity (White, Non-Hispanic)"
                                       , "Disease Duration"
                                       , "Follow-up Duration"
                                       , "Baseline, Total Healthcare Utilization"
                                       , "Baseline PDDS Risk"
                                       , "Year 1 PDDS Risk"
               )
               )
      ) %>% 
      ggplot(aes(x=DMT, y=Variable, color=Value*100)) + 
      geom_point(color="black", shape=15, size=6) + 
      geom_point(shape=15, size=5) + 
      theme_minimal() + 
      scale_color_distiller(type="seq", palette="RdBu"
                            , guide = guide_colorbar(frame.colour = "black"
                                                     , ticks.colour = "black")
                            , limits=c(-35, 35)
                            ) + 
      labs(color="Deviation % from\nOverall Mean") + 
      scale_size(guide="none") + 
      theme(legend.position="bottom", legend.key.width = unit(2, "cm")
            , axis.text.x = element_text(angle = 30, face="bold")
            ) + 
      scale_y_discrete(limits = rev)
      
# By Outcome #### 
  ## 2-year #### 
    plot.df.2yr = data %>% 
        mutate(PDDSOutcome_2yr = ifelse(is.na(PDDSOutcome_2yr), "Missing", as.character(PDDSOutcome_2yr))) %>% 
        union_all(data %>% mutate(PDDSOutcome_2yr = "Overall")) %>%
        select(-c("SI_3yr", "SD_3yr")) %>% 
        group_by(PDDSOutcome_2yr) %>%
        # group_by(Registry) %>%
        summarize(n = n()
                  , Age = mean(Age_at_DMT_Initiation)
                  , Disease_Duration = mean(Years_First_MS_PheCode_to_Study_DMT)
                  , Followup_Duration = mean(Years_FirstPheCode_to_Study_DMT)
                  , Female = mean(1-subject_sex_ch)
                  , White_NonHispanic = mean(White_NonHispanic)
                  , Utilization_Total = mean(Utilization_Total)
                  , pdds.baseline.All = mean(pdds.baseline.All)
                  , pdds.baseline.12mo = mean(pdds.baseline.12mo)
                  # , pdds.LF.All.3yr = mean(pdds.LF.All.3yr)
                  , pdds.LF.12mo.3yr = mean(pdds.LF.12mo.3yr)
        )
      
      overall = plot.df.2yr %>% filter(PDDSOutcome_2yr=="Overall")
      plot.df.2yr[,3:ncol(plot.df.2yr)] = (plot.df.2yr[,3:ncol(plot.df.2yr)] / slice(overall[,3:ncol(overall)], rep(1, 5))) - 1
      
      plot.2yr = plot.df.2yr %>% 
        mutate(PDDSOutcome_2yr = 
                 case_when(str_detect(PDDSOutcome_2yr, "Improvement") ~ "Improvement"
                           , str_detect(PDDSOutcome_2yr, "Progression") ~ "Progression"
                           , str_detect(PDDSOutcome_2yr, "Unlabel") ~ "Unlabeled"
                           , str_detect(PDDSOutcome_2yr, "No Sustained") ~ "No Change"
                 ) 
        ) %>% 
        filter(PDDSOutcome_2yr!="Overall") %>%
        mutate(PDDSOutcome_2yr = paste0(PDDSOutcome_2yr, "\n(n=", format(n, big.mark=',', trim=T), ")")) %>% 
        rename(`Race and Ethnicity (White, Non-Hispanic)` = White_NonHispanic
               , `Age at Target Treatment Initiation` = Age
               , `Self-Reported Gender (Women)` = Female
               , `Utilization` = Utilization_Total
               , `Follow-up Duration` = Followup_Duration
               , `Disease Duration` = Disease_Duration
               , `Baseline PDDS Risk` = pdds.baseline.All
               , `Year 1 PDDS Risk` = pdds.LF.12mo.3yr
               , `Baseline, Total Healthcare Utilization` = Utilization_Total
               ) %>% 
        pivot_longer(names_to = "Variable", values_to = "Value"
                     , cols = `Age at Target Treatment Initiation`:`Year 1 PDDS Risk`) %>% 
        mutate(Variable = factor(Variable, levels = sort(unique(Variable)))) %>% 
        # inner_join(coef.df.wide, by="Variable") %>% 
        mutate(Variable = 
                 fct_relevel(Variable, c("Age at Target Treatment Initiation"
                                         , "Self-Reported Gender (Women)"
                                         , "Race and Ethnicity (White, Non-Hispanic)"
                                         , "Disease Duration"
                                         , "Follow-up Duration"
                                         , "Baseline, Total Healthcare Utilization"
                                         , "Baseline PDDS Risk"
                                         , "Year 1 PDDS Risk"
                 )
                 )
        ) %>% 
        ggplot(aes(x=PDDSOutcome_2yr, y=Variable, color=Value*100)) + 
        geom_point(color="black", shape=15, size=7) + 
        geom_point(shape=15, size=6) + 
        theme_minimal() + 
        scale_color_distiller(type="seq", palette="RdBu", limits = c(-30, 30), breaks = seq(-30, 30, 10)) + 
        labs(color="Deviation % from\nOverall Mean") + 
        scale_size(guide="none") +  
        scale_y_discrete(limits = rev)

      
    ## 3-year #### 
      plot.df.3yr = data %>% 
        mutate(PDDSOutcome_3yr = ifelse(is.na(PDDSOutcome_3yr), "Missing", as.character(PDDSOutcome_3yr))) %>% 
        union_all(data %>% mutate(PDDSOutcome_3yr = "Overall")) %>%
        select(-c("SI_2yr", "SD_2yr")) %>% 
        group_by(PDDSOutcome_3yr) %>%
        # group_by(Registry) %>%
        summarize(n = n()
                  , Age = mean(Age_at_DMT_Initiation)
                  , Disease_Duration = mean(Years_First_MS_PheCode_to_Study_DMT)
                  , Followup_Duration = mean(Years_FirstPheCode_to_Study_DMT)
                  , Female = mean(1-subject_sex_ch)
                  , White_NonHispanic = mean(White_NonHispanic)
                  , Utilization_Total = mean(Utilization_Total)
                  , pdds.baseline.All = mean(pdds.baseline.All)
                  , pdds.baseline.12mo = mean(pdds.baseline.12mo)
                  # , pdds.LF.All.3yr = mean(pdds.LF.All.3yr)
                  , pdds.LF.12mo.3yr = mean(pdds.LF.12mo.3yr)
        )
      
      overall = plot.df.3yr %>% filter(PDDSOutcome_3yr=="Overall")
      plot.df.3yr[,3:ncol(plot.df.3yr)] = (plot.df.3yr[,3:ncol(plot.df.3yr)] / slice(overall[,3:ncol(plot.df.3yr)], rep(1, 5))) - 1
      
      plot.3yr = plot.df.3yr %>% 
        filter(PDDSOutcome_3yr!="Overall") %>%
        mutate(PDDSOutcome_3yr = 
                 case_when(str_detect(PDDSOutcome_3yr, "Improvement") ~ "Sustained\nImprovement"
                           , str_detect(PDDSOutcome_3yr, "Progression") ~ "Sustained\nWorsening"
                           , str_detect(PDDSOutcome_3yr, "Unlabel") ~ "Unlabelled"
                           , str_detect(PDDSOutcome_3yr, "No Sustained") ~ "No Sustained\nChange"
                           ) 
               ) %>% 
        mutate(PDDSOutcome_3yr = paste0(PDDSOutcome_3yr, "\n(n=", format(n, big.mark=',', trim=T), ")")) %>% 
        select(-'pdds.baseline.12mo') %>% 
        rename(`Race and Ethnicity (White, Non-Hispanic)` = White_NonHispanic
               , `Age at Target Treatment Initiation` = Age
               , `Self-Reported Gender (Women)` = Female
               , `Utilization` = Utilization_Total
               , `Follow-up Duration` = Followup_Duration
               , `Disease Duration` = Disease_Duration
               , `Baseline PDDS Risk` = pdds.baseline.All
               , `Year 1 PDDS Risk` = pdds.LF.12mo.3yr
               , `Baseline, Total Healthcare Utilization` = Utilization_Total
        ) %>% 
        pivot_longer(names_to = "Variable", values_to = "Value"
                     , cols = `Age at Target Treatment Initiation`:`Year 1 PDDS Risk`) %>% 
        mutate(Variable = 
                 fct_relevel(Variable, c("Age at Target Treatment Initiation"
                                              , "Self-Reported Gender (Women)"
                                              , "Race and Ethnicity (White, Non-Hispanic)"
                                              , "Disease Duration"
                                              , "Follow-up Duration"
                                              , "Baseline, Total Healthcare Utilization"
                                              , "Baseline PDDS Risk"
                                              , "Year 1 PDDS Risk"
                                              )
                             )
                 ) %>% 
        ggplot(aes(x=PDDSOutcome_3yr, y=Variable, color=Value*100)) + 
        geom_point(color="black", shape=15, size=7) + 
        geom_point(shape=15, size=6) + 
        theme_minimal() + 
        scale_color_distiller(type="seq", palette="RdBu", limits = c(-30, 30), breaks = seq(-30, 30, 10)) + 
        labs(color="Deviation % from\nOverall Mean") + 
        scale_size(guide="none") + 
        scale_y_discrete(limits = rev)
      
# (Fig 2)  Combined Plots #### 
  all.plts = plot_grid(plt.dmt + theme(legend.position = "none" # plot.3yr
                                   , axis.text.x = element_text(face="bold")
                                   ) +
                  xlab("") + ylab("") + 
                    ggtitle("Treatment") + 
                  theme(axis.text.y = element_text(face="bold")
                        , axis.text.x = element_text(vjust=0.5, face="bold")
                        , plot.title = element_text(face="bold", hjust=0.5)
                        )
                  
                , plot.3yr + ggtitle("Outcome") + # plt.dmt
                  theme(legend.position = "right"
                        , legend.key.width = unit(0.5, "cm")
                        , legend.key.height = unit(1.75, "cm")
                        , legend.text = element_text(face="bold")
                        , legend.title = element_text(face="bold")
                        , axis.text.x = element_text(face="bold")
                        , plot.title = element_text(face="bold", hjust=0.5)
                        ) + xlab("") + ylab("") +
                  theme(axis.text.y = element_blank()
                        , axis.text.x = element_text(angle = 30, vjust=0.5, face="bold")
                        , legend.title= element_text(hjust=0.7))
                , ncol = 2
                , nrow = 1
                , rel_widths =c(0.85, 1)
                , align = "h"
                , axis = "btl"
                )
      
  if(.save) cowplot::save_plot(all.plts
                               , filename = here("Results", "BCD_NTZ", "Fig2_AllPlots_version.png")
                               , base_height=6
                               , base_asp=1.5
                               )    
  
  

# Supplemental tables ####

  ## (4a) Year-2 ####
      .tbl.s2 = data %>% 
        mutate(DMT = ifelse(str_detect(DMT, "BCD"), "BCD", "NTZ")) %>% 
        group_by(PDDSOutcome_2yr) %>% 
        summarize(n = as.character(n())
                  , `BCD: n (%)` = npct(DMT, "BCD")
                  , `NTZ: n (%)` = npct(DMT, "NTZ")
                  , `Age at Target Treatment Initiation: years, mean (SD)` = meanSD(Age_at_DMT_Initiation)
                  , `Disease Duration: years, mean (SD)` = meanSD(Years_First_MS_PheCode_to_Study_DMT)
                  , `Follow-up Duration` = meanSD(Years_FirstPheCode_to_Study_DMT)
                  , `Self-Reported Gender (Male): n (%)` = npct(subject_sex_ch)
                  , `Self-Reported Gender (Female): n (%)` = npct(1-subject_sex_ch)
                  , `Race/Ethnicity - White, Non-Hispanic: n (%)` = npct(White_NonHispanic)
                  , `Race/Ethnicity - Non-white or Hispanic: n (%)` = npct(1-White_NonHispanic)
                  , `Utilization: number of codes, mean (SD) ` = meanSD(Utilization_Total)
                  # , `Baseline PDDS Risk (12 mo.)` = meanSD(pdds.baseline.12mo)
                  , `Baseline PDDS Risk (All): mean (SD)` = meanSD(pdds.baseline.All)
                  , `Year 1 PDDS Risk: mean (SD)` = meanSD(pdds.LF.12mo.3yr)

                  # n = as.character(n())
                  # , `Age at target treatment initiation, mean (SD)` = meanSD(Age_at_DMT_Initiation)
                  # , `BCD Treatment Class` = npct(DMT, "BCD")
                  # , `Disease Duration` = meanSD(Years_First_MS_PheCode_to_Study_DMT)
                  # , `Follow-up Duration` = meanSD(Years_FirstPheCode_to_Study_DMT)
                  # , `Female Sex` = percent(mean(1-subject_sex_ch), accuracy=0.01)
                  # , `White, Non-Hispanic ` = meanSD(White_NonHispanic)
                  # , Utilization = meanSD(Utilization_Total)
                  # # , `Baseline PDDS Risk (12 mo.)` = meanSD(pdds.baseline.12mo)
                  # , `Baseline PDDS Risk (All)` = meanSD(pdds.baseline.All)
                  # , `Year 1 PDDS Risk` = meanSD(pdds.LF.12mo.2yr)
        ) %>% 
        mutate(PDDSOutcome_2yr =
                 case_when(PDDSOutcome_2yr == "Decrease/Improvement" ~ "Sustained Improvement"
                           , PDDSOutcome_2yr == "Increase/Progression" ~ "Sustained Worsening"
                           , str_detect(PDDSOutcome_2yr, "Insufficient") ~ "Unlabeled"
                           , T ~ PDDSOutcome_2yr
                           )
               , PDDSOutcome_2yr = paste0(PDDSOutcome_2yr, " (n=", n, ")")
               ) %>%
        pivot_longer(values_to = "Value", cols = -PDDSOutcome_2yr) %>%
        pivot_wider(values_from = "Value", names_from = "PDDSOutcome_2yr") %>% 
        # mutate() %>% 
        arrange(name!="n", name) %>% 
        select(Variable = name, starts_with("No"), starts_with("Sustain"), starts_with("Unlabel")) %>% 
        filter(Variable!="n") %>% 
        mutate(order = 
                 case_when(str_detect(Variable, 'Age at') ~ 1
                           , str_detect(Variable, 'Male') ~ 2
                           , str_detect(Variable, 'Female') ~ 3
                           , str_detect(Variable, 'White,') ~ 4
                           , str_detect(Variable, 'Non-white') ~ 5
                           , str_detect(Variable, 'Disease Duration') ~ 6
                           , str_detect(Variable, 'Follow-up') ~ 7
                           , str_detect(Variable, "BCD") ~ 7.5
                           , str_detect(Variable, "NTZ") ~ 7.75
                           , str_detect(Variable, "Utilization") ~ 8
                           , str_detect(Variable, 'Baseline') ~ 9
                           , str_detect(Variable, 'Year 1') ~ 10
                           , T ~ Inf) 
        ) %>% 
        arrange(order)
      
      if(.save) write.csv(.tbl.s2, here("Results", "BCD_NTZ", "SuppTable_Yr2.csv"))
    
  
  ## (4b) Year-3 ####
      .tbl.s3 = data %>% 
        mutate(DMT = ifelse(str_detect(DMT, "BCD"), "BCD", "NTZ")) %>% 
        group_by(PDDSOutcome_3yr) %>% 
        summarize(n = as.character(n())
                  , `BCD: n (%)` = npct(DMT, "BCD")
                  , `NTZ: n (%)` = npct(DMT, "NTZ")
                  , `Age at Target Treatment Initiation: years, mean (SD)` = meanSD(Age_at_DMT_Initiation)
                  , `Disease Duration: years, mean (SD)` = meanSD(Years_First_MS_PheCode_to_Study_DMT)
                  , `Follow-up Duration` = meanSD(Years_FirstPheCode_to_Study_DMT)
                  , `Self-Reported Gender (Male): n (%)` = npct(subject_sex_ch)
                  , `Self-Reported Gender (Female): n (%)` = npct(1-subject_sex_ch)
                  , `Race/Ethnicity - White, Non-Hispanic: n (%)` = npct(White_NonHispanic)
                  , `Race/Ethnicity - Non-white or Hispanic: n (%)` = npct(1-White_NonHispanic)
                  , `Utilization: number of codes, mean (SD) ` = meanSD(Utilization_Total)
                  # , `Baseline PDDS Risk (12 mo.)` = meanSD(pdds.baseline.12mo)
                  , `Baseline PDDS Risk (All): mean (SD)` = meanSD(pdds.baseline.All)
                  , `Year 1 PDDS Risk: mean (SD)` = meanSD(pdds.LF.12mo.3yr)

                  # n = as.character(n())
                  # , Age = meanSD(Age_at_DMT_Initiation)
                  # , `BCD Treatment Class` = npct(DMT, "BCD")
                  # , `Disease Duration` = meanSD(Years_First_MS_PheCode_to_Study_DMT)
                  # , `Follow-up Duration` = meanSD(Years_FirstPheCode_to_Study_DMT)
                  # , `Female Sex` = percent(mean(1-subject_sex_ch), accuracy=0.01)
                  # , `White, Non-Hispanic ` = meanSD(White_NonHispanic)
                  # , Utilization = meanSD(Utilization_Total)
                  # # , `Baseline PDDS Risk (12 mo.)` = meanSD(pdds.baseline.12mo)
                  # , `Baseline PDDS Risk (All)` = meanSD(pdds.baseline.All)
                  # , `Year 1 PDDS Risk` = meanSD(pdds.LF.12mo.3yr)
        ) %>% 
        mutate(PDDSOutcome_3yr =
                 case_when(PDDSOutcome_3yr == "Decrease/Improvement" ~ "Sustained Improvement"
                           , PDDSOutcome_3yr == "Increase/Progression" ~ "Sustained Worsening"
                           , str_detect(PDDSOutcome_3yr, "Insufficient") ~ "Unlabeled"
                           , T ~ PDDSOutcome_3yr
                 )
               , PDDSOutcome_3yr = paste0(PDDSOutcome_3yr, " (n=", n, ")")
        ) %>%
        pivot_longer(values_to = "Value", cols = -PDDSOutcome_3yr) %>%
        pivot_wider(values_from = "Value", names_from = "PDDSOutcome_3yr") %>% 
        # mutate() %>% 
        arrange(name!="n", name) %>% 
        select(Variable = name, starts_with("No"), starts_with("Sustain"), starts_with("Unlabel")) %>% 
        filter(Variable!="n") %>% 
        mutate(order = 
                 case_when(str_detect(Variable, 'Age at') ~ 1
                           , str_detect(Variable, 'Male') ~ 2
                           , str_detect(Variable, 'Female') ~ 3
                           , str_detect(Variable, 'White,') ~ 4
                           , str_detect(Variable, 'Non-white') ~ 5
                           , str_detect(Variable, 'Disease Duration') ~ 6
                           , str_detect(Variable, 'Follow-up') ~ 7
                           , str_detect(Variable, "BCD") ~ 7.5
                           , str_detect(Variable, "NTZ") ~ 7.75
                           , str_detect(Variable, "Utilization") ~ 8
                           , str_detect(Variable, 'Baseline') ~ 9
                           , str_detect(Variable, 'Year 1') ~ 10
                           , T ~ Inf) 
        ) %>% 
        arrange(order)
      
      
      if(.save) write.csv(.tbl.s3, here("Results", "BCD_NTZ", "SuppTable_Yr3.csv"))
      
      
    ## (4c; 5ab) Treatment ####
        # Results ofr 5ab also generated in SuppTbl_BalancedCovariates.R
      
      .tbl.trt = data %>% 
        mutate(DMT = case_when(str_detect(DMT, "BCD") ~ "BCD"
                               , T ~ "NTZ")
               ) %>% 
        group_by(DMT) %>% 
        summarize(n = as.character(n())
                  , `Age at Target Treatment Initiation: years, mean (SD)` = meanSD(Age_at_DMT_Initiation)
                  , `Disease Duration: years, mean (SD)` = meanSD(Years_First_MS_PheCode_to_Study_DMT)
                  , `Follow-up Duration` = meanSD(Years_FirstPheCode_to_Study_DMT)
                  , `Self-Reported Gender (Male): n (%)` = npct(subject_sex_ch)
                  , `Self-Reported Gender (Female): n (%)` = npct(1-subject_sex_ch)
                  , `Race/Ethnicity - White, Non-Hispanic: n (%)` = npct(White_NonHispanic)
                  , `Race/Ethnicity - Non-white or Hispanic: n (%)` = npct(1-White_NonHispanic)
                  , `Utilization: number of codes, mean (SD) ` = meanSD(Utilization_Total)
                  # , `Baseline PDDS Risk (12 mo.)` = meanSD(pdds.baseline.12mo)
                  , `Baseline PDDS Risk (All): mean (SD)` = meanSD(pdds.baseline.All)
                  , `Year 1 PDDS Risk: mean (SD)` = meanSD(pdds.LF.12mo.3yr)
                  , `Sustained Worsening: n (%)` = npct(PDDSOutcome_3yr, "Increase/Progression")
                  , `Sustained Improvement: n (%)` = npct(PDDSOutcome_3yr, "Decrease/Improvement")
                  , `No Sustained Change: n (%)` = npct(PDDSOutcome_3yr, "No Sustained Change")
                  , `Unlabeled: n (%)` = npct(PDDSOutcome_3yr, "Insufficient Data (Unlabeled)")
                  
                  
                  # n = as.character(n())
                  # , Age = meanSD(Age_at_DMT_Initiation)
                  # , `BCD Treatment Class` = npct(DMT, "BCD")
                  # , `Disease Duration` = meanSD(Years_First_MS_PheCode_to_Study_DMT)
                  # , `Follow-up Duration` = meanSD(Years_FirstPheCode_to_Study_DMT)
                  # , `Female Sex` = percent(mean(1-subject_sex_ch), accuracy=0.01)
                  # , `White, Non-Hispanic ` = meanSD(White_NonHispanic)
                  # , Utilization = meanSD(Utilization_Total)
                  # # , `Baseline PDDS Risk (12 mo.)` = meanSD(pdds.baseline.12mo)
                  # , `Baseline PDDS Risk (All)` = meanSD(pdds.baseline.All)
                  # , `Year 1 PDDS Risk` = meanSD(pdds.LF.12mo.2yr)
        ) %>% 
        mutate(DMT = paste0(DMT, " (n=", n, ")")) %>% 
        pivot_longer(values_to = "Value", cols = -DMT) %>%
        pivot_wider(values_from = "Value", names_from = "DMT") %>% 
        # mutate() %>% 
        arrange(name!="n", name) %>% 
        select(Variable = name, starts_with("BCD"), starts_with("NTZ")) %>% 
        filter(Variable!="n") %>% 
        mutate(order = 
                 case_when(str_detect(Variable, 'Age at') ~ 1
                           , str_detect(Variable, 'Male') ~ 2
                           , str_detect(Variable, 'Female') ~ 3
                           , str_detect(Variable, 'White,') ~ 4
                           , str_detect(Variable, 'Non-white') ~ 5
                           , str_detect(Variable, 'Disease Duration') ~ 6
                           , str_detect(Variable, 'Follow-up') ~ 7
                           , str_detect(Variable, "Utilization") ~ 8
                           , str_detect(Variable, 'Baseline') ~ 9
                           , str_detect(Variable, 'Year 1') ~ 10
                           , T ~ Inf) 
        ) %>% 
        arrange(order)
      
      
        if(.save) write.csv(.tbl.trt, here("Results", "BCD_NTZ", "SuppTable_Trt.csv"))
          
        ### Inference for Balancing Comparison ####
       
          data.test = data %>% 
            mutate(DMT = case_when(str_detect(DMT, "BCD") ~ "BCD", T ~ "NTZ"))
          
          wilcox.vars = c('Age_at_DMT_Initiation', 'Years_First_MS_PheCode_to_Study_DMT'
                          , 'Years_FirstPheCode_to_Study_DMT', 'Utilization_Total'
                          , 'pdds.baseline.All', "pdds.LF.All.3yr"
                          )
          
          chisq.vars = c('subject_sex_ch', 'White_NonHispanic')
          
          
          .wilcox.df = Reduce(`rbind`
                              , lapply(wilcox.vars, .wilcox.fn, data=data.test, group.var="DMT")
                              )
          .chisq.df = Reduce(`rbind`
                             , lapply(chisq.vars, .chisq.fn, data=data.test, group.var="DMT")
                             )
          
          .pval.df = as_tibble(union(.wilcox.df, .chisq.df)) %>% 
            mutate(var = case_when(
              var == 'Age_at_DMT_Initiation' ~ 'Age at Target Treatment Initiation: years, mean (SD)'
              , var == 'Years_First_MS_PheCode_to_Study_DMT' ~ 'Disease Duration: years, mean (SD)'
              , var == 'Years_FirstPheCode_to_Study_DMT' ~ 'Follow-up Duration'
              , var == 'Utilization_Total' ~ 'Utilization: number of codes, mean (SD) '
              , var == 'pdds.baseline.All' ~ 'Baseline PDDS Risk (All): mean (SD)'
              , var == 'pdds.LF.All.3yr' ~ 'Year 1 PDDS Risk: mean (SD)'
              , var == 'subject_sex_ch' ~ 'Self-Reported Gender (Male): n (%)'
              , var=='White_NonHispanic' ~ 'Race/Ethnicity - White, Non-Hispanic: n (%)'
            )) %>% 
            rename(name=var) %>% 
            mutate(p.val = formatC(p.val, format = 'f', flag='0', digits = 3)) %>% 
            mutate(p.val = case_when(p.val=="0.000" ~ "<0.001"
                                     , T ~ p.val)
                   )
          
          
          .tbl.trt.pval = .tbl.trt %>% 
            left_join(.pval.df, by=c('Variable'='name')) %>% 
            rename(name=Variable) %>% 
            mutate(order = 
                     case_when(str_detect(name, 'Age at') ~ 1
                               , str_detect(name, 'Male') ~ 2
                               , str_detect(name, 'Female') ~ 3
                               , str_detect(name, 'White,') ~ 4
                               , str_detect(name, 'Non-white') ~ 5
                               , str_detect(name, 'Disease Duration') ~ 6
                               , str_detect(name, 'Follow-up') ~ 7
                               , str_detect(name, "Utilization") ~ 8
                               , str_detect(name, 'Baseline') ~ 9
                               , str_detect(name, 'Year 1') ~ 10
                               , T ~ Inf) 
            ) %>% 
            arrange(order)
          
          if(.save) write.csv(.tbl.trt.pval, here("Results", "BCD_NTZ", "SuppTable_Trt_Pvalues.csv"))
        
      
# utilization and duration check ####  
  all.equal(bcd_ntz_2yr_joint$Utilization, bcd_ntz_3yr_joint$Utilization)
      
  all.equal(bcd_ntz_2yr_joint$Disease_Duration, bcd_ntz_3yr_joint$Disease_Duration)

  
  
  tmp2yr = data %>% 
    mutate(Utilization_TotalMean = mean(Utilization_Total)
           , DD_TotalMean = mean(Disease_Duration)) %>% 
    group_by(PDDSOutcome_2yr) %>% 
    summarize(Util_TotalMean = first(Utilization_TotalMean)
              , Util_Mean = mean(Utilization_Total)
              # , Utilization_SD = sd(Utilization)
              , DD_TotalMean = first(DD_TotalMean)
              , DD_Mean = mean(Disease_Duration)
              # , DD_SD = mean(Disease_Duration)
              ) %>% 
    mutate(Util_Diff = (Util_Mean - Util_TotalMean) / Util_TotalMean
           , DD_Diff = (DD_Mean - DD_TotalMean) / DD_TotalMean) %>% 
    select(PDDSOutcome_2yr, Util_Diff, DD_Diff)
  
  tmp2yr
  plot.df.2yr %>% select(PDDSOutcome_2yr, Utilization_Total, Disease_Duration)
  
  

  tmp3yr = data %>% 
    mutate(Utilization_TotalMean = mean(Utilization_Total)
           , DD_TotalMean = mean(Disease_Duration)) %>% 
    group_by(PDDSOutcome_3yr) %>% 
    summarize(Util_TotalMean = first(Utilization_TotalMean)
              , Util_Mean = mean(Utilization_Total)
              # , Utilization_SD = sd(Utilization)
              , DD_TotalMean = first(DD_TotalMean)
              , DD_Mean = mean(Disease_Duration)
              # , DD_SD = mean(Disease_Duration)
    ) %>% 
    mutate(Util_Diff = (Util_Mean - Util_TotalMean) / Util_TotalMean
           , DD_Diff = (DD_Mean - DD_TotalMean) / DD_TotalMean) %>% 
    select(PDDSOutcome_3yr, Util_Diff, DD_Diff)
  
  tmp3yr
  plot.df.3yr %>% select(PDDSOutcome_3yr, Utilization_Total, Disease_Duration)

    
  ## change percentage ####  
    data %>% 
      count(PDDSOutcome_2yr, PDDSOutcome_3yr) %>% 
      mutate(`%` = 100*round(n / sum(n), 4)) %>% 
      filter(PDDSOutcome_2yr!=PDDSOutcome_3yr) %>% 
      arrange(desc(n))
  
  
  
  
  
  
  # Numbers for biorender graph and derivation ####
    .df = bcd_ntz_3yr_joint
    bcd.deriv = read_csv(here("Data", "created_data", "BCD"
                              , gsub(".csv", "_Derivation.csv", .bcd_filename_3yr)
    ))
    ntz.deriv = read_csv(here("Data", "created_data", "NTZ"
                              , gsub(".csv", "_Derivation.csv", .ntz_filename_3yr)
    ))
    
    
  ## Final Iteration (Jan 2025) #####
    ### EHR - Total  ####
      bcd.deriv %>% filter(str_detect(Step, 'Total Patients'))
      # ntz.deriv %>% filter(str_detect(Step, 'Total Patients')) # equivalent
    
    ### EHR -  MS  ####
      bcd.deriv %>% filter(str_detect(Step, 'MS Patients'))
      # ntz.deriv %>% filter(str_detect(Step, 'MS Patients')) # equivalent
    
    ### EHR - BCD/NTZ (meeting inclusion criteria)  ####
      .df %>% 
        filter(is.na(id_participant)) %>% 
        count(DMT)
        # count(!is.na(PDDS_Change)) %>% mutate(`%` = 100*n/sum(n))
      
    ### Registry - BCD/NTZ (meeting inclusion criteria) ####
      .df %>% 
        filter(!is.na(id_participant)) %>% 
        count(DMT)
        # group_by(DMT) %>% count(!is.na(PDDS_Change)) %>% mutate(`%` = 100*n/sum(n))
    
    ### Labeled/Unlabeled   
      .df %>% 
        count(is.na(PDDS_Change))
  
  ### Registry ####
  # Labelled/Unlabelled Totals 
    .df %>% 
      filter(!is.na(id_participant)) %>% 
      count(!is.na(PDDS_Change)) %>% mutate(`%` = 100*n/sum(n))
    .df %>% 
      filter(!is.na(id_participant)) %>% 
      group_by(DMT) %>% count(!is.na(PDDS_Change)) %>% mutate(`%` = 100*n/sum(n))
    
    
      
  ## Prior Iteration
    # DMT boxes 
    
  
    sum(!is.na(bcd_ntz_3yr_joint$id_participant)) # registry 
    sum(is.na(bcd_ntz_3yr_joint$id_participant)) # EHR

    
    analytic_object_3yr$X %>% 
      filter(N_MS_PheCodes==0) %>% 
      count(N_MS_CUIs)
    
    
    bcd.deriv = read_csv(here("Data", "created_data", "BCD"
                              , gsub(".csv", "_Derivation.csv", .bcd_filename_3yr)
                              ))
    ntz.deriv = read_csv(here("Data", "created_data", "NTZ"
                              , gsub(".csv", "_Derivation.csv", .ntz_filename_3yr)
                              ))
    
    # starting no 
    (bcd.deriv %>% slice(4) %>% pull(N)) + (ntz.deriv %>% slice(4) %>% pull(N))
    
    # ever chemotherapy 
    bcd.deriv %>% slice(4:10)
    ntz.deriv %>% slice(4:10)
    
    bcd.deriv %>% filter(str_detect(Step, "chemo") | str_detect(lead(Step), "chemo"))
    ntz.deriv %>% filter(str_detect(Step, "chemo") | str_detect(lead(Step), "chemo"))
    
    
    # previous high efficacy 
    (bcd.deriv %>% slice(4) %>% pull(N)) - (bcd.deriv %>% slice(6) %>% pull(N)) + 
    (ntz.deriv %>% slice(4) %>% pull(N)) - (ntz.deriv %>% slice(6) %>% pull(N))
    
    # drug switching 
    (bcd.deriv %>% slice(6) %>% pull(N)) - (bcd.deriv %>% slice(7) %>% pull(N)) + 
    (ntz.deriv %>% slice(6) %>% pull(N)) - (ntz.deriv %>% slice(7) %>% pull(N))
    
    
    # chemotherapy 
    (bcd.deriv %>% slice(7) %>% pull(N)) - (bcd.deriv %>% slice(8) %>% pull(N)) + 
    (ntz.deriv %>% slice(7) %>% pull(N)) - (ntz.deriv %>% slice(8) %>% pull(N))
    
    # missing covariates 
    (bcd.deriv %>% slice(8) %>% pull(N)) - (bcd.deriv %>% slice(10) %>% pull(N)) + 
    (ntz.deriv %>% slice(8) %>% pull(N)) - (ntz.deriv %>% slice(10) %>% pull(N))
    
    
    
    # label proportion 
    data %>% group_by(DMT) %>% count(PDDSOutcome_2yr) %>% mutate(`%` = 100*n/sum(n))
    data %>% group_by(DMT) %>% count(PDDSOutcome_3yr) %>% mutate(`%` = 100*n/sum(n))