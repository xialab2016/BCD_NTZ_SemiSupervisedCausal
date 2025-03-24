
# Set-Up #### 
library(pacman)
p_load(here, ggplot2, Matrix, dplyr, RColorBrewer, cowplot, scales
       , forcats
       , rlang
       , weights
       , sjstats # "weights" option in  mann_whitney_test fn
       )
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

coef.df %>% 
  group_by(Variable) %>%
  filter(max(abs(Coefficient))>0) %>% 
  arrange(Variable) %>% 
  ggplot(aes(x = Model, y = Variable, color=Coefficient.Std)) + labs(color="Coefficient\n(max-normalized)") +
  # ggplot(aes(x = Model, y = Variable, color=Coefficient)) + theme() + labs(color="Coefficient") +
  geom_point() + 
  scale_color_distiller(type="div", palette = "RdBu") + 
  theme_classic()


## Data Import ####

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

ehr_ft_list = rbind(once_disability_ehr, once_ms_ehr)

# Model Objects (loaded for weights) 

si.3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Increase_3yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
si.2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Increase_2yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))

sd.3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Decrease_3yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
sd.2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Decrease_2yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))

# all.equal(si.2yr$DiPS, sd.2yr$DiPS)
  # both 2 year outcome models only have intercept/treatment status as non-zero coefficients
  # thus only the underlying treatment model contributes to DiPS, and this is the same treatment model in both instances

data = data.frame(PATIENT_NUM = bcd_ntz_2yr_joint$PATIENT_NUM
                  , Registry = ifelse(!is.na(bcd_ntz_2yr_joint$id_participant), 1, 0)
                  , DMT = ifelse(analytic_object_2yr$A, "BCD", "NTZ")
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
                  , DiPS.2yr = sd.2yr$DiPS 
                  , DiPS.3yr.SI = si.3yr$DiPS
                  , DiPS.3yr.SD = sd.3yr$DiPS 
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


# DiPS are propensity scores of receiving BCD over NTZ, invert weights for NTZ observations 
data.w = data %>% mutate(across(starts_with("DiPS")
                              , ~ case_when(DMT=="BCD" ~ 1/.x
                                            , T ~ 1/(1-.x) )
                              )
                       )



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



# using calculation which also includes (N-1)/N correction for number of non-zero weights 
  # equivalent to npreg::wtd.var, compare to commentd functions below/Hmisc::wtd.var 
  # weighted.var = function(x, w) weighted.mean( ( x-weighted.mean(x, w) )^2, w)*sum(w)/(sum(w)-1)
  # weighted.sd = function(x, w) sqrt(weighted.var(x, w))

weighted.var = function(x, w){
  N = sum(w!=0)
  (N/(N-1) )*sum( w*(x-weighted.mean(x, w) )^2 )/sum(w)
}
weighted.sd = function(x, w) sqrt(weighted.var(x, w))  

w.meanSD = function(x, w){
  .wm = round(weighted.mean(x, w), 3)
  .wsd = round(weighted.sd(x, w), 3)
  
  m = formatC(.wm, format = 'f', flag='0', digits = 2)
  sd = formatC(.wsd, format = 'f', flag='0', digits = 2)
  
  paste0(m, " (", sd, ")")
}

w.npct = function(x, w, cat=1){
  .p = weighted.mean(x==cat, w)
  n = round(sum( (x==cat)*w ), 0)
  
  paste0(format(n, big.mark=",")
         , " ("
         , scales::percent(.p, accuracy = 0.01)
         , ")")
}


balance.tbl = function(data, w.var){
  
  if(is.null(w.var) | !is.character(w.var)){
    message("Null or invalid weight var, unweighted statistics calculated")
    data = data %>% mutate(weight=1) 
    w.var="weight"
  }
  
  .wvar = rlang::ensym(w.var)
  
  .bal.tab = data %>% 
    # mutate(across(starts_with("DiPS"), ~ pmin(pmax(., 0.1), 0.9))) %>%
    group_by(DMT) %>% 
    summarise(n = paste0("(n=", n(), ")")
              , `Weighted n` = paste0("(weighted n=", round(sum(!!.wvar), 3), ")")
              , `Age at Target Treatment Initiation: years, mean (SD)` = w.meanSD(Age_at_DMT_Initiation, w=!!.wvar)
              
              , `Disease Duration: years, mean (SD)` = w.meanSD(Years_First_MS_PheCode_to_Study_DMT, w=!!.wvar)
              , `Follow-up Duration` = w.meanSD(Years_FirstPheCode_to_Study_DMT, w=!!.wvar)
              , `Self-Reported Gender (Male): n (%)` = w.npct(subject_sex_ch, !!.wvar)
              , `Self-Reported Gender (Female): n (%)` = w.npct(1-subject_sex_ch, !!.wvar)
              , `Race/Ethnicity - White, Non-Hispanic: n (%)` = w.npct(White_NonHispanic, !!.wvar)
              , `Race/Ethnicity - Non-white or Hispanic: n (%)` = w.npct(1-White_NonHispanic, !!.wvar)
              , `Utilization: number of codes, mean (SD) ` = w.meanSD(Utilization_Total, w = !!.wvar)
              , `Baseline PDDS Risk (All): mean (SD)` = w.meanSD(pdds.baseline.All, w = !!.wvar)
              , `Year 1 PDDS Risk: mean (SD)` = w.meanSD(pdds.LF.12mo.3yr, w = !!.wvar)
              # , `Sustained Worsening: n (%)` = w.npct(PDDSOutcome_3yr, cat="Increase/Progression", w=!!.wvar)
              # , `Sustained Improvement: n (%)` = w.npct(PDDSOutcome_3yr, cat="Decrease/Improvement", w=!!.wvar)
              # , `No Sustained Change: n (%)` = w.npct(PDDSOutcome_3yr, cat="No Sustained Change", w=!!.wvar)
              # , `Unlabeled: n (%)` = w.npct(PDDSOutcome_3yr, cat="Insufficient Data (Unlabeled)", w=!!.wvar)
              
              , Weight = w.var
              ) 
  
  return(.bal.tab)
}


weighted.tests = function(data, w, wilcox.vars, chisq.vars){
  
  # mann_whitney_test via sjstats package (with survey package's svyranktest)
  # wtd.chi.sq via weights package 
  
  wilcox.df = Reduce(`rbind`
                     , lapply(wilcox.vars
                              , function(x) data.frame(var=x, p.val = mann_whitney_test(data=data, select=x, by='DMT', weight=w)$p)
                     )
  )
  
  chisq.df = Reduce(`rbind`
                    , lapply(chisq.vars
                             , function(x) {
                               data.frame(var=x, p.val = wtd.chi.sq(data[[x]], data[['DMT']], weight=data[[w]])['p.value'])
                             }
                             )
                    )
  
  .op = as_tibble(union(wilcox.df, chisq.df)) %>% 
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
    mutate(p.val = formatC(p.val, format = 'f', flag='0', digits = 3))
  
  if(anyNA(.op)) stop("Naming issue in output dataframe, merging not possible")
  
  return(.op)  
}

wilcox.vars = c('Age_at_DMT_Initiation', 'Years_First_MS_PheCode_to_Study_DMT'
                , 'Years_FirstPheCode_to_Study_DMT', 'Utilization_Total'
                , 'pdds.baseline.All', "pdds.LF.All.3yr"
                )
chisq.vars = c('subject_sex_ch', 'White_NonHispanic')

weighted.tests(data.w, 'DiPS.2yr', wilcox.vars, chisq.vars)
weighted.tests(data.w, 'DiPS.3yr.SD', wilcox.vars, chisq.vars)
weighted.tests(data.w, 'DiPS.3yr.SI', wilcox.vars, chisq.vars)


## 2 year weights ####
  bal.tab.2yr = balance.tbl(data=data.w, w.var="DiPS.2yr") %>% 
    pivot_longer(values_to = "Value"
                 , cols=n:Weight
                 ) %>% 
    pivot_wider(id_cols=name, values_from=Value, names_from=DMT) %>% 
    left_join(weighted.tests(data.w, 'DiPS.2yr', wilcox.vars, chisq.vars)
              , 'name') %>% 
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
  
  
  bal.tab.2yr


## 3 year - SD weights ####
  bal.tab.3yr.sd = balance.tbl(data=data.w, w.var="DiPS.3yr.SD") %>% 
    pivot_longer(values_to = "Value"
                 , cols=n:Weight
    ) %>% 
    pivot_wider(id_cols=name, values_from=Value, names_from=DMT) %>% 
    left_join(weighted.tests(data.w, 'DiPS.3yr.SD', wilcox.vars, chisq.vars)
              , 'name') %>% 
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
  
  
  bal.tab.3yr.sd

## 3 year - SI weights ####
  bal.tab.3yr.si = balance.tbl(data=data.w, w.var="DiPS.3yr.SI") %>% 
    pivot_longer(values_to = "Value"
                 , cols=n:Weight
    ) %>% 
    pivot_wider(id_cols=name, values_from=Value, names_from=DMT) %>% 
    left_join(weighted.tests(data.w, 'DiPS.3yr.SI', wilcox.vars, chisq.vars)
              , 'name') %>% 
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
  
  bal.tab.3yr.si

## Export ####
  .save = F
  if(.save) write.csv(bal.tab.2yr, here("Results", "BCD_NTZ", "SuppTbl_CovarBalancing_DiPS2yr.csv"))
  if(.save) write.csv(bal.tab.3yr.si, here("Results", "BCD_NTZ", "SuppTbl_CovarBalancing_DiPS3yr_SI.csv"))
  if(.save) write.csv(bal.tab.3yr.sd, here("Results", "BCD_NTZ", "SuppTbl_CovarBalancing_DiPS3yr_SD.csv"))



  
# Visualization/Loveplot style? ####


# bal.tab.no.weight = balance.tbl(data=data.w, w.var=NULL) %>% 
#   pivot_longer(values_to = "Value"
#                , cols=n:Weight
#   ) %>% 
#   pivot_wider(id_cols=name, values_from=Value, names_from=DMT)
# 
# 
# 
# .tmp %>% 
#   filter(name %nin% c("n", "Weighted n", "Weight")) %>% 
#   rowwise() %>% 
#   mutate(BCD.mean = as.numeric(gsub(',', '', str_split(BCD, " ")[[1]][1]))
#          , NTZ.mean = as.numeric(gsub(',', '', str_split(NTZ, " ")[[1]][1]))
#          
#          , BCD.sd = as.numeric(gsub('%|\\(|\\)', '', str_split(BCD, " ")[[1]][2]))
#          , NTZ.sd = as.numeric(gsub('%|\\(|\\)', '', str_split(NTZ, " ")[[1]][2]))
#          
#          , .mad = abs(BCD.mean - NTZ.mean)
#          , .mad.std = abs(BCD.mean - NTZ.mean) / (BCD.sd + NTZ.sd)
#          , .mad.std.cont = case_when(str_detect('%', name) ~ .mad
#                                      , T ~ .mad.std)
#   ) 
# 
