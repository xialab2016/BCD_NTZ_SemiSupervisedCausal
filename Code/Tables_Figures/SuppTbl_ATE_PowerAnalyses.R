# ~post-hoc power analysis 

library(pacman)
p_load(here, readr, tibble, dplyr)

.bcd_filename_3yr <- "bcd_allMS_3yr_2025_02_05.csv" #"bcd_all_3yr_2023_12_20.csv"
.bcd_filename_2yr <- "bcd_allMS_2yr_2025_02_05.csv" #"bcd_all_2yr_2023_12_20.csv"

.ntz_filename_3yr <- "ntz_allMS_3yr_2025_02_05.csv" #"ntz_all_3yr_2023_12_20.csv"  
.ntz_filename_2yr <- "ntz_allMS_2yr_2025_02_05.csv" #"ntz_all_2yr_2023_12_20.csv"  

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


  ## Creating Shared Cohorts and Import Objects ####
  bcd_ntz_3yr_joint <- bcd_ntz_3yr %>% 
    inner_join(bcd_ntz_2yr[,"PATIENT_NUM"], "PATIENT_NUM")
  
  ### Semi-supervised ####
  perturb_increase_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Increase_3yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
  perturb_increase_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Increase_2yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
  perturb_decrease_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Decrease_3yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
  perturb_decrease_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Decrease_2yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
  
  
  ### Complete Case DiPS ####
  perturb_CC_increase_3yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "CompleteCase_2025_Increase_3yr_ATE_BCD_NTZ_2.RDS"))
  perturb_CC_increase_2yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "CompleteCase_2025_Increase_2yr_ATE_BCD_NTZ_2.RDS"))
  perturb_CC_decrease_3yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "CompleteCase_2025_Decrease_3yr_ATE_BCD_NTZ_2.RDS"))
  perturb_CC_decrease_2yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "CompleteCase_2025_Decrease_2yr_ATE_BCD_NTZ_2.RDS"))

  # perturb_CC_increase_3yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "CompleteCase_3mo_2025_Increase_3yr_ATE_BCD_NTZ_2.RDS"))
  # perturb_CC_increase_2yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "CompleteCase_3mo_2025_Increase_2yr_ATE_BCD_NTZ_2.RDS"))
  # perturb_CC_decrease_3yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "CompleteCase_3mo_2025_Decrease_3yr_ATE_BCD_NTZ_2.RDS"))
  # perturb_CC_decrease_2yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "CompleteCase_3mo_2025_Decrease_2yr_ATE_BCD_NTZ_2.RDS"))


  ### 3-month Window Semi-supervised #### 
    perturb_3mo_increase_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity", "Increase_2yr_3monthWindow_SensitivityAnalysis_2025_ATE_BCD_NTZ.RDS"))
    perturb_3mo_increase_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity", "Increase_3yr_3monthWindow_SensitivityAnalysis_2025_ATE_BCD_NTZ.RDS"))
    perturb_3mo_decrease_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity", "Decrease_2yr_3monthWindow_SensitivityAnalysis_2025_ATE_BCD_NTZ.RDS"))
    perturb_3mo_decrease_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity", "Decrease_3yr_3monthWindow_SensitivityAnalysis_2025_ATE_BCD_NTZ.RDS"))
    
  
# Utils/Functions #### 
  .mad = function(x, na.rm=T) mean(abs(x-mean(x, na.rm=na.rm)), na.rm=na.rm)
  
  normal.twoside.power <- function(beta.MD
                                   , alpha
                                   , .std.err
                                   , model
                                   , root){
    
    if(missing(.std.err)) {
      if(missing(model)) stop("Must include model object or specify .std.err")
      
      # higher threshold for complete-case - lower total perturbations
      # if( mean(is.na(model$ATE_Perturbed)) > 0.05) stop(">5% missing perturbations, review prior to analysis")
      .p = na.omit(model$ATE_Perturbed)
      # .std.err = sd(.p)
      .std.err = .mad(.p)
    }
    
    power = pnorm(-beta.MD/.std.err - qnorm(1-alpha/2)) + pnorm(beta.MD/.std.err - qnorm(1-alpha/2))
    
    if(!missing(root)) return(power - root)
    
    
    op.list = list(power = power
                   , ate.md = beta.MD)
    class(op.list) = "power.list"
    return(op.list)
  
  }
  
  print.power.list = function(obj) print(paste0("Power = ", round(obj$power, 4), "; ATE (Min Detect) = ", round(obj$ate.md, 4)))
  
  pwr.root.find = function(model, alpha = 0.05, root = 0.8, interval = c(0.001, 0.4)){
    root.find = uniroot(normal.twoside.power, alpha = alpha, model = model, root = root
                        , interval = interval
    )
    
    normal.twoside.power(root.find$root, alpha, model=model)
  }

# MDE Root Finding #### 
  # pwr.root.find(perturb_CC_average_2yr, root = 0.8, interval = c(0, 1)) # example single calc 
  
  .perturb.ss = setdiff(grep("perturb_", ls(), value=T)
                        , c(grep("perturb_CC", ls(), value=T), grep("perturb_3mo", ls(), value=T))
                        )
  .perturb.ss.pwr = unlist(
    lapply(.perturb.ss, function(x){
      .x = get(x)
      pwr.root.find(.x, root=0.8, interval=c(0, 10))$ate.md
      })
    )
  names(.perturb.ss.pwr) = .perturb.ss
  
  .perturb.cc = grep("perturb_CC", ls(), value=T)
  .perturb.cc.pwr = unlist(
    lapply(.perturb.cc, function(x){
      .x = get(x)
      pwr.root.find(.x, root=0.8, interval=c(0.1, 10))$ate.md
    })
  )
  names(.perturb.cc.pwr) = .perturb.cc
  

  .perturb.3mo = grep("perturb_3mo", ls(), value=T)
  .perturb.3mo.pwr = unlist(
    lapply(.perturb.3mo, function(x){
      .x = get(x)
      pwr.root.find(.x, root=0.8, interval=c(0, 10))$ate.md
    })
  )
  names(.perturb.3mo.pwr) = .perturb.3mo


  # cbind(.perturb.ss.pwr, .perturb.cc.pwr, .perturb.3mo.pwr) %>% 
  cbind(.perturb.ss.pwr, .perturb.cc.pwr) %>% round(4) # MDEs 


# Supp Table 
  .df.fn = function(x){
    .x = get(x)
    .df = data.frame(ATE = .x$ATE_full
                     , MeanAD = .mad(.x$ATE_Perturbed, na.rm=T)
                     , SD = sd(.x$ATE_Perturbed, na.rm=T)
                     , CI.L = quantile(.x$ATE_Perturbed, 0.025, na.rm=T)
                     , CI.H = quantile(.x$ATE_Perturbed, 0.975, na.rm=T)
                     , Model = x) %>% 
      mutate(ATE.L.Diff = ATE-CI.L 
             , ATE.H.Diff = CI.H - ATE) %>% 
      mutate(across(ATE:CI.H, round, 4))
    rownames(.df) = NULL
    return(.df)
  }
  
  .cc.df = lapply(.perturb.cc, .df.fn) %>% bind_rows()
  .ss.df = lapply(.perturb.ss, .df.fn) %>% bind_rows() 
  .ss.3mo.df = lapply(.perturb.3mo, .df.fn) %>% bind_rows()
  
  .pwr.df = data.frame(Model = "Semi-Supervised (6-month PDDS separation)"
                       , MDE = round(.perturb.ss.pwr, 3)
                       ) %>% 
    union_all(data.frame(Model = "Semi-Supervised (3-month PDDS separation)"
                     , MDE = round(.perturb.3mo.pwr, 3)
                     )
          ) %>% 
    union_all(data.frame(Model = "Labelled Data Only"
                     , MDE = round(.perturb.cc.pwr, 3)
                     )
          ) %>% 
    rownames_to_column("tmp") %>% 
    mutate(Outcome = case_when(str_detect(tmp, "average") ~ "Average PDDS"
                               , str_detect(tmp, "increase") ~ "Sustained Worsening"
                               , str_detect(tmp, "decrease") ~ "Sustained Improvement"
                               )
           , Endpoint = ifelse(str_detect(tmp, "3yr"), "3-year", "2-year")
           ) %>% 
    select(-tmp)
    
  
  
  .output.tbl = .ss.df %>% 
    mutate_if(is.numeric, round, 3) %>% 
    mutate(CI = paste0("(", CI.L, ", ", CI.H, ")")) %>% 
    select(Model, ATE, CI, `Std. Err` = SD) %>% 
    union(  .cc.df %>% 
              mutate_if(is.numeric, round, 3) %>% 
              mutate(CI = paste0("(", CI.L, ", ", CI.H, ")")) %>% 
              select(Model, ATE, CI, `Std. Err` = SD)
    ) %>% 
    union(  .ss.3mo.df %>% 
              mutate_if(is.numeric, round, 3) %>% 
              mutate(CI = paste0("(", CI.L, ", ", CI.H, ")")) %>% 
              select(Model, ATE, CI, `Std. Err` = SD)
    ) %>% 
    # filter(str_detect(Model, "3yr")) %>% # 3-year results only 
    mutate(Outcome = case_when(str_detect(Model, "average") ~ "Average PDDS"
                               , str_detect(Model, "increase") ~ "Sustained Worsening"
                               , str_detect(Model, "decrease") ~ "Sustained Improvement"
                               )
           , Endpoint = ifelse(str_detect(Model, "3yr"), "3-year", "2-year")
           , Model = case_when(str_detect(Model, "CC") ~ "Labelled Data Only"
                               , str_detect(Model, "3mo") ~ "Semi-Supervised (3-month PDDS separation)"
                               , T ~ "Semi-Supervised (6-month PDDS separation)")
           ) %>% 
    select(Outcome, Model, everything()) %>% 
    arrange(Outcome, Model) %>% 
    left_join(.pwr.df, by = c("Model", "Outcome", "Endpoint"))
  
# Supp Table Export #### 
  .save = T
  
  ## 3-year ####
    .output.tbl %>% filter(Endpoint=="3-year"
                           & !str_detect(Outcome, "Average")
                           ) 
  
    if(.save) write.csv(.output.tbl %>% filter(Endpoint=="3-year" & 
                                                 !str_detect(Outcome, "Average")
                                               ) 
                        , here("Results", "BCD_NTZ", "SuppTbl_Power3yr.csv")
                        )
      # write.csv(here("Results", "BCD_NTZ", "Archive", "compare_3mo_6mo.csv"))
      

  ## 2-year ####
    .output.tbl %>% filter(Endpoint=="2-year"
                           & !str_detect(Outcome, "Average")
                           )   
  
  if(.save) write.csv(.output.tbl %>% filter(Endpoint=="2-year" & 
                                               !str_detect(Outcome, "Average")
                                             ) 
                      , here("Results", "BCD_NTZ", "SuppTbl_Power2yr.csv")
                      )


# ## Comparing label/unlabel proportion #### 
#     # ~post-hoc power analysis 
#     
#     library(pacman)
#     p_load(here, readr, tibble, dplyr)
#     
#     .bcd_filename_3yr <- "bcd_allMS_3yr_2024_03_07.csv" #"bcd_all_3yr_2023_12_20.csv"
#     .bcd_filename_2yr <- "bcd_allMS_2yr_2024_03_07.csv" #"bcd_all_2yr_2023_12_20.csv"
#     .bcd_filename_3yr.3mo <- "bcd_allMS_3monthWindow_3yr_2024_06_05.csv" #"bcd_all_3yr_2023_12_20.csv"
#     .bcd_filename_2yr.3mo <- "bcd_allMS_3monthWindow_2yr_2024_06_05.csv" #"bcd_all_2yr_2023_12_20.csv"
#     
#     .ntz_filename_3yr <- "ntz_allMS_3yr_2024_03_07.csv" #"ntz_all_3yr_2023_12_20.csv"  
#     .ntz_filename_2yr <- "ntz_allMS_2yr_2024_03_07.csv" #"ntz_all_2yr_2023_12_20.csv"  
#     .ntz_filename_3yr.3mo <- "ntz_allMS_3monthWindow_3yr_2024_06_05.csv" #"ntz_all_3yr_2023_12_20.csv"  
#     .ntz_filename_2yr.3mo <- "ntz_allMS_3monthWindow_2yr_2024_06_05.csv" #"ntz_all_2yr_2023_12_20.csv"  
#     
#     # Data Creation/Import #### 
#     .skip_build <- T  # skips cohort_building.R in DMT script 
#     source(here("Code", "driver_scripts", "bcd.R"))
#     bcd_3yr <- read_csv(here("Data", "created_data", "BCD", .bcd_filename_3yr))
#     bcd_2yr <- read_csv(here("Data", "created_data", "BCD", .bcd_filename_2yr))
#     bcd_3yr.3mo <- read_csv(here("Data", "created_data", "BCD", .bcd_filename_3yr.3mo))
#     bcd_2yr.3mo <- read_csv(here("Data", "created_data", "BCD", .bcd_filename_2yr.3mo))
#     
#     .skip_build <- T  # skips cohort_building.R in DMT script 
#     source(here("Code", "driver_scripts", "ntz.R"))
#     ntz_3yr <- read_csv(here("Data", "created_data", "NTZ", .ntz_filename_3yr))
#     ntz_2yr <- read_csv(here("Data", "created_data", "NTZ", .ntz_filename_2yr))
#     ntz_3yr.3mo <- read_csv(here("Data", "created_data", "NTZ", .ntz_filename_3yr.3mo))
#     ntz_2yr.3mo <- read_csv(here("Data", "created_data", "NTZ", .ntz_filename_2yr.3mo))
#     
#     bcd_ntz_3yr <- rbind(bcd_3yr, ntz_3yr)
#     bcd_ntz_2yr <- rbind(bcd_2yr, ntz_2yr)
#     
#     bcd_ntz_3yr.3mo <- rbind(bcd_3yr.3mo, ntz_3yr.3mo)
#     bcd_ntz_2yr.3mo <- rbind(bcd_2yr.3mo, ntz_2yr.3mo)
#     
#     table(is.na(bcd_ntz_3yr$PDDS_Change))
#     table(is.na(bcd_ntz_3yr.3mo$PDDS_Change))
#     
#     table(is.na(bcd_ntz_2yr$PDDS_Change))
#     table(is.na(bcd_ntz_2yr.3mo$PDDS_Change))
#     
#     table(bcd_ntz_3yr$PDDS_Change)
#     table(bcd_ntz_3yr$PDDS_Change
#           , bcd_ntz_3yr.3mo$PDDS_Change
#           , dnn = c("6-month", "3-month")
#           , useNA='always')
#     
#     
