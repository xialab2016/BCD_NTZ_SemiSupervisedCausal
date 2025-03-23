
source(here("Code", "setup.R"))

if("ehr_2004_2022" %in% ls()){
  ehr <- ehr_2004_2022
  ehr$patient_num <- as.character(ehr$patient_num)
}else if ("ehr" %nin% ls()){
  ehr <- read_csv(here("Data", "Latest"
                       , ehr_codified_rolled)
                  , show_col_types = FALSE)
  # ehr$patient_num <- as.character(ehr$patient_num)
}

# cleaning phecodes using prescribed rules 
  # creating function, used later within analytic_wrangling
  #### Just keep parent code if one digit does not appear; remove two digits ####
  PheCodeCleaning = function(features){
    feat.set = unique(features)
    phecode.all = feat.set[grep("PheCode", feat.set)]
    
    tmp = matrix(unlist(lapply(strsplit(phecode.all, split = "\\.")
                               , function(x){cbind(len=length(x), numdigit=nchar(x[2]))})
                        )
    , ncol=2, byrow = TRUE)
    
    colnames(tmp)= c("len", "numdigit")
    
    phecode.nonedigit = phecode.all[tmp[,"len"] == 1]
    phecode.onedigit = phecode.all[which(tmp[,"numdigit"] == 1)]
    phecode.onedigit.convert.none = unlist(lapply(strsplit(phecode.onedigit, split="\\."), function(x){x[1]}))
    # phecode.twodigit = phecode.all[which(tmp[,"numdigit"] == 2)] # not retained, do not need to track 
    
    feat.set.red = c(phecode.onedigit # all one digit codes 
                 , setdiff(phecode.nonedigit, phecode.onedigit.convert.none) # all parent codes without a one digit child code 
                 , setdiff(feat.set, phecode.all) # all non PheCode features 
    )
    
  }

  # ehr.analytic = ehr[ehr$feature_id %in% feat.set,]
  # .excluded_fts = setdiff(ehr$feature_id, ehr.analytic$feature_id)

# minimal cost to load these lazily 

ms_dmt = read_excel(here("Data", "Latest", "MS Drugs Mechanism Citations Efficacy 20230107.xlsx")) 
ms_rxnorms = paste0("RXNORM:", ms_dmt$`RxNorm Ingredient id`)


once_ms_ehr = read_csv(here("Data", "Latest", "ONCE_Multiple sclerosis_Codified.csv")) %>% 
  filter(expanded_features==1) %>%
  select(term = Variable
         , expanded_features
         , phenotyping_features) %>% 
  mutate(set = "MS") %>% 
  filter(term != 'PheCode:335' 
         & !str_detect(tolower(term), pattern = "local|Local|lab|Lab")
         # & term %nin% ms_rxnorms
         ) 

once_ms_nlp = read_csv(here("Data", "Latest", "ONCE_Multiple sclerosis_CUI.csv")) %>% 
  # filter(expanded_features==1) %>%
  select(cui
         , expanded_features
         , phenotyping_features) %>% 
  mutate(set = "MS") 


once_disability_ehr = read_csv(here("Data", "Latest", "ONCE_disability_PheCode296.2_cos0.165.csv")) %>%  
  filter(phenotyping_features == 1) %>%
  select(term = Variable
         , expanded_features
         , phenotyping_features) %>% 
  mutate(set = "disability") %>% 
  filter(term != 'PheCode:335' 
         & !str_detect(tolower(term), pattern = "local|Local|lab|Lab")
         # & term %nin% ms_rxnorms
  ) 

once_disability_nlp = read_csv(here("Data", "Latest", "ONCE_disability_C0231170_titlecos0.5_titlecut0.3_exactFALSE.csv")) %>%  
  # filter(phenotyping_features == 1) %>%
  select(cui
         , expanded_features
         , phenotyping_features) %>% 
  mutate(set = "disability")

CUI = rbind(once_disability_nlp, once_ms_nlp)
ehr_ft_list = rbind(once_disability_ehr, once_ms_ehr)


analytic_wrangling <- function(.analytic_df
                               , pre_process = T
                               , .lookback_window = dmonths(6)
                               , .lookforward_window = dmonths(0)
                               , .EHR_trim = 0.9
                               , .NLP = T
                               , .EHR_all = F
                               , .NLP_all = F
                               , .import_NLP = F # only T if need to re-import NLP data (slow process, large file)
){
  
  # Set-Up #### 
  .analytic_df$subject_sex_ch = as.factor(.analytic_df$subject_sex_ch)
  .analytic_df$patient_num = as.numeric(.analytic_df$PATIENT_NUM)
  
  if(!.NLP_all){
    CUI = CUI %>% filter( (expanded_features==1 & set=="MS") | (phenotyping_features==1 & set=="disability") )
  }
  
  if("NLP" %nin% ls(envir = .GlobalEnv) | .import_NLP){
    print("Importing NLP data to Global environment (bottleneck step, ideally set .import_NLP to F and previously import NLP file to avoid redundant, intensive import step")
    NLP <- read_csv(here("Data", "Latest", paste0(nlp_filename, ".csv")))
    NLP = setDT(NLP) # fread function not completing import, read_csv preferred for now 
    assign("NLP", NLP, envir=.GlobalEnv)
  }
  
  
  
  # NLP_pts = NLP[as.data.table(.analytic_df), on = .(patient_num), nomatch = NULL, allow.cartesian = T] # inner_join
  NLP_cui = NLP[as.data.table(CUI), on = .(feature_id==cui), nomatch = NULL] # inner_join
  NLP_pts = NLP_cui[as.data.table(.analytic_df), on = .(patient_num), nomatch = NULL, allow.cartesian = T] # inner_join
  NLP_sub = NLP_pts[start_date < date + .lookforward_window
                    & start_date >= date - .lookback_window
                    , .(patient_num, start_date, feature_id) ,]
    
  rm(NLP_pts)
  
  # EHR features ####
    if (!(exists("ehr", where=.GlobalEnv))){
      stop("EHR data not present in R environment. Import EHR features using `ehr_codified_rolled`")
    }
  
  ## Codified/Correlated Subset ####
    if(!.EHR_all){
      ehr.set <- ehr %>% filter(feature_id %in% ehr_ft_list$term)
  
      ehr_valid <- ehr.set %>% 
        inner_join(.analytic_df[,c("patient_num", "date")], by="patient_num"
                   , relationship="many-to-many"
                   ) %>% 
        filter(start_date < date + .lookforward_window
               & start_date >= date - .lookback_window
               )
  
      ehr_valid$DMT_Study_Start <- NULL
      
      ehr_count_long <- ehr_valid %>%
        filter(substr(feature_id, 1, 4) %in% c("PheC", "RXNO", "LOIN", "CCS-")) %>%
        group_by(patient_num, feature_id) %>% 
        summarize(count=n(),.groups = "drop")

      ehr_count = tidyr::spread(ehr_count_long, key = feature_id, value = count, fill = 0) 
      
      .cols_non0_bool <- apply(ehr_count %>% select(-patient_num), MARGIN = 2, function(x) sum(x==0) / length(x)) < .EHR_trim
      .cols_non0 = names(.cols_non0_bool)[.cols_non0_bool]
      .cols_non0_roll = PheCodeCleaning(.cols_non0)
      # ehr_count <- ehr_count[, c(T, .cols_non0)] # T keeps patients_num
      ehr_count <- ehr_count[, c("patient_num", .cols_non0_roll)] # T keeps patients_num
      
      .incl_fts <- ehr_count %>% select(-patient_num) %>% colnames
      
      .ehr_ft_summary = list(Fts_Considered = unique(ehr_ft_list$term)
                             , Fts_Obs = unique(ehr_count_long$feature_id)
                             , Fts_Obs_Trim = unique(.cols_non0)
                             , Ft_Obs_Trim_Roll = unique(.cols_non0_roll) # colnames(ehr_count)[2:ncol(ehr_count)]
                             )
      
    }else## All Observed %>% Fts #### 
    {
      ehr.pts <- ehr %>% 
        inner_join(.analytic_df[,c("patient_num", "date")], by="patient_num"
                   , relationship="many-to-many"
                   ) %>% 
        filter(start_date < date + .lookforward_window
               & start_date >= date - .lookback_window
               )    
      
      ehr_count_all_long <- ehr.pts %>%
        filter(substr(feature_id, 1, 4) %in% c("PheC", "RXNO", "LOIN", "CCS-")
        ) %>%
        group_by(patient_num, feature_id) %>% 
        summarize(count=n(),.groups = "drop") 
      
      # ehr_count_all_long_rollup = ehr_count_all_long[ehr_count_all_long$feature_id %in% PheCodeCleaning(ehr_count_all_long$feature_id),]
      
      ehr_count_all = tidyr::spread(ehr_count_all_long, key = feature_id, value = count, fill = 0) 
      
      .cols_non0_bool <- apply(ehr_count_all %>% select(-patient_num), MARGIN = 2, function(x) sum(x==0) / length(x)) < .EHR_trim
      .cols_non0 = names(.cols_non0_bool)[.cols_non0_bool]
      .cols_non0_roll = PheCodeCleaning(.cols_non0)
      ehr_count_all <- ehr_count_all[, c("patient_num", .cols_non0_roll)] # T keeps patients_num
      
      .incl_fts <- ehr_count_all %>% select(-patient_num) %>% colnames
      .ehr_ft_summary = list(Fts_Considered = unique(ehr_count_all_long$feature_id)
                             , Ft_Obs = unique(ehr_count_all_long$feature_id)
                             , Ft_Obs_Trim = unique(.cols_non0)
                             , Ft_Obs_Trim_Roll = unique(.cols_non0_roll)
      )
    }
  
  ## NLP ####
    .incl_nlp_fts = .nlp_ft_summary =  NULL
    
    NLP_Utilization <- NLP_sub %>% count(patient_num, name = 'Utilization_NLP')
    
    N_MS_CUIs <- NLP_sub %>% 
      filter(feature_id=="C0026769") %>%
      count(patient_num, name = "N_MS_CUIs")
    # N_MS_CUIs$PATIENT_NUM <- N_MS_CUIs$patient_num
    # N_MS_CUIs$patient_num <- NULL
    
    if(.NLP){
      
      if(.NLP_all){
        NLP_High_Count <- NLP_sub %>% 
          group_by(patient_num, feature_id) %>% 
          summarize(count=n(),.groups = "drop")
        
        nlp_high_count = tidyr::spread(NLP_High_Count, key = feature_id, value = count, fill = 0) 
        
        .cols_non0 <- apply(nlp_high_count %>% select(-patient_num), MARGIN = 2, function(x) sum(x==0) / length(x)) < .EHR_trim
        nlp_high_count <- nlp_high_count[, c(T, .cols_non0)] # T keeps patients_num
  
        .incl_nlp_fts <- nlp_high_count %>% select(-patient_num) %>% colnames
        .nlp_ft_summary = list(Fts_Considered = unique(NLP_High_Count$feature_id)
                               , Fts_Observed = unique(NLP_High_Count$feature_id)
                               , Fts_Included_PostTrim = colnames(nlp_high_count)[2:ncol(nlp_high_count)]
                               )
                               
      }else{
        NLP_High_Count <- NLP_sub %>% 
          filter(feature_id %in% CUI$cui) %>% 
          group_by(patient_num, feature_id) %>% 
          summarize(count=n(),.groups = "drop")
        
        nlp_high_count = tidyr::spread(NLP_High_Count, key = feature_id, value = count, fill = 0) 
        
        .cols_non0 <- apply(nlp_high_count %>% select(-patient_num), MARGIN = 2, function(x) sum(x==0) / length(x)) < .EHR_trim
        nlp_high_count <- nlp_high_count[, c(T, .cols_non0)] # T keeps patients_num
  
        .incl_nlp_fts <- nlp_high_count %>% select(-patient_num) %>% colnames
        .nlp_ft_summary = list(Fts_Considered = unique(CUI$cui)
                               , Fts_Observed = unique(NLP_High_Count$feature_id)
                               , Fts_Included_PostTrim = colnames(nlp_high_count)[2:ncol(nlp_high_count)]
        )
        
      }
    }
    
  
  # Merging into df and normalizing #### 
  
  df_tot <- .analytic_df %>% left_join(NLP_Utilization, by="patient_num") %>% 
    mutate(Utilization_Total = N_AllCodes_PrePDDS + coalesce(Utilization_NLP, 0)) %>% 
    select(-c(N_AllCodes_PrePDDS, Utilization_NLP ))
  
  if(.EHR_all){
    ehr_fts <- .incl_fts
    if(.NLP){
      
      .tmp <- left_join(df_tot, ehr_count_all, by = "patient_num") %>%
        left_join(nlp_high_count, by="patient_num")
      
      analytic_df <- .tmp %>% 
        mutate(across(.cols = !!.incl_fts, ~ .x / Utilization_Total )) %>% 
        mutate(across(.cols = !!.incl_nlp_fts, ~ .x / Utilization_Total )) 
      
      
    }else{
      .tmp <- left_join(df_tot, ehr_count_all, by="patient_num")
      
      analytic_df <- .tmp %>% 
        mutate(across(.cols = !!.incl_fts, ~ .x / Utilization_Total)) 
      
    }
  }else{ # so EHR subset only
    ehr_fts <- .incl_fts
    
    if(.NLP){
      
      .tmp <- left_join(df_tot, ehr_count, by="patient_num") %>%
        left_join(nlp_high_count, by="patient_num")
      
      analytic_df <- .tmp %>% 
        mutate(across(.cols = !!.incl_fts, ~ .x / Utilization_Total )) %>% 
        mutate(across(.cols = !!.incl_nlp_fts, ~ .x / Utilization_Total )) 
      
    }else{
      .tmp <- left_join(df_tot, ehr_count, by="patient_num")
      
      analytic_df <- .tmp %>% 
        mutate(across(.cols = !!.incl_fts, ~ .x / Utilization_Total)) 
    }
    
  }
  
  if(.NLP){
    analytic_df[.incl_nlp_fts][is.na(analytic_df[.incl_nlp_fts])] <- 0
  }
    analytic_df[.incl_fts][is.na(analytic_df[.incl_fts])] <- 0

  ## Cleaning up X dataframe #### 
    
    X = analytic_df[,c("patient_num"
                       , "White_NonHispanic"
                       , "subject_sex_ch"
                       , "Age_at_PDDS"
                       , "Disease_Subtype"
                       , "N_MS_PheCodes"
                       , "Days_FirstPheCode_to_PDDS"
                       , "Days_First_MS_PheCode_to_PDDS"
                       , "Utilization_Total"
                       , .incl_fts
                       , .incl_nlp_fts
    )] %>% 
    left_join(N_MS_CUIs, "patient_num") %>% 
      select(-patient_num)

  X[is.na(X$N_MS_CUIs), "N_MS_CUIs"] <- 0
  
  X <- X %>% mutate(N_MS_CUIs = N_MS_CUIs / Utilization_Total
                    , N_MS_PheCodes = N_MS_PheCodes / Utilization_Total)
  
  X$subject_sex_ch = as.factor(X$subject_sex_ch=="M")
  X$Disease_Subtype_1 = as.numeric(X$Disease_Subtype==1)
  X$Disease_Subtype_2 = as.numeric(X$Disease_Subtype==2)
  X$Disease_Subtype = NULL
  
  
  
  # Removing redundancies, these are present as N_MS_PheCodes and N_MS_CUIs respectively
  if("PheCode:335" %in% colnames(X)) X$`PheCode:335` <- NULL
  if("C0026769" %in% colnames(X)) X$C0026769 <- NULL
  
  
  # Pre-processing #### 
  
  if(pre_process){
    X_preproc <- X %>% mutate(across(-c(White_NonHispanic
                                     , subject_sex_ch
                                     , Disease_Subtype_1
                                     , Disease_Subtype_2
    )
    , .fns = ~ (.x-mean(.x)) / sd(.x)
    )
    )
  }else{
    X_preproc <- X
  }
  
  # Output ####
  message("Output data.frame = analytic_df")
  message("Design matrix = X; Outcome vector Y; Dataframe object of [X,Y] in analytic_df")
  
  # if(.NLP){
    return(list(X = X_preproc, Y = analytic_df$score
                , analytic_df = analytic_df
                , EHR_fts = ehr_fts
                , NLP_fts = .incl_nlp_fts
                , .ehr_ft_summary = .ehr_ft_summary
                , .nlp_ft_summary = .nlp_ft_summary
                , input_df = .analytic_df))
  # }else{
  #   return(list(X = X_pdds, A = A, Y_avg = Y_avg, Y_change = Y_change, Y_NegChange = Y_NegChange
  #               , EHR_fts = ehr_fts, NLP_fts = NULL
  #               , input_df = .analytic_df))
  # }
  
}

