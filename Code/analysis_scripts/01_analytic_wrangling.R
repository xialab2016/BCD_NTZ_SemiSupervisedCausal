
ehr <- read_csv(here("Data", "Latest"
                     , ehr_codified_rolled)
                , show_col_types = FALSE)

ehr$patient_num <- as.character(ehr$patient_num)


# cleaning phecodes using prescribed rules 
# easiest to do here and not track throughout analytic wrangling 
#### Just keep parent code if one digit does not appear; remove two digits ####
# feat.set = unique(ehr$feature_id)
# phecode.all = feat.set[grep("PheCode", feat.set)]
# tmp = matrix(unlist(lapply(strsplit(phecode.all, split = "\\.")
#                            , function(x){cbind(len=length(x), numdigit=nchar(x[2]))})
#                     )
#              , ncol=2, byrow = TRUE)
# 
# 
# colnames(tmp)= c("len", "numdigit")
# 
# phecode.nonedigit = phecode.all[tmp[,"len"] == 1]
# phecode.onedigit = phecode.all[which(tmp[,"numdigit"] == 1)]
# phecode.onedigit.convert.none = unlist(lapply(strsplit(phecode.onedigit, split="\\."), function(x){x[1]}))
# # phecode.twodigit = phecode.all[which(tmp[,"numdigit"] == 2)] # not retained, do not need to track 
# 
# feat.set = c(phecode.onedigit # all one digit codes 
#              , setdiff(phecode.nonedigit, phecode.onedigit.convert.none) # all parent codes without a one digit child code 
#              , setdiff(feat.set, phecode.all) # all non PheCode features 
#              )

  # creating as custom function  
  phecode.retain = function(ft.set){
    # function to apply PheCode roll-up procedure, for a given input list of all observed PheCodes
    # to be applied to codified features AFTER limiting to observed, non-sparse set of codified features 
    
    phecode.all = ft.set[grep("PheCode", ft.set)]
    
    
    tmp = matrix(unlist(lapply(strsplit(phecode.all, split = "\\.")
                               , function(x){cbind(len=length(x), numdigit=nchar(x[2]))})
    )
    , ncol=2, byrow = TRUE)
    
    colnames(tmp)= c("len", "numdigit")
    
    phecode.nonedigit = phecode.all[tmp[,"len"] == 1]
    phecode.onedigit = phecode.all[which(tmp[,"numdigit"] == 1)]
    phecode.onedigit.convert.none = unlist(lapply(strsplit(phecode.onedigit, split="\\."), function(x){x[1]}))
    phecode.twodigit = phecode.all[which(tmp[,"numdigit"] == 2)] # not retained, do not need to track
    phecode.twodigit.convert.none = unlist(lapply(strsplit(phecode.twodigit, split="\\."), function(x){x[1]}))
    
    phecode.none.keep = setdiff(phecode.nonedigit, phecode.onedigit.convert.none) # all parent codes without a one digit child code 
    # lastly keeping any "orphan" two digit codes, which don't hve a one or no digit parent 
    phecode.twodigit.keep = phecode.twodigit[phecode.twodigit.convert.none %nin% c(phecode.onedigit.convert.none, phecode.nonedigit)]
    
    phecode.retain = c(phecode.onedigit # all one digit codes, keep by default 
                       , phecode.none.keep # parent, no-digit codes with no children one-digit codes 
                       , phecode.twodigit.keep
                       , setdiff(ft.set, phecode.all) # all nonPheCode fts 
                       )
    
    return(phecode.retain)
    # return(ft.set)
                                  
  }

  # uncomment if applying roll-up BEFORE feature selection/sparsity check 
  # ehr = ehr %>% filter(feature_id %in% phecode.retain(ft.set = unique(ehr$feature_id)))
  
  # length(unique(ehr$feature_id[substr(ehr$feature_id, 1, 4)=="PheC"]))
  # length(unique(ehr.analytic$feature_id[substr(ehr.analytic$feature_id, 1, 4)=="PheC"]))
  # removes ~1/3 of PheCodes, seems reasonable 


source(here("Code", "00_shared_setup.R")) # for nlp_filename

ms_dmt = read_excel(here("Data", "Latest", "MS Drugs Mechanism Citations Efficacy 20230107.xlsx")) 
ms_rxnorms = paste0("RXNORM:", ms_dmt$`RxNorm Ingredient id`)

once_ms_ehr = once_ms_ehr_full %>% # read_csv(here("Data", "Latest", "ONCE_Multiple sclerosis_Codified.csv")) %>% 
  filter(expanded_features==1) %>%
  select(term = Variable
         , expanded_features
         , phenotyping_features) %>% 
  mutate(set = "MS") %>% 
  filter(term != 'PheCode:335' 
         & !str_detect(tolower(term), pattern = "local|Local|lab|Lab")
         & term %nin% ms_rxnorms
  ) 

once_ms_nlp = once_ms_nlp_full %>% # read_csv(here("Data", "Latest", "ONCE_Multiple sclerosis_CUI.csv")) %>% 
  filter(expanded_features==1) %>%
  select(cui
         , expanded_features
         , phenotyping_features) %>% 
  mutate(set = "MS") %>% 
  filter(cui %nin% c("C0528175", "C1529600"))

once_disability_ehr = once_disability_ehr_full %>%  #read_csv(here("Data", "Latest", "ONCE_disability_PheCode296.2_cos0.165.csv")) %>%  
  filter(phenotyping_features == 1) %>%
  select(term = Variable
         , expanded_features
         , phenotyping_features) %>% 
  mutate(set = "disability") %>% 
  filter(term != 'PheCode:335' 
         & !str_detect(tolower(term), pattern = "local|Local|lab|Lab")
         & term %nin% ms_rxnorms
  ) 

once_disability_nlp = once_disability_nlp_full %>% #read_csv(here("Data", "Latest", "ONCE_disability_C0231170_titlecos0.5_titlecut0.3_exactFALSE.csv")) %>%  
  filter(phenotyping_features == 1) %>%
  select(cui
         , expanded_features
         , phenotyping_features) %>% 
  mutate(set = "disability") %>% 
  filter(cui %nin% c("C0528175", "C1529600"))


CUI = rbind(once_disability_nlp, once_ms_nlp)
ehr_ft_list = rbind(once_disability_ehr, once_ms_ehr)
# setdiff(ehr_ft_list$term, feat.set) # should be non-empty, filtering by roll-up procedure applied later after sparsity trimming 

analytic_wrangling <- function(.analytic_df
                               , pre_process = T
                               , efficacy = NULL
                               , .EHR_trim = 0.8
                               , .NLP = T
                               , .EHR_all = F
                               , .NLP_all = F
                               , .RXNORM_Exclude_All = T
                               , .DMT_reference = NULL
                               , .diff_window = dmonths(6)
                               , .lookback_window = Inf
){
  
  # Set-Up #### 
  if(is.null(.DMT_reference)){
    stop("Enter DMT reference  group (character object)")  
  }
  
  if(is.null(efficacy)){
    stop("No efficacy argument specified. Please enter high/H or standard/S/std (case insensitive)")
  }
  
  .analytic_df$subject_sex_ch = as.factor(.analytic_df$subject_sex_ch)
  .analytic_df$patient_num = as.numeric(.analytic_df$PATIENT_NUM)
  
  
  df = .analytic_df[,c("PATIENT_NUM"
                       , "id_participant"
                       , "DMT"
                       , "White_NonHispanic"
                       , "subject_sex_ch"
                       , "Age_at_DMT_Initiation"
                       , "Disease_Subtype"
                       , "Utilization"
                       , "Days_Earliest_DMT_to_Study_DMT"
                       , "N_MS_PheCodes"
                       , "Followup_Duration", "Disease_Duration"
                       , "Average_PDDS"
                       , "PDDS_Diff", "PDDS_Diff_Date"
                       , "PDDS_Change", "PDDS_NegChange", "DMT_Study_Start")] #%>% 
    # filter(!is.na(PDDS_Change))
  
  if(nrow(df[is.na(df$Days_Earliest_DMT_to_Study_DMT) | df$Days_Earliest_DMT_to_Study_DMT<0,])>0){
    stop("Non-positive alues of Prior DMT Duration detected, check code or make this a warning to clean")
    df$Days_Earliest_DMT_to_Study_DMT[is.na(df$Days_Earliest_DMT_to_Study_DMT)] <- 0
    df$Days_Earliest_DMT_to_Study_DMT[df$Days_Earliest_DMT_to_Study_DMT < 0 ] <- 0
  }
  
  if(nrow(df[is.na(df$Disease_Duration) | df$Disease_Duration<0,])>0){
    stop("Non-positive alues of Disease Duration detected, check code or make this a warning to clean")
    df$Disease_Duration[df$Disease_Duration < 0 ] <- 0
  }
  
  df$DMT = as.factor(df$DMT)
  df$subject_sex_ch = as.factor(df$subject_sex_ch)
  # df$PATIENT_NUM = as.character(df$PATIENT_NUM)
  df$patient_num = as.character(df$PATIENT_NUM)
  
  
  # EHR features ####
  
  if (!(exists("ehr", where=.GlobalEnv))){
    stop("EHR data not present in R environment. Import EHR features using `ehr_codified_rolled`")
  }
  
  ## Codified/Correlated Subset ####
  if(!.EHR_all){
    feat.set.thresholded = ehr_ft_list 
    ehr.set <- ehr %>% filter(feature_id %in% feat.set.thresholded$term)

    ehr_valid <- ehr.set %>% 
      inner_join(df[,c("patient_num", "DMT_Study_Start")], by="patient_num") %>% 
      filter(start_date < DMT_Study_Start
              & start_date >= DMT_Study_Start - .lookback_window)

    if(.RXNORM_Exclude_All){
      .ehr_ft_colnames <- setdiff(unique(ehr_valid$feature_id), .rxnorm_exclude)
    }else{  
      .ehr_ft_colnames <- unique(ehr_valid$feature_id)
    }
    
    ehr_valid$DMT_Study_Start <- NULL
    
    ehr_count_long <- ehr_valid %>%
      filter(feature_id %in% .ehr_ft_colnames
             & substr(feature_id, 1, 4) %in% c("PheC", "RXNO", "LOIN", "CCS-")
      ) %>%
      group_by(patient_num, feature_id) %>% 
      summarize(count=n(),.groups = "drop")
    
    ehr_count = tidyr::spread(ehr_count_long, key = feature_id, value = count, fill = 0) 
    
    .cols_non0 <- apply(ehr_count %>% select(-patient_num), MARGIN = 2, function(x) sum(x==0) / length(x)) < .EHR_trim
    # phecode.retain(names(which(.cols_non0)))
    # setdiff(names(which(.cols_non0)), phecode.retain(names(which(.cols_non0))))
    .cols_non0_keep = as.logical(.cols_non0 * (names(.cols_non0) %in% phecode.retain(names(which(.cols_non0))) ))
    
    ehr_count <- ehr_count[, c(T, .cols_non0_keep)] # T keeps patients_num
    
    .incl_fts <- ehr_count %>% select(-patient_num) %>% colnames
    .ehr_ft_summary = list(Fts_Considered = feat.set.thresholded$term
                           , Fts_Observed = unique(ehr_count_long$feature_id)
                           , Fts_Included_PostTrim = colnames(ehr_count)[2:ncol(ehr_count)]
                           )
    # apply(ehr_count %>% select(-patient_num), MARGIN = 2, function(x) sum(x==0) / length(x)) < 0.9
  }else{
    ## All Observed EHR Fts #### 
    
    ehr.pts <- ehr %>% 
      inner_join(df[,c("patient_num", "DMT_Study_Start")], by="patient_num") %>% 
      filter(start_date < DMT_Study_Start
              & start_date >= DMT_Study_Start - .lookback_window)    # %>% 
    # filter(feature_id %nin% c(paste0("RXNORM:", HighEff_Exclude_RxNormIDs)
    #                           , paste0("RXNORM:", StdEff_Exclude_RxNormIDs))
    #        )
    
    if(.RXNORM_Exclude_All){
      .ehr_ft_colnames <- setdiff(unique(ehr.pts$feature_id), .rxnorm_exclude)
    }else{  
      .ehr_ft_colnames <- unique(ehr.pts$feature_id)
    }
    
    ehr.pts$DMT_Study_Start <- NULL
    
    ehr_count_all_long <- ehr.pts %>%
      filter(feature_id %in% .ehr_ft_colnames
             & substr(feature_id, 1, 4) %in% c("PheC", "RXNO", "LOIN", "CCS-")
      ) %>%
      group_by(patient_num, feature_id) %>% 
      summarize(count=n(),.groups = "drop") 
    
    ehr_count_all = tidyr::spread(ehr_count_all_long, key = feature_id, value = count, fill = 0) 
    
    .cols_non0 <- apply(ehr_count_all %>% select(-patient_num), MARGIN = 2, function(x) sum(x==0) / length(x)) < .EHR_trim
    .cols_non0_keep = as.logical(.cols_non0 * (names(.cols_non0) %in% phecode.retain(names(which(.cols_non0))) ))
    
    ehr_count_all <- ehr_count_all[, c(T, .cols_non0_keep)] # T keeps patients_num
    
    .incl_fts_all <- ehr_count_all %>% select(-patient_num) %>% colnames
    
    .ehr_ft_summary = list(Fts_Considered = unique(ehr_count_all_long$feature_id)
                           , Fts_Observed = unique(ehr_count_all_long$feature_id)
                           , Fts_Included_PostTrim = colnames(ehr_count_all)[2:ncol(ehr_count_all)]
    )
    
    # apply(ehr_count %>% select(-patient_num), MARGIN = 2, function(x) sum(x==0) / length(x)) < 0.9
  }
  
  # NLP Prep ####
  
  if("NLP" %nin% ls(envir = .GlobalEnv)){
    print("Importing NLP data to Global environment (bottleneck step, ideally set .import_NLP to F and previously import NLP file to avoid redundant, intensive import step")
    NLP <- read_csv(here("Data", "Latest", nlp_filename))
    NLP = setDT(NLP) # fread function not completing import, read_csv preferred for now 
    assign("NLP", NLP, envir=.GlobalEnv)
  }
  
  NLP_cui = NLP[as.data.table(CUI), on = .(feature_id==cui), nomatch = NULL] # inner_join
  NLP_pts = NLP_cui[as.data.table(.analytic_df), on = .(patient_num), nomatch = NULL, allow.cartesian = T] # inner_join
  NLP_sub = NLP_pts[start_date < DMT_Study_Start 
                    & start_date >= DMT_Study_Start - .lookback_window
                    , .(patient_num, start_date, feature_id) , ]    
  
  rm(NLP_cui, NLP_pts)
  NLP_sub$patient_num = as.character(NLP_sub$patient_num)
  
  NLP_Utilization <- NLP_sub %>% count(patient_num, name = 'Utilization_NLP')
  
  N_MS_CUIs <- NLP_sub %>% 
    filter(feature_id=="C0026769") %>%
    count(patient_num)
  # N_MS_CUIs$PATIENT_NUM <- N_MS_CUIs$patient_num
  # N_MS_CUIs$patient_num <- NULL
  
  if(.NLP){
    
    if(.NLP_all){
      NLP_High_Count_long <- NLP_sub %>% 
        group_by(patient_num, feature_id) %>% 
        summarize(count=n(),.groups = "drop")
      
      nlp_high_count = tidyr::spread(NLP_High_Count_long, key = feature_id, value = count, fill = 0) 
      
      .cols_non0 <- apply(nlp_high_count %>% select(-patient_num), MARGIN = 2, function(x) sum(x==0) / length(x)) < .EHR_trim
      nlp_high_count <- nlp_high_count[, c(T, .cols_non0)] # T keeps patients_num
      nlp_high_count$patient_num <- as.character(nlp_high_count$patient_num)
      
      .incl_nlp_fts <- nlp_high_count %>% select(-patient_num) %>% colnames
      .nlp_ft_summary = list(Fts_Considered = unique(NLP_High_Count_long$feature_id)
                          , Fts_Observed = unique(NLP_High_Count_long$feature_id)
                          , Fts_Included = colnames(nlp_high_count)[2:ncol(nlp_high_count)]
                          )
    }else{
      .selected_cuis = unique(CUI$cui)
      
      NLP_High_Count_long <- NLP_sub %>% 
        filter(feature_id %in% .selected_cuis) %>% 
        group_by(patient_num, feature_id) %>% 
        summarize(count=n(),.groups = "drop")
      
      nlp_high_count = tidyr::spread(NLP_High_Count_long, key = feature_id, value = count, fill = 0) 
      
      .cols_non0 <- apply(nlp_high_count %>% select(-patient_num), MARGIN = 2, function(x) sum(x==0) / length(x)) < .EHR_trim
      nlp_high_count <- nlp_high_count[, c(T, .cols_non0)] # T keeps patients_num
      nlp_high_count$patient_num <- as.character(nlp_high_count$patient_num)
      
      .incl_nlp_fts <- nlp_high_count %>% select(-patient_num) %>% colnames
      .nlp_ft_summary = list(Fts_Considered = unique(.selected_cuis)
                          , Fts_Observed = unique(NLP_High_Count_long$feature_id)
                          , Fts_Included = colnames(nlp_high_count)[2:ncol(nlp_high_count)]
                          )
    }
  }
  
  
  # Merging into df and normalizing #### 
  
  df_tot <- df %>% left_join(NLP_Utilization, by="patient_num") %>% 
    mutate(Utilization_Total = Utilization + coalesce(Utilization_NLP, 0)) %>% 
    select(-c(Utilization, Utilization_NLP ))
  
  if(.EHR_all){
    ehr_fts <- .incl_fts_all
    if(.NLP){
      
      .analytic_df <- left_join(df_tot, ehr_count_all, by = "patient_num") %>%
        left_join(nlp_high_count, by="patient_num")
      
      analytic_df <- .analytic_df %>% 
        mutate(across(.cols = !!.incl_fts_all, ~ .x / Utilization_Total )) %>% 
        mutate(across(.cols = !!.incl_nlp_fts, ~ .x / Utilization_Total )) 
      
      
    }else{
      .analytic_df <- left_join(df_tot, ehr_count_all, by="patient_num")
      
      analytic_df <- .analytic_df %>% 
        mutate(across(.cols = !!.incl_fts_all, ~ .x / Utilization_Total)) 
      
    }
  }else{ # so EHR subset only
    ehr_fts <- .incl_fts
    
    if(.NLP){
      
      .analytic_df <- left_join(df_tot, ehr_count, by="patient_num") %>%
        left_join(nlp_high_count, by="patient_num")
      
      analytic_df <- .analytic_df %>% 
        mutate(across(.cols = !!.incl_fts, ~ .x / Utilization_Total )) %>% 
        mutate(across(.cols = !!.incl_nlp_fts, ~ .x / Utilization_Total )) 
      
    }else{
      .analytic_df <- left_join(df_tot, ehr_count, by="patient_num")
      
      analytic_df <- .analytic_df %>% 
        mutate(across(.cols = !!.incl_fts, ~ .x / Utilization_Total)) 
    }
    
  }
  
  if(.NLP){
    analytic_df[.incl_nlp_fts][is.na(.analytic_df[.incl_nlp_fts])] <- 0
  }
  if(.EHR_all){
    analytic_df[.incl_fts_all][is.na(.analytic_df[.incl_fts_all])] <- 0
  } else{
    analytic_df[.incl_fts][is.na(.analytic_df[.incl_fts])] <- 0
  }
  
  A = as.numeric(analytic_df$DMT != .DMT_reference)
  Y_avg = as.numeric(analytic_df$Average_PDDS)
  Y_change = as.numeric(analytic_df$PDDS_Change)
  Y_NegChange = as.numeric(analytic_df$PDDS_NegChange)
  Y_diff = as.numeric(analytic_df$PDDS_Diff)  
  Y_diff[analytic_df$PDDS_Diff_Date - analytic_df$DMT_Study_Start > .diff_window] = NA
  Y_AnyIncrease = as.numeric(analytic_df$PDDS_Diff>=1)  
  Y_AnyIncrease[analytic_df$PDDS_Diff_Date - analytic_df$DMT_Study_Start > .diff_window] = NA
  Y_AnyDecrease = as.numeric(analytic_df$PDDS_Diff<=-1)  
  Y_AnyDecrease[analytic_df$PDDS_Diff_Date - analytic_df$DMT_Study_Start > .diff_window] = NA
  
  ## Cleaning up X dataframe #### 
  X = analytic_df %>%
    left_join(N_MS_CUIs, "patient_num") %>% 
    select(-c("patient_num", "PATIENT_NUM","id_participant","DMT",
              "Average_PDDS", "PDDS_Change", "PDDS_NegChange", "DMT_Study_Start"
              , "PDDS_Diff", "PDDS_Diff_Date")
           , N_MS_CUIs = n)
  
  X[is.na(X$N_MS_CUIs), "N_MS_CUIs"] <- 0
  
  
  X <- X %>% mutate(N_MS_CUIs = N_MS_CUIs / Utilization_Total
                    , N_MS_PheCodes = N_MS_PheCodes / Utilization_Total)
  
  # X[X$Utilization_Total==0,]$N_MS_CUIs = 0
  # X[X$Utilization_Total==0,]$N_MS_PheCodes = 0

  X$subject_sex_ch = as.numeric(X$subject_sex_ch=="M")
  X$Disease_Subtype_1 = as.numeric(X$Disease_Subtype==1)
  X$Disease_Subtype_2 = as.numeric(X$Disease_Subtype==2)
  X$Disease_Subtype = NULL
  
  
  
  # Removing redundancies, these are present as N_MS_PheCodes and N_MS_CUIs respectively
  if("PheCode:335" %in% colnames(X)) X$`PheCode:335` <- NULL
  if("C0026769" %in% colnames(X)) X$C0026769 <- NULL
  
  
  
  
  # Pre-processing #### 
  
  if(pre_process){
    X_pdds <- X %>% mutate(across(-c(White_NonHispanic
                                     , subject_sex_ch
                                     , Disease_Subtype_1
                                     , Disease_Subtype_2
    )
    , .fns = ~ (.x-mean(.x)) / sd(.x)
    )
    )
  }else{
    X_pdds <- X
  }
  
  # for searching subtype disease subtype ms subtype
  X_pdds$Disease_Subtype_1 = NULL
  X_pdds$Disease_Subtype_2 = NULL
  
  # Output ####
  message("Output data.frame = analytic_df")
  message("Design matrix = X; Treatment vector  = A; Outcome vector = Y_avg, Y_change")
  
  if(.NLP){
    return(list(X = X_pdds, A = A
                , Y_avg = Y_avg
                , Y_change = Y_change
                , Y_NegChange = Y_NegChange
                , Y_diff = Y_diff
                , Y_AnyIncrease = Y_AnyIncrease
                , Y_AnyDecrease = Y_AnyDecrease
                , EHR_fts = ehr_fts, NLP_fts = .incl_nlp_fts
                , .ehr_ft_summary = .ehr_ft_summary
                , .nlp_ft_summary = .nlp_ft_summary
                , input_df = .analytic_df))
  }else{
    return(list(X = X_pdds, A = A
                , Y_avg = Y_avg
                , Y_change = Y_change
                , Y_NegChange = Y_NegChange
                , Y_diff = Y_diff
                , Y_AnyIncrease = Y_AnyIncrease
                , Y_AnyDecrease = Y_AnyDecrease
                , EHR_fts = ehr_fts, NLP_fts = NULL
                , .ehr_ft_summary = .ehr_ft_summary
                , .nlp_ft_summary = NULL
                , input_df = .analytic_df))
  }
  
}

