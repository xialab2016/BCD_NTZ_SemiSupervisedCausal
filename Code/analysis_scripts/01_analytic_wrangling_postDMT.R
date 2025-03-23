
if("ehr_2004_2022" %in% ls()){
  ehr <- ehr_2004_2022
  ehr$patient_num <- as.character(ehr$patient_num)
}else if ("ehr" %nin% ls(all.names=T)){
  ehr <- read_csv(here("Data", "Latest"
                       , ehr_codified_rolled)
                  , show_col_types = FALSE)
  ehr$patient_num <- as.character(ehr$patient_num)
}

# cleaning phecodes using prescribed rules 
# easiest to do here and not track throughout analytic wrangling 
#### Just keep parent code if one digit does not appear; remove two digits ####
feat.set = unique(ehr$feature_id)
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

feat.set = c(phecode.onedigit # all one digit codes 
             , setdiff(phecode.nonedigit, phecode.onedigit.convert.none) # all parent codes without a one digit child code 
             , setdiff(feat.set, phecode.all) # all non PheCode features 
)

# ehr.analytic = ehr[ehr$feature_id %in% feat.set,]
# ehr = ehr[ehr$feature_id %in% feat.set,]

# length(unique(ehr$feature_id[substr(ehr$feature_id, 1, 4)=="PheC"]))
# length(unique(ehr.analytic$feature_id[substr(ehr.analytic$feature_id, 1, 4)=="PheC"]))
# removes ~1/3 of PheCodes, seems reasonable 

source(here("Code", "00_shared_setup.R")) # for nlp_filename


analytic_wrangling_postDMT <- function(object
                                       , ehr_ft_list
                                       , nlp_ft_list
                                       , pre_process = T
                                       , .lag = lubridate::dyears(1)
                                       , .lookback_window = Inf
){

  .n_nlp = length(nlp_ft_list) 
  
  df = object[,c("PATIENT_NUM", "id_participant", "DMT_Study_Start", "Utilization"
                 , "White_NonHispanic", "subject_sex_ch", "Age_at_DMT_Initiation"
                 , "Days_Earliest_DMT_to_Study_DMT", "Followup_Duration", "Disease_Duration")]  
  df$patient_num = as.character(df$PATIENT_NUM)
  
  # EHR features ####
  
  if (!(exists("ehr", where=.GlobalEnv))){
    stop("EHR data not present in R environment. Import EHR features using `ehr_codified_rolled`")
  }
  
  ## Codified/Correlated Subset ####
  
    # extracting updated N_MS_PheCode count 
    .n_ms_phecode_df = ehr %>%
      inner_join(df[,c("patient_num", "DMT_Study_Start")], by="patient_num") %>% 
      filter(feature_id == "PheCode:335" 
             & start_date < DMT_Study_Start + .lag
             & start_date >= DMT_Study_Start - .lookback_window
             ) %>% 
      distinct(patient_num, start_date) %>% 
      group_by(patient_num) %>% 
      count()

    ehr.set <- ehr %>% filter(feature_id %in% ehr_ft_list)

    ehr_valid <- ehr.set %>% 
      inner_join(df[,c("patient_num", "DMT_Study_Start")], by="patient_num") %>% 
      filter(start_date < DMT_Study_Start + .lag
             & start_date >= DMT_Study_Start - .lookback_window)
    
    ehr_valid$DMT_Study_Start <- NULL
    
    ehr_count <- ehr_valid %>%
      filter(feature_id %in% ehr_ft_list
             & substr(feature_id, 1, 4) %in% c("PheC", "RXNO", "LOIN", "CCS-")
      ) %>%
      group_by(patient_num, feature_id) %>% 
      summarize(count=n(),.groups = "drop")
    
    ehr_count = tidyr::spread(ehr_count, key = feature_id, value = count, fill = 0) 
    
    
  # NLP Prep ####
  if(.n_nlp){
    if("NLP" %nin% ls(envir = .GlobalEnv)){
      print("Importing NLP data to Global environment (bottleneck step, ideally set .import_NLP to F and previously import NLP file to avoid redundant, intensive import step")
      NLP <- read_csv(here("Data", "Latest", nlp_filename))
      NLP = setDT(NLP) # fread function not completing import, read_csv preferred for now 
      assign("NLP", NLP, envir=.GlobalEnv)
    }
    
    dflag = df
    dflag$DMT_Study_Start = dflag$DMT_Study_Start + .lag
    dflag$patient_num = as.numeric(dflag$patient_num)
    
    NLP_cui = NLP[data.table(cui=nlp_ft_list), on = .(feature_id==cui), nomatch = NULL] # inner_join
    NLP_pts = NLP_cui[as.data.table(dflag), on = .(patient_num), nomatch = NULL, allow.cartesian = T] # inner_join
    NLP_sub = NLP_pts[start_date < DMT_Study_Start + .lag
                      & start_date >= DMT_Study_Start - .lookback_window
                      , .(patient_num, start_date, feature_id) ,]    
    
    rm(NLP_cui, NLP_pts)
    NLP_sub$patient_num = as.character(NLP_sub$patient_num)
    
    NLP_Utilization <- NLP_sub %>% count(patient_num, name = 'Utilization_NLP')
    
    N_MS_CUIs <- NLP_sub %>% 
      filter(feature_id=="C0026769") %>%
      count(patient_num)
    # N_MS_CUIs$PATIENT_NUM <- N_MS_CUIs$patient_num
    # N_MS_CUIs$patient_num <- NULL
  
      
      NLP_High_Count <- NLP_sub %>% 
        filter(feature_id %in% nlp_ft_list) %>% 
        group_by(patient_num, feature_id) %>% 
        summarize(count=n(),.groups = "drop")
      
      nlp_high_count = tidyr::spread(NLP_High_Count, key = feature_id, value = count, fill = 0) 
      
      nlp_high_count$patient_num <- as.character(nlp_high_count$patient_num)
  }
  # Merging into df and normalizing #### 
  
  df_tot <- df %>% left_join(NLP_Utilization, by="patient_num") %>% 
    mutate(Utilization_Total = Utilization + coalesce(Utilization_NLP, 0)) %>% 
    select(-c(Utilization, Utilization_NLP ))
  
  .analytic_df <- left_join(df_tot, ehr_count, by = "patient_num") %>%
    left_join(nlp_high_count, by="patient_num")
      
  analytic_df <- .analytic_df %>% 
    mutate(across(.cols = !!ehr_ft_list, ~ .x / Utilization_Total )) %>% 
    mutate(across(.cols = !!nlp_ft_list, ~ .x / Utilization_Total )) 
      
    analytic_df[nlp_ft_list][is.na(analytic_df[nlp_ft_list])] <- 0
    analytic_df[ehr_ft_list][is.na(analytic_df[ehr_ft_list])] <- 0

  
  ## Cleaning up X dataframe #### 
  X = analytic_df %>%
    left_join(N_MS_CUIs, "patient_num") %>% 
    left_join(.n_ms_phecode_df, "patient_num") %>% 
    select(-c("patient_num", "id_participant")
           , N_MS_CUIs = n.x
           , N_MS_PheCodes = n.y)
  
    X[is.na(X$N_MS_CUIs), "N_MS_CUIs"] <- 0
    X[is.na(X$N_MS_PheCodes), "N_MS_PheCodes"] <- 0
    
  
  X <- X %>% mutate(N_MS_CUIs = N_MS_CUIs / Utilization_Total
                    , N_MS_PheCodes = N_MS_PheCodes / Utilization_Total)
  
  # X[X$Utilization_Total==0,]$N_MS_CUIs = 0
  # X[X$Utilization_Total==0,]$N_MS_PheCodes = 0
  
  
  # Removing redundancies, these are present as N_MS_PheCodes and N_MS_CUIs respectively
  if("PheCode:335" %in% colnames(X)) X$`PheCode:335` <- NULL
  if("C0026769" %in% colnames(X)) X$C0026769 <- NULL
  
  
  # Pre-processing ####
  X$subject_sex_ch = as.integer(X$subject_sex_ch=="F")
  
  if(pre_process){
    X_pdds <- X %>% mutate(across(-c(PATIENT_NUM, DMT_Study_Start, White_NonHispanic, subject_sex_ch)
                                  , .fns = ~ (.x-mean(.x)) / sd(.x)))
  }else{
    X_pdds <- X
  }
  

  # Output ####
  message("Output data.frame = analytic_df")
  message("Design matrix = X; Treatment vector  = A; Outcome vector = Y_avg, Y_change")
  warning("Data frame contains identifiers, necessary for PDDS merging but must be removed prior to analysis: PATIENT_NUM, DMT_Study_Start")
  return(X_pdds)

}

