# PDDS Imputation Cohort 
## Dominic DiSanto

if("pacman" %in% installed.packages()[,1]){
  library(pacman)
}else{
  install.packages("pacman")
  library(pacman)
}

p_load(here)


source(here("Code", "setup.R"))

# Data Import #### 

    ehr <- read_csv(here("Data", "Latest"
                         , ehr_codified_rolled )
                    , col_names = T
    ) %>% 
      mutate(PATIENT_NUM = as.character(patient_num))
  
  
  dup_rm <- read_excel(here("Data", "ID Linkage", 
                            "Duplicate_IDs_Reconciled.xlsx")
                       , sheet = "Remove") %>% 
    mutate(RemoveTag=1
           , PATIENT_NUM=as.character(PATIENT_NUM)) %>% 
    select(-c(n_id_part, n_pt_num))
  
  id_link_nodups <- id_link %>% 
    left_join(dup_rm, by=c("PATIENT_NUM", "id_participant")) %>%
    filter(RemoveTag!=1 | is.na(RemoveTag)) %>% 
    filter(!is.na(PATIENT_NUM)) %>% # 6 id_participant values with no PATIENT_NUM
    select(-RemoveTag) %>% 
    distinct()
  
  ## PDDS Import #### 
  
  .pdds <- read_csv(here("Data", "Latest", "pdds_all.csv"), col_select = -`...1`) %>% 
    filter(!is.na(id_participant)) %>% # 3 obs with no identifier 
    left_join(id_link, by="id_participant") %>% 
    select(PATIENT_NUM, id_participant, date, score)
  
  pdds_MayPull <- read_excel(here("Data", "Latest", "MS_Annotation_PDDS_Final_20230502.xlsx")) %>%
    mutate(id_participant = NA
           , PATIENT_NUM = PERSON_ID) %>%
    select(PATIENT_NUM 
           , id_participant
           , date = `PDDS Date (if NA, no \r\nPDDS scores available\r\n in Epic)`
           , score = `PDDS Score (if NA, no \r\nPDDS scores available\r\n in Epic)`
    ) %>%
    mutate(date = as_date(date)) %>% 
    filter(!is.na(PATIENT_NUM))
  
  pdds_pull1 <- read_excel(here("Data", "data_pulls", "fromUPMC", "1. PDDS Treatment of Interest (BCD vs NTZ).xlsx")
                           , guess_max=1e6) %>% 
    select(PATIENT_NUM
           , id_participant
           , date = pdds_date
           , score = pdds_score ) %>% 
    mutate(score = as.double(score)
           , date = as.Date(as.integer(date), origin="1899/12/30")
    )
  
  # if import fails, check date column for "NA" values (forces to string w/ irrecoverable dates)
  
  pdds_pull2 <- read_csv(here("Data", "data_pulls", "fromUPMC", "2. PRT PDDS 2023-06-23.csv")) %>% 
    left_join(id_link, by = "id_participant") %>% 
    select(PATIENT_NUM
           , id_participant
           , date 
           , score 
    ) %>% 
    filter(id_participant!=1) %>% # Not valid ID, removing participant (removed for missing demographics as well, this just matters for accurate derivation summary)
    mutate(id_participant = gsub("_x", "", id_participant))
  
  pdds_pull4 <- read_excel(here("Data", "data_pulls", "fromUPMC", "4. MS Annotation PDDS Final.xlsx")
                           , guess_max = 1e5) %>% 
    select(PATIENT_NUM = PERSON_ID
           , everything()
    ) %>% 
    left_join(id_link, by = "PATIENT_NUM") %>%
    select(PATIENT_NUM
           , id_participant
           , date = `PDDS Date (if NA, no \r\nPDDS scores available\r\n in Epic)`
           , score = `PDDS Score (if NA, no \r\nPDDS scores available\r\n in Epic)`
    )  %>% 
    mutate(score = as.double(score)
           , date = as.Date(as.integer(date), origin="1899/12/30")
    ) %>% 
    filter(!is.na(date))
  
  pdds_all <- rbind(.pdds
                , pdds_MayPull
                , pdds_pull1 
                , pdds_pull2
                , pdds_pull4
                ) %>% 
    filter(year(date)<2024) %>% # some missing dates, two erroneous 2028 dates
    left_join(dup_rm, c("PATIENT_NUM", "id_participant")) %>% 
    filter(is.na(RemoveTag)) %>%
    select(-RemoveTag) %>%
    # filter(!is.na(PATIENT_NUM)) %>%
    distinct() # pdds files contain duplicated records 

  
  # all(pdds_all$id_participant[is.na(pdds_all$PATIENT_NUM)] %in% id_link_nodups$id_participant)
  # sum(pdds_all$id_participant[is.na(pdds_all$PATIENT_NUM)] %nin% id_link_nodups$id_participant)
  
  
# Cohort Derivation ####
  
  
  ## 1 - Initial Cohort #### 
      # Any patient with a valid PDDS score 
      # Specifically month-windows for unique PATIENT_NUM-date observations
      # Allows for multiple observations for a given patient, as long as they are 1+ months apart
  
  cohort1 = pdds_all %>%
    filter(!is.na(PATIENT_NUM)) %>% 
    group_by(PATIENT_NUM, id_participant) %>% 
    arrange(PATIENT_NUM, id_participant, date) %>% 
    mutate(.ddiff = date - lag(date)) %>% 
    mutate(.calc_date = case_when(.ddiff>.pdds_window | is.na(.ddiff) ~ date
                                  , T ~ NA_Date_
                                  )) %>% 
    fill(.calc_date, .direction = "down") %>% 
    mutate(.ddiff2 = date - lag(.calc_date)) %>% 
    filter(.ddiff2 > .pdds_window | is.na(.ddiff2)) %>% 
    mutate(.calc_date = case_when(row_number()>1 ~ .calc_date)) %>% 
    group_by(PATIENT_NUM, id_participant, .calc_date) %>% 
    filter(row_number()==1) %>% 
    ungroup() %>% 
    select(PATIENT_NUM, id_participant, date ,score)
  
  
  cohort_steps <- data.frame(list(Step = c("(0a) Total Patients with any PDDS score"
                                           , "(0b) Total Patients/Obs with PATIENT_NUM identifier"
                                           , "(1) PDDS Scores Observed 1+ Months Apart")
                                  , N_pts = c(nrow(unique(pdds_all[,c("PATIENT_NUM", "id_participant")]))
                                              , nrow(unique(pdds_all[!is.na(pdds_all$PATIENT_NUM),c("PATIENT_NUM", "id_participant")]))
                                              , nrow(unique(cohort1[,c("PATIENT_NUM", "id_participant")]))
                                              )
                                  , N_obs = c(nrow(pdds_all)
                                              , nrow(pdds_all[!is.na(pdds_all$PATIENT_NUM),])
                                              , nrow(cohort1))
                                  )
                             )
  
    
  # random checker   
    # .rand_id <- sample(pdds_all$PATIENT_NUM, 1)
    # pdds_all %>% filter(PATIENT_NUM==.rand_id) %>% arrange(date) %>% print(n=300)
    # cohort1 %>% filter(PATIENT_NUM==.rand_id) %>% arrange(date) %>% print(n=300)

  
## 2 - Demographics Covariates ####
    # Some records excluded for missing a subset of demographic variables
    # Oddly some variables do not have any record in our registry demographics file 
    # Registry DOB "imputed" (day not provided, year/month only) (see Notes.MD, 8-21 to 8-25 section, for detail bullet)
  
    source(here("Code", "demographic_covariates.R"))
    
    cohort_demo_inclNA <- cohort1 %>% 
      left_join(demo_covars, by=c("PATIENT_NUM", "id_participant")) %>% 
      mutate(Age_at_PDDS = (date - dob) / dyears(1) ) 
    
    cohort_demo <- cohort_demo_inclNA %>% 
                    filter(if_all(.cols = c(white_race:subject_sex_ch)
                                  , .fns = ~ !is.na(.)
                    ))

    
    .demo_missing <- cohort_demo_inclNA %>% 
      filter(if_any(.cols = c(Age_at_PDDS, White_NonHispanic, subject_sex_ch) 
                    , is.na
      )
      ) %>% 
      distinct(PATIENT_NUM, Age_at_PDDS, White_NonHispanic, subject_sex_ch) %>% 
      mutate(across(.cols = c(Age_at_PDDS, White_NonHispanic, subject_sex_ch)
                    , is.na
      )) 
    
      write_csv(x = count(x = .demo_missing, Age_at_PDDS, White_NonHispanic, subject_sex_ch)
                , file = here("Data", "created_data", paste0(gsub(pattern = ".csv", "", .output_file)
                                                                       , "_demo_missing.csv"))
                )
        
    
    
    cohort_steps <- rbind(cohort_steps
                          , data.frame(
                              Step = "(2) Missing any demographic covariates"
                              , N_pts = nrow(unique(cohort_demo[,c("PATIENT_NUM", "id_participant")]))
                              , N_obs = nrow(cohort_demo)
                            )
                          )    
    
    
  ## 3 - Clinical Covariates ####
      # Preforms many:many matches which produces (expected) warnings
  .clin_covar_input <- cohort_demo
  source(here("Code", "clinical_covariates.R"))
  
  cohort_clinical_inclNA <- cohort_demo %>% 
    left_join(clinical_covars, by=c("PATIENT_NUM", "id_participant", "date")) %>% 
    mutate(N_MS_PheCodes = coalesce(N_MS_PheCodes, 0)
           , N_AllCodes_PrePDDS = coalesce(N_AllCodes_PrePDDS, 0)
    ) %>% 
    group_by(PATIENT_NUM, id_participant, date) %>% 
    mutate(Days_FirstPheCode_to_PDDS = pmax(date - First_PheCode, ddays(0), na.rm = F)
           , Days_First_MS_PheCode_to_PDDS = pmax(date - First_MS_PheCode, ddays(0), na.rm = F)
    ) %>% 
    ungroup()
  
  cohort_clinical_inclNA %>% 
    filter(if_any(.cols = c(N_AllCodes_PrePDDS
                            # , Disease_Subtype
                            , Days_FirstPheCode_to_PDDS
                            , Days_First_MS_PheCode_to_PDDS
                            , N_MS_PheCodes) 
                  , .fns = function(x) is.na(x) | is.infinite(x)
    )
    ) %>% 
    distinct(PATIENT_NUM
             , id_participant
             , N_AllCodes_PrePDDS
             , Disease_Subtype
             , Days_FirstPheCode_to_PDDS
             , Days_First_MS_PheCode_to_PDDS
             , N_MS_PheCodes) %>% 
    mutate(across(.cols = c(Disease_Subtype
                            , N_AllCodes_PrePDDS
                            , Days_FirstPheCode_to_PDDS
                            , Days_First_MS_PheCode_to_PDDS
                            , N_MS_PheCodes)
                  , is.na
    )) %>% 
    group_by(Disease_Subtype
          , N_AllCodes_PrePDDS
          , Days_FirstPheCode_to_PDDS
          , Days_First_MS_PheCode_to_PDDS
          , N_MS_PheCodes) %>%
    summarize(N_obs = n()
           , N_pts = n_distinct(id_participant, PATIENT_NUM)) %>% 
    write_csv(file = here("Data", "created_data", paste0(gsub(pattern = ".csv", "", .output_file)
                                                         , "_clinical_missing.csv"))
    )
    
  
  
  cohort_clinical <- cohort_clinical_inclNA %>% 
    tidyr::drop_na(N_AllCodes_PrePDDS
                   # , Disease_Subtype
                   , Days_FirstPheCode_to_PDDS
                   , Days_First_MS_PheCode_to_PDDS
                   , N_MS_PheCodes
                   ) %>% 
    filter(!if_any(.cols = c(N_AllCodes_PrePDDS
                            , Days_FirstPheCode_to_PDDS
                            , Days_First_MS_PheCode_to_PDDS
                            , N_MS_PheCodes) 
                  , is.infinite
    )
    ) 
    
  ### Cohort Steps #### 
  
  cohort_steps <- rbind(cohort_steps
                        , data.frame(Step = "(3) Complete Clinical Covariates (excl. Disease Subtype)"
                                     , N_pts = nrow(unique(cohort_clinical[,c("PATIENT_NUM", "id_participant")])) 
                                     , N_obs = nrow(cohort_clinical)
                                     )
                        , data.frame(Step = "(3b) Complete Clinical Covariates AND Non-Missing Disease Subtype"
                                     , N_pts = nrow(unique(cohort_clinical[!is.na(cohort_clinical$Disease_Subtype),c("PATIENT_NUM", "id_participant")])) 
                                     , N_obs = nrow(cohort_clinical[!is.na(cohort_clinical$Disease_Subtype),])
                                     )
  )
  
  
  
  
  
# Export #### 
  cohort_output <- cohort_clinical %>% 
    select(PATIENT_NUM
           , id_participant
           , date
           , score
           , White_NonHispanic
           , subject_sex_ch
           , Age_at_PDDS
           , N_AllCodes_PrePDDS
           , Disease_Subtype
           , N_MS_PheCodes 
           , Days_FirstPheCode_to_PDDS
           , Days_First_MS_PheCode_to_PDDS) 
  
    cohort_output %>% write_csv(here("Data", "created_data"
                                     , paste0(.output_file, ".csv")
                                     )
                                )

    cohort_steps %>% write_csv(here("Data", "created_data"
                                    , paste0(.output_file, "_Derivation", ".csv")
                                    )
                               )
    
  message("Output dataframe from cohort construction created (object titled cohort_output)")
  message(paste0("Data file exported to "
                 , here("Data", "created_data", paste0(.output_file, ".csv"))
                 )
        )
  