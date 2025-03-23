
# Data Processing (EHR & CLIMB)
## Dominic DiSanto
## Updated: 7.4.2023 


# Preface
  # Print 
    print("Running cohort_bulding.R script (including sub-called programs)...")
    if(.high_eff){
          print("High-efficacy drug specified by `.high_eff` object. . This affects inclusion/exclusion criteria related to pre-study DMT and DMT switching")
      } else{
          print("Standard-efficacy drug specified by `.high_eff` object. This affects inclusion/exclusion criteria related to pre-study DMT and DMT switching")
      }

# Data Import #### 

  if("ehr_2004_2022" %nin% ls()) { # bit lazy, don't re-import if I already have from upstream script
    ehr_2004_2022 <- read_csv(here("Data", "Latest"
                                   , ehr_codified_rolled )
                              , col_names = T
    ) %>% 
      mutate(PATIENT_NUM = as.character(patient_num)) # step 0/1
  }
  
  registry_dmt <- read_csv(here("Data", "Latest", "DMTHistory.csv")
                           , col_names=T) %>% select(-`...1`) # step 0/1 
  
  id_link <- read_excel(here("Data", "ID Linkage"
                                     , id_link_filename)) %>% select(-LINE_NUM) # step 0/1
  
  reg_demo <- read_csv(here("Data", "Latest", "ClinicalDemographics_230425.csv")) %>% 
  # reg_demo <- read_csv(here("Data", "Latest", "ClinicalDemographics.csv")) %>% 
    select(-`...1`)# step 2 
  
  
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
    )
  
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
    mutate(score = suppressWarnings(as.double(score))
           , date = suppressWarnings(as.Date(as.integer(date), origin="1899/12/30"))
    ) %>% 
    filter(!is.na(date))
  
  pdds_ids <- rbind(.pdds
                , pdds_MayPull
                , pdds_pull1 
                , pdds_pull2
                , pdds_pull4
  ) %>% distinct(PATIENT_NUM, id_participant) 


# Cohort Derivation ####
  
  ## 1 -  Participants must be in EHR data *and* registry ####
  # UPMC provided file with identifiers to be removed/kept, as some 
    # PATIENT_NUM values had multiple id_participant ID's (so multiple registry ID's)
    # per UPMC team, these were individuals consented twice and errant assigned a second ID
    # UPMC team provided files with guidance on retaining/removing ID's
  # 
  
  dup_rm <- read_excel(here("Data", "ID Linkage", 
                            "Duplicate_IDs_Reconciled.xlsx")
                       , sheet = "Remove") %>% 
    mutate(RemoveTag=1
           , PATIENT_NUM=as.character(PATIENT_NUM)) %>% 
    select(-c(n_id_part, n_pt_num))
  
  id_link_nodups <- id_link %>% 
    union(pdds_ids) %>% 
    left_join(dup_rm, by=c("PATIENT_NUM", "id_participant")) %>%
    filter(RemoveTag!=1 | is.na(RemoveTag)) %>% 
    filter(!is.na(PATIENT_NUM)) %>% # 6 id_participant values with no PATIENT_NUM
    # select(-RemoveTag) %>% 
    distinct()
  
  ms_notes = read_csv(here("Data", "created_data", "UPMC_MS_note3_patient_status_90_specificity.csv")) %>%
    select(patient_num, pred_y)
  
  tmp = read_excel(here("Data", "Gold Standard Labels", "PRT Dx 2024-02-09.xlsx"))
  
  # combining phenotyping algorithm with gold standard labels 
  ms.keep = c(ms_notes[ms_notes$pred_y==1,][["patient_num"]]
              , tmp[tmp$subject_group=="MS",][["PATIENT_NUM"]]
              )
  
  # Removing non-MS ID's (predicted MS but not by gold standard label, trust gold standard over algorithm)
  ms.keep = ms.keep[ms.keep %nin% tmp[tmp$subject_group!="MS",][["PATIENT_NUM"]]]
  
  cohort1 <- ehr_2004_2022 %>% 
    left_join(id_link_nodups, by = "PATIENT_NUM") %>% 
    filter(RemoveTag != 1 | is.na(RemoveTag)) %>% 
    select(-RemoveTag) 
  
  cohort1 = cohort1[which(cohort1$patient_num %in% ms.keep),]
  
  cohort_steps <- data.frame(list(Step = c("(0a) Total Patients in EHR"
                                           , "(0b) MS Patients in EHR"
                                           )
                                  , N=c(n_distinct(ehr_2004_2022$PATIENT_NUM)
                                        , length(unique(cohort1$patient_num))
                                        )
                                  )
                             )

    
## 2 - Excluding if EHR features 7/8/12 present (non-MS patients included in registry) #### 
  # excluding if any record with these exclusion diagnoses occurs
  # only one participant has multiple enrollment diagnoses, both meet exclusion criteria (see appendix)
  cohort2 <- cohort1 %>% 
    left_join(select(.data = reg_demo, id_participant, enrollment_diagnosis)
              , by="id_participant")  %>% 
    group_by(PATIENT_NUM) %>% 
    mutate(ehr_exclude_ever = enrollment_diagnosis %in% c(7, 8, 11, 12)) %>% 
    filter(ehr_exclude_ever==0) %>% 
    select(-ehr_exclude_ever) %>% 
    ungroup() #%>% 
    # filter(!is.na(enrollment_diagnosis))
  
  cohort_steps <- rbind(cohort_steps
                        , data.frame(list(Step = "(2) Removing non-MS enrollment diagnoses via Registry enrollment_diagnosis"
                                          , N = n_distinct(cohort2$patient_num)
                        )
                        )
  )
  
  ## 3 - Received DMT of Interest #### 
  .dmt_IDs <- cohort2 %>% distinct(PATIENT_NUM, id_participant)
  source(here("Code", "01_dmt_cleaning.R"))  
  # currently outputs a warning for no non-missing arguments
  # this is related to Chemotherapy drugs in registry Rx data (no chemo drugs present)
  # currently not an error, the registry simply does not include any chemotherapy drugs 
  #  for exclusion criteria in step 5
  
  cohort3a <- cohort2 %>% 
    # filter(!is.na(enrollment_diagnosis)) %>% 
    # left_join(dmt_total, by="PATIENT_NUM") %>%
    left_join(dmt_total, by=c("PATIENT_NUM", "id_participant")) %>%
    filter(!is.na(EHR_DMT_Study_Start) | !is.na(Registry_DMT_Study_Start)) %>%
    group_by(PATIENT_NUM, id_participant) %>% 
    mutate(DMT_Study_Start = pmin(EHR_DMT_Study_Start, Registry_DMT_Study_Start, na.rm=T)) %>% 
    ungroup() %>% 
    mutate(Eval_End_Date = date(ymd(DMT_Study_Start)+.eval_pd))
  
  cohort3b <- cohort3a %>% 
    filter(DMT_Earliest_Start>='2014-01-01')
  
  cohort_steps <- rbind(cohort_steps
                        , data.frame(list(Step = paste0("(3a) 1+ Record of DMT (", dmt_text, ") on or after 1-1-2014")
                                          , N = n_distinct(cohort3a$patient_num)
                        )
                        )
                        , data.frame(list(Step = paste("(3b) Removing records with pre-study record of", dmt_text ,"1-1-2014")
                                          , N = n_distinct(cohort3b$patient_num)
                        )
                        )
  )
  
  
  
  ## 4 - Receiving any high-efficacy drug prior to or during evaluation period (1-1-2014) ####  
  # other than drug of interest (right?) 
  # Necessary data are generated from 01_dmt_cleaning.R (which was run in step 3)
  
  # source(here("Code", "01_dmt_cleaning.R")) # exceuted previously, in step 3
  
  if(.high_eff){
    cohort4 <- cohort3b %>% 
      left_join(dmt_exclude, by=c("PATIENT_NUM", "id_participant")) %>% 
      filter(is.na(DMT_HighEff_Recent) |  # equivalent to no other high efficacy DMT given 
               DMT_HighEff_Earliest>=DMT_Study_Start) # or any DMT was given after first BCD 
    # switching taken care of in step 5
    
    
    cohort_steps <- rbind(cohort_steps
                          , data.frame(Step = paste0("(4) No Prior High Efficacy Treatment (other than ", dmt_text, ")")
                                       , N = n_distinct(cohort4$PATIENT_NUM)
                          )
    )
    
  } else{ # standard efficacy
    
    cohort4 <- cohort3b %>% 
      left_join(dmt_exclude, by=c("PATIENT_NUM", "id_participant")) %>% 
      filter(
        (is.na(DMT_HighEff_Recent) & is.na(DMT_StdEff_Recent)) |  # equivalent to no other high efficacy DMT given
          (DMT_StdEff_Earliest>=DMT_Study_Start & DMT_HighEff_Earliest>=DMT_Study_Start)
      )  
    # switching taken care of in step 5
    
    
    cohort_steps <- rbind(cohort_steps
                          , data.frame(Step = paste0("(4) No Prior High or Standard Efficacy Treatment (other than ", dmt_text, ")")
                                       , N = n_distinct(cohort4$PATIENT_NUM)
                          )
    )
    
  }
  
  
  ## 5 - No DMT switching #### 
  # (i.e. no treatment after DMT initiation for specific DMT of interest, high or standard efficacy )
  # with exception OF our DMT of interest
  # e.g. if we care about NTZ, a patient can have multiple records of NTZ prescriptions during the 
  # evaluation period and not be excluded, as long as no OTHER drug was prescribed during the eval period 
  
  cohort5 <- cohort4 %>% 
    filter(
      (is.na(DMT_StdEff_Recent) | DMT_StdEff_Earliest>Eval_End_Date | DMT_StdEff_Recent<DMT_Study_Start
      )
      &
        (is.na(DMT_HighEff_Recent) | DMT_HighEff_Earliest>Eval_End_Date | DMT_HighEff_Recent<DMT_Study_Start)
    )
  # a LOT of folks receive a standard efficacy drug within 3 years 
  # reduces cohort to 145 by including lines 152-3
  # if we only consider high efficacy 
  
  cohort_steps <- rbind(cohort_steps
                        , data.frame(Step = "(5) No switching to other DMT (high-efficacy) during study period (user specified time period in `.eval_pd`)"
                                     , N = n_distinct(cohort5$PATIENT_NUM)
                        )
  )
  
  
  ## 6 - No chemotherapy drugs before outcome assessment (including total lookback period, before DMT initiation/study start date) ####
  
  cohort6 <- cohort5 %>% 
    filter(DMT_Chemo_Earliest > Eval_End_Date | # earliest record of chemo happened after our evaluation period
             is.na(DMT_Chemo_Earliest)) # or no chemo at all
  
  
  cohort_steps <- rbind(cohort_steps
                        , data.frame(Step = paste0("(6) No chemotherapy drugs anytime pre-study or during eval period, ", .eval_pd, " after DMT")
                                     , N = n_distinct(cohort6$PATIENT_NUM)
                        )
  )
  
  
  
  
  # 7/8 - Covariates #### 
  # Not completely filtering for missing covariates yet, but including the code in order to count 
  
  ## 7 - Demographic ####
  source(here("Code", "02_demographic_covariates.R"))
  
  cohort_demo_inclNA <- cohort6 %>% 
    left_join(demo_covars, by=c("PATIENT_NUM", "id_participant")) %>% 
    mutate(Age_at_DMT_Initiation = (DMT_Study_Start - dob_ehr) / dyears(1) ) 
  
  .demo_missing <- cohort_demo_inclNA %>% 
    filter(if_any(.cols = c(Age_at_DMT_Initiation, White_NonHispanic, subject_sex_ch) 
                  , is.na
    )
    ) %>% 
    distinct(PATIENT_NUM, Age_at_DMT_Initiation, White_NonHispanic, subject_sex_ch) %>% 
    mutate(across(.cols = c(Age_at_DMT_Initiation, White_NonHispanic, subject_sex_ch)
                  , is.na
    )) 
    
  demo_wb <- openxlsx::createWorkbook()
  
    openxlsx::addWorksheet(demo_wb, "Count")
    openxlsx::addWorksheet(demo_wb, "Data")
    
    openxlsx::writeData(wb = demo_wb, sheet = "Count"
                        , x = count(x = .demo_missing, Age_at_DMT_Initiation, White_NonHispanic, subject_sex_ch)
    )
    
    openxlsx::writeData(wb = demo_wb, sheet = "Data"
                        , x = .demo_missing
    )    
    
  openxlsx::saveWorkbook(demo_wb
                         ,file = here("Data", "created_data", dmt_text, paste0(.export_filename, "_demo_missing.xlsx"))
                         , overwrite = T)
    
    
    

  
  cohort_demo <- cohort_demo_inclNA %>% 
    tidyr::drop_na(Age_at_DMT_Initiation, White_NonHispanic, subject_sex_ch)
  
  cohort_steps <- rbind(cohort_steps
                        , data.frame(Step = "(7) Complete Demographic Covariates"
                                     , N = n_distinct(cohort_demo$PATIENT_NUM)
                        )
  )
  
  ## 8 - Clinical #### 
  .clin_covar_input <- cohort_demo
  source(here("Code", "03_clinical_covariates.R"))
  
  cohort_clinical_inclNA <- cohort_demo %>% 
    left_join(clinical_covars, by=c("PATIENT_NUM", "id_participant")) %>% 
    mutate(N_MS_PheCodes = coalesce(N_MS_PheCodes, 0)
           , N_AllCodes_PreStudyDMT = coalesce(N_AllCodes_PreStudyDMT, 0)
    ) %>% 
    group_by(PATIENT_NUM, id_participant) %>% 
    mutate(earliest_dmt = min(DMT_StdEff_Earliest, DMT_Earliest_Start, na.rm=T)
           , Days_Earliest_DMT_to_Study_DMT = coalesce(DMT_Study_Start - earliest_dmt, ddays(0) )
           # , recent_dmt = max(DMT_StdEff_Recent, DMT_Recent_Prestudy, na.rm=T)
           # , Days_Recent_DMT_to_Study_DMT = coalesce(DMT_Study_Start - recent_dmt, ddays(0) )
           , Days_FirstPheCode_to_Study_DMT = coalesce(DMT_Study_Start - First_PheCode, ddays(0) )
           , Days_First_MS_PheCode_to_Study_DMT = coalesce(DMT_Study_Start - First_MS_PheCode, ddays(0) )
    ) %>% 
    ungroup()
  
  cohort_clinical_inclNA %>% 
    filter(if_any(.cols = c(N_AllCodes_PreStudyDMT
                            # , Disease_Subtype
                            , Days_FirstPheCode_to_Study_DMT
                            , Days_First_MS_PheCode_to_Study_DMT
                            , N_MS_PheCodes) 
                  , is.na
    )
    ) %>% 
    distinct(PATIENT_NUM
             , N_AllCodes_PreStudyDMT
             , Disease_Subtype
             , Days_FirstPheCode_to_Study_DMT
             , Days_First_MS_PheCode_to_Study_DMT
             , N_MS_PheCodes) %>% 
    mutate(across(.cols = c(Disease_Subtype
                            , N_AllCodes_PreStudyDMT
                            , Days_FirstPheCode_to_Study_DMT
                            , Days_First_MS_PheCode_to_Study_DMT
                            , N_MS_PheCodes)
                  , is.na
    )) %>% 
    count(Disease_Subtype
          , N_AllCodes_PreStudyDMT
          , Days_FirstPheCode_to_Study_DMT
          , Days_First_MS_PheCode_to_Study_DMT
          , N_MS_PheCodes) %>%
    write_csv(file = here("Data", "created_data", dmt_text, paste0(.export_filename, "_clinicalcovar_missing.csv")))
  
  
  cohort_clinical <- cohort_clinical_inclNA %>% 
    tidyr::drop_na(N_AllCodes_PreStudyDMT
                   # , Disease_Subtype
                   , Days_FirstPheCode_to_Study_DMT
                   , Days_First_MS_PheCode_to_Study_DMT
                   , N_MS_PheCodes) %>% 
    filter(N_AllCodes_PreStudyDMT>0)
  
  ### Cohort Steps #### 
  
  cohort_steps <- rbind(cohort_steps
                        , data.frame(Step = "(8) Complete Clinical Covariates (and non-zero utilization)"
                                     , N = n_distinct(cohort_clinical$PATIENT_NUM) 
                        )
  )
  
  
  
  
  ## 9 - Missing/Unavailable PDDS Outcome/Values #### 
  .pdds_input <- cohort_clinical
  source(here("Code", "04_pdds.R"))
  
  cohort_outcome <- cohort_clinical %>% 
    left_join(pdds_outcome, by="PATIENT_NUM") %>% 
    mutate(PDDS_NegChange = case_when(PDDS_Change == "Improvement/Decrease" ~ 1
                                                  , PDDS_Change %in% c("No sustainment", "Progression/Increase") ~ 0
                                                  , T ~ NA_integer_
                                      )
           , PDDS_Change = case_when(PDDS_Change == "Progression/Increase" ~ 1
                                    , PDDS_Change %in% c("No sustainment", "Improvement/Decrease") ~ 0
                                    , T ~ NA_integer_
                                    )
           )
  

  cohort_steps <- rbind(cohort_steps
                        , data.frame(Step = c("(9a) Valid PDDS Difference Outcome"
                                              , "(9b) Valid PDDS Sustained Change Outcome (3+ obs)"
                        )
                        , N = c(n_distinct(cohort_outcome[!is.na(cohort_outcome$PDDS_Diff), "PATIENT_NUM"])
                                , n_distinct(cohort_outcome[!is.na(cohort_outcome$PDDS_Change), "PATIENT_NUM"])
                        )
                        )
  )
  
  
  
  
  
  
  # Export #### 
  ## Some cohort cleaning ####
  cohort_output <- cohort_outcome %>% 
    # select(-c(start_date, feature_id, RemoveTag)) %>% 
    # distinct() %>% 
    select(PATIENT_NUM
           , id_participant
           , DMT_Study_Start # date of DMT start/initiation
           , Eval_End_Date # date of end of study look-up period 
           , Age_at_DMT_Initiation
           , White_NonHispanic
           , white_race # used to create White_NonHispanic 
           , nonhispanic_ethnicity # used to create White_NonHispanic
           , subject_sex_ch 
           , Disease_Subtype
           , enrollment_diagnosis # used to create Disease_Subtype
           , Disease_Duration = Days_First_MS_PheCode_to_Study_DMT # Disease_Duration
           , N_MS_PheCodes 
           , Followup_Duration = Days_FirstPheCode_to_Study_DMT # Follow-up Duration
           , Utilization = N_AllCodes_PreStudyDMT # healthcare utilization metric
           , Days_Earliest_DMT_to_Study_DMT # prior DMT use 
           # , PriorPDDS_Merged
           # , PDDS_Imputed_flag
           # , # elixhauser
           # , # charlson
           , DMT 
           , Average_PDDS
           , N_PDDS 
           , PDDS_Diff
           , PDDS_Diff_Date 
           , PDDS_Change 
           , PDDS_NegChange
    ) %>% distinct() %>% 
    mutate(Disease_Duration = 
             case_when(Disease_Duration>=0 ~ time_length(Disease_Duration, unit="days")
                       , T ~ time_length(0, "days")
             )
           , Followup_Duration = time_length(Followup_Duration, unit="days")
           , Days_Earliest_DMT_to_Study_DMT = time_length(Days_Earliest_DMT_to_Study_DMT, unit="days")
    ) %>% 
    mutate(Disease_Duration_Yrs = duration(Disease_Duration, "days") / dyears(1)
           , Followup_Duration_Yrs = duration(Followup_Duration, "days") / dyears(1)
           , Days_Earliest_DMT_to_Study_DMT_Yrs = duration(Days_Earliest_DMT_to_Study_DMT, "days") / dyears(1)
    )
  
  ## Data ####
  .skip_build = ifelse(".skip_build" %in% ls(all.names=T), .skip_build, F)
  
  if(.skip_build){
    warning("Not exporting final data set, .skip_build=T found")
    cohort_output %>% write_csv(here("Data", "created_data", dmt_text
                                     , paste0(.export_filename, ".csv")
    )
    )
  }
  
  # Consort/Derivation Numbers 
  if(.export_cohort_steps){
    cohort_steps %>% write_csv(here("Data", "created_data", dmt_text
                                    , paste0(.export_filename, "_Derivation", ".csv")
    )
    )
  }
  
  assign(.final_cohort_name, cohort_output, envir = .GlobalEnv) 
  
  print(paste0("Output dataframe from cohort construction:", .final_cohort_name))
  