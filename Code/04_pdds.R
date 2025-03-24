


# Preface ####  
print("Running 04_pdds.R ....")

## Notes ####

if(".pdds_input" %nin% ls(all.names = T)){
  stop("No input dataframe specified. Please assign a dataframe to the .pdds_input object to procede. This dataframe must contain both identifiers (PATIENT_NUM, id_participant) and DMT_Study_Start variable in order to assign some covariates")
}


# Post-DMT PDDS as Outcome #### 
if(".unittest" %nin% ls(all.names = T)){
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
           , date = as.Date(as.integer(date), origin = "1899/12/30")
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
    mutate(score = as.double(score)
           , date = as.Date(as.integer(date), origin="1899/12/30")
           )
  
  pdds <- rbind(.pdds
                , pdds_MayPull
                , pdds_pull1 
                , pdds_pull2
                , pdds_pull4
                ) %>% distinct() 

  # pdds <- .pdds
}
  
  pdds_dmt <- .pdds_input %>% 
    group_by(PATIENT_NUM) %>% 
    # mutate(DMT_Start_Date = min(c(EHR_DMT_Earliest_Start, Registry_DMT_Earliest_Start), na.rm=T)) %>% 
    distinct(PATIENT_NUM, id_participant, DMT_Study_Start) %>% 
    right_join(pdds, by="PATIENT_NUM")
    # right_join(pdds, by="id_participant")
  
# subset to after and within three years of DMT initiation 
  .pdds_outcome <- pdds_dmt %>%
    filter(date >= DMT_Study_Start 
           & date <= (DMT_Study_Start + .eval_pd)
           )
  source(here("Code", "PDDS_SustainFunction.R"))
  
  ## Average over  window #### 
    pdds_avg <- .pdds_outcome %>% 
      group_by(PATIENT_NUM) %>%
    arrange(PATIENT_NUM, date) %>%
    mutate(score_diff = lead(score) - score
           , score_diff_dt = lead(date)) %>%
      summarize(Average_PDDS = mean(score)
                , N_PDDS = n()
                , First_PostDMT_PDDS = score[which.min(date)]
                , test = first(score)
                , First_PostDMT_PDDS_dt = min(date)
                # , Recent_PostDMT_PDDS = max(date)
                , Recent_PostDMT_PDDS_dt = max(date)
                , PDDS_Diff = first(score_diff)
                , PDDS_Diff_Date = first(score_diff_dt)
      ) %>%
    ungroup()

    # pdds_avg = pdds_average_fn(data = .pdds_outcome
    #                            , date = "date", score="score", id="PATIENT_NUM"
    #                            , .cleanup = T, .indx = T
    #                            , window = dmonths(.change_window)
    #                            )

  ## Change by 1+ (requiring two measurements) #### 
    # measurements must be separated by at user-specified time period (typically 6 months)
  
  pdds_change = pdds_sustain(data = .pdds_outcome
                             , date = "date", score="score", id="PATIENT_NUM"
                             , .cleanup = T, .indx = T
                             , window = dmonths(.change_window)
                             ) %>% 
    distinct(PATIENT_NUM, PDDS_Change=Sustainment)
  
  pdds_outcome <- full_join(pdds_avg, pdds_change, by="PATIENT_NUM") 
  
# Conclusion ####
  gc()
  
  print("Script 04_pdds_outcome.R complete")  
  print("Output data.frame: pdds_outcome")
  print("Includes both 'change' and 'average' PDDS outcome measures")
  cat("\n\n")
  
  