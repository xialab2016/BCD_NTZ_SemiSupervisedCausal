
# Demographic Covariates 

# Data Import ####
  if("reg_demo" %nin% ls()) {
    reg_demo <- read_csv(here("Data", "Latest", "ClinicalDemographics_230425.csv")) %>% 
      select(-`...1`) 
  }

# EHR Demographics ####

  ehr_demo <- read_delim(here("Data", "Latest", "ms_demographics.dsv"), delim = "\t"
                         , guess_max = 17000) %>% 
    mutate(PATIENT_NUM = as.character(PATIENT_STUDY_ID)
           , female_gender = 
             as.factor(
               case_when(GENDER_CODE %in% c("1", "F") ~ "F"
                       , GENDER_CODE %in% c("2", "M") ~ "M"
                       , T ~ NA_character_)
               )
           , white_race = 
             case_when(RACE_CODE == 26 ~ 1
                       , !is.na(RACE_CODE) ~ 0)
           , nonhispanic_ethnicity  = 
             case_when(
               ETHNIC_CODE == 3 ~ 0
               , !is.na(ETHNIC_CODE) ~ 1
             )
           , dob_ehr = lubridate::dmy(BIRTH_DATE)
    ) %>% 
    rowwise() %>% 
    mutate(White_NonHispanic = white_race*nonhispanic_ethnicity) %>% 
    select(PATIENT_NUM, dob_ehr, white_race, nonhispanic_ethnicity, White_NonHispanic, subject_sex_ch = female_gender) %>% 
    left_join(id_link, "PATIENT_NUM")
    # select(PATIENT_NUM, dob_ehr)

# Registry Demographics #### 
  reg_demo_covars <- reg_demo %>% 
    mutate(white_race = 
             case_when(
               race == 4 ~ 1
               , !is.na(race) ~ 0
             )
           , nonhispanic_ethnicity = 
             case_when(
               ethnicity == 2 ~ 1
               , !is.na(ethnicity) ~ 0
             )
           , subject_sex_ch = 
             case_when(
               subject_sex == 1 ~ "M"
               , subject_sex == 2 ~ "F"
             )
           , dob_reg = lubridate::my(dob)
           ) %>% 
    rowwise() %>% 
    mutate(White_NonHispanic = white_race*nonhispanic_ethnicity) %>% 
    select(id_participant, dob_reg, white_race, nonhispanic_ethnicity, White_NonHispanic, subject_sex_ch) %>% 
    left_join(id_link, "id_participant")
    
  # lubridate passes some dates to 2000's that are 1900's (e.g. Jan-50 under lubridate::my becomes January 1st, 2050). Must clean manually   
  year(reg_demo_covars$dob_reg[year(reg_demo_covars$dob_reg)>2030 & !is.na(reg_demo_covars$dob_reg)]) = year(reg_demo_covars$dob_reg[year(reg_demo_covars$dob_reg)>2030 & !is.na(reg_demo_covars$dob_reg)]) - 100
    
# Merging ####
  # Updated to prioritize registry in join (but using both sources)
  demo_covars <- ehr_demo %>% 
    # left_join(id_link, "PATIENT_NUM") %>% 
    full_join(reg_demo_covars, c("PATIENT_NUM", "id_participant")) %>% 
    mutate(
      white_race = coalesce(white_race.y, white_race.x)
      , nonhispanic_ethnicity = coalesce(nonhispanic_ethnicity.y, nonhispanic_ethnicity.x)
      , White_NonHispanic = coalesce(White_NonHispanic.y, White_NonHispanic.x)
      , subject_sex_ch = coalesce(subject_sex_ch.y, subject_sex_ch.x)
      , dob = coalesce(dob_ehr, dob_reg)
      ) %>% 
    select(id_participant, PATIENT_NUM, white_race, nonhispanic_ethnicity, White_NonHispanic, subject_sex_ch, dob) %>% 
    distinct()
  
  
# Conclusion #### 
  gc()
  
  print("Script demographic_covariates.R complete")
  print("Output data.frame: demo_covars")
  print("Contains EHR, Registry, and combined demographic covariates. Registry covariates prioritized over EHR")
  cat("\n\n")
  
