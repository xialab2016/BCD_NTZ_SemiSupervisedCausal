## Clinical Covariates

  # Note: Prior DMT Use Duration is created in cohort_building.R under the clinical covariates section
  # This requires the dmt_cleaning program and resulting variables, which are already merged into 
  # the cohort at this point in the construction.



# Preface #### 
  print("Running 03_clinical_covariates.R ....")

  .ls_init <- ls(all.names = T)

  if(".clin_covar_input" %nin% ls(all.names = T)){
    stop("No clinical covariate input dataframe specified. Please assign a dataframe to the .clin_covar_input object to procede. This dataframe must contain both identifiers (PATIENT_NUM, id_participant) and DMT_Study_Start variable in order to assign some covariates")
  }else if(any(c("DMT_Study_Start", "PATIENT_NUM", "id_participant") %nin% colnames(.clin_covar_input))){
    stop("Dataframe specified for .clin_covar_input does NOT contain one of the following required variables: PATIENT_NUM, id_participant, DMT_Study_Start")
  }
  
# Set-Up #### 
  
  
  
  
if("ehr_2004_2022" %in% ls()){
  ehr <- ehr_2004_2022
  ehr$patient_num <- as.character(ehr$patient_num)
}else if ("ehr" %nin% ls(all.names=T)){
  ehr <- read_csv(here("Data", "Latest"
                       , ehr_codified_rolled )
                  , col_names = T
  ) %>% 
    mutate(PATIENT_NUM = as.character(patient_num)) # step 0/1
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
  
  ehr_2004_2022 = ehr[ehr$feature_id %in% feat.set,]
  # length(unique(ehr$feature_id[substr(ehr$feature_id, 1, 4)=="PheC"]))
  # length(unique(ehr.analytic$feature_id[substr(ehr.analytic$feature_id, 1, 4)=="PheC"]))
  # removes ~1/3 of PheCodes, seems reasonable 
  
  

if("registry_dmt" %nin% ls()) {
  registry_dmt <- read_csv(here("Data", "Latest", "DMTHistory.csv") 
                           , col_names=T) %>% select(-`...1`)
}

if("reg_demo" %nin% ls()) {
  reg_demo <- read_csv(here("Data", "Latest", "ClinicalDemographics_230425.csv")) %>% 
    select(-`...1`)# step 2 
}
  
  if("node" %nin% ls()){
    node <- read_csv(here("Data", "Latest", "node.csv"))
  }
  
  # Disease Subtype ####
  # groupings based on enrollment_diagnosis in registry demo dat (defined in Weijing's Notes)
  # (1-3) - Cat1 / (4-5) - Cat2 / (6, 9-10) - Cat3
  .clin_covar.subtype <- reg_demo %>% 
    mutate( Disease_Subtype = 
              case_when(
                enrollment_diagnosis %in% c(1, 2, 3) ~ 1
                , enrollment_diagnosis %in% c(4, 5) ~ 2 
                , enrollment_diagnosis %in% c(6, 9, 10) ~ 3
              )
    ) %>% 
    distinct(id_participant, Disease_Subtype)
  

  
  # Disease Duration(s) ####
  # Time between first MS PheCode and treatment initiation
  # Time Elapsed from first PheCode (any) to treatment initiation
  
  # I'm only pulling the dates of these first PheCodes, the true calculation is then complete
  # in cohort_building.R 
  
  # .clintime_RDS_exists <- file.exists(here("Data", "created_data", "_intermediate_RDS", ".clin_covar.time_elapsed.RDS"))
  # 
  # if(.clintime_RDS_exists){
  #   .clin_covar.time_elapsed <- readRDS(here("Data", "created_data", "_intermediate_RDS", ".clin_covar.time_elapsed.RDS"))
  # }else{
  
  .clin_covar.time_elapsed <- ehr_2004_2022 %>% 
    # filter(substr(feature_id, 1, 7)=="PheCode") %>% 
    filter(str_detect(feature_id, "PheCode")) %>%
    group_by(PATIENT_NUM) %>% 
    mutate(First_PheCode = min(start_date)
           , First_MS_PheCode = 
             min(
               case_when(
                 feature_id == "PheCode:335" ~ start_date
               )
               , na.rm=T
             )
    ) %>% 
    distinct(PATIENT_NUM, First_PheCode, First_MS_PheCode)
  
  .clin_covar.count_all_codes <- .clin_covar_input %>% 
    filter(start_date < DMT_Study_Start) %>% 
    distinct(PATIENT_NUM, id_participant, start_date, feature_id
    ) %>% 
    filter(str_detect(feature_id, "PheCode|LOINC|RXNORM|CCS")
           & substr(feature_id, 1, 5)!="LOCAL") %>%
  count(PATIENT_NUM, id_participant)
  
  # Count of MS PheCodes (before DMT initiation) ####
  .clin_covar.count_MS_phecodes <- .clin_covar_input %>% 
    filter(start_date <= DMT_Study_Start & feature_id=="PheCode:335") %>% 
    distinct(PATIENT_NUM, id_participant, start_date, feature_id
             # , corvalue
    ) %>% 
    count(PATIENT_NUM, id_participant, name = "N_MS_PheCodes")   
  

  # Merging Output File ####
  clinical_covars <- .clin_covar.subtype %>% 
    full_join(.clin_covar.count_all_codes, by="id_participant") %>% 
    full_join(.clin_covar.time_elapsed, by="PATIENT_NUM") %>%  
    full_join(.clin_covar.count_MS_phecodes, by=c("PATIENT_NUM", "id_participant")) %>% 
    select(PATIENT_NUM, id_participant, N_AllCodes_PreStudyDMT = n
             , Disease_Subtype, First_PheCode, First_MS_PheCode
             , N_MS_PheCodes) %>% 
    distinct()
  
  # Conclusion ####
  .ls_new <- setdiff(ls(all.names = T), c(.ls_init, "clinical_covars"))
  rm(.ls_new)
  
  gc()
  
  print("Script 03_clinical_covariates.R complete")
  print("Output data.frame: clinical_covars")
  cat("\n\n")