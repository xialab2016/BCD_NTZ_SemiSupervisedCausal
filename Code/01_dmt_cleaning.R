
# Preface #### 
  print("Running 01_dmt_cleaning.R ....")

## Notes ####
# note this code does not include dmt_switching, which requires PDDS outcome date (which in turn requires DMT initiation)
  # start with DMT initiation (this program)
  # then move to pdds_outcome cleaning (02_pdds_outcome)
  # then exclude for switching in subsequent programs 



# Set-Up ####
if("dmt_update" %nin% ls()){
  dmt_update = read_excel(here("Data", "Latest", "Additional MS DMT updated 20240229.xlsx")) %>% 
    mutate(patient_num = PATIENT_NUM
           , PATIENT_NUM = as.character(PATIENT_NUM)
           , start_date = lubridate::date(date)) %>%
    select(patient_num, start_date, feature_id, PATIENT_NUM)
}

if("ehr_2004_2022" %nin% ls()) { # large file, manually skipping if present in environment
  ehr_2004_2022 <- read_csv(here("Data", "Latest"
                                 , ehr_codified_rolled )
                            , col_names = T
  ) %>% 
    mutate(PATIENT_NUM = as.character(patient_num)) # step 0/1
}

if("registry_dmt" %nin% ls()) {
registry_dmt <- read_csv(here("Data", "Latest", "DMTHistory.csv") 
                         , col_names=T) %>% select(-`...1`)
}

.ls_init <- ls(all.names = T)

# (EHR) DMT Start Date #### 

ehr_rxnorm <- ehr_2004_2022 %>% 
  union_all(dmt_update) %>% 
  inner_join(.dmt_IDs, by="PATIENT_NUM") %>% 
  filter(str_detect(feature_id, "RXNORM")) %>% 
  mutate(RxNorm_ID = gsub("\\D*", "", feature_id)
         , treatment =
           toupper(
             case_when(RxNorm_ID == "1876366" ~ "BCD:Ocrelizumab"
               , RxNorm_ID == "121191" ~ "BCD:Rituximab"
               , RxNorm_ID == "712566" ~ "BCD:Ofatumumab" 
               , RxNorm_ID == "354770" ~ "Natalizumab" 
               , RxNorm_ID %in% c("84375", "214582") ~ "GA"
               , RxNorm_ID %in% c("1012892", "2121085", "2288236", "2532300") ~ "S1P Receptor"
               , RxNorm_ID %in% c("1373478", "1546433", "2261783") ~ "Fumarate"
               , RxNorm_ID %in% c("117055", "3002", "7005", "44157", "6851", "265323") ~ "Exclude"
               , T ~ "Other"
             )
           )
         , Source = "EHR"
         , end_date = NA 
         , High_Eff_Exclude = tolower(RxNorm_ID) %in% tolower(HighEff_Exclude_RxNormIDs)
         , Std_Eff_Exclude = tolower(RxNorm_ID) %in% tolower(StdEff_Exclude_RxNormIDs)
         , Chemo_Exclude = tolower(RxNorm_ID) %in% tolower(Chemo_Exclude_RxNormIDs)
         )


  # limit to drug of interest
    ehr_rxnorm_sub <- ehr_rxnorm %>% 
      # filter(treatment=="FUMARATE") 
    filter(RxNorm_ID %in% drugs$RxNormID)
    
  # earliest and most recent dates of drug of interest
    .dts <- ehr_rxnorm_sub %>% 
      group_by(PATIENT_NUM) %>% 
      mutate(DMT_Earliest_Date = min(start_date, na.rm=T)
             , DMT_Most_Recent_Date = max(start_date, na.rm=T)) %>% 
      distinct(PATIENT_NUM, id_participant, DMT_Earliest_Date, DMT_Most_Recent_Date)

  # earliest date within study period (after 1-1-2014)
      # (and merging earliest/latest dates created above in .dts)
    dmt_dts_ehr <- ehr_rxnorm_sub %>% 
      filter(start_date >= '2014-01-01') %>% 
      group_by(PATIENT_NUM) %>% 
      mutate(DMT_Earliest_Study_Date = min(start_date, na.rm=T)) %>% 
      distinct(PATIENT_NUM, id_participant, DMT_Earliest_Study_Date) %>% 
      full_join(.dts, by=c("PATIENT_NUM", "id_participant")) %>% 
      select(PATIENT_NUM
             , id_participant
             , EHR_DMT_Study_Start = DMT_Earliest_Study_Date
             , EHR_DMT_Earliest_Start = DMT_Earliest_Date
             , EHR_DMT_MostRecent = DMT_Most_Recent_Date
             ) %>%
      ungroup()

    
    # Jan 2025 edits to include sub-drugs in Table 1 #f    
    # dmt_dts_ehr = dmt_dts_ehr %>%
    #   left_join(ehr_rxnorm_sub %>% select(PATIENT_NUM, start_date, treatment) %>% distinct()
    #             , by = c("PATIENT_NUM", "EHR_DMT_Study_Start"="start_date")) 
    
    
    # dmt_dts_ehr %>% left_join(ehr_rxnorm_sub %>% select(PATIENT_NUM, start_date, treatment) %>% distinct()
    #                           , by = c("PATIENT_NUM", "EHR_DMT_Study_Start"="start_date")) %>% 
    #   arrange(PATIENT_NUM) %>% group_by(PATIENT_NUM) %>% filter(n()>1)

    # dmt_dts_ehr %>% left_join(ehr_rxnorm_sub %>% select(PATIENT_NUM, start_date, treatment) %>% distinct()
    #                           , by = c("PATIENT_NUM", "EHR_DMT_Study_Start"="start_date")) %>% dim()
    # dmt_dts_ehr %>% dim()
  
    # note the above dataset contains a small number of PATIENT_NUMs with multiple id_participants 
      # these are excluded at downstream steps, and are currently ignored 
      # as of 1.26.2023, we are unsure of how to deal with these observations so continuing to exclude
    
    
# (Registry) DMT Start Date #### 
    
    registry_dmt_cln <- registry_dmt %>% 
      filter(!is.na(start)) %>% 
      mutate(
        treatment = 
          toupper(
            case_when(
              type_brand == "Tecfidera" ~ "Fumarate"
              , type_brand %in% c("Aubagio", "Other") ~ "Other" 
              , type_brand %in% c("Copaxone") ~ "GA"
              , type_brand %in% c("Novantrone") ~ "CPB" 
              , type_brand == "Tysabri" ~ "Natalizumab" # why is this separate and no in "EXCLUDE"?
              , type_brand %in% c("Avonex", "Rebif", "Betaseron", "Plegridy", "Extavia")  ~ "Interferon-Beta"
              , type_brand == "Rituxan" ~ "BCD:Rituximab"
              , type_brand == "Ocrevus" ~ "BCD:Ocrelizumab"
              , type_brand == "Kesimpta" ~ "BCD:Kesimpta" # not present in data now but still valid DMT (generic: Ofatumab)
              , type_brand == "Gilenya" ~ "S1P Receptor"
            )
          )
        , High_Eff_Exclude = (tolower(type_brand) %in% tolower(HighEff_Exclude_BrandNames) ) | 
                             (tolower(type_generic) %in% tolower(HighEff_Exclude_GenericNames))
        
        , Std_Eff_Exclude = (tolower(type_brand) %in% tolower(StdEff_Exclude_BrandNames) ) | 
                             (tolower(type_generic) %in% tolower(StdEff_Exclude_GenericNames))
        
        , Chemo_Exclude = (tolower(type_brand) %in% tolower(Chemo_Exclude_BrandNames)) |
                           (tolower(type_generic) %in% tolower(Chemo_Exclude_GenericNames))
        , Source="Registry"
      )
    
    
    # limit to drug of interest
      registry_dmt_sub <- registry_dmt_cln %>%
        filter(tolower(type_generic) %in% tolower(drugs$Type_Generic) |
                 tolower(type_brand) %in% tolower(drugs$Type_BrandName)) 
        
    # earliest date of drug of interest
      .dts <- registry_dmt_sub %>% 
        group_by(id_participant) %>% 
        mutate(DMT_Earliest_Date = min(start, na.rm=T)
               , DMT_Most_Recent_Date = max(start, na.rm=T)) %>% 
        distinct(id_participant, DMT_Earliest_Date, DMT_Most_Recent_Date)
      
      # registry_dmt_sub %>% 
      #   group_by(id_participant) %>% 
      #   mutate(DMT_Earliest_Date = min(start, na.rm=T)
      #          , DMT_Most_Recent_Date = max(start, na.rm=T)) %>% 
      #   distinct(id_participant, DMT_Earliest_Date, DMT_Most_Recent_Date, type_generic) %>% 
      #   arrange(id_participant) %>% group_by(id_participant)
      

    # earliest date within study period (after 1-1-2014)
    # (and merging earliest/latest dates created above in .dts)
      dmt_dts_reg <- registry_dmt_sub %>% 
        filter(start >= '2014-01-01') %>% 
        group_by(id_participant) %>% 
        mutate(DMT_Earliest_Study_Date = min(start, na.rm=T)) %>% 
        distinct(id_participant, DMT_Earliest_Study_Date) %>% 
        full_join(.dts, by="id_participant") %>% 
        select(id_participant
               , Registry_DMT_Study_Start = DMT_Earliest_Study_Date
               , Registry_DMT_Earliest_Start = DMT_Earliest_Date
               , Registry_DMT_MostRecent = DMT_Most_Recent_Date) %>% 
        ungroup()


      

# Merging ####
  # dmt_dts_ehr; dmt_dts_reg
  dmt_total <- full_join(dmt_dts_ehr, dmt_dts_reg, by="id_participant") %>% 
    group_by(PATIENT_NUM) %>% 
    mutate(DMT = dmt_text
           , DMT_Earliest_Start = min(EHR_DMT_Earliest_Start, Registry_DMT_Earliest_Start, na.rm=T)
           # , DMT_Recent_Prestudy = min(EHR_DMT_MostRecent, Registry_DMT_MostRecent, na.rm=T)
  ) %>% 
    ungroup()
  
    # some id_participants may not have a corresponding PATIENT_NUM
      # this indicates no EHR record of the DMT of interest, but does NOT 
      # mean a patient is only in the registry (registry is subset of EHR)
      # this can be verified by merging back in the ID link as a sanity check that the PATIENT_NUM exists


# Exclusion criteria ####
  
  ## EHR ####
    .ehr_HighEff_exclude <- ehr_rxnorm %>%
      filter(High_Eff_Exclude == 1) %>% 
      group_by(PATIENT_NUM, id_participant) %>% 
      summarize(EHR_HighEff_Earliest = min(start_date, na.rm=T)
                , EHR_HighEff_Recent = max(start_date, na.rm=T))
    
    .ehr_StdEff_exclude <- ehr_rxnorm %>%
      filter(Std_Eff_Exclude == 1) %>% 
      group_by(PATIENT_NUM, id_participant) %>% 
      summarize(EHR_StdEff_Earliest = min(start_date, na.rm=T)
                , EHR_StdEff_Recent = max(start_date, na.rm=T))
    
    .ehr_Chemo_exclude <- ehr_rxnorm %>%
      filter(Chemo_Exclude == 1) %>% 
      group_by(PATIENT_NUM, id_participant) %>% 
      summarize(EHR_Chemo_Earliest = min(start_date, na.rm=T)
                , EHR_Chemo_Recent = max(start_date, na.rm=T))
  
    .ehr_dmt_exclude <- .ehr_HighEff_exclude %>% 
      full_join(.ehr_StdEff_exclude, by=c("PATIENT_NUM", "id_participant")) %>% 
      full_join(.ehr_Chemo_exclude, by=c("PATIENT_NUM", "id_participant")) 
      
    
  ## Registry #### 
    .reg_HighEff_exclude <- registry_dmt_cln %>%
      filter(High_Eff_Exclude == 1) %>% 
      group_by(id_participant) %>% 
      summarize(Reg_HighEff_Earliest = min(start, na.rm=T)
                , Reg_HighEff_Recent = max(start, na.rm=T))
    
    .reg_StdEff_exclude <- registry_dmt_cln %>%
      filter(Std_Eff_Exclude == 1) %>% 
      group_by(id_participant) %>% 
      summarize(Reg_StdEff_Earliest = min(start, na.rm=T)
                , Reg_StdEff_Recent = max(start, na.rm=T))
    
    .reg_Chemo_exclude <- registry_dmt_cln %>%
      filter(Chemo_Exclude == 1) %>% 
      group_by(id_participant) %>% 
      summarize(Reg_Chemo_Earliest = min(start, na.rm=T)
                , Reg_Chemo_Recent = max(start, na.rm=T))
    
    .reg_dmt_exclude <- .reg_HighEff_exclude %>% 
      full_join(.reg_StdEff_exclude, by="id_participant") %>% 
      full_join(.reg_Chemo_exclude, by="id_participant") 
    
  
  ## Merging Together 
    dmt_exclude <- .ehr_dmt_exclude %>% 
      full_join(.reg_dmt_exclude, by="id_participant") %>% 
      # rowwise() %>% 
      mutate(DMT_HighEff_Earliest = pmin(EHR_HighEff_Earliest, Reg_HighEff_Earliest, na.rm=T)
             , DMT_HighEff_Recent = pmax(EHR_HighEff_Earliest, Reg_HighEff_Earliest, na.rm=T)
             
             , DMT_StdEff_Earliest = pmin(EHR_StdEff_Earliest, Reg_StdEff_Earliest, NA, na.rm=T)
             , DMT_StdEff_Recent = pmax(EHR_StdEff_Earliest, Reg_StdEff_Earliest, NA, na.rm=T)
             , DMT_Chemo_Earliest = pmin(EHR_Chemo_Earliest, Reg_Chemo_Earliest, na.rm=T)
             , DMT_Chemo_Recent = pmax(EHR_Chemo_Earliest, Reg_Chemo_Earliest, na.rm=T)             
      ) %>%
      select(PATIENT_NUM, id_participant, DMT_HighEff_Earliest:DMT_Chemo_Recent)
        # comment the select statement above if we end up interested in EHR and Registry specific
          # exclusion for drug classes
    
# Conclusion ####
  ## Cleaning #### 
    .ls_new <- setdiff(ls(all.names = T), c(.ls_init, "dmt_exclude"))
    rm(.ls_new)
  
  ## Printing #### 
    
    gc()

    print("Script 01_dmt_cleaning.R complete")
    print("Primary outputs: `dmt_total` and `dmt_exclude` data.frames")
    print("'`dmt_total` contains records of earliest dates for DMT of interest (both total and only within study timeperiod)")
    print("`dmt_exclude` data frame contains earliest and most recent date of receiging non-DMT of interest high- and standard efficacy drugs as well as chemotherapy excluding drugs")
    cat("\n\n")
    
      
# Appendix ####
  # Pre-study DMT of interest exclusion
    dmt_dts_ehr %>% filter(EHR_DMT_Earliest_Start < EHR_DMT_Study_Start)
      # n=15 for BCD would be exxcluded via EHR 
    dmt_dts_reg %>% filter(Registry_DMT_Earliest_Start < Registry_DMT_Study_Start)
      # n=1 would be excluded for BCD via registry) (and would already be excluded by EHR)
    
  # comparison table tbd 
  # dmt_total %>% count(EHR_DMT_After2014 = !is.na(EHR_DMT_Study_Start)
  #                     , Registry_DMT_After2014 = !is.na(Registry_DMT_Study_Start)
  #                     )
