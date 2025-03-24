# Shared-Setup 
## Dominic DiSanto
## Updated 3.15.2023


print("Running 00_shared_setup.R...")

# Library Loading ####
  # including installation if not installed 
  if("pacman" %in% installed.packages()[,1]){
      library(pacman)
    }else{
      install.packages("pacman")
      library(pacman)
    }
    
    pacman::p_load(dplyr, magrittr, stringr, lubridate, tidyr # data cleaning
                   , readr, readxl, here, openxlsx #import/export
                   , broom, table1, purrr, flextable # table 1 script 
                   , grid, gridExtra, ggplot2 # viz, table/PDF outputs 
                   , glmnet # analysis, 00_causal_core.R
                   , ordinalNet, rms # ordinal modelling 
                   , DT, data.table # datatable for large NLP dataframe handling
                   , boot # bootstrap estimates in analytic scripts 
                   , ranger, randomForest # RF
                   , np , npreg # smoothing spline regression (np only used)
                   , update=F)    
  
    
# `%nin%` function #### 
  # useful "not in" function (negation of R's %in% function)
  `%nin%` <- Negate(`%in%`)
  # other necessary cleaning 
    summarize = dplyr::summarize 
    
# File Names nad ID Import #### 
    
  id_link_filename <- "R3_2468_Xia_PROMOTE_20220223_studyid_NoPHI.xls"
  
  id_link <- read_xls(here("Data", "ID Linkage", id_link_filename)) %>% 
    select(-LINE_NUM)
              
  ehr_codified_rolled = "UPMC_MS_2004_to_2022_Codified_processed_data_2023-08-25.csv"

  # pt_level_icd <- "UPMC_MS_2004_to_2022_Codified_patient_level_ICD_count_2023-02-17.csv"
    # used for Charlson/Elixhauser data 
  
  registry_dmt_filename <- "DMTHistory.csv"

  nlp_filename = "UPMC_MS_2011_to_2021_NLP_processed_data_2023-10-21.csv"
  
  
  
  
# Knowledge Graph Import for Selection #### 
  
  once_ms_ehr_full = read_csv(here("Data", "Latest", "ONCE_Multiple sclerosis_Codified.csv")) 
  once_ms_nlp_full = read_csv(here("Data", "Latest", "ONCE_Multiple sclerosis_CUI.csv"))
  
  once_disability_ehr_full = read_csv(here("Data", "Latest", "ONCE_disability_PheCode296.2_cos0.165.csv")) 
  once_disability_nlp_full = read_csv(here("Data", "Latest", "ONCE_disability_C0231170_titlecos0.5_titlecut0.3_exactFALSE.csv"))

# Drug Exclusion Lists ####  
  # stdeff RxNorm ID's added manually via RxNav
    # https://mor.nlm.nih.gov/RxNav/
  # Ingredient ID (IN/MIN in RxNav) used (confirmed via Weijing/UPMC team)
  
  ## EHR Codes ####
  
  ## High Eff 
    .HighEff_Exclude_RxNormIDs <- c("117055" # Alemtuzumab/Lemtrada  
                          , "44157" # Mavenclad/cladribine
                          , "7005" #novantrone/mitoxantrone
                          , "354770" # ntz/tysabri
                          , "121191" # rituximab/rituxan (BCD)
                          , "1876366" # Ocrelizumab/ocrevus (BCD)
                          # , "3002" # Cytoxan # non-FDA approved, not used in exclusion list
                          , "712566" # Ofatumumab/kesimpta (BCD)
                          )
    
    .HighEff_Exclude_GenericNames <- c("alemtuzumab"
                                      , "cladribine"
                                      , "mitoxantrone"
                                      , "natalizumab" 
                                      , "ocrelizumab" 
                                      , "ofatumumab"
                                      , "rituximab"
                                      # , "cyclophosphamide" # non-FDA approved, not used in exclusion list
                                      )
    
    .HighEff_Exclude_BrandNames <- c("Lemtrada" # alemtuzumab
                                     , "Mavenclad" # cladribine
                                     , "Novantrone" # mitoxatrone
                                     , "Tysabri" # NTZ 
                                     , "Ocrevus" # ocrelizumab 
                                     , "Kesimpta" # ofatumumab 
                                     , "Rituxan" # rituximab
                                     # , "Cytoxan" # , "cyclophosphamide" # non-FDA approved, not used in exclusion list
                                     )
  
  ## Standard Eff
    .StdEff_Exclude_RxNormIDs <- c("2261783" # diroximel fumarate
                                  , "2532300" # Ponesimod
                                  , "1546433" # monomethyl fumarate 
                                  , "1012892" # Fingolimod: 
                                  , "2121085" # Siponimod: 
                                  , "2288236" # Ozanimod: 
                                  , "1373478" # dimethyl fumarate
                                  # below here I've added (Dominic 2.1.2022)
                                  , "190353" # , "daclizumab"
                                  , "214582" # , "glatiramer acetate"
                                  , "75917" # , "interferon beta-1a"
                                  , "72257" # , "interferon beta-1b"
                                  , "1546168" # , "peginterferon beta-1a"
                                  , "1310520" # , "teriflunomide"
    )
    
    
    .StdEff_Exclude_GenericNames <- c("daclizumab" #?
                                      , "dimethyl fumarate" # 
                                      , "dioroximel fumarate" # 
                                      , "fingolimod" # 
                                      , "glatiramer acetate" #?
                                      , "interferon beta-1a" #?
                                      , "interferon beta-1b" #?
                                      , "monomethyl fumarate" # 
                                      , "ozanimod" # 
                                      , "peginterferon beta-1a" #? 
                                      , "ponesimod" # 
                                      , "siponimod" # 
                                      , "teriflunomide" #?
    )
    
    
    .StdEff_Exclude_BrandNames <- c("Zenapax"
                                    , "Tecfidera"
                                    , "Vumerity"
                                    , "Gilenya"
                                    , "Copaxone"
                                    , "Avonex"
                                    , "Rebif"
                                    , "Betaseron"
                                    , "Extavia"
                                    , "Bafiertam"
                                    , "Zeposia"
                                    , "Plegridy"
                                    , "Ponvory" 
                                    , "Mayzent"
                                    , "Aubagio"
    )
    
    ## Chemo
    
    .Chemo_Exclude_list <- c("117055" # Alemtuzumab/Lemtrada
                             , "3002" # Cytoxan/cyclophosphoamide 
                             , "6851" # methotrexate
    )
    
    .Chemo_Exclude_GenericNames <- c("alemtuzumab"
                                     , "cyclophosphamide"
                                     , "methotrexate"
    )
    
    .Chemo_Exclude_BrandNames <- c("Lemtrada"
                                   , "Cytoxan"
                                   , "Methotrexate"
    )
    
  
    if ("drugs" %in% ls()){
      HighEff_Exclude_RxNormIDs <- setdiff(tolower(.HighEff_Exclude_RxNormIDs), tolower(drugs$RxNormID))
      HighEff_Exclude_GenericNames <- setdiff(tolower(.HighEff_Exclude_GenericNames), tolower(drugs$Type_Generic))
      HighEff_Exclude_BrandNames <- setdiff(tolower(.HighEff_Exclude_BrandNames), tolower(drugs$Type_BrandName))
      
      StdEff_Exclude_RxNormIDs <- setdiff(tolower(.StdEff_Exclude_RxNormIDs), tolower(drugs$RxNormID))
      StdEff_Exclude_GenericNames <- setdiff(tolower(.StdEff_Exclude_GenericNames), tolower(drugs$Type_Generic))
      StdEff_Exclude_BrandNames <- setdiff(tolower(.StdEff_Exclude_BrandNames), tolower(drugs$Type_BrandName))
      
      Chemo_Exclude_RxNormIDs <- setdiff(tolower(.Chemo_Exclude_list), tolower(drugs$RxNormID))
      Chemo_Exclude_GenericNames <- setdiff(tolower(.Chemo_Exclude_GenericNames), tolower(drugs$RxNormID))
      Chemo_Exclude_BrandNames <- setdiff(tolower(.Chemo_Exclude_BrandNames), tolower(drugs$RxNormID))
    }

    cui.pull = function(.drugs){
      .op = unlist(sapply(
        tolower(.drugs)
        , function(x) once_ms_nlp_full %>% 
          filter(str_detect(tolower(term), {{x}})) %>% 
          pull(term)
      )
      )
      names(.op) = NULL 
      return(.op)
    }    

    .StdEff_Exclude_CUIs = cui.pull(c(.StdEff_Exclude_BrandNames, .StdEff_Exclude_GenericNames)) 
    .HighEff_Exclude_CUIs = cui.pull(c(.HighEff_Exclude_BrandNames, .HighEff_Exclude_GenericNames)) 
    
    
    # All drug RXNORM list 
    .rxnorm_exclude <- paste0("RXNORM:", c(.StdEff_Exclude_RxNormIDs, .HighEff_Exclude_RxNormIDs))
    
    # Conclusion ####
    gc() # cleaning unused memory
    
    print("Script 00_shared_setup.R complete\n")
    print("Output includes vectors of drug exclusion (e.g. HighEff_Exclude_RxNormIDs) and dataframe `id_link`")
    cat("\n\n")
    