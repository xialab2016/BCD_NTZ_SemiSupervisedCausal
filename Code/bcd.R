# BCD 
## Dominic DiSanto
## Updated 5.22.2023


# Loading `here` and `pacman` packages ####
  # Full library loading contained in 00_shared_setup.R 
  if("pacman" %in% installed.packages()[,1]){
    library(pacman)
  }else{
    install.packages("pacman")
    library(pacman)
  }
  
  p_load(here)


# Parameters ####

  # DMT of interest (see appendix for example)
  drugs <- data.frame(RxNormID = c("1876366", "121191", "712566")
                      , Type_Generic = c("Ocrelizumab", "Rituximab", "Ofatumumab")
                      , Type_BrandName = c("Ocrevus", "Kesimpta", "Rituxan"))
  
  dmt_text <- "BCD"
  
  .high_eff <- T
  
  .eval_yrs = 2
  .eval_pd <- lubridate::dyears(.eval_yrs)   # Post-DMT evaluation period of interest (must be lubridate object in days, months, or years)
   
  .change_window <- 6 # specified in months, can take fractions, specified 1/4.35 returns one week     
  
  .export_cohort_steps <- T # defaults to (.export_filename)_derivation.csv under "Results"
  
  .export_data <- T
  
  if(.change_window==6){
    .export_filename <- .final_cohort_name <- paste0("bcd_allMS_", .eval_yrs, "yr_2024_10_07")  
  }else{
    .export_filename <- .final_cohort_name <- paste0("bcd_allMS_3monthWindow_", .eval_yrs, "yr_2024_10_07")
  }
  
  
# Cohort Creation #### 
  cat("\nDrug(s) of Interest (stored in data.frame `.drugs`):\n\n")
  print(drugs)
  
  dir.create(here("Data", "created_data", dmt_text), showWarnings=F)
  
  source(here("Code", "00_shared_setup.R"))
  
  if(".skip_build" %nin% ls(all.names = T)){ # used to skip creation when calling driver scripts from analytic files
    # .ls_init <- ls(all.names = T)
    source(here("Code", "cohort_building.R"))
    # rm(list = setdiff(ls(all.names = T), .ls_init))
    gc() 
  }
