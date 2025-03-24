
.mad = function(x, na.rm=F) mean(abs(x-mean(x, na.rm=na.rm)), na.rm=na.rm)

library(pacman)
p_load(here, dplyr, ggplot2, RColorBrewer, stringr)

source(here("Code", "analysis_scripts", "00_causal_core_Robust_Perturbed.R"))

.save = F

perturb_increase_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Increase_2yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
perturb_increase_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Increase_3yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
perturb_decrease_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Decrease_2yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
perturb_decrease_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Decrease_3yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
perturb_average_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Average_2yr_BCD_NTZ_ATE_BCD_NTZ.RDS"))
perturb_average_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Average_3yr_BCD_NTZ_ATE_BCD_NTZ.RDS"))

.p = setdiff(grep("perturb_", ls(), value=T), grep("perturb_CC", ls(), value=T))
lapply(.p, function(x){
  .x = get(x)
  .op =  quantile(.x$ATE_Perturbed
            , c(0.025, 0.975)
            , na.rm=T
            )
  # .op = sum(is.na(.x$ATE_Perturbed))
  names(.op) = x
  return(.op)
} )

# # check for non-zero treatment coefficient in oucome models 
  # lapply(list(perturb_increase_2yr, perturb_increase_3yr
  #             , perturb_decrease_2yr, perturb_decrease_2yr
  #             , perturb_average_2yr, perturb_average_3yr)
  #        , function(x) x[["Outcome_ImputeModel_Coefs"]][["OC"]][["Trt"]])


.results = function(perturb
                    , .alpha=0.05
                    , round = 4){

  # .cl = lapply(lapply(perturb[["ATE_Perturbed"]], unlist), class) %>% unlist()
  # table(.cl) %>% prop.table()
  # # <1% error rate (~6 of 1300)
  # perturb[["ATE_Perturbed"]][which(.cl=="list")]
  
  n_labelled = sum(!is.na(perturb$Outcome_ImputeModel$outcome))
  perturbATE = unlist(perturb[["ATE_Perturbed"]])#[which(.cl=="numeric")])#[1:1000]
  summary(perturbATE)
  
  
  c(ATE = perturb$ATE_full
    , `Std. Err` = sd(perturbATE)# / sqrt(n_labelled)
    , `P-Value` = perturb.pvalue(perturbATE)
    , Normal.CI.L = perturb$ATE_full + qnorm(.alpha/2)*sd(perturbATE) #/ sqrt(n_labelled)
    , Normal.CI.H = perturb$ATE_full + qnorm(1-.alpha/2)*sd(perturbATE) #/ sqrt(n_labelled)
    # , MAD = mean(abs(perturbATE - mean(perturbATE)))
    # , Normal.CI.L = perturb$ATE_full + qnorm(.alpha/2)*mean(abs(perturbATE - mean(perturbATE))) / sqrt(n_labelled)
    # , Normal.CI.H = perturb$ATE_full + qnorm(1-.alpha/2)*mean(abs(perturbATE - mean(perturbATE))) / sqrt(n_labelled)
    , Empirical.CI.L = quantile(perturbATE, .alpha/2, names=F)
    , Empirical.CI.H = quantile(perturbATE, 1-.alpha/2, names = F)
    ) %>% round(round)
  
}  


ATE.forestplot = function(main.df, .xlim=c(-0.3, 0.3), .subtitle=NULL){
  
  main.ATE = main.df %>%
    mutate(across(-c(Outcome, `Eval. Window`), as.numeric)) %>% 
    rename(`Evaluation Window` = `Eval. Window`) %>% 
    ggplot(aes(x = ATE, y=Outcome, color=`Evaluation Window`)) + 
    geom_point(size=3, stroke = 0.5,position=position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin=Empirical.CI.L, xmax=Empirical.CI.H)
                   , height=.1, lwd=0.4, position=position_dodge(width=0.5)) +
    # labs(title=paste0("BCD (vs. NTZ) ATE")
    #      , subtitle = .subtitle) + 
    xlab("ATE (95% CI)") + 
    ylab("Outcome") + 
    geom_vline(xintercept=0, color='black') +
    theme_minimal() +
    # scale_linetype_manual(name = "95% CI", values = c("Normal Approx"=1, "Perturbation"=2))
    theme(legend.position=c(0.85, 0.8)
          , axis.text=element_text(size=10)
          , title = element_text(face="bold")) +
    # scale_color_manual(values = brewer.pal(n=3, name="Dark2")) +
    scale_color_manual(values = brewer.pal(9, "RdBu")[c(1, 9)]) +
    # scale_color_discrete(type = c("#EF8A62", "#67A9CF"))
    xlim(.xlim)
  
  return(main.ATE)
}

ATE.forestplot.3yr = function(main.df, .xlim=c(-0.3, 0.3), .subtitle=NULL){
  
  main.ATE = main.df %>%
    mutate(across(-c(Outcome, `Eval. Window`), as.numeric)) %>% 
    rename(`Evaluation Window` = `Eval. Window`) %>% 
    ggplot(aes(x = ATE, y=Outcome, color=`Evaluation Window`)) + 
    geom_point(size=3, stroke = 0.5,position=position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin=Empirical.CI.L, xmax=Empirical.CI.H)
                   , height=.1, lwd=0.4, position=position_dodge(width=0.5)) +
    # labs(title=paste0("BCD (vs. NTZ) ATE")
    #      , subtitle = .subtitle) + 
    xlab("ATE (95% CI)") + 
    ylab("Outcome") + 
    geom_vline(xintercept=0, color='black') +
    theme_minimal() +
    # scale_linetype_manual(name = "95% CI", values = c("Normal Approx"=1, "Perturbation"=2))
    theme(legend.position=c(0.85, 0.8)
          , axis.text=element_text(size=10)
          , title = element_text(face="bold")) +
    # scale_color_manual(values = brewer.pal(n=3, name="Dark2")) +
    scale_color_manual(values = brewer.pal(9, "RdBu")[c(1, 9)]) +
    # scale_color_discrete(type = c("#EF8A62", "#67A9CF"))
    xlim(.xlim)
  
  return(main.ATE)
}


SI.df = bind_rows(c(Model = "3-year", .results(perturb_increase_3yr))
                  , c(Model = "2-year", .results(perturb_increase_2yr))
                  ) 

SD.df = bind_rows(c(Model = "3-year", .results(perturb=perturb_decrease_3yr))
                  , c(Model = "2-year", .results(perturb=perturb_decrease_2yr))
                  )

Avg.df = bind_rows(c(Model = "3-year", .results(perturb_average_3yr))
                   , c(Model = "2-year", .results(perturb_average_2yr))
                   )

main.df = union(SI.df %>% mutate(Outcome = "Sustained Worsening")
                , SD.df %>% mutate(Outcome = "Sustained Improvement")) %>%
          # union(Avg.df %>% mutate(Outcome = "Average PDDS")) %>%
  mutate(`Eval. Window` = Model) %>% select(-Model)

main.df %>% select(Outcome, Yr = `Eval. Window`, ATE, Empirical.CI.L, Empirical.CI.H) %>%
  mutate(across(ATE:Empirical.CI.H, as.numeric)) %>% 
  mutate(LowDiff = ATE - Empirical.CI.L
         , HighDiff = abs(ATE - Empirical.CI.H)) 

ATE.forestplot(main.df = main.df
               # , .xlim=c(-1, 1)
               , .xlim=c(-0.45, 0.45)
               )

if(.save) ggsave(here("Results", "BCD_NTZ", "Fig3_Main_ATE_Forest.png")
                 , height = 3.5, width=7, units="in")


  ## Storing Table 1 #### 
    SI.df = SI.df %>% mutate(across(ATE:Empirical.CI.H, ~ round(as.numeric(.x), 3) )) %>% 
      mutate(across(ATE:Empirical.CI.H, ~ sprintf("%.3f", .x)))
    
    SD.df = SD.df %>% mutate(across(ATE:Empirical.CI.H, ~ round(as.numeric(.x), 3) )) %>% 
      mutate(across(ATE:Empirical.CI.H, ~ sprintf("%.3f", .x)))
    
    tbl1 = SI.df %>%
      mutate(Outcome = "Sustained Worsening") %>% 
      union(SD.df %>% mutate(Outcome = "Sustained Improvement")) %>% 
      mutate(`Adjusted P-Value` = sprintf("%.3f", p.adjust(`P-Value`, "holm"))) %>% 
      mutate(`Normal CI` = paste0("(", Normal.CI.L, ", ", Normal.CI.H, ")")
             , `Perturbed CI` = paste0("(", Empirical.CI.L, ", ", Empirical.CI.H, ")")
      ) %>% select(Model, Outcome, ATE, `Std. Err`
                   , `Perturbed CI`, `P-Value`, `Adjusted P-Value`) %>% 
      arrange(Outcome, Model)
    
    # tbl1 %>% filter(Model=="3-year") %>% mutate(`Adjusted P-Value` = sprintf("%.3f", p.adjust(`P-Value`, "holm")))
    if(.save) write.csv(tbl1, here("Results", "BCD_NTZ", "Table1_ATEResults.csv"))

