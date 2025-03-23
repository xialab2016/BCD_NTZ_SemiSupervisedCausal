
.mad = function(x, na.rm=F) mean(abs(x-mean(x, na.rm=na.rm)), na.rm=na.rm)

library(pacman)
p_load(here, dplyr, ggplot2, RColorBrewer, stringr)

source(here("Code", "analysis_scripts", "00_causal_core_Robust_Perturbed.R"))

.save = F

perturb_increase_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Increase_2yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
perturb_increase_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Increase_3yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
perturb_decrease_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Decrease_2yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
perturb_decrease_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Decrease_3yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
perturb_average_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Average_3yr_BCD_NTZ_ATE_BCD_NTZ.RDS"))
perturb_average_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Average_2yr_BCD_NTZ_ATE_BCD_NTZ.RDS"))

.results = function(perturb
                    , .alpha=0.05
                    , round = 4){

  # .cl = lapply(lapply(perturb[["ATE_Perturbed"]], unlist), class) %>% unlist()
  # table(.cl) %>% prop.table()
  # # <1% error rate (~6 of 1300)
  # perturb[["ATE_Perturbed"]][which(.cl=="list")]
  
  n_labelled = sum(!is.na(perturb$Outcome_ImputeModel$outcome))
  perturbATE = unlist(perturb[["ATE_Perturbed"]])#[which(.cl=="numeric")])#[1:1000]
  if( mean(is.na(perturbATE)) > 0.01) stop(">1% missing perturbations, review prior to analysis")
  perturbATE = na.omit(perturbATE)
  
  
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


# Main Analyses #### 
  # Needed to include on plots 

SI.df = bind_rows(c(Model = "3-year", .results(perturb_increase_3yr))
                  , c(Model = "2-year", .results(perturb_increase_2yr))
                  ) 

SD.df = bind_rows(c(Model = "3-year", .results(perturb_decrease_3yr))
                  , c(Model = "2-year", .results(perturb_decrease_2yr))
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

# Sensitivity Analyses (Individual Plots) #### 

## 95% Trim ####
  SI_2yr_95trim = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                               # , "Archive"
                               , "Increase_2yr_95trim_SensitivityAnalysis_2025_95trim_SensitivityAnalysis.RDS"))
  SI_3yr_95trim = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                               , "Increase_3yr_95trim_SensitivityAnalysis_2025_95trim_SensitivityAnalysis.RDS"))
  SD_2yr_95trim = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                               , "Decrease_2yr_95trim_SensitivityAnalysis_2025_95trim_SensitivityAnalysis.RDS"))
  SD_3yr_95trim = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                               , "Decrease_3yr_95trim_SensitivityAnalysis_2025_95trim_SensitivityAnalysis.RDS"))
  
  SI.95.df = bind_rows(c(Model = "3-year", .results(SI_3yr_95trim))
                    , c(Model = "2-year", .results(SI_2yr_95trim))
                    )
  
  SD.95.df = bind_rows(c(Model = "3-year", .results(SD_3yr_95trim))
                    , c(Model = "2-year", .results(SD_2yr_95trim))
                    )
  
  
  df.95 = union(SI.95.df %>% mutate(Outcome = "Sustained Increase")
                  , SD.95.df %>% mutate(Outcome = "Sustained Decrease")) %>% 
    mutate(`Eval. Window` = Model) %>% select(-Model)
  
  sens.plot.95 = ATE.forestplot(main.df = df.95, .subtitle = "95% Feature Trim Sensitivity Analysis")
  if(.save) ggsave(here("Results", "BCD_NTZ", "Sensitivity", "Sensitivity_AllMS_ATE_Forest_95Trim.png"), sens.plot.95)



## All Features ####
  SI_2yr_AllFts = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                               , "Increase_2yr_AllFt_SensitivityAnalysis_2025_AllFt_SensitivityAnalysis.RDS"))
  SI_3yr_AllFts = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                               , "Increase_3yr_AllFt_SensitivityAnalysis_2025_AllFt_SensitivityAnalysis.RDS"))
  SD_2yr_AllFts = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                               , "Decrease_2yr_AllFt_SensitivityAnalysis_2025_AllFt_SensitivityAnalysis.RDS"))
  SD_3yr_AllFts = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                               , "Decrease_3yr_AllFt_SensitivityAnalysis_2025_AllFt_SensitivityAnalysis.RDS"))
  
  SI.AllFt.df = bind_rows(c(Model = "3-year", .results(SI_3yr_AllFts))
                       , c(Model = "2-year", .results(SI_2yr_AllFts))
                       )
  
  SD.AllFt.df = bind_rows(c(Model = "3-year", .results(SD_3yr_AllFts))
                       , c(Model = "2-year", .results(SD_2yr_AllFts))
                       )
  
  df.AllFt = union(SI.AllFt.df %>% mutate(Outcome = "Sustained Increase")
                , SD.AllFt.df %>% mutate(Outcome = "Sustained Decrease")) %>% 
    mutate(`Eval. Window` = Model) %>% select(-Model)

  
  sens.plot.AllFt = ATE.forestplot(main.df = df.AllFt, .subtitle="All Feature Inclusion Sensitivity Analysis")
  if(.save) ggsave(here("Results", "BCD_NTZ", "Sensitivity", "Sensitivity_AllMS_ATE_Forest_AllFt.png"), sens.plot.AllFt)

## 6 Month Lookback ####
  SI_2yr_6mo = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                               , "Increase_2yr_6month_SensitivityAnalysis_2025_6month_SensitivityAnalysis.RDS"))
  SI_3yr_6mo = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                               , "Increase_3yr_6month_SensitivityAnalysis_2025_6month_SensitivityAnalysis.RDS"))
  SD_2yr_6mo = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                               , "Decrease_2yr_6month_SensitivityAnalysis_2025_6month_SensitivityAnalysis.RDS"))
  SD_3yr_6mo = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                               , "Decrease_3yr_6month_SensitivityAnalysis_2025_6month_SensitivityAnalysis.RDS"))
  
  SI.6mo.df = bind_rows(c(Model = "3-year", .results(SI_3yr_6mo))
                          , c(Model = "2-year", .results(SI_2yr_6mo))
                          )
  
  SD.6mo.df = bind_rows(c(Model = "3-year", .results(SD_3yr_6mo))
                          , c(Model = "2-year", .results(SD_2yr_6mo))
                          )
  
  df.6mo = union(SI.6mo.df %>% mutate(Outcome = "Sustained Increase")
                   , SD.6mo.df %>% mutate(Outcome = "Sustained Decrease")) %>% 
    mutate(`Eval. Window` = Model) %>% select(-Model)
  
  
  sens.plot.6mo = ATE.forestplot(main.df = df.6mo, .subtitle = "6-month pre-DMT Lookback Window Sensitivity Analysis"
                                 , .xlim = c(-0.35, 0.35))
  if(.save) ggsave(here("Results", "BCD_NTZ", "Sensitivity", "Sensitivity_AllMS_ATE_Forest_6moLookback.png"), sens.plot.6mo)


## 12 Month Lookback ####
  SI_2yr_12mo = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                            , "Increase_2yr_12month_SensitivityAnalysis_2025_12month_SensitivityAnalysis.RDS"))
  SI_3yr_12mo = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                            , "Increase_3yr_12month_SensitivityAnalysis_2025_12month_SensitivityAnalysis.RDS"))
  SD_2yr_12mo = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                            , "Decrease_2yr_12month_SensitivityAnalysis_2025_12month_SensitivityAnalysis.RDS"))
  SD_3yr_12mo = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                            , "Decrease_3yr_12month_SensitivityAnalysis_2025_12month_SensitivityAnalysis.RDS"))
  
  SI.12mo.df = bind_rows(c(Model = "3-year", .results(SI_3yr_12mo))
                        , c(Model = "2-year", .results(SI_2yr_12mo))
                        )
  
  SD.12mo.df = bind_rows(c(Model = "3-year", .results(SD_3yr_12mo))
                        , c(Model = "2-year", .results(SD_2yr_12mo))
                        )
  
  
  df.12mo = union(SI.12mo.df %>% mutate(Outcome = "Sustained Increase")
                 , SD.12mo.df %>% mutate(Outcome = "Sustained Decrease")) %>% 
    mutate(`Eval. Window` = Model) %>% select(-Model)
  
  sens.plot.12mo = ATE.forestplot(main.df = df.12mo, .subtitle = "12-month pre-DMT Lookback Window Sensitivity Analysis"
                                  , .xlim=c(-1, 1)*0.35)
  if(.save) ggsave(here("Results", "BCD_NTZ", "Sensitivity", "Sensitivity_SD_AllMS_ATE_Forest_12moLookback.png"), sens.plot.12mo)


  ## Three-Month Window Outcome Definition ####
    SI_2yr_3moWindow = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                                    , "Increase_2yr_3monthWindow_SensitivityAnalysis_2025_ATE_BCD_NTZ.RDS"))
    SI_3yr_3moWindow = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                                    , "Increase_3yr_3monthWindow_SensitivityAnalysis_2025_ATE_BCD_NTZ.RDS"))
    SD_2yr_3moWindow = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                                    , "Decrease_2yr_3monthWindow_SensitivityAnalysis_2025_ATE_BCD_NTZ.RDS"))
    SD_3yr_3moWindow = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity"
                                    , "Decrease_3yr_3monthWindow_SensitivityAnalysis_2025_ATE_BCD_NTZ.RDS"))
    
    SI.3moWindow.df = bind_rows(c(Model = "3-year", .results(SI_3yr_3moWindow))
                           , c(Model = "2-year", .results(SI_2yr_3moWindow))
    )
    
    SD.3moWindow.df = bind_rows(c(Model = "3-year", .results(SD_3yr_3moWindow))
                           , c(Model = "2-year", .results(SD_2yr_3moWindow))
    )
    
    
    df.3moWindow = union(SI.3moWindow.df %>% mutate(Outcome = "Sustained Increase")
                    , SD.3moWindow.df %>% mutate(Outcome = "Sustained Decrease")) %>% 
      mutate(`Eval. Window` = Model) %>% select(-Model)
    
    sens.plot.3moWindow = ATE.forestplot(main.df = df.3moWindow, .subtitle = "12-month pre-DMT Lookback Window Sensitivity Analysis"
                                    , .xlim=c(-1, 1)*0.35)
    if(.save) ggsave(here("Results", "BCD_NTZ", "Sensitivity", "Sensitivity_SD_AllMS_ATE_Forest_12moLookback.png"), sens.plot.12mo)
  

# Combining #### 

  ## (Supp Fig) SI Sens Plot ####
    SI.AllATEs = bind_rows(SI.df %>% mutate(Analysis = "Main")
                           , SI.95.df %>% mutate(Analysis = "95% Trim")
                           , SI.AllFt.df %>% mutate(Analysis = "All Analytic Features")
                           , SI.6mo.df %>% mutate(Analysis = "6-month Lookback")
                           , SI.12mo.df %>% mutate(Analysis = "12-month Lookback")
                           , SI.3moWindow.df %>% mutate(Analysis = "3-month Outcome Window")
                           ) %>% 
      mutate(Sens = factor(Analysis, levels = rev(c("Main"
                                                    , "All Analytic Features"
                                                    , "95% Trim"
                                                    , "3-month Outcome Window"
                                                    , "6-month Lookback"
                                                    , "12-month Lookback")
                                                  )
                           )
             )
    
      
      .si.all.plot = SI.AllATEs %>% 
        mutate(Sensitivity = Analysis!="Main"
               , Year = paste0(Model, " Endpoint")
        ) %>% 
        mutate(across(-c(Year, Model, Sens, Analysis), as.numeric)) %>% 
        ggplot(aes(x = ATE, y=Sens, color=factor(Sensitivity))) + 
        geom_point(size=3) +
        # geom_errorbarh(aes(xmin=Normal.CI.L, xmax=Normal.CI.H, lty="Normal Approx")
        #                , height=.1, lwd=0.6, color="navy") +
        geom_errorbarh(aes(xmin=Empirical.CI.L, xmax=Empirical.CI.H, lty="Perturbation")
                       , height=.1, lwd=0.4) +
        # scale_y_continuous(breaks=1:nrow(df), labels=df$study) +
        # labs(title="BCD (vs. NTZ) ATE (Full Results, Sustained Increase)") + 
        xlab("ATE (95% CI)") + 
        ylab("Model") + 
        geom_vline(xintercept=0, color='black') +
        theme_minimal() +
        scale_color_manual(values=c("black", "gray")) + 
        theme(legend.position="none")  +
        xlim(c(-0.4, 0.4)) + 
        facet_wrap(~Year, nrow = 2, strip.position="left") + 
        theme(strip.placement = "outside"
              , strip.text.y = element_text(face="bold", vjust = 3)
        )
      
      
      if(.save) ggsave(plot = .si.all.plot
                       , here("Results", "BCD_NTZ", "Sensitivity", "SI_All_AllMS_ATE_Plot.png")
                       , height = 4, width=7, units="in")

  ## (Supp Fig) SD Sens Plot ####
    SD.AllATEs = bind_rows(SD.df %>% mutate(Analysis = "Main")
                           , SD.95.df %>% mutate(Analysis = "95% Trim")
                           , SD.AllFt.df %>% mutate(Analysis = "All Analytic Features")
                           , SD.6mo.df %>% mutate(Analysis = "6-month Lookback")
                           , SD.12mo.df %>% mutate(Analysis = "12-month Lookback")
                           , SD.3moWindow.df %>% mutate(Analysis = "3-month Outcome Window")
                           ) %>% 
        mutate(Sens = factor(Analysis, levels = rev(c("Main"
                                                  , "All Analytic Features"
                                                  , "95% Trim"
                                                  , "3-month Outcome Window"
                                                  , "6-month Lookback"
                                                  , "12-month Lookback")))
               )
  
    
    .sd.all.plot = SD.AllATEs %>%
      mutate(Sensitivity = Analysis!="Main"
             , Year = paste0(Model, " Endpoint")
             ) %>% 
      mutate(across(-c(Year, Model, Sens, Analysis), as.numeric)) %>% 
      ggplot(aes(x = ATE, y=Sens, color=factor(Sensitivity))) + 
      geom_point(size=3) +
      # geom_errorbarh(aes(xmin=Normal.CI.L, xmax=Normal.CI.H, lty="Normal Approx")
      #                , height=.1, lwd=0.6, color="navy") +
      geom_errorbarh(aes(xmin=Empirical.CI.L, xmax=Empirical.CI.H, lty="Perturbation")
                     , height=.1, lwd=0.4) +
      # scale_y_continuous(breaks=1:nrow(df), labels=df$study) +
      # labs(title="BCD (vs. NTZ) ATE (Full Results, Sustained Decrease)") + 
      xlab("ATE (95% CI)") + 
      ylab("") +
      geom_vline(xintercept=0, color='black') +
      theme_minimal() +
      scale_color_manual(values=c("black", "gray")) + 
      theme(legend.position="none")  +
      xlim(c(-0.4, 0.4)) + 
        theme(axis.title.x = element_text(face="bold")
              , axis.text.x =  element_text(face="bold")
              , axis.text.y =  element_text(face="bold")) + 
      facet_wrap(~Year, nrow = 2, strip.position="left") + 
      theme(strip.placement = "outside"
            , strip.text.y = element_text(face="bold", vjust = 3)
            )

    if(.save) ggsave(plot = .sd.all.plot
                     , here("Results", "BCD_NTZ", "Sensitivity", "SD_All_AllMS_ATE_Plot.png")
                     , height = 4, width=7, units="in")
    

    
  ## (Supp Tbl) SI Sens Table ####
    .si.supp.tbl = SI.AllATEs %>% 
      mutate(ATE = round(as.numeric(ATE), 3) 
             , `Std. Err` = round(as.numeric(`Std. Err`), 3)
             , Empirical.CI.L = round(as.numeric(Empirical.CI.L), 3)
             , Empirical.CI.H = round(as.numeric(Empirical.CI.H), 3)
      ) %>% 
      mutate(Outcome = "Sustained Worsening"
             , `95% CI` = paste0("(", Empirical.CI.L, ", ", Empirical.CI.H, ")") 
             ) %>% 
      select(Outcome, Endpoint = Model, Sens
             , ATE, `Std. Err`
             , `95% CI`) %>% 
      arrange(Endpoint, desc(Sens))
    
    if(.save) write.csv(.si.supp.tbl, here("Results", "BCD_NTZ", "Sensitivity", "SI_SensitivityTable.csv"))
    
  ## (Supp Tbl) SD Sens Table ####
    .sd.supp.tbl = SD.AllATEs %>% 
      mutate(ATE = round(as.numeric(ATE), 3) 
             , `Std. Err` = round(as.numeric(`Std. Err`), 3)
             , Empirical.CI.L = round(as.numeric(Empirical.CI.L), 3)
             , Empirical.CI.H = round(as.numeric(Empirical.CI.H), 3)
      ) %>% 
      mutate(Outcome = "Sustained Improvement"
             , `95% CI` = paste0("(", Empirical.CI.L, ", ", Empirical.CI.H, ")") 
      ) %>% 
      select(Outcome, Endpoint = Model, Sens
             , ATE, `Std. Err`
             , `95% CI`) %>% 
      arrange(Endpoint, desc(Sens))
    
    if(.save) write.csv(.sd.supp.tbl, here("Results", "BCD_NTZ", "Sensitivity", "SD_SensitivityTable.csv"))
    
# Deprecated Below #### 
# 
# # Sensitivity #### 
# 
# ### (old) Sensitivity Analysis Tables ####
# 
# 
# complete.tbls2 = union(SI.AllATEs %>% mutate(Outcome = "Sustained Increase")
#                        , SD.AllATEs %>% 
#                          mutate(Outcome = "Sustained Decrease")
# ) %>% 
#   mutate(across(ATE:Empirical.CI.H, ~ sprintf("%.3f", round(as.numeric(.x), 3)
#   ) )) %>% 
#   mutate(`95% Perturbation CI` = paste0("(", Empirical.CI.L, ", ", Empirical.CI.H, ")")) %>% 
#   select(Outcome, Model, Analysis, ATE, `Std. Err`, `95% Perturbation CI`) %>% 
#   arrange(Outcome, Model, Analysis==desc("Main"))
# 
# if(.save){
#   # 2yr SI 
#   complete.tbls2 %>% 
#     filter(Outcome=="Sustained Increase" & Model=="2-year") %>% 
#     write.csv(here("Results", "BCD_NTZ", "Sensitivity", "Supp_Tbl3a_SI_2yr_SensResults.csv"))
#   
#   # 3yr SI 
#   complete.tbls2 %>% 
#     filter(Outcome=="Sustained Increase" & Model=="3-year") %>% 
#     write.csv(here("Results", "BCD_NTZ", "Sensitivity", "Supp_Tbl3b_SI_3yr_SensResults.csv"))
#   
#   
#   # 2yr SI 
#   complete.tbls2 %>% 
#     filter(Outcome=="Sustained Decrease" & Model=="2-year") %>% 
#     write.csv(here("Results", "BCD_NTZ", "Sensitivity", "Supp_Tbl3c_SD_2yr_SensResults.csv"))
#   
#   # 2yr SI 
#   complete.tbls2 %>% 
#     filter(Outcome=="Sustained Decrease" & Model=="3-year") %>% 
#     write.csv(here("Results", "BCD_NTZ", "Sensitivity", "Supp_Tbl3d_SD_2yr_SensResults.csv"))
# }
# 
# 
# # No Crump Trim Comparison 
# # 
# #   ## SI 
# #     SI.df = bind_rows(c(Model = "3-year", .results(perturb_increase_3yr))
# #                       , c(Model = "2-year", .results(perturb_increase_2yr))
# #                       ) 
# #     
# #     notrimInc.df = bind_rows(c(Model = "3-year", .results(inc_3yr_notrim))
# #                              , c(Model = "2-year", .results(inc_2yr_notrim))
# #                              )
# #     
# #     union(SI.df %>% mutate(Trim = "Yes"), notrimInc.df %>% mutate(Trim="No")) %>% 
# #       select(Model, Trim, ATE, Empirical.CI.L, Empirical.CI.H) %>% 
# #       mutate(across(ATE:Empirical.CI.H, as.numeric)) %>% 
# #       mutate(LowDiff = ATE - Empirical.CI.L
# #              , HighDiff = abs(ATE - Empirical.CI.H))
# #     
# #     union(SI.df %>% mutate(Trim = "Yes")
# #           , notrimInc.df %>% mutate(Trim="No")) %>% 
# #       select(Model, Trim, ATE, Empirical.CI.L, Empirical.CI.H) %>%
# #       mutate(across(ATE:Empirical.CI.H, as.numeric)) %>%
# #       mutate(LowDiff = ATE - Empirical.CI.L
# #              , HighDiff = abs(ATE - Empirical.CI.H)) %>%
# #       mutate(Outcome = Trim, `Eval. Window` = Model) %>% 
# #       ATE.forestplot()
# #   
# #   ## SD 
# #     SD.df = bind_rows(c(Model = "3-year", .results(perturb_decrease_3yr))
# #                       , c(Model = "2-year", .results(perturb_decrease_2yr))
# #                       )
# #     
# #     notrimDec.df = bind_rows(c(Model = "3-year", .results(dec_3yr_notrim))
# #                              , c(Model = "2-year", .results(dec_2yr_notrim))
# #                              )
# #     
# #     main.df = union(SI.df %>% mutate(Outcome = "Sustained Increase")
# #                     , SD.df %>% mutate(Outcome = "Sustained Decrease")) %>% 
# #       union(Avg.df %>% mutate(Outcome = "Average PDDS")) %>% 
# #       mutate(`Eval. Window` = Model) %>% select(-Model)
# #     
# # 
# #     union(SD.df %>% mutate(Trim = "Yes"), notrimDec.df %>% mutate(Trim="No")) %>% 
# #       select(Model, Trim, ATE, Empirical.CI.L, Empirical.CI.H) %>% 
# #       mutate(across(ATE:Empirical.CI.H, as.numeric)) %>% 
# #       mutate(LowDiff = ATE - Empirical.CI.L
# #              , HighDiff = abs(ATE - Empirical.CI.H))
# #     
# #     union(SD.df %>% mutate(Trim = "Yes")
# #           , notrimDec.df %>% mutate(Trim="No")) %>% 
# #       select(Model, Trim, ATE, Empirical.CI.L, Empirical.CI.H) %>%
# #       mutate(across(ATE:Empirical.CI.H, as.numeric)) %>%
# #       mutate(LowDiff = ATE - Empirical.CI.L
# #              , HighDiff = abs(ATE - Empirical.CI.H)) %>%
# #       mutate(Outcome = Trim, `Eval. Window` = Model) %>% 
# #       ATE.forestplot()
# #   
# 
# ### Sensitivity (CC vs SS-3mo vs SS-6mo) #### 
#     
#   #### Quick Check #### 
#     
#     .perturb_3mo = grep("perturb_3mo", ls(), value=T)
#     .perturb = setdiff(grep("perturb_", ls(), value=T), .perturb_3mo)
#     
#     lapply(.perturb
#            , function(x) {
#              .x = get(x)
#              .output = sd(.x$ATE_Perturbed)
#              names(.output) = x
#              return(.output)
#            }
#     ) %>% unlist()
#     
#     lapply(.perturb_3mo
#            , function(x) {
#              .x = get(x)
#              .output = sd(.x$ATE_Perturbed)
#              names(.output) = x
#              return(.output)
#            }
#     ) %>% unlist()
#     
#     
#     
#   #### Forest Plot #### 
#     
#     perturb_CC_increase_3yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "perturb_ATE_3yrs_Increase_ATE_BCD_NTZ.RDS"))
#     perturb_CC_increase_2yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "perturb_ATE_2yrs_Increase_ATE_BCD_NTZ.RDS"))
#     perturb_CC_decrease_3yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "perturb_ATE_3yrs_Decrease_ATE_BCD_NTZ.RDS"))
#     perturb_CC_decrease_2yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "perturb_ATE_2yrs_Decrease_ATE_BCD_NTZ.RDS"))
#     perturb_CC_average_3yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "perturb_ATE_3yrs_Average_ATE_BCD_NTZ.RDS"))
#     perturb_CC_average_2yr = readRDS(here("Results", "BCD_NTZ", "Sensitivity", "CompleteCase", "perturb_ATE_2yrs_Average_ATE_BCD_NTZ.RDS"))
#     
#     perturb_3mo_increase_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity", "Increase_2yr_3monthWindow_SensitivityAnalysis_ATE_BCD_NTZ.RDS"))
#     perturb_3mo_increase_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity", "Increase_3yr_3monthWindow_SensitivityAnalysis_ATE_BCD_NTZ.RDS"))
#     perturb_3mo_decrease_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity", "Decrease_2yr_3monthWindow_SensitivityAnalysis_ATE_BCD_NTZ.RDS"))
#     perturb_3mo_decrease_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity", "Decrease_3yr_3monthWindow_SensitivityAnalysis_ATE_BCD_NTZ.RDS"))
#     perturb_3mo_average_3yr = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity", "Average_3yr_3monthWindow_SensitivityAnalysis_ATE_BCD_NTZ.RDS"))
#     perturb_3mo_average_2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "Sensitivity", "Average_2yr_3monthWindow_SensitivityAnalysis_ATE_BCD_NTZ.RDS"))
#     
#     
#     SI.df = bind_rows(c(Endpoint = "3-year", Model = "Semi-Supervised (6-month)", .results(perturb_increase_3yr))
#                       , c(Endpoint = "2-year", Model = "Semi-Supervised (6-month)", .results(perturb_increase_2yr))
#                       , c(Endpoint = "3-year", Model = "Semi-Supervised (3-month)", .results(perturb_3mo_increase_3yr))
#                       , c(Endpoint = "2-year", Model = "Semi-Supervised (3-month)", .results(perturb_3mo_increase_2yr))
#                       , c(Endpoint = "3-year", Model = "Complete Case", .results(perturb_CC_increase_3yr))
#                       , c(Endpoint = "2-year", Model = "Complete Case", .results(perturb_CC_increase_2yr))
#                       )
#     
#     SD.df = bind_rows(c(Endpoint = "3-year", Model = "Semi-Supervised (6-month)", .results(perturb_decrease_3yr))
#                       , c(Endpoint = "2-year", Model = "Semi-Supervised (6-month)", .results(perturb_decrease_2yr))
#                       , c(Endpoint = "3-year", Model = "Semi-Supervised (3-month)", .results(perturb_3mo_decrease_3yr))
#                       , c(Endpoint = "2-year", Model = "Semi-Supervised (3-month)", .results(perturb_3mo_decrease_2yr))
#                       , c(Endpoint = "3-year", Model = "Complete Case", .results(perturb_CC_decrease_3yr))
#                       , c(Endpoint = "2-year", Model = "Complete Case", .results(perturb_CC_decrease_2yr))
#                       )    
#     
#     Avg.df = bind_rows(c(Endpoint = "3-year", Model = "Semi-Supervised (6-month)", .results(perturb_average_3yr))
#                       , c(Endpoint = "2-year", Model = "Semi-Supervised (6-month)", .results(perturb_average_2yr))
#                       , c(Endpoint = "3-year", Model = "Semi-Supervised (3-month)", .results(perturb_3mo_average_3yr))
#                       , c(Endpoint = "2-year", Model = "Semi-Supervised (3-month)", .results(perturb_3mo_average_2yr))
#                       , c(Endpoint = "3-year", Model = "Complete Case", .results(perturb_CC_average_3yr))
#                       , c(Endpoint = "2-year", Model = "Complete Case", .results(perturb_CC_average_2yr))
#                       )
#     
#     main.df = union(SI.df %>% mutate(Outcome = "Sustained Increase")
#                     , SD.df %>% mutate(Outcome = "Sustained Decrease")) %>% 
#       union(Avg.df %>% mutate(Outcome = "Average PDDS"))
#     
#     # main.df %>% select(Outcome, Endpoint, ATE, Empirical.CI.L, Empirical.CI.H) %>%
#     #   mutate(across(ATE:Empirical.CI.H, as.numeric)) %>% 
#     #   mutate(LowDiff = ATE - Empirical.CI.L
#     #          , HighDiff = abs(ATE - Empirical.CI.H))
# 
#     main.df %>%
#       filter(!str_detect(Outcome, "Average")) %>% 
#       mutate(across(-c(Outcome, Endpoint, Model), as.numeric)) %>% 
#       ggplot(aes(x = ATE, y=Outcome, color=Model
#                  # , lty=Endpoint
#                  )
#              ) + 
#       geom_point(size=3, stroke = 0.5,position=position_dodge(width = 0.5)) +
#       geom_errorbarh(aes(xmin=Empirical.CI.L, xmax=Empirical.CI.H)
#                      , height=.1, lwd=0.4, position=position_dodge(width=0.5)) +
#       xlab("ATE (95% CI)") + 
#       ylab("Outcome") + 
#       geom_vline(xintercept=0, color='black') +
#       theme_minimal() +
#       theme(legend.position=c(0.85, 0.8)
#             , axis.text=element_text(size=10)
#             , title = element_text(face="bold")) + 
#       facet_wrap(~Endpoint) + 
#       theme(legend.position = c("bottom"))
# 
#     