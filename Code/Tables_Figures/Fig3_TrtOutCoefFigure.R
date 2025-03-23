
# Set-Up ####
  library(pacman)
  p_load(here, ggplot2, Matrix, dplyr, tibble, RColorBrewer, cowplot)
  source(here("Code", "00_shared_setup.R"))
  .save = F

## Upload KG Dictionaries ####
once_ms_ehr = read_csv(here("Data", "Latest", "ONCE_Multiple sclerosis_Codified.csv"))
once_ms_nlp = read_csv(here("Data", "Latest", "ONCE_Multiple sclerosis_CUI.csv"))
once_disability_ehr = read_csv(here("Data", "Latest", "ONCE_disability_PheCode296.2_cos0.165.csv")) 
once_disability_nlp = read_csv(here("Data", "Latest", "ONCE_disability_C0231170_titlecos0.5_titlecut0.3_exactFALSE.csv")) 

ft.dict = rbind(once_ms_ehr %>% select(Variable, term=Description)  
      , once_ms_nlp %>% select(Variable=cui, term)
      , once_disability_ehr %>% select(Variable, term=Description)  
      , once_disability_nlp %>% select(Variable=cui, term)
      )


# (Fig 4) SI-SD-Trt #### 
model.SI = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Increase_3yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
model.SD = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Decrease_3yr_BCD_NTZ_2025_ATE_BCD_NTZ.RDS"))
all.equal(model.SI$Outcome_ImputeModel_Coefs$Trt, model.SD$Outcome_ImputeModel_Coefs$Trt) # TRUE 


coef.df = bind_rows(data.frame(Coefficient = model.SI$Outcome_ImputeModel_Coefs$OC
                               , Model = "Sustained Worsening") %>% 
                      rownames_to_column("Variable")
                    , data.frame(Coefficient = model.SD$Outcome_ImputeModel_Coefs$OC, Model = "Sustained Improvement") %>% 
                      rownames_to_column("Variable")
                    , data.frame(Coefficient = model.SD$Outcome_ImputeModel_Coefs$Trt, Model = "Treatment") %>% 
                      rownames_to_column("Variable")
                    ) %>% 
  as_tibble() %>% 
  mutate(Variable = 
           case_when(Variable == "Trt" ~ "Treatment Class"
                     , Variable == "White_NonHispanic" ~ "White, Non-Hispanic"
                     , Variable == "subject_sex_ch" ~ "Male"
                     , Variable == "Age_at_DMT_Initiation" ~ "Age at Treatment" 
                     , Variable == "Days_Earliest_DMT_to_Study_DMT" ~ "Treatment Duration" 
                     , Variable == "Followup_Duration" ~ "Follow-up Duration" 
                     , Variable == "Disease_Duration" ~ "Disease Duration" 
                     , Variable == "Utilization_Total" ~ "Baseline Healthcare Utilization" 
                     , Variable == "N_MS_CUIs" ~ "MS CUI Frequency" 
                     , Variable == "N_MS_PheCodes" ~ "MS PheCode Frequency"
                     , Variable == "ridge_All_EHR_NLP_latent" ~ "Baseline PDDS Risk"
                     , T ~ Variable
           ) 
         ) %>% 
  filter(Variable!="(Intercept)") %>%
  mutate(Variable = gsub("RXNORM", "RxNorm", Variable)) %>% 
  mutate(Variable = case_when(Variable == "Male" ~ "Gender (Women)"
                              , T ~ Variable
                              )
         , Coefficient = case_when(Variable == "Gender (Women)" ~ -1*Coefficient
                                   , T ~ Coefficient)
         )


  ## Outcome ##### 
    plot.OC = coef.df %>%
  filter(abs(Coefficient)>0 & Model != "Treatment") %>%
      # group_by(Model) %>% 
      # mutate(Coefficient.Std = Coefficient / max(abs(Coefficient))) %>%
      group_by(Variable) %>% 
      mutate(Coefficient.Rank = max(Coefficient)) %>% 
      ungroup() %>%
             # , Y.order = ifelse(Model=="Sustained Worsening", rank(Coefficient), -Inf)
      mutate(Y.order = rank(Coefficient.Rank)) %>%
      arrange(desc(Coefficient)) %>% 
      arrange(Y.order) %>% 
      mutate(Model = relevel(factor(Model), ref="Sustained Worsening")) %>%
      ggplot(aes(x = Model
                 , y=reorder(Variable, Y.order, decreasing=F)
                 , fill = Coefficient
                 # , fill = Coefficient.Std
                 , size = abs(Coefficient)
      )
      ) + 
      labs(color="Coefficient") +
      # labs(color="Coefficient\n(max-normalized)") +
      geom_point(pch=21, color="black") +
      scale_fill_distiller(type="div", palette = "RdBu"
                            # , limits = c(-1, 1)
                            ) + 
      scale_size(name = "Coefficient\nMagnitude") + 
      theme_classic() + 
      ylab("Variable") + 
      scale_y_discrete(expand=expansion(mult=c(0.02,0.02)))

  ## Treatment ####
    plot.Trt = coef.df %>%
      filter(Model == "Treatment") %>%
      mutate(Coefficient.Std = Coefficient / max(abs(Coefficient))) %>%
      filter(abs(Coefficient)>0) %>% 
      arrange(Variable) %>% 
      ggplot(aes(x = Model
                 , y = reorder(Variable, Coefficient, decreasing=F)
                 , fill = Coefficient
                 # , fill = Coefficient.Std
                 , size = abs(Coefficient)
      )
      ) + 
      labs(color="Coefficient") +
      # labs(color="Coefficient\n(max-normalized)") +
      geom_point(pch=21, color="black") +
      scale_size(name = "Coefficient\nMagnitude") + 
      theme_classic() + 
      xlab("") + ylab("") + 
      scale_fill_distiller(type="div", palette="RdBu"
                          , guide = guide_colorbar(frame.colour = "black"
                                                   , ticks.colour = "black")
                          # , limits = c(-1, 1)
                          ) + 
  scale_y_discrete(expand=expansion(mult=c(0.02,0.02)))
    
  ## Combining ####
  coef.plts.combo = plot_grid(plot.OC + theme(legend.position="none"
                            , axis.title.x = element_text(face="bold")
                            , axis.title.y = element_text(face="bold")
                            , axis.text.y = element_text(face="bold")
                            , axis.text.x = element_text(angle=20
                                                         , hjust=0.65
                                                         , vjust=0.77
                                                         , face="bold"
                                                         )
                            )  + xlab("")
            , NULL
            , plot.Trt + theme(legend.position="none"
                               , axis.text.x = element_text(face="bold")
                               , axis.text.y = element_text(face="bold")
                               ) + 
              coord_cartesian(clip="off")
            , get_legend(plot.Trt)
            , rel_widths = c(1.5, 0.5, 1.25, 0.5)
            , nrow = 1
            , align = "h"
            , axis="btl"
            , labels = c("A", "", "B", "")
            )

    
    coef.plts.combo = plot_grid(plot.Trt + theme(legend.position="none"
                                                       , axis.text.x = element_text(face="bold")
                                                       , axis.text.y = element_text(face="bold")
                                                 , axis.title.y = element_text(face="bold")
                                                       ) + coord_cartesian(clip="off") + ylab("Variable")
    , NULL
    , plot.OC + theme(legend.position="none"
                                                , axis.title.x = element_text(face="bold")
                                                , axis.text.y = element_text(face="bold")
                                                , axis.text.x = element_text(angle=20
                                                                             , hjust=0.65
                                                                             , vjust=0.77
                                                                             , face="bold"
                                                )
    )  + xlab("") + ylab("") 
    , get_legend(plot.Trt)
    , rel_widths = c(1.25, 0.1, 1.5, 0.5)
    , nrow = 1
    , align = "h"
    , axis="btl"
    , labels = c("A", "", "B", "")
    )

    

  ## Exporting #### 
  if(.save)  cowplot::save_plot(plot = coef.plts.combo
                                , filename = here("Results", "BCD_NTZ", "Fig4_DiPS_Coefficients.png")
                                , base_height = 7
                                , base_asp = 1.3)

# OMIT (Sup Fig) 2 Year Results #### 

  model.SI.2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Increase_2yr_BCD_NTZ_ATE_BCD_NTZ.RDS"))
  model.SD.2yr = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Decrease_2yr_BCD_NTZ_ATE_BCD_NTZ.RDS"))
  all.equal(model.SI.2yr$Outcome_ImputeModel_Coefs$Trt
            , model.SD.2yr$Outcome_ImputeModel_Coefs$Trt) # TRUE 

  supp.coef.df = bind_rows(data.frame(Coefficient = model.SI.2yr$Outcome_ImputeModel_Coefs$OC
                                 , Model = "Sustained Worsening") %>% 
                        rownames_to_column("Variable")
                      , data.frame(Coefficient = model.SD.2yr$Outcome_ImputeModel_Coefs$OC, Model = "Sustained Improvement") %>% 
                        rownames_to_column("Variable")
                      , data.frame(Coefficient = model.SD.2yr$Outcome_ImputeModel_Coefs$Trt, Model = "Treatment") %>% 
                        rownames_to_column("Variable")
  ) %>% 
    as_tibble() %>% 
    mutate(Variable = 
             case_when(Variable == "Trt" ~ "Treatment Class"
                       , Variable == "White_NonHispanic" ~ "White, Non-Hispanic"
                       , Variable == "subject_sex_ch" ~ "Male Sex"
                       , Variable == "Age_at_DMT_Initiation" ~ "Age" 
                       , Variable == "Days_Earliest_DMT_to_Study_DMT" ~ "Treatment Duration" 
                       , Variable == "Followup_Duration" ~ "Follow-up Duration" 
                       , Variable == "Disease_Duration" ~ "Disease Duration" 
                       , Variable == "Utilization_Total" ~ "Baseline, Total Healthcare Utilization" 
                       , Variable == "N_MS_CUIs" ~ "MS CUI Frequency" 
                       , Variable == "N_MS_PheCodes" ~ "MS PheCode Frequency"
                       , Variable == "ridge_All_EHR_NLP_latent" ~ "Baseline PDDS Risk"
                       , T ~ Variable
             ) 
    ) %>% 
    filter(Variable!="(Intercept)") %>%
    mutate(Variable = gsub("RXNORM", "RxNorm", Variable))
  
  

  ## Outcome ####
  supp.plot.OC = supp.coef.df %>%
    filter(Model != "Treatment") %>%
    group_by(Model) %>% 
    mutate(Coefficient.Std = Coefficient / max(abs(Coefficient))
           , Y.order = rank(Coefficient) 
    ) %>%
    ungroup() %>% 
    filter(abs(Coefficient)>0) %>% 
    arrange(Variable) %>% 
    mutate(Model = relevel(factor(Model), ref="Sustained Worsening")) %>% 
    ggplot(aes(x = Model
               , y=reorder(Variable, Y.order, decreasing=F)
               , fill = Coefficient
               # , fill = Coefficient.Std
               , size = abs(Coefficient)
    )
    ) + 
    labs(color="Coefficient") +
    # labs(color="Coefficient\n(max-normalized)") +
    geom_point(pch=21, color="black") +
    scale_fill_distiller(type="div", palette = "RdBu"
                         # , limits = c(-1, 1)
    ) + 
    scale_size(name = "Coefficient\nMagnitude") + 
    theme_classic() + 
    ylab("Variable")
  
  ## Treatment ####
  supp.plot.Trt = supp.coef.df %>%
    filter(Model == "Treatment") %>%
    mutate(Coefficient.Std = Coefficient / max(abs(Coefficient))) %>%
    filter(abs(Coefficient)>0) %>% 
    arrange(Variable) %>% 
    ggplot(aes(x = Model
               , y=reorder(Variable, Coefficient, decreasing=T)
               , fill = Coefficient
               # , fill = Coefficient.Std
               , size = abs(Coefficient)
    )
    ) + 
    labs(color="Coefficient") +
    # labs(color="Coefficient\n(max-normalized)") +
    geom_point(pch=21, color="black") +
    scale_size(name = "Coefficient\nMagnitude") + 
    theme_classic() + 
    xlab("") + ylab("") + 
    scale_fill_distiller(type="div", palette="RdBu"
                         , guide = guide_colorbar(frame.colour = "black"
                                                  , ticks.colour = "black")
                         # , limits = c(-1, 1)
    )
  
  ## Combining ####
  supp.coef.plts.combo = plot_grid(supp.plot.OC + theme(legend.position="none"
                                              , axis.title.x = element_text(face="bold")
                                              , axis.title.y = element_text(face="bold")
                                              , axis.text.y = element_text(face="bold")
                                              , axis.text.x = element_text(angle=20
                                                                           , hjust=0.65
                                                                           , vjust=0.77
                                                                           , face="bold"
                                              )
  )  + xlab("")
  , NULL
  , supp.plot.Trt + theme(legend.position="none"
                     , axis.text.x = element_text(face="bold")
                     , axis.text.y = element_text(face="bold")
  ) + 
    coord_cartesian(clip="off")
  , get_legend(plot.Trt)
  , rel_widths = c(1.5, 0.5, 1.25, 0.5)
  , nrow = 1
  , align = "h"
  , axis="btl"
  , labels = c("A", "", "B", "")
  )

  
  ## Exporting #### 
  # if(.save)  cowplot::save_plot(plot = supp.coef.plts.combo
  #                               , filename = here("Results", "BCD_NTZ", "SuppFig_DiPS_Coefficients.png")
  #                               , base_height = 7
  #                               , base_asp = 1.3)
  
  
# Deprecated #### 
  
# Plotting 2 and 3-year results combined ####
  ### Deprecated, results not used # 
coef.plot = function(model.2yr, model.3yr){

model.2yr.oc = tibble(Coefficient = model.2yr$Outcome_ImputeModel_Coefs$OC#model.2yr$OC.coef
                      , Variable = names(model.2yr$Outcome_ImputeModel_Coefs$OC)) %>% #names(model.2yr$OC.coef)) %>% 
  mutate(Model = "Outcome", Year = "2-years")

model.3yr.oc = tibble(Coefficient = model.3yr$Outcome_ImputeModel_Coefs$OC #model.3yr$OC.coef
                      , Variable =  names(model.2yr$Outcome_ImputeModel_Coefs$OC)) %>% #names(model.3yr$OC.coef)) %>% 
  mutate(Model = "Outcome", Year = "3-years")


model.2yr.trt = tibble(Coefficient = model.2yr$Outcome_ImputeModel_Coefs$Trt #model.2yr$Trt.coef
                       , Variable = names(model.2yr$Outcome_ImputeModel_Coefs$Trt)) %>% #names(model.2yr$Trt.coef)) %>% 
  mutate(Model = "Treatment", Year = "2-years")

model.3yr.trt = tibble(Coefficient = model.3yr$Outcome_ImputeModel_Coefs$Trt #model.3yr$Trt.coef
                       , Variable = names(model.3yr$Outcome_ImputeModel_Coefs$Trt)) %>% #names(model.3yr$Trt.coef)) %>% 
  mutate(Model = "Treatment", Year = "3-years")


coef.df = bind_rows(model.2yr.oc, model.3yr.oc
                      , model.2yr.trt, model.3yr.trt) %>% 
  mutate(Variable = 
           case_when(Variable == "Trt" ~ "Treatment Class"
                     , Variable == "White_NonHispanic" ~ "White, Non-Hispanic"
                     , Variable == "subject_sex_ch" ~ "Male Sex"
                     , Variable == "Age_at_DMT_Initiation" ~ "Age" 
                     , Variable == "Days_Earliest_DMT_to_Study_DMT" ~ "Treatment Duration" 
                     , Variable == "Followup_Duration" ~ "Follow-up Duration" 
                     , Variable == "Disease_Duration" ~ "Disease Duration" 
                     , Variable == "Utilization_Total" ~ "Utilization" 
                     , Variable == "N_MS_CUIs" ~ "MS CUI Frequency" 
                     , Variable == "N_MS_PheCodes" ~ "MS PheCode Frequency"
                     , Variable == "ridge_All_EHR_NLP_latent" ~ "Baseline PDDS Risk"
                     , T ~ Variable
           ) 
         )

# coef.df = tibble(Coefficient = model$OC.coef
#        , Variable = names(model$OC.coef)) %>% mutate(Model = "Outcome") %>% 
#   union_all(
#         tibble(Coefficient = model$Trt.coef
#                , Variable = names(model$Trt.coef)
#                ) %>% mutate(Model = "Treatment")
#         ) %>% 
#   filter(Variable!="(Intercept)") %>% 
#   group_by(Model) %>% 
#   mutate(Coefficient.Std = Coefficient / max(abs(Coefficient)))


## All Combined #### 

coef.df %>% 
  mutate(Coefficient.Std = Coefficient / max(abs(Coefficient))) %>% 
  group_by(Variable) %>%
  filter(max(abs(Coefficient))>0) %>% 
  ungroup() %>% 
  arrange(Variable) %>% 
  ggplot(aes(x = Model, reorder(Variable, Coefficient, decreasing=F)
             # , color = Coefficient.Std
             , color = Coefficient
             , size = abs(Coefficient)
             )
         ) + 
  labs(color="Coefficient") +
  # labs(color="Coefficient\n(max-normalized)") +
  geom_point() +
  scale_color_distiller(type="div", palette = "RdBu"
                        # , limits = c(-1, 1)
                        ) + 
  theme_classic() + facet_wrap(~Year)


## Treat/Outcome Separate ####
  oc.plt = coef.df %>% 
    filter(Model == "Outcome" & Variable!="(Intercept)") %>% 
    mutate(Coefficient.Std = Coefficient / max(abs(Coefficient))) %>% 
    rename(`Eval. Window` = Year) %>% 
    group_by(Variable) %>%
    filter(max(abs(Coefficient))>0) %>% 
    ungroup() %>% 
    arrange(Variable) %>% 
    ggplot(aes(x = `Eval. Window`, y = reorder(Variable, Coefficient, decreasing=F)
               # , color = Coefficient.Std
               , color = Coefficient
               , size = abs(Coefficient)
    )
    ) + 
    labs(color="Coefficient") +
    # labs(color="Coefficient\n(max-normalized)") +
    geom_point() +
    scale_color_distiller(type="div", palette = "RdBu"
                          # , limits = c(-1, 1)
    ) + theme_classic() + ylab("") + 
  scale_size(range = c(0, 5))


  trt.plt = coef.df %>% 
    mutate(Variable = case_when(Variable == "Days_Earliest_DMT_to_Study_DMT" ~ "Prior_DMT_Use_Length"
                                , T ~ Variable)) %>% 
    filter(Model == "Treatment" & Variable!="(Intercept)") %>% 
    mutate(Coefficient.Std = Coefficient / max(abs(Coefficient))) %>% 
    rename(`Eval. Window` = Year) %>% 
    group_by(Variable) %>%
    filter(max(abs(Coefficient))>0) %>% 
    ungroup() %>% 
    arrange(Variable) %>% 
    ggplot(aes(x = `Eval. Window`, y = reorder(Variable, Coefficient, decreasing=F)
               # , color = Coefficient.Std
               , color = Coefficient
               , size = abs(Coefficient)
    )
    ) + 
    labs(color="Coefficient") +
    # labs(color="Coefficient\n(max-normalized)") +
    geom_point() +
    scale_color_distiller(type="div", palette = "RdBu"
                          # , limits = c(-1, 1)
    ) + theme_classic() + ylab("") + 
    scale_size(range = c(0, 5), guide = "none")
  
  
  plts.combined = plot_grid(plot_grid(oc.plt + theme(legend.position="none")
                                      , trt.plt + theme(legend.position="none") 
                                      , ncol=2
                                      , labels = c("Outcome", "Treatment")
                                      # , labels = c("A", "B")
                                      # , label_y = 1.02
                                      # , vjust = 1
                                      )
            , get_legend(trt.plt + theme(legend.position="bottom"
                                         , legend.key.width=unit(1.5, "cm")))
            , nrow=2
            , rel_heights = c(1, 0.1)
            )
            
  
  
  trt.plt.solo = coef.df %>% 
    mutate(Variable = case_when(Variable == "Days_Earliest_DMT_to_Study_DMT" ~ "Prior_DMT_Use_Length"
                                , T ~ Variable)) %>% 
    filter(Model == "Treatment" & Variable!="(Intercept)") %>% 
    mutate(Coefficient.Std = Coefficient / max(abs(Coefficient))
           , Model = "") %>% 
    group_by(Variable) %>%
    filter(max(abs(Coefficient))>0) %>% 
    ungroup() %>% 
    arrange(Variable) %>% 
    ggplot(aes(x = Model, y = reorder(Variable, Coefficient, decreasing=F)
               # , color = Coefficient.Std
               , color = Coefficient
               , size = abs(Coefficient)
    )
    ) + 
    labs(color="Coefficient") +
    # labs(color="Coefficient\n(max-normalized)") +
    geom_point() +
    scale_color_distiller(type="div", palette = "RdBu"
                          # , limits = c(-1, 1)
    ) + theme_classic() + ylab("") + 
    scale_size(range = c(0, 5), guide = "none") + 
    xlab("Treatment Model")
  
  
  return(list(oc.plt = oc.plt
              , trt.plt = trt.plt 
              , trt.plt.solo = trt.plt.solo
              , plts.combined = plts.combined
              , coef.df = coef.df))
  
}

model.2yr.SD = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Decrease_2yr_BCD_NTZ_ATE_BCD_NTZ.RDS"))
model.3yr.SD = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Decrease_3yr_BCD_NTZ_ATE_BCD_NTZ.RDS"))
.outcome = "SD"
SD.plts = coef.plot(model.2yr = model.2yr.SD, model.3yr = model.3yr.SD)

model.2yr.SI = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Increase_2yr_BCD_NTZ_ATE_BCD_NTZ.RDS"))
model.3yr.SI = readRDS(here("Results", "BCD_NTZ", "ATE", "ATE_Increase_3yr_BCD_NTZ_ATE_BCD_NTZ.RDS"))
.outcome = "SI"
SI.plts = coef.plot(model.2yr = model.2yr.SI, model.3yr = model.3yr.SI)

SD.plts$oc.plt
SI.plts$oc.plt
SI.plts$trt.plt.solo
  


plts.combined = plot_grid(plot_grid(get_legend(SI.plts$oc.plt + theme(legend.position = "none")))
                          , plot_grid(SI.plts$oc.plt + theme(legend.position="none")
                                    , SD.plts$oc.plt + theme(legend.position="none")
                                    , SI.plts$trt.plt.solo + theme(legend.position="none")
                                    , ncol=3
                                    , labels = c("SI", "SD", "Treatment")
                                    # , labels = c("A", "B")
                                    # , label_y = 1.02
                                    , vjust = 0
                            )
                            , get_legend(SI.plts$oc.plt + scale_size_continuous(guide="none") + theme(legend.position="bottom"
                                                         , legend.key.width=unit(2, "cm")))
                            , nrow=3
                            , rel_heights = c(0.1, 1, 0.1)
                            )


if(.save)  cowplot::save_plot(plot = plts.combined
                     , filename = here("Results", "BCD_NTZ", paste0("Fig4_DiPS_Coefficients.png"))
                     , base_height = 10
                     , base_asp = 1)



# common coefs  ####
  ## SI #### 
  SI.plts$coef.df %>% filter(Variable %in% (SI.plts$coef.df %>% 
                                      filter(abs(Coefficient)>0 & Year=="3-years") %>% 
                                      count(Variable) %>% 
                                      filter(n>1) %>% 
                                      pull(Variable)
                                      )
                             & Year=="3-years"
                     ) %>% 
  arrange(Variable) %>%  
  left_join(ft.dict, by="Variable", relationship="many-to-many") %>% 
  print(n=100)

  ## SD ####
  SD.plts$coef.df %>% filter(Variable %in% (SD.plts$coef.df %>% 
                                              filter(abs(Coefficient)>0 & Year=="3-years") %>% 
                                              count(Variable) %>% 
                                              filter(n>1) %>% 
                                              pull(Variable)
  )
  & Year=="3-years"
  ) %>% 
    arrange(Variable) %>% 
  left_join(ft.dict, by="Variable", relationship="many-to-many") %>% 
  print(n=100)


# data checks ####

all.equal(model.2yr.SD$Outcome_ImputeModel_Coefs$Trt
          , model.3yr.SD$Outcome_ImputeModel_Coefs$Trt)


all.equal(model.2yr.SI$Outcome_ImputeModel_Coefs$Trt
          , model.3yr.SI$Outcome_ImputeModel_Coefs$Trt)

all.equal(model.2yr.SI$Outcome_ImputeModel_Coefs$Trt
          , model.2yr.SD$Outcome_ImputeModel_Coefs$Trt)

# varname check 
model.2yr.SI$Outcome_ImputeModel_Coefs$Trt %>% 
  tibble(Variable = names(.), Coef=.) %>% 
  filter(abs(Coef)>0) %>% 
  filter(substr(Variable, 1, 1)!="C") %>% 
  pull(Variable) %>% unique()


tibble(Variable = names(model.2yr.SI$Outcome_ImputeModel_Coefs$OC)
       , OC2yr = model.2yr.SI$Outcome_ImputeModel_Coefs$OC
       , OC3yr = model.3yr.SI$Outcome_ImputeModel_Coefs$OC
       ) %>% 
  filter(abs(OC2yr)>0 | abs(OC3yr)>0) %>% 
  filter(substr(Variable, 1, 1)!="C"
         ) %>%
  pull(Variable) %>% unique() 


