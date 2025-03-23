
pdds_impute = function(.model, .analytic_object, .print=F, .pred_class=F){  
  
  if(all(class(.model)=="ordinalNetTune")){
    L = which.min(apply(.model$misclass, 1, mean))
    coefs = .model$fit$coefs[L,]
    .mod_names = names(coefs)
  }else{
    coefs = .model$coefficients
    .mod_names = names(coefs)
  }
  
  .onames = c("White_NonHispanic=1", "subject_sex_ch=TRUE", "Age_at_PDDS", "Days_FirstPheCode_to_PDDS", "Days_First_MS_PheCode_to_PDDS")
  .nnames = c("White_NonHispanic", "subject_sex_ch", "Age_at_DMT_Initiation", "Followup_Duration", "Disease_Duration")
  
  indx <- match(.onames,.mod_names)
  .mod_names[indx[!is.na(indx)]] <- .nnames[which(!is.na(indx))]
  
  .missing_vars = .mod_names[.mod_names %nin% colnames(.analytic_object$X)][-(1:7)]
  
  for(var in .missing_vars){
    if(exists(var, where=.analytic_object$X)) stop("variable found in dataframe")
    .analytic_object$X[,var] = 0   
  }
  # check 
  # if(sum(colnames(.analytic_object$X) %in% .mod_names)  != length(.mod_names)-7) stop("Not all columns correctly added")
  
  if(.print){
    print(dim(.analytic_object$X[,.mod_names[-(1:7)]]))
    print(length(.mod_names))
  }
  
  .predX = .analytic_object$X[,.mod_names[-(1:7)]] %>% as.matrix()
  
  if(all(class(.model) == "ordinalNetTune")){
    .type = ifelse(.pred_class, "class", "link")
    pred_vec = predict(object=.model$fit, .predX, whichLambda=L, type=.type)
    if(!.pred_class) pred_vec = pred_vec[,1] else pred_vec=pred_vec-1
    # equivalently for latent prediction 
    # pred_vec  = drop(cbind(rep(1, nrow(.predX)), .predX) %*% coefs[-(2:7)])
  }else{
    
    if(.pred_class){
      pred_vec = factor(apply(predict(.model, newdata = .predX, type = "fitted.ind"), 1, which.max)-1)
    }else{
      pred_vec = predict(object=.model, .predX, type="lp")
      # pred_vec2 = predictrms(fit=.model, newdata = .predX, type="lp")
    }
    
  }
  
  if(.print){
    if(.pred_class) print(table(pred_vec)) else print(summary(pred_vec))
  }      
  return(pred_vec)
}



# Copied from Weijing's code 02.22.2023
expit = function(x)
{
  return(1/(1+exp(-x)))
}

# A simple adaptive lasso
adalasso = function(..., lambda.min, .return_model=F)
{
  ridgefit = cv.glmnet(alpha = 0, ...)
  adawgt = 1/abs(coef(ridgefit, s = "lambda.min")[-1])
  if(missing(lambda.min))
  {
    adafit = cv.glmnet(penalty.factor = adawgt, ...)
    lambda.min = adafit$lambda.min
  }else
  {
    adafit = glmnet(penalty.factor = adawgt, ...)
  }
  if(.return_model){
    return(adafit)
  }else{
    return(drop(coef(adafit, s = lambda.min)))
  }
}


ordinalCoefs= function(...
                       , .link_fn="logit", .alpha=1, .nlambda=20, .nFolds=5
                       ){
  
  # ordinalNet fitting 
  .oNet <- ordinalNet(#... # x,y sufficient
                                    # ,
                                   x, factor(y)
                                   , alpha = .alpha
                                   , nLambda = .nlambda
                                   , family = "cumulative"
                                   , link = .link_fn
                                   # , nFolds = .nFolds # needed for ordinalNetTune
                                   )
  
  # coefficient extraction 
  .oNet_beta = coef(.oNet
                    # , whichLambda = which.min(apply(.oNet$misclass, 1, mean))
                    )
  
  return(.oNet_beta)
}


PS_kernel=function(Xi, Ti, Yi, .family, .outcome, .nfolds=10){
  if(is.data.frame(Xi)) Xi = as.matrix(Xi) else if(!is.matrix(Xi)) stop("Xi not matrix or dataframe cannot continue")
  
  # treatment ALASSO (Cheng eqn 9)
  Trt.coef = adalasso(x = Xi
                      , y = Ti
                      , family = .family 
                      , nfolds=.nfolds)

  
  # outcome ALASSO (Cheng eqn 10) 
  if(.outcome=="Y_diff"){
    OC.coef = ordinalCoefs(x = cbind(Xi[!is.na(Yi),], Trt=Ti[!is.na(Yi)])
                           , y = Yi[is.na(Yi)])
  }else{
    OC.coef = adalasso(x = cbind(Xi[!is.na(Yi),]
                                 , Trt=Ti[!is.na(Yi)]
                                 )
                       , y = Yi[!is.na(Yi)]
                       , family = .family
                       , nfolds = .nfolds)
  }
  # Kernel estimator based on coefs from 9/10 above 
  a = cbind(rep(1,nrow(Xi)), Xi) %*% t(t(Trt.coef))
  
  OC.coef.noTrt = OC.coef[-which(names(OC.coef)=="Trt")]
  b = cbind(rep(1,nrow(Xi)), Xi) %*% t(t(OC.coef.noTrt))
  S = cbind(a, b)
  
  S_preproc = apply(S, 2, function(x) (x-mean(x))/sd(x) )
  S_cdf = apply(S_preproc, 2, pnorm)
  
  .bw_preproc = np::npregbw(xdat = S_cdf, ydat=Ti, ckerorder=4, regtype="lc"
                            , nmulti=5)

  np_preproc_pred = predict(npreg(.bw_preproc))

  return(list(PS.Kernel = np_preproc_pred
              , KernelModel = .bw_preproc))
}

# Robust ATE (piecing together all other materials)
  AIPW_robust = function(Ti, Yi, PS, .IPW_stabilize, min.ipw=0.1, max.ipw=10, z_star){
    
    if(any(min(PS)<=0, max(PS)>=1) & !.IPW_stabilize ) warning("IPW values beyond (0,1) detected but stabilizing turned off!!!")
    IPW = Ti/PS + (1-Ti)/(1-PS)
    
    if(.IPW_stabilize){
      IPW = pmin(pmax(IPW, 0.1), 10)
    }
    
    mu_1 = sum(Ti*Yi*IPW) / sum(Ti*IPW)
    mu_0 = sum((1-Ti)*Yi*IPW) / sum((1-Ti)*IPW)
    
    return(mu_1 - mu_0)
  }


  
  
  
# attempt at catch-all function 
  
AIPW_robust_Total = function(.analytic_df
                             , .outcome # c("Y_change", "Y_AnyIncrease", "Y_diff")
                             , ...){
  
  analytic_object <- analytic_wrangling(.analytic_df = .analytic_df
                                        , ...
                                        # , pre_process = T
                                        # , efficacy="high"
                                        # , .EHR_trim = 0.90 # also NLP trim
                                        # , .NLP = T
                                        # , .NLP_all = F
                                        # , .EHR_all = F
                                        # , .RXNORM_Exclude_All = T
                                        # , .DMT_reference = DMT_reference
                                        # , .diff_window = dmonths(6)
                                        )

  X_1yr = analytic_wrangling_postDMT(object = .analytic_df
                                     , ehr_ft_list = analytic_object[["EHR_fts"]]
                                     , nlp_ft_list = analytic_object[["NLP_fts"]]
                                     , pre_process = T
                                     , .lag = lubridate::dyears(1) 
                                     ) %>% suppressWarnings() # my own identifier warning, taken care of in-house (in-function, 197ish)

  # setdiff(colnames(X_1yr), colnames(analytic_object$X)) # shoudl only diff in PATIENT_NUM, DMT_Study_Start in 1yr matrix 
  # setdiff(colnames(analytic_object$X), colnames(X_1yr)) # no diff
          
  
  ## Clean-up for latent scoring #### 
  lasso_logit_latent = pdds_impute(.model = lasso_logit, .analytic_object = analytic_object
                                   , .print=F, .pred_class=F)
  lasso_cauchit_latent = pdds_impute(.model = lasso_cauchit, .analytic_object = analytic_object
                                     , .print=F, .pred_class=F)
  ordinal_cauchit_latent = pdds_impute(.model = ordinal_KG_NLP$Model, .analytic_object = analytic_object
                                       , .print=F, .pred_class=F)
  
  lasso_logit_class = pdds_impute(.model = lasso_logit, .analytic_object = analytic_object
                                  , .print=T, .pred_class=T)
  
  .pdds_input = X_1yr #analytic_object$input_df
  .pdds_input$PATIENT_NUM = as.character(.pdds_input$PATIENT_NUM)
  .pdds_input$id_participant = 0 # lazily lets me source same pdds code from cohor building 
  assign(".pdds_input", .pdds_input, envir=.GlobalEnv)
  
  source(here("Code", "04_pdds.R")) # just for pdds dataframe 
  rm(.pdds_input) # probably bad practice
  pdds$PATIENT_NUM = suppressWarnings(as.double(pdds$PATIENT_NUM))
  # will provide NA warning due to missingness, expected behavior 
  
  .lag = lubridate::dyears(1)
  
  pdds_1yr = X_1yr %>% 
    left_join(pdds, by="PATIENT_NUM") %>% 
    distinct(PATIENT_NUM, date, score, DMT_Study_Start) %>% 
    filter(date < DMT_Study_Start + .lag) %>% 
    group_by(PATIENT_NUM) %>% 
    dplyr::summarize(PDDS_1yr_Avg = mean(score)) %>% 
    left_join(x = X_1yr[,"PATIENT_NUM"], y = ., "PATIENT_NUM") %>% 
    bind_cols(pred = lasso_logit_class) %>% 
    mutate(PDDS_1yr = coalesce(PDDS_1yr_Avg, pred)) %>% 
    # select(PATIENT_NUM, PDDS_1yr)
    pull(PDDS_1yr)
  
  # quick checks for ability to lazily cbind over merging 
  # all.equal(pdds_1yr$PATIENT_NUM, X_1yr$PATIENT_NUM) # T
  # all.equal(pdds_1yr$PATIENT_NUM, analytic_object$input_df$PATIENT_NUM) # T
  
  # Creating full design matrix for outcome imputation 
  # source(here("Code", "analysis_scripts", "00_causal_core_Robust.R"))
  
  ## Utility Weight Calculation ####
  # tmpX = base::subset(X_1yr, select = -c(PATIENT_NUM, DMT_Study_Start))
  tmpX = cbind(lasso_logit_latent, lasso_cauchit_latent, ordinal_cauchit_latent
               , analytic_object$X
               )
  
  .alpha = adalasso(x = as.matrix(tmpX)
                    , y = analytic_object$A
                    , family="binomial"
                    , .return_model = F)
  
  U_X = as.matrix(cbind(rep(1, nrow(tmpX)), tmpX))
  .pi = pmax(pmin(expit(drop(U_X %*% t(t(.alpha)))), 0.999), 0.001)
  .inv_pi = pmax(pmin(1/.pi, 10), 0.1)
  U_pi = pmax(pmin(analytic_object$A/.pi + (1-analytic_object$A)/(1-.pi), 10), 0.1)
  # equivalent to predict(.alpha_mod, tmpX, s="lambda.min") for alpha_mod = adalasso(.return_model=T)
  
  ## Combining 
  Xfull = cbind(
    # pre-DMT latent score 
    lasso_logit_latent, lasso_cauchit_latent, ordinal_cauchit_latent
    # design matrix through 1-yer post-DMT 
    , base::subset(X_1yr, select = -c(PATIENT_NUM, DMT_Study_Start))
    # 
    , pdds_1yr 
    , analytic_object$A # T 
    , U_pi 
  )
  
  Xfull$subject_sex_ch = as.integer(Xfull$subject_sex_ch=="F") #
  
  ## Outcome imputation #### 
  
  indx = !is.na(analytic_object[[.outcome]])
  
  if(.outcome == "Y_diff"){
    Xi_Outcome_Coefs = ordinalCoefs(x = as.matrix(Xfull[indx,])
                                    , y = analytic_object[[.outcome]][indx]
                                    )
  }else{
    Xi_Outcome_Coefs = adalasso(x = as.matrix(Xfull[indx,])
                                , y = analytic_object[[.outcome]][indx]
                                , family = "binomial"
                                , nfolds = 10
                                , .return_model=T)
  }
  
  Y_pred = predict(Xi_Outcome_Coefs, newx = as.matrix(Xfull), s="lambda.min", type = "class")                  
  
  Y_dag = coalesce(analytic_object[[.outcome]], as.integer(Y_pred))
  # check, should be TRUE 
  # all.equal(Y_dag[!is.na(analytic_object[[.outcome]])], analytic_object[[.outcome]][!is.na(analytic_object[[.outcome]])])
  
  
  ## Kernel Smoothing #### 
  
  # print(analytic_object[[.outcome]])
  propens = PS_kernel(Xi = as.matrix(analytic_object$X)
                      , Ti = analytic_object$A
                      , Yi = analytic_object[[.outcome]]
                      , .outcome = .outcome
                      , .family = "binomial"
                      , .nfolds=10
  )
  
  
  
  ## ATE Estimation ####
  ATE_robust = function(Ti, Yi, PS, .IPW_stabilize, min.ipw=0.1, max.ipw=10, z_star){
    
    if(any(min(PS)<=0, max(PS)>=1) & !.IPW_stabilize ) warning("IPW values beyond (0,1) detected but stabilizing turned off!!!")
    IPW = Ti/PS + (1-Ti)/(1-PS)
    
    if(.IPW_stabilize){
      IPW = pmin(pmax(IPW, min.ipw), max.ipw)
      IPW[is.na(IPW) & PS==1 & Ti==1] = 1
      IPW[is.na(IPW) & PS==0 & Ti==0] = 1
    }
    
    mu_1 = sum(Ti*Yi*IPW) / sum(Ti*IPW)
    mu_0 = sum((1-Ti)*Yi*IPW) / sum((1-Ti)*IPW)
    
    return(mu_1 - mu_0)
  }
  
  bcd_ntz_ate = ATE_robust(Ti = analytic_object$A
                            , Yi = Y_dag
                            , PS = propens$PS.Kernel
                            , min.ipw = 0.01 
                            , max.ipw = 10
                            , .IPW_stabilize = T)
  
  # 
  # B = 500
  # .ate = c()
  # for(iter in 1:500){
  #   indx = sample(1:nrow(analytic_object$X), size=nrow(analytic_object$X)*0.8, replace=T)
  #   .ate = c(.ate, ATE_robust(Ti = analytic_object$A[indx]
  #                              , Yi = Y_dag[indx]
  #                              , PS = propens$PS.Kernel[indx]
  #                              , min.ipw = 0.01
  #                              , max.ipw = 10
  #                              , .IPW_stabilize=T)
  #   )
  # }

  
  return(list(ATE = bcd_ntz_ate
              # , CI = quantile(.ate, c(0.025, 0.975))
              # , .boot = .ate
              # , EHR_fts = analytic_object$EHR_fts
              # , NLP_fts = analytic_object$NLP_fts
              , analytic_object = analytic_object 
              , X_1yr = X_1yr 
              )
         )
}



# Function to use in boot::boot bootstrap implementaion
boot_AIPW_robust_Total <- function(.analytic_df, indx, .outcome, ...){
  .ate = AIPW_robust_Total(.analytic_df = .analytic_df[indx,]
                            , .outcome = .outcome
                            , ...)
  
  return(.ate$ATE )
}



