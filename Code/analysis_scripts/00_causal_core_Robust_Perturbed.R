
gaus4 = function(x) 0.5*(3-x^2)*dnorm(x)

################## #
# PDDS Imputation 
################## #

pdds_impute = function(.model, .analytic_X, .print=F, type="response"){  
    
  ### Valid types (equivalent to glmnet args plus my custom expectation)
  # link; response (probs); class
  L = which.max(apply(.model$loglik, 1, mean))
  coefs = .model$fit$coefs[L,]
  .mod_names = names(coefs)

  .onames = c("White_NonHispanic=1", "subject_sex_ch=TRUE", "Age_at_PDDS", "Days_FirstPheCode_to_PDDS", "Days_First_MS_PheCode_to_PDDS")
  .nnames = c("White_NonHispanic", "subject_sex_ch", "Age_at_DMT_Initiation", "Followup_Duration", "Disease_Duration")
  
  indx <- match(.onames,.mod_names)
  .mod_names[indx[!is.na(indx)]] <- .nnames[which(!is.na(indx))]
  
  .missing_vars = .mod_names[.mod_names %nin% colnames(.analytic_X)][-(1:7)]
  
  for(var in .missing_vars){
    if(exists(var, where=.analytic_X)) stop("variable found in dataframe")
    .analytic_X[,var] = 0   
  }
  # check 
  # if(sum(colnames(.analytic_X) %in% .mod_names)  != length(.mod_names)-7) stop("Not all columns correctly added")
  
  if(.print){
    print(dim(.analytic_X[,.mod_names[-(1:7)]]))
    print(length(.mod_names))
  }
  
  .predX = .analytic_X[,.mod_names[-(1:7)]] %>% as.matrix()

  pred_vec = predict(object=.model$fit, .predX, whichLambda=L, type=type)
  
  if(type=="class"){
    pred_vec_cln = pred_vec - 1
  }else if (type=="link"){
    pred_vec_cln = pred_vec[,1] # choosing first intercept arbitrarily
    # equivalently for latent prediction 
    # pred_vec  = drop(cbind(rep(1, nrow(.predX)), .predX) %*% coefs[-(2:7)]) 
  }else if (type %in% c("response", "expectation", "E")){
    pred_vec_cln = pred_vec %*% 0:7 
  }

  if(.print){
    if(type=="class") print(table(pred_vec_cln)) else print(summary(pred_vec_cln))
  }      
  return(pred_vec_cln)
}



############################## # 
# expit and adalasso 
############################## #

# Copied from Weijing's code 02.22.2023
expit = function(x)
{
  return(1/(1+exp(-x)))
}

# A simple adaptive lasso
adalasso = function(..., lambda.min, .return_model=F, ridge.penalty, .penalize_Trt = T)
{
  if(!missing(ridge.penalty)){
    ridgefit = cv.glmnet(alpha = 0, ..., penalty.factor = ridge.penalty) 
  }
  else{
    ridgefit = cv.glmnet(alpha = 0, ...) 
  }
  
  adawgt = 1/abs(coef(ridgefit, s = "lambda.min"))
  if(!.penalize_Trt) adawgt[rownames(coef(ridgefit))=="Trt"] = 0
  adawgt = adawgt[-1]
  
  if(missing(lambda.min))
  {
    adafit = cv.glmnet(penalty.factor = adawgt, ...)
    lambda.min = adafit$lambda.min
  }else
  {
    adafit = glmnet(penalty.factor = adawgt, ...)
  }
  # print(length(adafit$lambda))
  if(.return_model){
    return(adafit)
  }else{
    return(drop(coef(adafit, s = lambda.min)))
  }
}


############################## #
# ordinal model coefficients
############################## #

ordinalCoefs= function(...
                       , .link_fn="logit", .alpha=1, .nlambda=20, .nFolds=5
                       ){
  
  # ordinalNet fitting 
  .oNet <- ordinalNet(#... # x,y sufficient
                                    # ,
                                   x
                                   , factor(y)
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


############################## #
# kernel smoothed pi 
############################## #

PS_kernel=function(Xi, Ti, Yi, .family, .outcome, .nfolds=5, weights=NULL){
  
  if(is.data.frame(Xi)) Xi = as.matrix(Xi) else if(!is.matrix(Xi)) stop("Xi not matrix or dataframe cannot continue")
  
  if(is.null(weights)) weights = rep(1, length(Ti))
  # treatment ALASSO (Cheng eqn 9)
  Trt.coef = adalasso(x = Xi
                      , y = Ti
                      , family = "binomial" 
                      , nfolds=.nfolds
                      , weight = weights
                      )

  
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
                       , nfolds = .nfolds
                       , weights = weights[!is.na(Yi)]
                       , ridge.penalty = c(rep(1, ncol(Xi)), 0)
                       , .penalize_Trt = F
                       # , .return_model = T
                       )
    
  }
  # Kernel estimator based on coefs from 9/10 above 
  a = cbind(rep(1,nrow(Xi)), Xi) %*% t(t(Trt.coef))
  
  OC.coef.noTrt = OC.coef[-which(names(OC.coef)=="Trt")]
  b = cbind(rep(1,nrow(Xi)), Xi) %*% t(t(OC.coef.noTrt))
  S = cbind(a, b)
  
  S_preproc = apply(S, 2, function(x) (x-mean(x))/ifelse(sd(x)==0, 1, sd(x)) )
  S_cdf = apply(S_preproc, 2, pnorm)
  
  plugin = ifelse(sd(S_cdf[,1])*nrow(Xi)^(-1/5)!=0
                  , sd(S_cdf[,1])*nrow(Xi)^(-1/5)
                  , sd(S_cdf[,2])*nrow(Xi)^(-1/5))
  # plugin = c(sd(S_cdf[,1])*nrow(Xi)^(-1/5)
             # , ifelse(sd(S_cdf[,2])*sum(!is.na(Yi))^(-1/5)==0, 1, sd(S_cdf[,2])*sum(!is.na(Yi))^(-1/5))
             # )
  
  np_preproc_pred = c()
  for(i in 1:nrow(S_cdf)){
    num1 = gaus4((-S_cdf[,1] + S_cdf[i,1]) / plugin)
    num2 = gaus4((-S_cdf[,2] + S_cdf[i,2]) / plugin)
    np_preproc_pred = c(np_preproc_pred, (t(num1*num2) %*% (Ti * weights)) / t(num1 * num2) %*% weights)
  }

  # np.model = np::npreg(xdat = S_cdf, ydat = Ti * weights, regtype="lc", ckerorder=4
  #                     , bws = plugin
  #                     )

  # np_preproc_pred = predict(np.model)
  # all.equal(predict(np.model), np_preproc_pred)
  
  return(list(PS.Kernel = np_preproc_pred
              , Trt_IntOnly = as.logical(sd(S[1,])==0)
              , Outcome_IntOnly = as.logical(sd(S[2,])==0)
              , Trt.coef = Trt.coef
              , OC.coef = OC.coef
              # , KernelModel = np.model
              )
         )
         
}

############## # 
# Robust ATE 
############## #

ATE_DiPS = function(Ti, Yi, PS, .IPW_stabilize, min.ipw=0.1, max.ipw=10
                    , z_star
                    , weights = NULL
                    , crump.trim = T
                    , labels){
  
    if(!is.null(weights)) {
      if(length(weights)!=length(Ti)) stop("Weights supplied but incorrect length, should match Ti vector length")
      # Ti = Ti * weights
    }
  
    trimmed = F 
    .labeled.removed = .unlabeled.removed = 0

    # if(crump.trim){
    #   .labeled.removed = sum((PS<0.05 | PS>0.95) & !is.na(labels) )
    #   .unlabeled.removed = sum((PS<0.05 | PS>0.95) & is.na(labels) )
    # if(.labeled.removed + .unlabeled.removed <= 0.4*length(labels)){
    #     trimmed = T
    #     .keep = which(PS>=0.1 & PS<=0.9)
    #     PS = PS[.keep]
    #     weights = weights[.keep]
    #     Ti = Ti[.keep]
    #     Yi = Yi[.keep]
    #   }
    # }

    if(crump.trim){
      .labeled.removed = sum((PS<0.1 | PS>0.9) & !is.na(labels) )
      .unlabeled.removed = sum((PS<0.1 | PS>0.9) & is.na(labels) )
      trimmed = T
      .keep = which(PS>=0.1 & PS<=0.9)
      PS = PS[.keep]
      weights = weights[.keep]
      Ti = Ti[.keep]
      Yi = Yi[.keep]
      }
    
    # if(any(min(PS)<=0, max(PS)>=1) & !.IPW_stabilize ) warning("IPW values beyond (0,1) detected but stabilizing turned off!!!")
    IPW = Ti/PS + (1-Ti)/(1-PS)
    
    # mu_1 = sum(Ti*Yi*IPW) / sum(Ti*IPW)
    # mu_0 = sum((1-Ti)*Yi*IPW) / sum((1-Ti)*IPW)
    mu_1 = sum(weights*Ti*Yi / PS) / sum(Ti*weights / PS)
    mu_0 = sum(weights * (1-Ti)* Yi / (1-PS)) / sum((1-Ti) * weights / (1-PS))
    
    .output = list(.ate = mu_1 - mu_0
                   , .labeled.removed = .labeled.removed
                   , .unlabeled.removed = .unlabeled.removed
                   , trimmed = trimmed)
    return(.output)
    
  }


############### #
# Full calculation
############### # 
AIPW_robust_analytic = function(analytic_object
                                , .analytic_object_1yr # misnomer, can take any lookforward, I just didn't bother to rename and ensure correctness
                                , .outcome # c("Y_change", "Y_NegChange", "Y_AnyIncrease", "Y_diff")
                                , pweights = NULL
                                , .lookforward_window = NULL 
                                , cv.outcome = T
                                ){
    
  if(is.null(.lookforward_window)) {
    warning("Setting lookforward to 1 year by default")
    .lookforward_window = dyears(1)
    }
  ## Kernel Smoothing #### 
  # ridge_12mo_EHR_NLP_latent = pdds_impute(.model = ridge_12mo_EHR_NLP, .analytic_X = analytic_object$X
  #                                         , .print=F, type="response")
  ridge_All_EHR_NLP_latent = pdds_impute(.model = ridge_All_EHR_NLP, .analytic_X = analytic_object$X
                                         , .print=F, type="response")
  
  X_inclPDDSlatent = cbind(analytic_object$X#, ridge_12mo_EHR_NLP_latent
                           , ridge_All_EHR_NLP_latent)
  # print(analytic_object[[.outcome]])
  # npseed(1) # for debugging
  .family = ifelse(.outcome %in% c("Y_change", "Y_NegChange"), "binomial", "gaussian")
  
  propens = PS_kernel(Xi = X_inclPDDSlatent # as.matrix(analytic_object$X)
                      , Ti = analytic_object$A
                      , Yi = analytic_object[[.outcome]]
                      , .outcome = .outcome
                      , .family = .family
                      , .nfolds=5
                      , weights = pweights
                      )
  
  
  ## Utility Covariate (truncating/cleaning)
  .pi = propens$PS.Kernel
  .pi_trim = pmax(pmin(.pi, 0.999), 0.001) # Stability truncation
  
  message(paste0("0s in PS = ", sum(.pi_trim==0)))
  U_pi = analytic_object$A/.pi_trim + (1-analytic_object$A)/(1-.pi_trim)

  
  ## Clean-up for latent scoring #### 
  ridge_12mo_EHR_NLP_lookforward_expect = pdds_impute(.model = ridge_12mo_EHR_NLP, .analytic_X = .analytic_object_1yr
                                            , .print=F, type="response")

  
  .pdds_input = .analytic_object_1yr #analytic_object$input_df
  .pdds_input$PATIENT_NUM = as.character(.pdds_input$PATIENT_NUM)
  .pdds_input$id_participant = 0 # lazily lets me source same pdds code from cohort building 
  assign(".pdds_input", .pdds_input, envir=.GlobalEnv)
  
  source(here("Code", "04_pdds.R")) # just for pdds dataframe 
  rm(.pdds_input) # probably bad practice
  pdds$PATIENT_NUM = suppressWarnings(as.double(pdds$PATIENT_NUM))
  # will provide NA warning due to missingness, expected behavior 

  
  pdds_lookforward = .analytic_object_1yr %>% 
    left_join(pdds, by="PATIENT_NUM") %>% 
    distinct(PATIENT_NUM, date, score, DMT_Study_Start) %>% 
    filter(date < DMT_Study_Start + .lookforward_window
           & date >= DMT_Study_Start) %>% 
    group_by(PATIENT_NUM) %>% 
    dplyr::summarize(PDDS_Lookforward_Avg = mean(score)) %>% 
    left_join(x = .analytic_object_1yr[,"PATIENT_NUM"], y = ., "PATIENT_NUM") %>% 
    bind_cols(pred = ridge_12mo_EHR_NLP_lookforward_expect) %>% 
    mutate(PDDS_Lookforward = coalesce(PDDS_Lookforward_Avg, pred)) %>% 
    # select(PATIENT_NUM, PDDS_1yr)
    pull(PDDS_Lookforward)
  
  # quick checks for ability to lazily cbind over merging 
  # all.equal(pdds_1yr$PATIENT_NUM, .analytic_object_1yr$PATIENT_NUM) # T
  # all.equal(pdds_1yr$PATIENT_NUM, analytic_object$input_df$PATIENT_NUM) # T
  
  # Creating full design matrix for outcome imputation 
  # source(here("Code", "analysis_scripts", "00_causal_core_Robust.R"))

  ## Combining 
  Xfull = cbind(
    # pre-DMT latent score 
    ridge_All_EHR_NLP_latent#, ridge_12mo_EHR_NLP_latent
    # design matrix through 1-yer post-DMT 
    , base::subset(.analytic_object_1yr, select = -c(PATIENT_NUM, DMT_Study_Start))
    , pdds_lookforward
    , analytic_object$A # T 
    , U_pi 
    )
  
  Xfull$subject_sex_ch = as.integer(Xfull$subject_sex_ch=="F") #
  
  ## Outcome imputation #### 
  if(cv.outcome){
    outcome_indx = which(!is.na(analytic_object[[.outcome]]))
    # simple sampling 
    cv_outcome_indx = sample(outcome_indx, replace=F, size = 0.8*length(outcome_indx))
    test_outcome_indx = setdiff(outcome_indx, cv_outcome_indx)
  # }else{
  #   outcome_indx = which(!is.na(analytic_object[[.outcome]]))
  #   cv_outcome_indx = outcome_indx
  # } 
    
    if(.outcome == "Y_diff"){
      # !!!! CV code not written for ordinal outcome!!!!
      Xi_Outcome_Coefs = ordinalCoefs(x = as.matrix(Xfull[outcome_indx,])
                                      , y = analytic_object[[.outcome]][outcome_indx]
      )
    }else{
      # adalasso 
      Xi_Outcome_Coefs = adalasso(x = as.matrix(Xfull[cv_outcome_indx,])
                                  , y = analytic_object[[.outcome]][cv_outcome_indx]
                                  , family = .family
                                  , type.measure = ifelse(.family=="binomial", "auc", "mse") # "auc"
                                  , nfolds = 5
                                  , .return_model=T
                                  , weights = pweights[cv_outcome_indx])
  
    }
    

    Y_test_class = predict(Xi_Outcome_Coefs, newx = as.matrix(Xfull[test_outcome_indx,]), s="lambda.min", type = "class")
    Y_test_prob = predict(Xi_Outcome_Coefs, newx = as.matrix(Xfull[test_outcome_indx,]), s="lambda.min", type = "response")
    test_acc = sum(Y_test_class==analytic_object[[.outcome]][test_outcome_indx]) / length(test_outcome_indx)
    test_auc = survival::concordancefit(y = analytic_object[[.outcome]][test_outcome_indx]
                                        , x = Y_test_prob)$concordance
    
    # rf 
    # Y_test_class <- predict(pdds_rf_binary, newdata=rfdf[test_outcome_indx,]) 
    # Y_test_prob = predict(pdds_rf_binary, newdata=rfdf[test_outcome_indx,], type = "prob")[,2]
    # test_acc = sum(Y_test_class==analytic_object[[.outcome]][test_outcome_indx]) / length(test_outcome_indx)
    # test_auc = survival::concordancefit(y = analytic_object[[.outcome]][test_outcome_indx]
    #                                     , x = Y_test_class)$concordance    
    
  }else{
    test_acc = test_auc = NULL
  }
  
  # Fitting models on full data 
    if(.outcome == "Y_diff"){
      # !!!! CV code not written for ordinal outcome!!!!
      Xi_Outcome_Coefs = ordinalCoefs(x = as.matrix(Xfull[which(!is.na(analytic_object[[.outcome]])),])
                                      , y = analytic_object[[.outcome]][which(!is.na(analytic_object[[.outcome]]))]
      )
    }else{
      # adalasso 
      Xi_Outcome_Coefs = adalasso(x = as.matrix(Xfull[which(!is.na(analytic_object[[.outcome]])),])
                                  , y = analytic_object[[.outcome]][which(!is.na(analytic_object[[.outcome]]))]
                                  , family = .family
                                  , type.measure = ifelse(.family=="binomial", "auc", "mse") # "auc"
                                  , nfolds = 5
                                  , .return_model=T
                                  , weights = pweights[which(!is.na(analytic_object[[.outcome]]))])
      
    }
    
  Y_pred = predict(Xi_Outcome_Coefs, newx = as.matrix(Xfull), s="lambda.min", type = "class")                  
  Y_dag = coalesce(analytic_object[[.outcome]], as.integer(Y_pred))
  # check, should be TRUE 
  # all.equal(Y_dag[!is.na(analytic_object[[.outcome]])], analytic_object[[.outcome]][!is.na(analytic_object[[.outcome]])])
  

  ## ATE Estimation ####
  
  .dips.output = ATE_DiPS(Ti = analytic_object$A
                         , Yi = Y_dag
                         , PS = .pi_trim
                         , min.ipw = 0.01 
                         , max.ipw = 10
                         , .IPW_stabilize = T
                         , weights = pweights
                         , labels = analytic_object[[.outcome]]
                         , crump.trim = T
                         )
  
  
  bcd_ntz_ate = .dips.output$.ate 
  
  
  if(cv.outcome){
    return(list(ATE = bcd_ntz_ate
                , Impute_Acc = test_acc
                , Impute_AUC = test_auc
                , Crump.Labeled.Removed = .dips.output$.labeled.removed
                , Crump.Unlabeled.Removed = .dips.output$.unlabeled.removed
                , Crump.Trimmed = .dips.output$trimmed
                , Trt_KS_IntOnly = propens$Trt_IntOnly
                , Outcome_KS_IntOnly = propens$Outcome_IntOnly
                , OC.coef = propens$OC.coef
                , Trt.coef = propens$Trt.coef
                , DiPS.notrim = propens$PS.Kernel
                , Y_dag = Y_dag
                , ImputeModel = list(Model = Xi_Outcome_Coefs
                                     , data = Xfull 
                                     , outcome = analytic_object[[.outcome]]
                                     , cv_outcome_indx = cv_outcome_indx
                                     , test_outcome_indx = test_outcome_indx
                                     )
    )
    )
  }else{
    return(list(ATE = bcd_ntz_ate
                , Trt_KS_IntOnly = propens$Trt_IntOnly
                , Outcome_KS_IntOnly = propens$Outcome_IntOnly
                , Crump.Labeled.Removed = .dips.output$.labeled.removed
                , Crump.Unlabeled.Removed = .dips.output$.unlabeled.removed
                , Crump.Trimmed = .dips.output$trimmed
                , OC.coef = propens$OC.coef
                , Trt.coef = propens$Trt.coef
                , DiPS.notrim = propens$PS.Kernel
                )
           )
  }

}





# Inverting CI's to P-values 
perturb.pvalue = function(.ATE){
  CI = quantile(.ATE, c(0.025, 0.975))
  
  mu = mean(.ATE)
  sgn = sign(mu)
  sigma = sd(.ATE)
  
  .md = (1 + (.ATE - mu)^2 / sigma)^-1
  
  .small = which(sgn*.ATE>0)[ which.min( sgn*.ATE[which(sgn*.ATE>0)] ) ]
  
  .small = which(sgn*.ATE>0)[which.min( sgn*.ATE[which(sgn*.ATE>0)])]
  return(1 - rank(-.md)[.small] / 1000)
}
