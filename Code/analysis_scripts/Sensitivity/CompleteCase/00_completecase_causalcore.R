
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
    OC.coef = ordinalCoefs(x = cbind(Xi, Trt=Ti)
                           , y = Yi)
  }else{
    OC.coef = adalasso(x = cbind(Xi
                                 , Trt=Ti
                                 )
                       , y = Yi
                       , family = .family
                       , nfolds = .nfolds
                       , weights = weights
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
  
  # if(crump.trim){
  #   .labeled.removed = sum((PS<0.05 | PS>0.95) & !is.na(labels) )
  #   .unlabeled.removed = sum((PS<0.05 | PS>0.95) & is.na(labels) )
  #   if(.labeled.removed + .unlabeled.removed <= 0.4*length(labels)){
  #     trimmed = T
  #     .keep = which(PS>=0.1 & PS<=0.9)
  #     PS = PS[.keep]
  #     weights = weights[.keep]
  #     Ti = Ti[.keep]
  #     Yi = Yi[.keep]
  #   }
  # }else{
  #   .n.removed = 0
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
                   , trimmed = trimmed
                   )
    return(.output)
    
  }


############### #
# Full calculation
############### # 
AIPW_robust_analytic = function(analytic_object
                                , .outcome # c("Y_change", "Y_NegChange", "Y_AnyIncrease", "Y_diff")
                                , pweights = NULL
                                ){
    

  ## Kernel Smoothing #### 
  # ridge_12mo_EHR_NLP_latent = pdds_impute(.model = ridge_12mo_EHR_NLP, .analytic_X = analytic_object$X
  #                                         , .print=F, type="response")
  ridge_All_EHR_NLP_latent = pdds_impute(.model = ridge_All_EHR_NLP, .analytic_X = analytic_object$X
                                         , .print=F, type="response")
  
  X_inclPDDSlatent = cbind(analytic_object$X
                           # , ridge_12mo_EHR_NLP_latent
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
  summary(propens$PS.Kernel)
  
  ## Utility Covariate (truncating/cleaning)
  .pi = propens$PS.Kernel 
  .pi_trim = pmax(pmin(.pi, 0.999), 0.001) # stability trim
  
  ## ATE Estimation ####
  
  .dips.output = ATE_DiPS(Ti = analytic_object$A
                         , Yi = analytic_object[[.outcome]]
                         , PS = .pi_trim
                         , min.ipw = 0.01 
                         , max.ipw = 10
                         , .IPW_stabilize = T
                         , weights = pweights
                         , labels = analytic_object[[.outcome]]
                         , crump.trim = T
                         )

  bcd_ntz_ate = .dips.output$.ate 
  
  return(list(ATE = bcd_ntz_ate
              , Trt_KS_IntOnly = propens$Trt_IntOnly
              , Outcome_KS_IntOnly = propens$Outcome_IntOnly
              , OC.coef = propens$OC.coef
              , Trt.coef = propens$Trt.coef
              , DiPS.notrim = propens$PS.Kernel
              , Crump.Labeled.Removed = .dips.output$.labeled.removed
              , Crump.Unlabeled.Removed = .dips.output$.unlabeled.removed
              )
         )

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
