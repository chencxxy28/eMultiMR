#'@title Calibration scores
#'@description This function can be used to obtain the calibration scores for the MultiR estimate
#'@usage multir.propensity.or(prob_fit1, prob_fit2,or_fit11, or_fit12,
#'or_fit01,or_fit02, trt_ind)
#'@param prob_fit1 The fitted probability in the treatment assignment model1
#'@param prob_fit2 The fitted probability in the treatment assignment model2
#'@param or_fit11 The estimating function for the CM model 1 for the treatment group
#'@param or_fit12 The estimating function for the CM model 2 for the treatment group
#'@param or_fit01 The estimating function for the CM model 1 for the control group
#'@param or_fit02 The estimating function for the CM model 2 for the control group
#'@param trt_ind Treatment assignment (eg., 0,1)
#'@return A vector of calibration scores
#'@export
#'@import rootSolve MASS

multir.propensity.or<-function(prob_fit1,prob_fit2,
                               or_fit11,or_fit12,
                               or_fit01,or_fit02,
                               trt_ind)
{
  trt_a<-trt_ind
  #iptw for trt
  prob_fitted_right_mean<-mean(prob_fit1)
  prob_fitted_wrong_mean<-mean(prob_fit2)
  #or for trt
  or_fit11.mean<-apply(or_fit11,1,mean)
  or_fit12.mean<-apply(or_fit12,1,mean)

  ZZ_trt_propensity<-rbind((as.vector(prob_fit1)-as.vector(prob_fitted_right_mean))[trt_a==1],
                           (as.vector(prob_fit2)-as.vector(prob_fitted_wrong_mean))[trt_a==1],
                           ((or_fit11)-(or_fit11.mean))[,trt_a==1],
                           ((or_fit12)-(or_fit12.mean))[,trt_a==1]
  )
  lambda_prop_trt<-lambda_propensity(ZZ_trt_propensity)
  num_trt<-sum(trt_a==1)
  Prop_trt<-apply(ZZ_trt_propensity,2,function(xx){1/(1+t(matrix(lambda_prop_trt,ncol=1))%*%xx)/num_trt})

  #iptw for control
  prob_fitted_right_control<-1-prob_fit1
  prob_fitted_wrong_control<-1-prob_fit2
  prob_fitted_right_mean_control<-mean(prob_fitted_right_control)
  prob_fitted_wrong_mean_control<-mean(prob_fitted_wrong_control)
  #or for control
  or_fit01.mean<-apply(or_fit01,1,mean)
  or_fit02.mean<-apply(or_fit02,1,mean)

  ZZ_control_propensity<-rbind((as.vector(prob_fitted_right_control)-as.vector(prob_fitted_right_mean_control))[trt_a==0],
                               (as.vector(prob_fitted_wrong_control)-as.vector(prob_fitted_wrong_mean_control))[trt_a==0],
                               ((or_fit01)-(or_fit01.mean))[,trt_a==0],
                               ((or_fit02)-(or_fit02.mean))[,trt_a==0]
  )
  lambda_prop_control<-lambda_propensity(ZZ_control_propensity)
  num_control<-sum(trt_a==0)
  Prop_control<-apply(ZZ_control_propensity,2,function(xx){1/(1+t(matrix(lambda_prop_control,ncol=1))%*%xx)/num_control})

  #combine propensity
  Prop_multir<-rep()
  Prop_multir[trt_a==1]<-Prop_trt
  Prop_multir[trt_a==0]<-Prop_control

  Prop_multir
}


#'@title The MultiR estimate
#'@description This function can be used to obtain the MultiR estimate.
#'@usage multir.est(Prop_multir, trt_ind, main_outcome)
#'@param Prop_multir The calibration score
#'@param trt_ind Treatment assignment (eg., 0,1)
#'@param main_outcome The main outcome in binary
#'@return A vector of estimates
#'@export
#'@import rootSolve MASS

multir.est<-function(Prop_multir,trt_ind,main_outcome)
{
  yy<-main_outcome
  n=length(main_outcome)
  trt_a<-trt_ind
  score_function_en<-function (eta)
  {
    sf=rep()
    for (i in 1:n)
    {
      x_i=as.matrix(c(1,trt_a[i]),ncol=1)
      mu_i=as.vector(1/(1+exp(-t(x_i)%*%eta)))
      sf_i=x_i*(yy[i]-mu_i)*Prop_multir[i]
      sf=cbind(sf,sf_i)
    }
    sf
  }
  nf_en<-function (eta)
  {
    score_function_en(eta)%*%rep(1,n)
  }
  fit<-glm(yy~trt_a,family = binomial(link = "logit"))
  eta_initial<-fit$coefficients
  eta_multir<-multiroot(f = nf_en, start = as.vector(eta_initial))$root

  eta_multir
}
