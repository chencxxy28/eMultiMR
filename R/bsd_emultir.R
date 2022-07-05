
#'@title The bootstrapping standard deviation for the eMultiR estimate
#'@description This function can be used to obtain the bootstrapping standard deviation for the eMultiR estimate.
#'@usage bootstrap.multir.plus.en.se(y,x,trt_a,yy,x_reg,r,time,proportion,id,dist,num.boost=50)
#'@param y The vector of secondary outcomes
#'@param x The covariate matrix for the secondary model
#'@param trt_a Treatment assignment (eg., 0,1)
#'@param yy The main outcome in binary
#'@param x_reg The covariate matrix used for both PS and CM models
#'@param r Missing indicator for subjects
#'@param time The number of visits
#'@param proportion The proportion of the subjects having the secondary outcome
#'@param id The id of subjects
#'@param dist The distribution of the secondary outcomes (so far only allow "gaussian")
#'@param num.boost The number of bootstrapping. The default is 50
#'@return The bootstrapping standard deviation
#'@export
#'@import rootSolve MASS
#'@importFrom stats fitted gaussian predict sd time

bootstrap.multir.plus.en.se<-function(y,x,trt_a,yy,x_reg,r,time,proportion,id,dist,num.boost=50)
{
  n=length(yy)
  est.boost<-list(eta_multir_rr=NULL,eta_multir_rw=NULL,
                  eta_multir_wr=NULL,eta_multir_ww=NULL,
                  eta_en1_multir_rr=NULL,eta_en1_multir_rw=NULL,
                  eta_en1_multir_wr=NULL,eta_en1_multir_ww=NULL)
  for(i in 1:num.boost)
  {
    # resample_ind_trt<-sample(which(trt_a==1),replace=TRUE)
    # resample_ind_ctl<-sample(which(trt_a==0),replace=TRUE)
    # resample_ind<-c(resample_ind_trt,resample_ind_ctl)
    resample_ind<-sample(1:n,replace=TRUE)
    yy_b<-yy[resample_ind]
    x_reg_b<-x_reg[resample_ind,]
    trt_a_b<-trt_a[resample_ind]
    resample_ind_long<-rep()
    for (i in 1:n)
    {
      aa<-resample_ind[i]
      resample_ind_long_i<-c((time*(aa-1)+1):(time*aa))
      resample_ind_long<-c(resample_ind_long,resample_ind_long_i)
    }
    x_b<-x[resample_ind_long,]
    y_b<-y[resample_ind_long]

    est.boost.single<-bootstrap.multir.plus.en(y=y_b,x=x_b,trt_a=trt_a_b,yy=yy_b,x_reg=x_reg_b,r=r,proportion=proportion,id=id,dist=dist)
    est.boost=Map(cbind, est.boost, est.boost.single)
  }
  lapply(est.boost,function (x)
  {
    apply(x,1,sd)[2]
  }
  )
}



#'@title The bootstrapping results for the eMultiR estimate based on one sampling
#'@description This function can be used to obtain the bootstrapping result for the eMultiR estimate
#'based on one sampling.
#'@usage bootstrap.multir.plus.en(y,x,trt_a,yy,x_reg,r,time,proportion,id,dist)
#'@param y The vector of secondary outcomes
#'@param x The covariate matrix for the secondary model
#'@param trt_a Treatment assignment (eg., 0,1)
#'@param yy The main outcome in binary
#'@param x_reg The covariate matrix used for both PS and CM models
#'@param r Missing indicator for subjects
#'@param time The number of visits
#'@param proportion The proportion of the subjects having the secondary outcome
#'@param id The id of subjects
#'@param dist The distribution of the secondary outcomes (so far only allow "gaussian")
#'@return The list of bootstrapping estimates
#'@export
#'@import rootSolve MASS
#'@importFrom stats fitted gaussian predict sd time

bootstrap.multir.plus.en<-function(y,x,trt_a,yy,x_reg,r,time,proportion,id,dist)
{
  n<-length(yy)
  Prop_scores<-MinBo.one(s.outcome=y,cov.mat=x,n=n,proportion=proportion,time=time,r=r,id=id,dist="gaussian")
  summary(Prop_scores)

  ######propensity score model
  #fit logistic regression under right model
  fit<-glm(trt_a~x_reg-1,family = binomial(link = "logit"))
  eta_pi_initial<-fit$coefficients
  prob_fitted_right<-fitted(fit)
  #prob_fitted_right[trt_a==0]<-(1-prob_fitted_right)[trt_a==0]

  #fit logistic regression under wrong model 1
  fit<-glm(trt_a~x_reg[,-3]-1,family = binomial(link = "logit"))
  eta_pi_initial<-fit$coefficients
  prob_fitted_wrong<-fitted(fit)
  #prob_fitted_wrong[trt_a==0]<-(1-prob_fitted_wrong)[trt_a==0]

  #fit logistic regression under wrong model 2
  fit<-glm(trt_a~x_reg[,-c(2)]-1,family = binomial(link = "logit"))
  eta_pi_initial<-fit$coefficients
  prob_fitted_wrong2<-fitted(fit)
  #prob_fitted_wrong[trt_a==0]<-(1-prob_fitted_wrong)[trt_a==0]

  #######regression model
  #fit correct model
  data.used.reg<-data.frame(x_reg[,-1],trt_a,x_reg[,-1]*trt_a)
  names(data.used.reg)<-c("x1","x2","x3","trt","x1trt","x2trt","x3trt")
  fit<-glm(yy~x1+x2+x3+trt+x1trt+x2trt+x3trt,family = binomial(link = "logit"),data=data.used.reg)

  newdata.trt1<-data.frame(x_reg[,-1],1,x_reg[,-1]*1)
  names(newdata.trt1)<-c("x1","x2","x3","trt","x1trt","x2trt","x3trt")
  pseudo_y.trt1.right<-predict(fit,newdata.trt1,type ="response")

  newdata.trt0<-data.frame(x_reg[,-1],0,x_reg[,-1]*0)
  names(newdata.trt0)<-c("x1","x2","x3","trt","x1trt","x2trt","x3trt")
  pseudo_y.trt0.right<-predict(fit,newdata.trt0,type ="response")

  data.impute<-data.frame(pseudo_y=c(pseudo_y.trt1.right,pseudo_y.trt0.right),
                          trt=c(rep(1,length(pseudo_y.trt1.right)),
                                rep(0,length(pseudo_y.trt0.right))))

  fit<-glm(pseudo_y~trt,family = gaussian(link = "logit"),data=data.impute)
  est.or.right.all<-fit$coefficients

  ####################  ####################  ####################  ####################

  #fit wrong model 1
  data.used.reg<-data.frame(x_reg[,c(2,4)],trt_a,x_reg[,c(2,4)]*trt_a)
  names(data.used.reg)<-c("x1","x3","trt","x1trt","x3trt")
  fit<-glm(yy~x1+x3+trt+x1trt+x3trt,family = binomial(link = "logit"),data=data.used.reg)

  newdata.trt1<-data.frame(x_reg[,c(2,4)],1,x_reg[,c(2,4)]*1)
  names(newdata.trt1)<-c("x1","x3","trt","x1trt","x3trt")
  pseudo_y.trt1.wrong1<-predict(fit,newdata.trt1,type ="response")

  newdata.trt0<-data.frame(x_reg[,c(2,4)],0,x_reg[,c(2,4)]*0)
  names(newdata.trt0)<-c("x1","x3","trt","x1trt","x3trt")
  pseudo_y.trt0.wrong1<-predict(fit,newdata.trt0,type ="response")

  data.impute<-data.frame(pseudo_y=c(pseudo_y.trt1.wrong1,pseudo_y.trt0.wrong1),
                          trt=c(rep(1,length(pseudo_y.trt1.wrong1)),
                                rep(0,length(pseudo_y.trt0.wrong1))))

  fit<-glm(pseudo_y~trt,family = gaussian(link = "logit"),data=data.impute)
  est.or.wrong1.all<-fit$coefficients

  ####################  ####################  ####################  ####################

  #fit wrong model 2
  data.used.reg<-data.frame(x_reg[,c(3,4)],trt_a,x_reg[,c(3,4)]*trt_a)
  names(data.used.reg)<-c("x2","x3","trt","x2trt","x3trt")
  fit<-glm(yy~x2+x3+trt+x2trt+x3trt,family = binomial(link = "logit"),data=data.used.reg)

  newdata.trt1<-data.frame(x_reg[,c(3,4)],1,x_reg[,c(3,4)]*1)
  names(newdata.trt1)<-c("x2","x3","trt","x2trt","x3trt")
  pseudo_y.trt1.wrong2<-predict(fit,newdata.trt1,type ="response")

  newdata.trt0<-data.frame(x_reg[,c(3,4)],0,x_reg[,c(3,4)]*0)
  names(newdata.trt0)<-c("x2","x3","trt","x2trt","x3trt")
  pseudo_y.trt0.wrong2<-predict(fit,newdata.trt0,type ="response")

  data.impute<-data.frame(pseudo_y=c(pseudo_y.trt1.wrong2,pseudo_y.trt0.wrong2),
                          trt=c(rep(1,length(pseudo_y.trt1.wrong2)),
                                rep(0,length(pseudo_y.trt0.wrong2))))

  fit<-glm(pseudo_y~trt,family = gaussian(link = "logit"),data=data.impute)
  est.or.wrong2.all<-fit$coefficients

  ####################  ####################  ####################  ####################
  ###all wrong models
  #estimating fct for or: trt group
  est.or.fct.trt11<-est.fct.or(eta=est.or.wrong2.all,pseudo_y=pseudo_y.trt1.wrong2,trt=1)
  est.or.fct.trt12<-est.fct.or(eta=est.or.wrong1.all,pseudo_y=pseudo_y.trt1.wrong1,trt=1)
  #control grp
  est.or.fct.trt01<-est.fct.or(eta=est.or.wrong2.all,pseudo_y=pseudo_y.trt0.wrong2,trt=0)
  est.or.fct.trt02<-est.fct.or(eta=est.or.wrong1.all,pseudo_y=pseudo_y.trt0.wrong1,trt=0)

  #caclulate multiple robust propensity score with or given two wrong models for both iptw and or
  Prop_multir_ww<-multir.propensity.or(prob_fit1=prob_fitted_wrong2,prob_fit2=prob_fitted_wrong,
                                       or_fit11=est.or.fct.trt11,or_fit12=est.or.fct.trt12,
                                       or_fit01=est.or.fct.trt01,or_fit02=est.or.fct.trt02,
                                       trt_ind=trt_a)
  #Prop_multir<-multir.propensity(prob_fit1=prob_fitted_wrong2,prob_fit2=prob_fitted_wrong,trt_ind=trt_a)

  #calculate multiple robust en1 casual effect
  eta_en1_multir_ww_pool<-multir.ib.est(Prop_multir=Prop_multir_ww,Prop_scores=Prop_scores,trt_ind=trt_a,main_outcome=yy)

  #calculate multir casual effect
  eta_multir_ww_pool<-multir.est(Prop_multir=Prop_multir_ww,trt_ind=trt_a,main_outcome=yy)

  ####################  ####################  ####################  ####################
  ###all wrong models for or, but one right for iptw

  #estimating fct for or: trt group
  est.or.fct.trt11<-est.fct.or(eta=est.or.wrong2.all,pseudo_y=pseudo_y.trt1.wrong2,trt=1)
  est.or.fct.trt12<-est.fct.or(eta=est.or.wrong1.all,pseudo_y=pseudo_y.trt1.wrong1,trt=1)
  #control grp
  est.or.fct.trt01<-est.fct.or(eta=est.or.wrong2.all,pseudo_y=pseudo_y.trt0.wrong2,trt=0)
  est.or.fct.trt02<-est.fct.or(eta=est.or.wrong1.all,pseudo_y=pseudo_y.trt0.wrong1,trt=0)

  #caclulate multiple robust propensity score with or given two wrong models for both iptw and or
  Prop_multir_rw<-multir.propensity.or(prob_fit1=prob_fitted_right,prob_fit2=prob_fitted_wrong,
                                       or_fit11=est.or.fct.trt11,or_fit12=est.or.fct.trt12,
                                       or_fit01=est.or.fct.trt01,or_fit02=est.or.fct.trt02,
                                       trt_ind=trt_a)
  #Prop_multir<-multir.propensity(prob_fit1=prob_fitted_wrong2,prob_fit2=prob_fitted_wrong,trt_ind=trt_a)

  #calculate multiple robust en1 casual effect
  eta_en1_multir_rw_pool<-multir.ib.est(Prop_multir=Prop_multir_rw,Prop_scores=Prop_scores,trt_ind=trt_a,main_outcome=yy)

  #calculate multir casual effect
  eta_multir_rw_pool<-multir.est(Prop_multir=Prop_multir_rw,trt_ind=trt_a,main_outcome=yy)

  ####################  ####################  ####################  ####################
  ###all wrong models for iptw, but one right for or

  #estimating fct for or: trt group
  est.or.fct.trt11<-est.fct.or(eta=est.or.wrong1.all,pseudo_y=pseudo_y.trt1.wrong1,trt=1)
  est.or.fct.trt12<-est.fct.or(eta=est.or.right.all,pseudo_y=pseudo_y.trt1.right,trt=1)
  #control grp
  est.or.fct.trt01<-est.fct.or(eta=est.or.wrong1.all,pseudo_y=pseudo_y.trt0.wrong1,trt=0)
  est.or.fct.trt02<-est.fct.or(eta=est.or.right.all,pseudo_y=pseudo_y.trt0.right,trt=0)

  #caclulate multiple robust propensity score with or given one right models for both iptw and or
  Prop_multir_wr<-multir.propensity.or(prob_fit1=prob_fitted_wrong2,prob_fit2=prob_fitted_wrong,
                                       or_fit11=est.or.fct.trt11,or_fit12=est.or.fct.trt12,
                                       or_fit01=est.or.fct.trt01,or_fit02=est.or.fct.trt02,
                                       trt_ind=trt_a)
  #Prop_multir<-multir.propensity(prob_fit1=prob_fitted_wrong2,prob_fit2=prob_fitted_wrong,trt_ind=trt_a)

  #calculate multiple robust en1 casual effect
  eta_en1_multir_wr_pool<-multir.ib.est(Prop_multir=Prop_multir_wr,Prop_scores=Prop_scores,trt_ind=trt_a,main_outcome=yy)

  #calculate multir casual effect
  eta_multir_wr_pool<-multir.est(Prop_multir=Prop_multir_wr,trt_ind=trt_a,main_outcome=yy)

  ####################  ####################  ####################  ####################
  ###one right model for iptw, and one right for or

  #estimating fct for or: trt group
  est.or.fct.trt11<-est.fct.or(eta=est.or.wrong1.all,pseudo_y=pseudo_y.trt1.wrong1,trt=1)
  est.or.fct.trt12<-est.fct.or(eta=est.or.right.all,pseudo_y=pseudo_y.trt1.right,trt=1)
  #control grp
  est.or.fct.trt01<-est.fct.or(eta=est.or.wrong1.all,pseudo_y=pseudo_y.trt0.wrong1,trt=0)
  est.or.fct.trt02<-est.fct.or(eta=est.or.right.all,pseudo_y=pseudo_y.trt0.right,trt=0)

  #caclulate multiple robust propensity score with or given one right models for both iptw and or
  Prop_multir_rr<-multir.propensity.or(prob_fit1=prob_fitted_right,prob_fit2=prob_fitted_wrong,
                                       or_fit11=est.or.fct.trt11,or_fit12=est.or.fct.trt12,
                                       or_fit01=est.or.fct.trt01,or_fit02=est.or.fct.trt02,
                                       trt_ind=trt_a)
  #Prop_multir<-multir.propensity(prob_fit1=prob_fitted_wrong2,prob_fit2=prob_fitted_wrong,trt_ind=trt_a)

  #calculate multiple robust en1 casual effect
  eta_en1_multir_rr_pool<-multir.ib.est(Prop_multir=Prop_multir_rr,Prop_scores=Prop_scores,trt_ind=trt_a,main_outcome=yy)

  #calculate multir casual effect
  eta_multir_rr_pool<-multir.est(Prop_multir=Prop_multir_rr,trt_ind=trt_a,main_outcome=yy)

  return(list(eta_multir_rr=eta_multir_rr_pool,eta_multir_rw=eta_multir_rw_pool,
              eta_multir_wr=eta_multir_wr_pool,eta_multir_ww=eta_multir_ww_pool,
              eta_en1_multir_rr=eta_en1_multir_rr_pool,eta_en1_multir_rw=eta_en1_multir_rw_pool,
              eta_en1_multir_wr=eta_en1_multir_wr_pool,eta_en1_multir_ww=eta_en1_multir_ww_pool))
}
