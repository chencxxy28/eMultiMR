#'@title Information borrow scores
#'@description This function can be used to obtain the information borrow scores for the eMultiR estimate
#'@usage MinBo.one(s.outcome,cov.mat,n,proportion,time,r,id,dist="gaussian")
#'@param s.outcome The vector of secondary outcomes
#'@param cov.mat The covariate matrix
#'@param n The sample size
#'@param proportion The proportion of the subjects having the secondary outcome
#'(so far only allow the value one)
#'@param time The number of visits
#'@param r Missing indicator for subjects
#'@param id The id of subjects
#'@param dist The distribution of the secondary outcomes (so far only allow "gaussian")
#'@return A vector of information borrow scores
#'@export
#'@import rootSolve MASS geepack

MinBo.one<-function(s.outcome,cov.mat,n,proportion,time,r,id,dist="gaussian")
{
  y<-s.outcome
  x<-cov.mat
  adata<-cbind(y,x)[1:(proportion*n*time),]
  m_1<-n*proportion
  dim(adata)

  ##get initial values for the algorithm
  fit<-geese(y~x-1,id=id)
  beta_initial<-fit$beta

  #function to find lambda, given beta and an
  lambda_find<-function(beta)
  {
    ZZ<-wgeef(beta=beta,adata,r=r,id=id,dist=dist,time=time,m_1=m_1)[,1:(m_1+1)]
    dim(ZZ)
    apply(ZZ,1,mean)

    gamma<-1
    c<-0
    lambda<-rep(0,nrow(ZZ))
    tol<-10e-8
    Delta_old<-0

    repeat{
      rl<-R1der(lambda,ZZ)
      rll<-R2der(lambda,ZZ)
      Delta<--ginv(rll)%*%rl
      if(mean(abs(Delta))<tol | mean(Delta-Delta_old)==0)
      {break}else{
        repeat{
          mm<-0
          repeat{
            delta<-gamma*Delta
            index_1<-apply(ZZ,2,function (xx)
            {ifelse(1+t(lambda+delta)%*%as.matrix(xx,ncol=1)<=0,1,0)}
            )
            if (sum(index_1)>0)
            {gamma<-gamma/2
            mm<-mm+1}else{break}}
          index_2<-ifelse(R0der(lambda+delta,ZZ)-R0der(lambda,ZZ)<0,1,0)
          if (index_2==1)
          {gamma<-gamma/2}else{break}
        }
        Delta_old<-Delta
      }
      lambda<-lambda+delta
      c<-c+1
      gamma<-(c)^(-0.5)
    }
    lambda
  }


  ####multiroot_method
  beta_ee<-function(beta)
  {
    total<-wgeef(beta=beta,adata,r=r,id=id,dist=dist,time=time,m_1=m_1)
    ZZ<-total[,1:(m_1+1)]
    ZZ_d<-total[,(m_1+2):((ncol(x)+1)*m_1+1+ncol(x))]
    beta_ee<-0
    lambda=lambda_find(beta)
    for (i in 1:(m_1+1))
    {
      scaler<-(1/(1+t(matrix(lambda,ncol=1))%*%ZZ[,i]))
      ee_i<-matrix((t(ZZ_d[,((i-1)*ncol(x)+1):(i*ncol(x))])%*%matrix(lambda,ncol=1)),nrow=ncol(x))*as.vector(scaler)
      beta_ee<-beta_ee+ee_i
    }
    beta_ee
  }

  ##get estimated beta
  beta<-multiroot(f = beta_ee, start = as.vector(beta_initial))$root
  lambda<-lambda_find(beta)

  #calculate P based on oracle
  total<-wgeef_oracle(beta=beta,adata,r=r,id=id,dist=dist,time=time,m_1=m_1)
  ZZ<-total[,1:(m_1)]
  ZZ_der<-total[,(m_1+1):((ncol(x)+1)*m_1)]
  Prop_scores<-apply(ZZ,2,function(xx){1/(1+t(matrix(lambda,ncol=1))%*%xx)/m_1})
  Prop_scores<-c(Prop_scores*m_1,rep(1,(1-proportion)*n))
  Prop_scores
}

#'@title The eMultiR estimate with information borrow
#'@description This function can be used to obtain the eMultiR estimate with information borrow.
#'@usage multir.ib.est(Prop_multir, Prop_scores, trt_ind, main_outcome)
#'@param Prop_multir The calibration score
#'@param Prop_scores The information borrow score
#'@param trt_ind Treatment assignment (eg., 0,1)
#'@param main_outcome The main outcome in binary
#'@return A vector of estimates
#'@export
#'@import rootSolve MASS

multir.ib.est<-function(Prop_multir,Prop_scores,trt_ind,main_outcome)
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
      sf_i=x_i*(yy[i]-mu_i)*Prop_scores[i]*Prop_multir[i]
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
  eta_en1_multir<-multiroot(f = nf_en, start = as.vector(eta_initial))$root

  eta_en1_multir
}
