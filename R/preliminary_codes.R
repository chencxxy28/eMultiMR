#'@title Construct a AR1 correlation structure
#'@description This function can be used to create a AR1 correlation structure.
#'@usage ar1_cor(n, rho)
#'@param n The number of time
#'@param rho The correlation coefficient
#'
#'@return A correlation matrix
#'@export

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                    (1:n - 1))
  rho^exponent
}


#'@title Estimating function for the CM model
#'@description This function can be used to create a estimating function for the CM model.
#'@usage est.fct.or(eta,pseudo_y,trt=1,n)
#'@param eta The CM estimate
#'@param pseudo_y The pseudo outcome for treated (1) or control group (0)
#'@param trt Treatment assignment corresponding to `pseudo_y` (either 0 or 1)
#'@param n The sample size
#'@return A correlation matrix
#'@export
#
est.fct.or<-function (eta,pseudo_y,trt=1,n)
{
  yy<-pseudo_y
  sf=rep()
  for (i in 1:n)
  {
    x_i=as.matrix(c(1,trt),ncol=1)
    mu_i=as.vector(1/(1+exp(-t(x_i)%*%eta)))
    sf_i=x_i*(yy[i]-mu_i)
    sf=cbind(sf,sf_i)
  }
  sf
}


##first derivative of -log EL
R1der<-function(lambda,ZZ)
{
  apply(ZZ,2,function(xx)
  {as.matrix(xx,ncol=1)/as.vector((1+t(lambda)%*%as.matrix(xx,ncol=1)))})%*%rep(1,ncol(ZZ))
}

##second derivative of -log EL
R2der<-function(lambda,ZZ)
{
  r2der<-0
  for(i in 1:ncol(ZZ))
  {
    r2der_i<--as.matrix(ZZ[,i],ncol=1)%*%t(as.matrix(ZZ[,i],ncol=1))/as.vector(1+t(lambda)%*%as.matrix(ZZ[,i],ncol=1))^2
    r2der<-r2der+r2der_i
  }
  r2der
}

##-log EL
R0der<-function(lambda,ZZ)
{
  apply(ZZ,2, function (xx) {log(as.vector(1+t(lambda)%*%as.matrix(xx,ncol=1)))})%*%rep(1,ncol(ZZ))
}

#lambda for propensity scores
lambda_propensity<-function(ZZ)
{
  ZZ<-ZZ
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
    if(mean(abs(Delta))<tol | sum(Delta-Delta_old)==0)
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



v<-function(x,beta,dist) #function used to calculate variance, mean, and derivative of random variables followed by binomial, poisson, and gaussian
{
  n<-length(x[,1])
  reg<-x%*%beta
  pi<-exp(reg)/(1+exp(reg))
  if (dist=="gaussian")
  {
    v<-rep(1,n)
    der<-x
    mu<-reg
  }
  else if (dist=="binomial")
  {
    v<-pi*(1-pi)
    der<-x*as.vector(v)
    mu<-as.vector(pi)
  }
  else if (dist=="poisson")
  {
    v<-exp(reg)
    der<-x*as.vector(exp(reg))
    mu<-as.vector(exp(reg))
  }
  list(v=v,der=der,mu=mu)
}

##based on adjusted empirical likelihood
wgeef<-function(beta,adata,r,id,dist,time,m_1)
{ #full wgee
  y<-adata[,1]
  x<-adata[,-1]
  n<-length(unique(id))
  A<-diag(1,time,time)
  # R1<-diag(1,time,time)
  # R2<-matrix(1,nrow=time,ncol=time)-diag(1,time,time)
  # R3<-R2
  # R3[1,4]<- R3[4,1]<- 0
  # R4<-diag(0,time,time)
  # R4[time,time]<-R4[1,1]<-1

  R1<-R2<-R3<-R4<-diag(0,time)
  R1[1,1]<-1
  R2[2,2]<-1
  R3[3,3]<-1
  R4[4,4]<-1

  W<-diag(1,time,time)
  V<-v(x,beta,dist)
  wgeef<-rep()
  dwgeef<-rep()
  sum_dwgee<-0

  #z.col<-ncol(z)

  y[which(is.na(y))]<-0
  x[which(is.na(x))]<-0
  for (i in 1:m_1)
  {
    index_obs<-which(r[((i-1)*time+1):(i*time)]==1)
    AA<-(A*V$v[((i-1)*time+1):(i*time)]^(-0.5))[index_obs,index_obs]
    WW<-W*r[((i-1)*time+1):(i*time)]
    WW_obs<-WW[index_obs,index_obs]
    R_obs1<-R1[index_obs,index_obs]
    R_obs2<-R2[index_obs,index_obs]
    R_obs3<-R3[index_obs,index_obs]
    R_obs4<-R4[index_obs,index_obs]
    e<-(y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)])
    e_obs<-e[index_obs]
    if (length(e_obs)==1)
    {
      wgeei1<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs1)*AA*WW_obs*e_obs)
      wgeei2<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs2)*AA*WW_obs*e_obs)
      wgeei3<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs3)*AA*WW_obs*e_obs)
      wgeei4<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs4)*AA*WW_obs*e_obs)
      dwgeei1<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs1)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei2<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs2)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei3<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs3)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei4<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs4)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))

      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
    }else
    {
      wgeei1<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs1)%*%AA)%*%WW_obs%*%e_obs
      wgeei2<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs2)%*%AA)%*%WW_obs%*%e_obs
      wgeei3<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs3)%*%AA)%*%WW_obs%*%e_obs
      wgeei4<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs4)%*%AA)%*%WW_obs%*%e_obs

      dwgeei1<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs1)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei2<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs2)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei3<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs3)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei4<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs4)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])

      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
    }
    wgeef<-cbind(wgeef,wgeei134)
    dwgeef<-cbind(dwgeef,dwgeei134)
    sum_dwgee<-sum_dwgee+dwgeei134
  }
  wgeef_adjusted<-cbind(wgeef,-max(log(m_1)/2,1)*apply(wgeef,1,mean))
  #l_index<-rep(1:time,time=n)
  dwgeef_adjusted<-cbind(dwgeef,-max(log(m_1)/2,1)*sum_dwgee/m_1)
  return(cbind(wgeef_adjusted,dwgeef_adjusted))
}

##based on empirical likelihood
wgeef_oracle<-function(beta,adata,r,id,dist,time,m_1)
{ #full wgee
  y<-adata[,1]
  x<-adata[,-1]
  n<-length(unique(id))
  A<-diag(1,time,time)
  # R1<-diag(1,time,time)
  # R2<-matrix(1,nrow=time,ncol=time)-diag(1,time,time)
  # R3<-R2
  # R3[1,4]<-R3[4,1]<-0
  # R4<-diag(0,time,time)
  # R4[time,time]<-R4[1,1]<-1

  R1<-R2<-R3<-R4<-diag(0,time)
  R1[1,1]<-1
  R2[2,2]<-1
  R3[3,3]<-1
  R4[4,4]<-1

  W<-diag(1,time,time)
  V<-v(x,beta,dist)
  wgeef<-rep()
  dwgeef<-rep()
  sum_dwgee<-0

  #z.col<-ncol(z)

  y[which(is.na(y))]<-0
  x[which(is.na(x))]<-0
  for (i in 1:m_1)
  {
    index_obs<-which(r[((i-1)*time+1):(i*time)]==1)
    AA<-(A*V$v[((i-1)*time+1):(i*time)]^(-0.5))[index_obs,index_obs]
    WW<-W*r[((i-1)*time+1):(i*time)]
    WW_obs<-WW[index_obs,index_obs]
    R_obs1<-R1[index_obs,index_obs]
    R_obs2<-R2[index_obs,index_obs]
    R_obs3<-R3[index_obs,index_obs]
    R_obs4<-R4[index_obs,index_obs]
    e<-(y[((i-1)*time+1):(i*time)]-V$mu[((i-1)*time+1):(i*time)])
    e_obs<-e[index_obs]
    if (length(e_obs)==1)
    {
      wgeei1<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs1)*AA*WW_obs*e_obs)
      wgeei2<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs2)*AA*WW_obs*e_obs)
      wgeei3<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs3)*AA*WW_obs*e_obs)
      wgeei4<-(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs4)*AA*WW_obs*e_obs)
      dwgeei1<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs1)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei2<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs2)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei3<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs3)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))
      dwgeei4<--(matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))%*%(AA*(R_obs4)*AA*WW_obs*t((matrix(V$der[((i-1)*time+1):(i*time),][index_obs,],ncol=1))))

      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
    }else
    {
      wgeei1<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs1)%*%AA)%*%WW_obs%*%e_obs
      wgeei2<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs2)%*%AA)%*%WW_obs%*%e_obs
      wgeei3<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs3)%*%AA)%*%WW_obs%*%e_obs
      wgeei4<-t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs4)%*%AA)%*%WW_obs%*%e_obs

      dwgeei1<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs1)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei2<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs2)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei3<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs3)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])
      dwgeei4<--t(V$der[((i-1)*time+1):(i*time),][index_obs,])%*%(AA%*%(R_obs4)%*%AA)%*%WW_obs%*%(V$der[((i-1)*time+1):(i*time),][index_obs,])

      wgeei134<-rbind(wgeei1,wgeei2,wgeei3,wgeei4)
      dwgeei134<-rbind(dwgeei1,dwgeei2,dwgeei3,dwgeei4)
    }
    wgeef<-cbind(wgeef,wgeei134)
    dwgeef<-cbind(dwgeef,dwgeei134)
    #sum_dwgee<-sum_dwgee+dwgeei12
  }
  #wgeef_adjusted<-cbind(wgeef,-max(log(n)/2,1)*apply(wgeef,1,mean))
  #l_index<-rep(1:time,time=n)
  #dwgeef_adjusted<-cbind(dwgeef,-max(log(n)/2,1)*sum_dwgee/n)
  return(cbind(wgeef,dwgeef))
}




