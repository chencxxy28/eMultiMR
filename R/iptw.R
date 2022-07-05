#'@title IPTW estimate
#'@description This function can be used to obtain the IPTW estimate.
#'@usage iptw.binary(prob_fit,trt_ind,main_outcome)
#'@param prob_fit The fitted probability in the treatment assignment model
#'@param trt_ind Treatment assignment
#'@param main_outcome The main outcome in binary
#'
#'@return A vector of estimates
#'@export
#'@import rootSolve MASS
#'@importFrom stats binomial glm


iptw.binary<-function(prob_fit,trt_ind,main_outcome)
{
  trt_a<-trt_ind
  yy<-main_outcome
  n<-length(yy)
  iptw_weight<-prob_fit
  iptw_weight[trt_a==0]<-(1-prob_fit)[trt_a==0]

  score_function<-function (eta)
  {
    sf=rep()
    for (i in 1:n)
    {
      x_i=as.matrix(c(1,trt_a[i]),ncol=1)
      mu_i=as.vector(1/(1+exp(-t(x_i)%*%eta)))
      sf_i=x_i*(yy[i]-mu_i)/iptw_weight[i]
      sf=cbind(sf,sf_i)
    }
    sf
  }
  nf<-function (eta)
  {
    score_function(eta)%*%rep(1,n)
  }
  fit<-glm(yy~trt_a,family = binomial(link = "logit"))
  eta_initial<-fit$coefficients
  eta<-multiroot(f = nf, start = as.vector(eta_initial))$root
  eta
}
