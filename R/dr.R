#'@title DR estimate
#'@description This function can be used to obtain the doubly robust estimate.
#'@usage dr.est(prob_fit,trt_ind,main_outcome,pseudo_outcome_obs,
#'pseudo_outcome_comb,pseudo_trt)
#'@param prob_fit The fitted probability in the treatment assignment model
#'@param trt_ind Treatment assignment (eg., 0,1)
#'@param main_outcome The main outcome in binary
#'@param pseudo_outcome_obs Pseudo outcome fitted from the first step of the CM model
#'@param pseudo_outcome_comb Pseudo outcome combined from treated and control groups
#' in the CM model
#'@param pseudo_trt Pseudo treatment assignment combined
#'@return A vector of estimates
#'@export
#'@import rootSolve MASS

dr.est<-function(prob_fit,trt_ind,main_outcome,
                 pseudo_outcome_obs,
                 pseudo_outcome_comb,
                 pseudo_trt
)
{
  trt_a<-trt_ind
  yy<-main_outcome
  n<-length(yy)
  iptw_weight<-prob_fit
  iptw_weight[trt_a==0]<-(1-prob_fit)[trt_a==0]

  sf=rep()
  for (i in 1:n)
  {
    x_i=as.matrix(c(1,trt_a[i]),ncol=1)
    mu_i=pseudo_outcome_obs[i]
    sf_i=x_i*(yy[i]-mu_i)/iptw_weight[i]
    sf=cbind(sf,sf_i)
  }
  ee1<-apply(sf,1,sum)

  ee2_function<-function (eta)
  {
    sf=rep()
    for (i in 1:(2*n))
    {
      x_i=as.matrix(c(1,pseudo_trt[i]),ncol=1)
      mu_i=as.vector(1/(1+exp(-t(x_i)%*%eta)))
      sf_i=x_i*(pseudo_outcome_comb[i]-mu_i)
      sf=cbind(sf,sf_i)
    }
    apply(sf,1,sum)+ee1
  }

  fit<-glm(yy~trt_a,family = binomial(link = "logit"))
  eta_initial<-fit$coefficients
  eta<-multiroot(f = ee2_function, start = as.vector(eta_initial))$root
  eta
}
