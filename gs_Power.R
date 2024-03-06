## this code computes power in a group sequential design
rm(list=ls())
library(mnormt)
library(rootSolve)

## Step 1. setup default parameter set

p_t<-0.58 #0.58 - to compute mid diff use H0 value for control
p_c<-0.6
delta<- -0.1
prop_psi<-0.6

# function to compute Information - 1/variance - of natural estimator = P_t - P_c - delta

I_k<-function(p_t_l,p_c_l,n_k_l){
  # P(success in test arm), P(success in comparator arm), n in each arm 
  # uses simple wald variance estimate p(1-p)/n
  thevar<-(p_t_l*(1-p_t_l))/n_k_l + (p_c_l*(1-p_c_l))/n_k_l
  return(thevar^-1); #information is reciprocal variance
}

cov_z1z2<-function(p_t_l,p_c_l,n_z1_l,n_z2_l){
  #compute covariance matrix between test statistics at interim and end of study
  covar_m<-matrix(rep(1,4),ncol=2)
  # covar is sqrt ratio of information
  covar_m[1,2]<-covar_m[2,1]<-sqrt(I_k(p_t_l,p_c_l,n_z1_l)/I_k(p_t_l,p_c_l,n_z2_l))
  return(covar_m)
  
}

mean_H0<-function(p_t_l,p_c_l,n_z1_l,n_z2_l,d_l){return(c(0,0))}

mean_H1<-function(p_t_l,p_c_l,n_z1_l,n_z2_l,d_l){
    return(c( (p_t-p_c-d_l)*sqrt(I_k(p_t,p_c,n_z1_l)),
            (p_t-p_c-d_l)*sqrt(I_k(p_t,p_c,n_z2_l)))
  )
}



myF<-function(x,pars){
  
  ## sort out extra parameters
  p_t_l=pars["p_t"]
  p_c_l=pars["p_c"]
  psi_l=pars["psi"]
  d_l=pars["d"]
  rho1_l=pars["rho1"]
  rho2_l=pars["rho2"]
  alpha_l=pars["alpha"]
  m=pars["N"]

  #print(p_t_l)
  #print(psi_l)
  
  c11=x[1]
  c12=x[2]
  c2=x[3]
  beta_l=x[4]
  
  n_z1_l=psi_l*m*m
  n_z2_l=m*m
  
  ## compute Type I and Type II errors
  
  ## type I at IA
  
  A1<-sadmvn(lower=c(c12,-Inf), upper=c(Inf,Inf), mean=mean_H0(p_t_l,p_c_l,n_z1_l,n_z2_l,d_l), 
                                                  varcov=cov_z1z2(p_t_l,p_c_l,n_z1_l,n_z2_l))
  ## type I at EoS
  
  A2<-sadmvn(lower=c(-Inf,c2), upper=c(c12,Inf), mean=mean_H0(p_t_l,p_c_l,n_z1_l,n_z2_l,d_l), 
                                                 varcov=cov_z1z2(p_t_l,p_c_l,n_z1_l,n_z2_l))
  ## type II at IA
  
  B1<-sadmvn(lower=c(-Inf,-Inf), upper=c(c11,Inf), mean=mean_H1(p_t_l,p_c_l,n_z1_l,n_z2_l,d_l), 
                                                   varcov=cov_z1z2(p_t_l,p_c_l,n_z1_l,n_z2_l))
  ## type II at EoS
  B2<-sadmvn(lower=c(c11,-Inf), upper=c(c12,c2), mean=mean_H1(p_t_l,p_c_l,n_z1_l,n_z2_l,d_l), 
                                                 varcov=cov_z1z2(p_t_l,p_c_l,n_z1_l,n_z2_l))
  
  ## compute error spends same order as above
  
  A1_S<- ((I_k(p_t_l,p_c_l,n_z1_l)/I_k(p_t_l,p_c_l,n_z2_l))^rho1_l)*alpha_l
  
  A2_S<- alpha_l-((I_k(p_t_l,p_c_l,n_z1_l)/I_k(p_t_l,p_c_l,n_z2_l))^rho1_l)*alpha_l
  
  B1_S<- ((I_k(p_t_l,p_c_l,n_z1_l)/I_k(p_t_l,p_c_l,n_z2_l))^rho2_l)*beta_l
  
  B2_S<- beta_l-((I_k(p_t_l,p_c_l,n_z1_l)/I_k(p_t_l,p_c_l,n_z2_l))^rho2_l)*beta_l
  
  ## crossing probability
  CR1<-sadmvn(lower=c(c12,-Inf), upper=c(Inf,Inf), mean=mean_H1(p_t_l,p_c_l,n_z1_l,n_z2_l,d_l), 
              varcov=cov_z1z2(p_t_l,p_c_l,n_z1_l,n_z2_l))
  
  ## min diff
  DF1<- d_l+c12/(sqrt(I_k(p_t_l,p_c_l,n_z1_l)))
  
  cat("cross and min diff and A1=",CR1," ",DF1," ",A1,"\n")
  
  res_l<-c(f1=A1-A1_S,f2=A2-A2_S,f3=B1-B1_S,f4=B2-B2_S)
  names(res_l)<-c("f1.alpha1","f2.alpha2","f3.beta1","f4.beta2")
  
  return(res_l)
  
  
  
}

## Step 1. setup the joint probability distributions

## sample size fixed apriori = 28.1162728^2

## setup additional parameters - fixed apriori
mypars=c(p_t=0.58,
         p_c=0.6,
         psi=0.5,
         d=-0.1,
         rho1=1.,
         rho2=3,
         alpha=0.025,
         N=28.8372205)

# test function call
cat("### Test objective function call ###\n")
print(myF(c(0.5,2.6,2.03,0.1),pars=mypars))

# run root solver
cat("### Run root finder ###\n")
# note - under-estimate n as a starting condition
# integral accuracy is sensitive to higher n as tails becoming vanishingly small 
print(res<-multiroot(f=myF,start=c(0.5,2.,2.0,0.01),parms=mypars))

# check solution
cat("### Check solution - should be close to zero###\n")
print(myF(res$root,pars=mypars))

