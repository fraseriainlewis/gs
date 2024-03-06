## generates data for use in contour plots - see figure 2 and figure 3

rm(list=ls())

library(rootSolve)

library(mnormt)


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
  return(c( (p_t_l-p_c_l-d_l)*sqrt(I_k(p_t_l,p_c_l,n_z1_l)),
            (p_t_l-p_c_l-d_l)*sqrt(I_k(p_t_l,p_c_l,n_z2_l)))
  )
}

gridN<-40; ## grid size to evaluate psi and rho1
#################################################################
################## PART 1 - compute sample size for fixed power
#################################################################

myF<-function(x,pars){
  
  ## sort out extra parameters
  p_t_l=pars["p_t"]
  p_c_l=pars["p_c"]
  psi_l=pars["psi"]
  d_l=pars["d"]
  rho1_l=pars["rho1"]
  rho2_l=pars["rho2"]
  alpha_l=pars["alpha"]
  beta_l=pars["beta"]
  
  #print(p_t_l)
  #print(psi_l)
  
  c11=x[1]
  c12=x[2]
  c2=x[3]
  m=x[4]
  
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
  
  #cat("cross and min diff=",CR1," ",DF1,"\n")
  
  res_l<-c(f1=A1-A1_S,f2=A2-A2_S,f3=B1-B1_S,f4=B2-B2_S)
  names(res_l)<-c("f1.alpha1","f2.alpha2","f3.beta1","f4.beta2")
  
  return(res_l)
  

}


##################

## compute min diff and cross prob
myF2<-function(list_l){
  
  x=list_l[[1]]$root
  pars=list_l[[2]]
  ## sort out extra parameters
  p_t_l=pars["p_t"]
  p_c_l=pars["p_c"]
  psi_l=pars["psi"]
  d_l=pars["d"]
  rho1_l=pars["rho1"]
  rho2_l=pars["rho2"]
  alpha_l=pars["alpha"]
  beta_l=pars["beta"]
  
  c11=x[1]
  c12=x[2]
  c2=x[3]
  m=x[4]
  
  n_z1_l=psi_l*m*m
  n_z2_l=m*m
  
  A1_S<- ((I_k(p_t_l,p_c_l,n_z1_l)/I_k(p_t_l,p_c_l,n_z2_l))^rho1_l)*alpha_l
  
  ## crossing probability
  CR1<-sadmvn(lower=c(c12,-Inf), upper=c(Inf,Inf), mean=mean_H1(p_t_l,p_c_l,n_z1_l,n_z2_l,d_l), 
              varcov=cov_z1z2(p_t_l,p_c_l,n_z1_l,n_z2_l))
  
  ## min diff
  DF1<- d_l+c12/(sqrt(I_k(p_t_l,p_c_l,n_z1_l)))
  
  return(c(CR1,DF1,m*m,A1_S))
  
  #return(c(psi_l,rho1_l,"mindiff"=DF1,"crossprob"=CR1,"A1_spend"=A1_S,"n"=m*m))
  
  
  
}

## Step 1. setup the joint probability distributions

## setup additional parameters - fixed apriori
mypars=c(p_t=0.58,
         p_c=0.6,
         psi=0.6,
         d=-0.1,
         rho1=2,
         rho2=3,
         alpha=0.025,
         beta=0.1)# 90% power

mylevels<-expand.grid(psi=seq(0.2,0.8,len=gridN),rho1=seq(0.01,4,len=gridN))

myres<-data.frame(mylevels,CR=NA,MD=NA,N=NA,alpha=NA)

for(i in 1:nrow(mylevels)){
  cat("i=",i,"\n")
  mypars[c("psi","rho1")]<-c(mylevels[i,1],mylevels[i,2])
  #print(res<-multiroot(f=myF,start=c(0.5,2.3,1.0,10),parms=mypars))
  res<-multiroot(f=myF,start=c(0.5,2.3,1.0,10),parms=mypars)
  restmp<-myF2(list(res,mypars))
  myres[i,3:6]<-restmp;
  #myres[[i]]<-list(res,mypars);
  #startme<-res$root
}

print(head(myres))

write.csv(myres[,c("psi","rho1","CR")],file="crP.csv",row.names=FALSE)
write.csv(myres[,c("psi","rho1","MD")],file="mdP.csv",row.names=FALSE)
write.csv(myres[,c("psi","rho1","N")],file="nP.csv",row.names=FALSE)
write.csv(myres[,c("psi","rho1","alpha")],file="alphaP.csv",row.names=FALSE)

#################################################################
################## PART 2 - compute power for fixed sample size 
#################################################################

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
  
  #cat("cross and min diff=",CR1," ",DF1,"\n")
  
  res_l<-c(f1=A1-A1_S,f2=A2-A2_S,f3=B1-B1_S,f4=B2-B2_S)
  names(res_l)<-c("f1.alpha1","f2.alpha2","f3.beta1","f4.beta2")
  
  return(res_l)
  
  
}


##################

## compute min diff and cross prob
myF2<-function(list_l){
  
  x=list_l[[1]]$root
  pars=list_l[[2]]
  ## sort out extra parameters
  p_t_l=pars["p_t"]
  p_c_l=pars["p_c"]
  psi_l=pars["psi"]
  d_l=pars["d"]
  rho1_l=pars["rho1"]
  rho2_l=pars["rho2"]
  alpha_l=pars["alpha"]
  m=pars["N"]
  
  c11=x[1]
  c12=x[2]
  c2=x[3]
  beta_l=x[4]
  
  n_z1_l=psi_l*m*m
  n_z2_l=m*m
  
  A1_S<- ((I_k(p_t_l,p_c_l,n_z1_l)/I_k(p_t_l,p_c_l,n_z2_l))^rho1_l)*alpha_l
  
  ## crossing probability
  CR1<-sadmvn(lower=c(c12,-Inf), upper=c(Inf,Inf), mean=mean_H1(p_t_l,p_c_l,n_z1_l,n_z2_l,d_l), 
              varcov=cov_z1z2(p_t_l,p_c_l,n_z1_l,n_z2_l))
  
  ## min diff
  DF1<- d_l+c12/(sqrt(I_k(p_t_l,p_c_l,n_z1_l)))
  
  return(c(CR1,DF1,1-beta_l,A1_S))
  
  #return(c(psi_l,rho1_l,"mindiff"=DF1,"crossprob"=CR1,"A1_spend"=A1_S,"n"=m*m))
  
  
  
}

## Step 1. setup the joint probability distributions

## setup additional parameters - fixed apriori
mypars=c(p_t=0.58,
         p_c=0.6,
         psi=0.6,
         d=-0.1,
         rho1=2,
         rho2=3,
         alpha=0.025,
         N=28.8372205)

mylevels<-expand.grid(psi=seq(0.2,0.8,len=gridN),rho1=seq(0.01,4,len=gridN))

myres<-data.frame(mylevels,CR=NA,MD=NA,powerB=NA,alpha=NA)

for(i in 1:nrow(mylevels)){
  cat("i=",i,"\n")
  mypars[c("psi","rho1")]<-c(mylevels[i,1],mylevels[i,2])
  #print(res<-multiroot(f=myF,start=c(0.5,2.3,1.0,10),parms=mypars))
  res<-multiroot(f=myF,start=c(0.5,2.3,1.0,0.1),parms=mypars)
  restmp<-myF2(list(res,mypars))
  myres[i,3:6]<-restmp;
  #myres[[i]]<-list(res,mypars);
  #startme<-res$root
}

print(head(myres))

write.csv(myres[,c("psi","rho1","CR")],file="crN.csv",row.names=FALSE)
write.csv(myres[,c("psi","rho1","MD")],file="mdN.csv",row.names=FALSE)
write.csv(myres[,c("psi","rho1","powerB")],file="powerN.csv",row.names=FALSE)
write.csv(myres[,c("psi","rho1","alpha")],file="alphaN.csv",row.names=FALSE)


