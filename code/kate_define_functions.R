#2.0 DEFINING THE FUNCTIONS

#homogeneous
pred.Exp<-function(theta, dose) {
  f<- 1-exp(-theta*dose)
  return(f)
}

###approximate beta-poisson
pred.betaPoisson<-function(alpha, betap, dose){
  f<-1-(1+(dose/betap))^-alpha
  return(f)
}


#gamma - in terms of shape and scale
pred.inf.gamma<-function(theta,k,omega,dose) {
  ((theta^(k-1))*(exp(-theta/omega)))/((omega^k)*(gamma(k)))*(exp(-theta*dose))
}


#integration function for gamma - shape and scale
int2<-function(k,omega,dose){
  n <- length(dose)
  ig=NA;id=NA
  for (i in 1:n){
    ig=integrate(pred.inf.gamma,0,Inf,k=k,omega=omega,dose=dose[i],stop.on.error = FALSE)$value
    id=rbind(id,ig)
    pi_pred=1-as.vector(id)
    pi_pred=pi_pred[!is.na(pi_pred)]
  } 
  return(pi_pred)
}

####LIKELIHOOD FUNCTIONS#######

#likelihood function for gamma (no laplace)

Lik.fun2=function(k,omega,dose) {# function for calculating the negative log likelihood
  print(c(k,omega)) #print current values for debugging
  pi_pred=int2(k,omega,dose)
  #run the rest  
  pi_obs<- pos.c/tot.c
  #deviance calcs
  #Y1<- pos.c*log((pi_pred+1e-15)/(pi_obs+1e-15))
  #Y2<- (tot.c-pos.c)*log((1-pi_pred+1e-15)/(1-pi_obs+1e-15))
  #Y<- (-2)*(sum(Y1)+sum(Y2))
  #return(Y)
  #Lik calcs 
  V.LL <- dbinom(x=pos.c,size=tot.c,prob=pi_pred,log=TRUE)
  tot.LL=V.LL
  print(-sum(tot.LL))
  return(-sum(tot.LL))
  
}

Lik.fun2(0.5,1,dose)


###LIKELIHOOD FUNCTION FOR BETA ALONE
Lik.fun.beta=function(alpha,betap,dose) {# function for calculating the negative log likelihood
  pi_pred.v=pred.betaPoisson(alpha=alpha,betap=betap,dose=dose)
  #run the rest  
  pi_obs.v<-pos.c/tot.c
  #Deviance calcs - vax
  Y1<- pos.c*log((pi_pred.v+1e-15)/(pi_obs.v+1e-15))
  Y2<- (tot.c-pos.c)*log((1-pi_pred.v+1e-15)/(1-pi_obs.v+1e-15))
  Y.v<- (-2)*(sum(Y1)+sum(Y2))
  Y<-Y.v
  return(Y)
  #likelihood calcs
  #V.LL <- dbinom(x=pos.c,size=tot.c,prob=pi_pred.v,log=TRUE)
  #tot.LL=V.LL
  #print(-sum(tot.LL))
  #return(-sum(tot.LL))
}

###LIKELIHOOD FOR HOMOGENEOUS MODEL ALONE###
Lik.fun.hom=function(theta,dose) {# function for calculating the negative log likelihood
  pi_pred.c=pred.Exp(theta=theta,dose=dose)
  #run the rest  
  pi_obs.c<- pos.c/tot.c
  #Deviance calcs - control
  Y1<- pos.c*log((pi_pred.c+1e-15)/(pi_obs.c+1e-15))
  Y2<- (tot.c-pos.c)*log((1-pi_pred.c+1e-15)/(1-pi_obs.c+1e-15))
  Y.c<- (-2)*(sum(Y1)+sum(Y2))
  Y<-Y.c
  return(Y)
  #V.LL <- dbinom(x=pos.c,size=tot.c,prob=pi_pred.c,log=TRUE)
  #tot.LL=V.LL
  #print(-sum(tot.LL))
  #return(-sum(tot.LL))
}


#ARRANGING VARIABLES FOR OPTIM STEP

##gamma model with shape and scale
nll.gamma.full <- function (par) {
  Lik.fun2(k=par[1],omega=par[2], dose=dose)#
}

nll.gamma.full(c(1,1))

#beta model
nll.bet <- function (par) {
  Lik.fun.beta(alpha=par[1],betap=par[2], dose=dose)#
}

#hom model
nll.hom <- function (par) {
  Lik.fun.hom(theta=par[1], dose=dose)#
}

