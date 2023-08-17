##### Plots SIR model comparisons from bootstrap parameters #####

require(deSolve)
library(lme4)
library(ggplot2)
library(gridExtra)


rm(list=ls()) # clear workspace
# dev.off() # close graphical device (if necessary)
cat("\014") # clear console

############################################################################
################################## MODELS ##################################
############################################################################


######### heterogeneous SIR model #########
sir.het.model <- function (t, x, params) {
  S <- x[1:bins]
  I <- x[bins+1]
  R <- x[bins+2]
  ## extract the parameters
  beta <- params["beta"]
  gamma <- params["gamma"]
  mu <- params["mu"]
  ## model equations
  dS.dt <- -(y*beta*I*S) # susceptible
  dI.dt <- sum(y*beta*I*S)-gamma*I-mu*I # infectious
  dR.dt <- gamma*I # recovered
  dC.dt <- sum(y*beta*I*S) # cumulative cases
  dxdt <- c(dS.dt,dI.dt,dR.dt,dC.dt)
  ## return result as a list
  list(dxdt)
}

####### homogeneous SIR model #########
sir.hom.model <- function (t, x, params) {
  ## first extract the state variables
  S <- x[1]
  I <- x[2]
  R <- x[3]
  ## extract the parameters
  beta <- params["beta"]
  VE <- params["VE"]
  gamma <- params["gamma"]
  mu <- params["mu"]
  ## model equations
  dS.dt <- -(VE*beta*I*S) # susceptible
  dI.dt <- (VE*beta*I*S)-gamma*I-mu*I # infectious
  dR.dt <- gamma*I # recovered
  dC.dt <- (VE*beta*I*S) # cumulative cases
  ## combine results into a single vector
  dxdt <- c(dS.dt,dI.dt,dR.dt,dC.dt)
  ## return result as a list
  list(dxdt)
}

####### Set Up Simulations #########

input_parms <- read.csv("../all_pseudosets.csv")

# Parameters
beta = 0.00275 # infectivity (= transmission) Williams 2014
gamma = 0.03 # recovery Williams 2014


##look at beta distribution for this range of values
bins=300 #number to discretize
x_s=seq(0,1,length=bins+1) #x value of rectangle

# Set the midpoints of the discretizing kernel. Using midpoint rule for evaluation
# these are the midpoints - the y's for susceptibility
y = 0.5 * (x_s[1:(bins)] + x_s[2:(bins+1)])

# Width of cells
h = x_s[2] - x_s[1]

############################################################################
############################## FUNCTION TO RUN MODELS ######################
############################################################################


################### Function to run Models ##################
outSIR <- function (betap,alpha,params){
  
  VE = alpha/(alpha+betap) # mean heterogeneity (true distribution)
  
  ic1 = pbeta(x_s,alpha,betap)
  ic2 = ic1[2:length(x_s)]-ic1[1:length(x_s)-1]
  VE1 = sum(ic2*y) # mean heterogeneity (discretized distribution)
  
  #set the frequency of the population in each group
  init_susc=ic2
  sum(init_susc)
  # plot(y,init_susc*100)
  
  # time
  tend = 300
  times<- seq(from=1, to=tend, by=1)
  
  # initial conditions
  N=100
  I_0=1
  R_0=0
  C_0=1
  # S_0 defined later
  
  ### heterogeneous model ###
  
  # initial condition
  S_0=init_susc*(N-I_0) #number of individuals in each susceptibility class
  xstart <-c(S=S_0,I=I_0,R=R_0,C=C_0)
  
  # parameters
  params_het <- params #c(beta=beta,gamma=gamma,mu=mu_high_prior) 
  
  # run heterogeneous model
  out <- ode(
    func=sir.het.model,
    y=xstart,
    times=times,
    parms=params_het
  )
  
  outtmp = as.data.frame(out)
  outtmp$S = rowSums(outtmp[,2:301])
  outtmp$D = 100-(rowSums(outtmp[,2:303]))# deaths
  keeps <- c("S","I","R","C","D")
  outtmp = outtmp[keeps]
  
  ### homogeneous model ###
  
  # initial condition
  S_0=(N-I_0) #number of individuals in each susceptibility class
  xstart <-c(S=S_0,I=I_0,R=R_0,C=C_0)
  
  params1 <- c(params,VE=VE1)
  params_hom <- params1
  
  # run homogeneous model
  out <- ode(
    func=sir.hom.model,
    y=xstart,
    times=times,
    parms=params_hom
  )
  
  outtmpH = as.data.frame(out)
  outtmpH$D = 100-(rowSums(outtmpH[,2:4]))# deaths
  keeps <- c("S","I","R","C","D")
  outtmpH = outtmpH[keeps]
  colnames(outtmpH) <-c('Sh','Ih','Rh','Ch','Dh')
  
  list(cbind(outtmp[tend,],outtmpH[tend,]))
  
}

############################################################################
############################## SET-UP SIMULATIONS ##########################
############################################################################

# Mortality parameters
mu_no_prior = 5.892373e-2 # no prior exposure mortality
mu_low_prior = 2.811438e-2 # low prior exposure mortality
mu_high_prior = 2.726825e-17 # high prior exposure mortality

# 750 VA

params <- c(beta=beta,gamma=gamma,mu=mu_low_prior) 
par.alpha = input_parms$par1.beta[input_parms$bird.groups=="750 va"]
par.beta = input_parms$par2.beta[input_parms$bird.groups=="750 va"]

mylist <- list()
for (j in 1:length(par.beta)){
  betap = par.beta[j]
  alpha = par.alpha[j]
  
  outtmp1 <- outSIR(betap,alpha,params)
  mylist[[j]] <- as.data.frame(outtmp1)
}
df750 <- do.call("rbind",mylist) #combine all vectors into a matrix
df750$Cdiff = df750$Ch-df750$C
df750$Ddiff = df750$Dh-df750$D

pC <- ggplot(df750, aes(x=Cdiff)) + 
  geom_density(fill="#00BFC4")+
  labs(title="Low dose (750)", x="Difference in epidemic size", y="Denisty") +
  coord_cartesian(xlim = c(-1, 60))+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=0),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position=c(.7,.65),
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"))
pD <- ggplot(df750, aes(x=Ddiff)) + 
  geom_density(fill="#00BFC4")+
  labs( x="Difference in deaths", y="Denisty") +
  coord_cartesian(xlim = c(-1, 20))+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=0),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position=c(.7,.65),
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"))



# 30000 VA
params <- c(beta=beta,gamma=gamma,mu=mu_high_prior)  
par.alpha = input_parms$par1.beta[input_parms$bird.groups=="30000 va"]
par.beta = input_parms$par2.beta[input_parms$bird.groups=="30000 va"]
 
mylist <- list()
for (j in 1:length(par.beta)){
  betap = par.beta[j]
  alpha = par.alpha[j]
  
  outtmp1 <- outSIR(betap,alpha,params)
  mylist[[j]] <- as.data.frame(outtmp1)
}
df <- do.call("rbind",mylist) #combine all vectors into a matrix
df$Cdiff = df$Ch-df$C
df$Ddiff = df$Dh-df$D
 
pC1 <- ggplot(df, aes(x=Cdiff)) + 
   geom_density(fill="#00BFC4")+
  labs(title="High dose (30000)", x="Difference in epidemic size", y="") +
  coord_cartesian(xlim = c(-1, 60))+
   theme_bw() +
   theme(axis.title=element_text(size=14),
         axis.text.y=element_text(size=0),
         axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
         panel.grid = element_blank(), 
         axis.line=element_line(),
         legend.position=c(.7,.65),
         legend.title=element_blank(),
         legend.text = element_text(size=20,face="italic"))
pD1 <- ggplot(df, aes(x=Ddiff)) + 
  geom_density(fill="#00BFC4")+
  labs(x="Difference in deaths ", y="") +
  coord_cartesian(xlim = c(-1, 20))+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=0),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position=c(.7,.65),
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"))
 
 
g1 <- arrangeGrob(pC, pC1, pD, pD1, nrow=2) #generates g
ggsave(file="Density_fittedmortality.png", g1) #saves g



############################################################################
######################### REPEAT WITH FIXED MORTALITY ######################
############################################################################


# Mortality parameters
mu_no_prior = 5.892373e-2 # no prior exposure mortality
mu_low_prior = mu_no_prior
mu_high_prior = mu_no_prior

# 750 VA

params <- c(beta=beta,gamma=gamma,mu=mu_low_prior) 
par.alpha = input_parms$par1.beta[input_parms$bird.groups=="750 va"]
par.beta = input_parms$par2.beta[input_parms$bird.groups=="750 va"]

mylist <- list()
for (j in 1:length(par.beta)){
  betap = par.beta[j]
  alpha = par.alpha[j]
  
  outtmp1 <- outSIR(betap,alpha,params)
  mylist[[j]] <- as.data.frame(outtmp1)
}
df750fix <- do.call("rbind",mylist) #combine all vectors into a matrix
df750fix$Cdiff = df750fix$Ch-df750fix$C
df750fix$Ddiff = df750fix$Dh-df750fix$D

pCfix <- ggplot(df750fix, aes(x=Cdiff)) + 
  geom_density(fill="#00BFC4")+
  labs(title="Low dose (750)", x="Difference in epidemic size", y="Denisty") +
  coord_cartesian(xlim = c(-1, 30))+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=0),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position=c(.7,.65),
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"))
pDfix <- ggplot(df750fix, aes(x=Ddiff)) + 
  geom_density(fill="#00BFC4")+
  labs( x="Difference in deaths", y="Denisty") +
  coord_cartesian(xlim = c(-1, 20))+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=0),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position=c(.7,.65),
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"))



# 30000 VA
params <- c(beta=beta,gamma=gamma,mu=mu_high_prior)  
par.alpha = input_parms$par1.beta[input_parms$bird.groups=="30000 va"]
par.beta = input_parms$par2.beta[input_parms$bird.groups=="30000 va"]

mylist <- list()
for (j in 1:length(par.beta)){
  betap = par.beta[j]
  alpha = par.alpha[j]
  
  outtmp1 <- outSIR(betap,alpha,params)
  mylist[[j]] <- as.data.frame(outtmp1)
}
dffix <- do.call("rbind",mylist) #combine all vectors into a matrix
dffix$Cdiff = dffix$Ch-dffix$C
dffix$Ddiff = dffix$Dh-dffix$D

pC1fix <- ggplot(dffix, aes(x=Cdiff)) + 
  geom_density(fill="#00BFC4")+
  labs(title="High dose (30000)", x="Difference in epidemic size", y="") +
  coord_cartesian(xlim = c(-1, 30))+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=0),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position=c(.7,.65),
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"))
pD1fix <- ggplot(dffix, aes(x=Ddiff)) + 
  geom_density(fill="#00BFC4")+
  labs(x="Difference in deaths ", y="") +
  coord_cartesian(xlim = c(-1, 20))+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=0),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position=c(.7,.65),
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"))


g2 <- arrangeGrob(pCfix, pC1fix, pDfix, pD1fix, nrow=2) #generates g2
ggsave(file="Density_fixedmortality.png", g2) #saves g2


########## Statistics on Distributions ########## 
quantile(df750$Cdiff,c(0.025, 0.975))
quantile(df750$Ddiff,c(0.025, 0.975))

quantile(df$Cdiff,c(0.025, 0.975))
quantile(df$Ddiff,c(0.025, 0.975))

quantile(df750fix$Cdiff,c(0.025, 0.975))
quantile(df750fix$Ddiff,c(0.025, 0.975))

quantile(dffix$Cdiff,c(0.025, 0.975))
quantile(dffix$Ddiff,c(0.025, 0.975))

