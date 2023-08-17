##### Plots SIR model trajectories and equilibrium comparisons from optimal parameters #####

# load required packages
require(deSolve)
library(lme4)
library(ggplot2)
library(gridExtra)
library(ggpattern)


rm(list=ls()) # clear workspace
#dev.off() # close graphical device (if necessary)
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

############################################################################
################################## INITIAL SETUP ###########################
############################################################################

##look at beta distribution for this range of values
bins=300 #number to discretize
x_s=seq(0,1,length=bins+1) #x value of rectangle

# Set the midpoints of the discretizing kernel. Using midpoint rule for evaluation
# these are the midpoints - the y's for susceptibility
y = 0.5 * (x_s[1:(bins)] + x_s[2:(bins+1)])

# Width of cells
h = x_s[2] - x_s[1]

# initial conditions
N=100
I_0=1
R_0=0
C_0=1
# S_0 defined later

# time
tend = 300
times<- seq(from=1, to=tend, by=1)

# Parameters
beta = 0.00275 # infectivity (= transmission) Williams 2014
gamma = 0.03 # recovery Williams 2014


#####################################################################
###################  FITTED MORTALITY  ##############################
#####################################################################

# Fitted mortality parameters
mu_no_prior = 5.892373e-2 # no prior exposure mortality
mu_low_prior = 2.811438e-2 # low prior exposure mortality
mu_high_prior = 2.726825e-17 # high prior exposure mortality

############  NO PRIOR EXPOSURE - HOMOGENEOUS ############ 

# no exposure
VE = 0.88646119

# initial conditions
S_0=(N-I_0) #number of individuals in each susceptibility class
xstart <-c(S=S_0, I=I_0, R=R_0, C=C_0)

# parameters
params_hom_no <- c(beta=beta,VE = VE, gamma=gamma, mu=mu_no_prior) 

# run homogeneous model
out <- ode(
  func=sir.hom.model,
  y=xstart,
  times=times,
  parms=params_hom_no
)

outhom = as.data.frame(out)
outhom$D = 100-(rowSums(outhom[,2:4]))# deaths
outhom$group = "no prior exposure (hom)"
tail(outhom)

############### LOW PRIOR EXPOSURE - HETEROGENEOUS ###############

# 750 VA
alpha= 0.20632603 
betap= 0.2316793

alpha/(alpha+betap) # mean heterogeneity

ic1 = pbeta(x_s,alpha,betap)
ic2 = ic1[2:length(x_s)]-ic1[1:length(x_s)-1]
VE1 = sum(ic2*y) # mean of discretized distribution

#set the frequency of the population in each group
init_susc=ic2
sum(init_susc)
# plot(y,init_susc*100)

# initial condition
S_0=init_susc*(N-I_0) #number of individuals in each susceptibility class
xstart <-c(S=S_0,I=I_0,R=R_0,C=C_0)

# parameters
params_het_low <- c(beta=beta,gamma=gamma,mu=mu_low_prior) 

# run heterogeneous model for low previous exposure
out <- ode(
  func=sir.het.model,
  y=xstart,
  times=times,
  parms=params_het_low
)

out750 = as.data.frame(out)
out750$D = 100-(rowSums(out750[,2:303]))# deaths
out750$S = rowSums(out750[,2:301])
keeps <- c("time","S","I","R","C","D")
out750 = out750[keeps]
out750$group = "low prior exposure (het)"
tail(out750)


############### LOW PRIOR EXPOSURE - HOMOGENEOUS  ###############

# 750 VA
alpha= 0.20632603 
betap= 0.2316793

# no exposure
VE = alpha/(alpha+betap) # homogeneous = mean heterogeneity (true distribution)

# initial conditions
S_0=(N-I_0) #number of individuals in each susceptibility class
xstart <-c(S=S_0,I=I_0,R=R_0,C=C_0)

# parameters
params_hom_low <- c(beta=beta, VE = VE1, gamma=gamma, mu=mu_low_prior)  

# run homogeneous model
out <- ode(
  func=sir.hom.model,
  y=xstart,
  times=times,
  parms=params_hom_low
)

outhom750 = as.data.frame(out) 
outhom750$D = 100-(rowSums(outhom750[,2:4]))# deaths
outhom750$group = "low prior exposure (hom)"
tail(outhom750)

############### HIGH PRIOR EXPOSURE - HETEROGENEOUS ###############

# 30000 VA
alpha = 0.10642501
betap = 0.16509577

alpha/(alpha+betap) #mean heterogeneity

ic1 = pbeta(x_s,alpha,betap)
ic2 = ic1[2:length(x_s)]-ic1[1:length(x_s)-1]
VE1 = sum(ic2*y)


#set the frequency of the population in each group
init_susc=ic2
sum(init_susc)
# plot(y,init_susc*100)

# initial condition
S_0=init_susc*(N-I_0) #number of individuals in each susceptibility class
xstart <-c(S=S_0,I=I_0,R=R_0,C=C_0)

# parameters
params_het_high <- c(beta=beta,gamma=gamma,mu=mu_high_prior) 

# run heterogeneous model
out <- ode(
  func=sir.het.model,
  y=xstart,
  times=times,
  parms=params_het_high
)

out30000 = as.data.frame(out)
out30000$D = 100-(rowSums(out30000[,2:303]))# deaths
out30000$S = rowSums(out30000[,2:301])
keeps <- c("time","S","I","R","C","D")
out30000 = out30000[keeps]
out30000$group = "high prior exposure (het)"
tail(out30000)

############### HIGH PRIOR EXPOSURE - HOMOGENEOUS ###############

# 30000 VA
alpha = 0.10642501
betap = 0.16509577

# no exposure
VE = alpha/(alpha+betap) # mean heterogeneity (true distribution)

# initial conditions
S_0=(N-I_0) #number of individuals in each susceptibility class
xstart <-c(S=S_0,I=I_0,R=R_0,C=C_0)

# parameters
params_hom_high <- c(beta=beta,VE = VE1,gamma=gamma,mu=mu_high_prior) 

# run homogeneous model
out <- ode(
  func=sir.hom.model,
  y=xstart,
  times=times,
  parms=params_hom_high
)

outhom30000 = as.data.frame(out)
outhom30000$D = 100-(rowSums(outhom30000[,2:4]))# deaths
outhom30000$group = "high prior exposure (hom)"
tail(outhom30000)


################################################################################
##################### COMBINING DATA - FITTED MORTALITY ########################
################################################################################

library(tidyverse)
library(viridis)
out.comb = bind_rows(outhom, out750, outhom750, out30000, outhom30000)
out.comb750 = bind_rows(out750, outhom750)
out.comb30000 = bind_rows(out30000, outhom30000)
head(out.comb)

out.comb750$exposure.group = "homogeneous"
out.comb750$exposure.group[out.comb750$group=="high prior exposure (het)"]="heterogeneous"
out.comb750$exposure.group[out.comb750$group=="low prior exposure (het)"]="heterogeneous"


out.comb30000$exposure.group = "homogeneous"
out.comb30000$exposure.group[out.comb30000$group=="high prior exposure (het)"]="heterogeneous"
out.comb30000$exposure.group[out.comb30000$group=="low prior exposure (het)"]="heterogeneous"

out.comb$exposure.group = "homogeneous"
out.comb$exposure.group[out.comb$group=="high prior exposure (het)"]="heterogeneous"
out.comb$exposure.group[out.comb$group=="low prior exposure (het)"]="heterogeneous"
out.comb$exposure.group[out.comb$group=="high prior exposure (hom)"]="homogeneous"
out.comb$exposure.group[out.comb$group=="low prior exposure (hom)"]="homogeneous"
out.comb$exposure.group[out.comb$group=="no prior exposure (hom)"]="homogeneous"


out.comb$prior.exposure = "homogeneous"
out.comb$prior.exposure[out.comb$group=="high prior exposure (het)"]="heterogeneous"
out.comb$prior.exposure[out.comb$group=="low prior exposure (het)"]="heterogeneous"
out.comb$prior.exposure[out.comb$group=="high prior exposure (hom)"]="homogeneous"
out.comb$prior.exposure[out.comb$group=="low prior exposure (hom)"]="homogeneous"
out.comb$prior.exposure[out.comb$group=="no prior exposure (hom)"]="homogeneous"


################################################################################
####################### TRAJECTORIES - FITTED MORTALITY ########################
################################################################################

p11 <- ggplot(outhom, aes(x=time)) + 
  geom_line(aes(y = C), color = "#7E1717",size = 0.8) +
  geom_line(aes(y = D), color="#7E1717", linetype="twodash",size = 0.8) +
  labs(title="No prior exposure",y="Cumulative",x="",legend ="")+
  ylim(0,100)+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position=c(.7,.65),
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"))

p12 = ggplot(data=out.comb750, aes(x=time)) +
  geom_line(aes( y=C, colour = exposure.group),size = 0.8) +
  geom_line(aes( y=D, colour = exposure.group), linetype="twodash",size = 0.8) +
  scale_color_manual(values = c("#068DA9", "#7E1717"))+
  labs(title="Low dose (750)",y="",x="",legend ="")+
  ylim(0,100)+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position="none")

p13 = ggplot(data=out.comb30000, aes(x=time)) +
  geom_line(aes( y=C, colour = exposure.group),size = 0.8) +
  geom_line(aes( y=D, colour = exposure.group), linetype="twodash",size = 0.8) +
  scale_color_manual(values = c("#068DA9", "#7E1717"))+
  labs(title="High dose (30000)",y="",x="",legend ="")+
  ylim(0,100)+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position="none")


p21 <- ggplot(outhom, aes(x=time))+
  geom_line(aes( y=I),colour="#7E1717",size = 0.8) +
  labs(y="Infectious population",x="Time",legend ="")+
  ylim(0,50)+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position=c(.7,.65),
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"))

p22 = ggplot(data=out.comb750, aes(x=time)) +
  geom_line(aes( y=I, colour = exposure.group),size = 0.8) +
  scale_color_manual(values = c("#068DA9", "#7E1717"))+
  labs(y="",x="Time",legend ="")+
  ylim(0,50)+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position="none")

p23 = ggplot(data=out.comb30000, aes(x=time)) +
  geom_line(aes( y=I, colour = exposure.group),size = 0.8) +
  scale_color_manual(values = c("#068DA9", "#7E1717"))+
  labs(y="",x="Time",legend ="")+
  ylim(0,50)+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position="none")


g1 <- arrangeGrob(p11, p12, p13, p21, p22, p23, nrow=2) #generates g1
ggsave(file="Trajectories_fittedmortality.png", g1,width=8,height=6,units="in",limitsize=FALSE) #saves g1


################################################################################
####################### BAR CHART - FITTED MORTALITY ###########################
################################################################################

## Data for bar chart
df_fitted <- data.frame(prior_exposure=rep(c("none", "low \n(750)", "high \n(30000)"),2),
                  Model = rep(c("homogeneous \n(exponential)", "heterogeneous \n(beta)"), each=3),
                  TotI2 = c(outhom$C[tend],outhom750$C[tend],outhom30000$C[tend],NA,out750$C[tend],out30000$C[tend]),
                  Died = c(outhom$D[tend],outhom750$D[tend],outhom30000$D[tend],NA,out750$D[tend],out30000$D[tend]))
df_fitted

## Bar chart
barfitted <- ggplot(data=df_fitted, aes(x=prior_exposure, y=TotI2, fill=Model,width=0.6)) +
  geom_bar(stat="identity", position=position_dodge(.7),alpha = 0.4)  + 
  scale_x_discrete(limits = c("none", "low \n(750)", "high \n(30000)")) +
  xlab("Prior exposure (dose)") +
  ylab("Total size of epidemic \n (darker color = died)") +
  geom_col(aes(x=prior_exposure,y=Died, fill = Model), position =position_dodge(.7)) + 
  ggtitle("Fitted Mortality Rate")+
  scale_fill_manual(values = c("#068DA9", "#7E1717"))+
  ylim(0,100)+
  theme_bw()+
  theme(strip.background = element_rect(fill="gray97"),
        strip.text.x = element_text (size = 15,hjust = 0.5, vjust = 0.5,face="italic"),
        axis.title=element_text(size=20),
        axis.text.y=element_text(size=15),
        plot.title = element_text(size=20),
        axis.text.x=element_text(size=15),
        #panel.grid = element_blank(), 
        panel.grid.major.y = element_line(color="lightgray", linetype="dotted"), 
        panel.grid.minor = element_blank(),
        axis.line=element_line(),
        legend.position="top",legend.title=element_blank(),
        legend.text = element_text(size=12))
barfitted

ggsave(file="Bar_fittedmortality.png",width=9,height=6,units="in",limitsize=FALSE)


############################################################################
#################### REPEAT WITH FIXED MORTALITY ###########################
############################################################################


# Fixed mortality parameters
mu_no_prior = 5.892373e-2 # no prior exposure mortality
mu_low_prior = mu_no_prior
mu_high_prior = mu_no_prior

############  NO PRIOR EXPOSURE - HOMOGENEOUS ############ 

# no exposure
VE = 0.88646119

# initial conditions
S_0=(N-I_0) #number of individuals in each susceptibility class
xstart <-c(S=S_0, I=I_0, R=R_0, C=C_0)

# parameters
params_hom_no <- c(beta=beta,VE = VE, gamma=gamma, mu=mu_no_prior) 

# run homogeneous model
out <- ode(
  func=sir.hom.model,
  y=xstart,
  times=times,
  parms=params_hom_no
)

outhomfix = as.data.frame(out)
outhomfix$D = 100-(rowSums(outhomfix[,2:4]))# deaths
outhomfix$group = "no prior exposure (hom)"
tail(outhomfix)



############### LOW PRIOR EXPOSURE - HETEROGENEOUS ###############

# 750 VA
alpha= 0.20632603 
betap= 0.2316793

alpha/(alpha+betap) # mean heterogeneity

ic1 = pbeta(x_s,alpha,betap)
ic2 = ic1[2:length(x_s)]-ic1[1:length(x_s)-1]
VE1 = sum(ic2*y)

#set the frequency of the population in each group
init_susc=ic2
sum(init_susc)
# plot(y,init_susc*100)

# initial condition
S_0=init_susc*(N-I_0) #number of individuals in each susceptibility class
xstart <-c(S=S_0,I=I_0,R=R_0,C=C_0)

# parameters
params_het_low <- c(beta=beta,gamma=gamma,mu=mu_low_prior) 

# run heterogeneous model for low previous exposure
out <- ode(
  func=sir.het.model,
  y=xstart,
  times=times,
  parms=params_het_low
)

out750fix = as.data.frame(out)
out750fix$D = 100-(rowSums(out750fix[,2:303]))# deaths
out750fix$S = rowSums(out750fix[,2:301])
keeps <- c("time","S","I","R","C","D")
out750fix = out750fix[keeps]
out750fix$group = "low prior exposure (het)"
tail(out750fix)

############### LOW PRIOR EXPOSURE - HOMOGENEOUS  ###############

# 750 VA
alpha= 0.20632603 
betap= 0.2316793

# no exposure
VE = alpha/(alpha+betap) # homogeneous = mean heterogeneity

# initial conditions
S_0=(N-I_0) #number of individuals in each susceptibility class
xstart <-c(S=S_0,I=I_0,R=R_0,C=C_0)

# parameters
params_hom_low <- c(beta=beta, VE = VE1, gamma=gamma, mu=mu_low_prior)  

# run homogeneous model
out <- ode(
  func=sir.hom.model,
  y=xstart,
  times=times,
  parms=params_hom_low
)

outhom750fix = as.data.frame(out) 
outhom750fix$D = 100-(rowSums(outhom750fix[,2:4]))# deaths
outhom750fix$group = "low prior exposure (hom)"
tail(outhom750fix)

############### HIGH PRIOR EXPOSURE - HETEROGENEOUS ###############

# 30000 VA
alpha = 0.10642501
betap = 0.16509577

alpha/(alpha+betap) #mean heterogeneity


ic1 = pbeta(x_s,alpha,betap)
ic2 = ic1[2:length(x_s)]-ic1[1:length(x_s)-1]
VE1 = sum(ic2*y)

#set the frequency of the population in each group
init_susc=ic2
sum(init_susc)
# plot(y,init_susc*100)

# initial condition
S_0=init_susc*(N-I_0) #number of individuals in each susceptibility class
xstart <-c(S=S_0,I=I_0,R=R_0,C=C_0)

# parameters
params_het_high <- c(beta=beta,gamma=gamma,mu=mu_high_prior) 

# run heterogeneous model
out <- ode(
  func=sir.het.model,
  y=xstart,
  times=times,
  parms=params_het_high
)

out30000fix = as.data.frame(out)
out30000fix$D = 100-(rowSums(out30000fix[,2:303]))# deaths
out30000fix$S = rowSums(out30000fix[,2:301])
keeps <- c("time","S","I","R","C","D")
out30000fix = out30000fix[keeps]
out30000fix$group = "high prior exposure (het)"
tail(out30000fix)

############### HIGH PRIOR EXPOSURE - HOMOGENEOUS ###############

# 30000 VA
alpha = 0.10642501
betap = 0.16509577

# no exposure
VE = alpha/(alpha+betap) #mean heterogeneity

# initial conditions
S_0=(N-I_0) #number of individuals in each susceptibility class
xstart <-c(S=S_0,I=I_0,R=R_0,C=C_0)

# parameters
params_hom_high <- c(beta=beta,VE = VE1,gamma=gamma,mu=mu_high_prior) 

# run homogeneous model
out <- ode(
  func=sir.hom.model,
  y=xstart,
  times=times,
  parms=params_hom_high
)

outhom30000fix = as.data.frame(out)
outhom30000fix$D = 100-(rowSums(outhom30000fix[,2:4]))# deaths
outhom30000fix$group = "high prior exposure (hom)"
tail(outhom30000fix)

################################################################################
################################ COMBINING DATA ################################
################################################################################

library(tidyverse)
library(viridis)
out.combfix = bind_rows(outhomfix, out750fix, outhom750fix, out30000fix, outhom30000fix)
out.comb750fix = bind_rows(out750fix, outhom750fix)
out.comb30000fix = bind_rows(out30000fix, outhom30000fix)
head(out.combfix)

out.comb750fix$exposure.group = "homogeneous"
out.comb750fix$exposure.group[out.comb750fix$group=="high prior exposure (het)"]="heterogeneous"
out.comb750fix$exposure.group[out.comb750fix$group=="low prior exposure (het)"]="heterogeneous"


out.comb30000fix$exposure.group = "homogeneous"
out.comb30000fix$exposure.group[out.comb30000fix$group=="high prior exposure (het)"]="heterogeneous"
out.comb30000fix$exposure.group[out.comb30000fix$group=="low prior exposure (het)"]="heterogeneous"

out.combfix$exposure.group = "homogeneous"
out.combfix$exposure.group[out.combfix$group=="high prior exposure (het)"]="heterogeneous"
out.combfix$exposure.group[out.combfix$group=="low prior exposure (het)"]="heterogeneous"
out.combfix$exposure.group[out.combfix$group=="high prior exposure (hom)"]="homogeneous"
out.combfix$exposure.group[out.combfix$group=="low prior exposure (hom)"]="homogeneous"
out.combfix$exposure.group[out.combfix$group=="no prior exposure (hom)"]="homogeneous"


out.combfix$prior.exposure = "homogeneous"
out.combfix$prior.exposure[out.combfix$group=="high prior exposure (het)"]="heterogeneous"
out.combfix$prior.exposure[out.combfix$group=="low prior exposure (het)"]="heterogeneous"
out.combfix$prior.exposure[out.combfix$group=="high prior exposure (hom)"]="homogeneous"
out.combfix$prior.exposure[out.combfix$group=="low prior exposure (hom)"]="homogeneous"
out.combfix$prior.exposure[out.combfix$group=="no prior exposure (hom)"]="homogeneous"


################################################################################
################################# TRAJECTORIES #################################
################################################################################



p11fix <- ggplot(outhomfix, aes(x=time)) + 
  geom_line(aes(y = C), color = "#7E1717",size = 0.8) +
  geom_line(aes(y = D), color="#7E1717", linetype="twodash",size = 0.8) +
  labs(title="No prior exposure",y="Cumulative",x="",legend ="")+
  #ylab("Cumulative")+
  #xlab("")+
  ylim(0,100)+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position=c(.7,.65),
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"))

p12fix = ggplot(data=out.comb750fix, aes(x=time)) +
  geom_line(aes( y=C, colour = exposure.group),size = 0.8) +
  geom_line(aes( y=D, colour = exposure.group), linetype="twodash",size = 0.8) +
  scale_color_manual(values = c("#068DA9", "#7E1717"))+
  labs(title="Low dose (750)",y="",x="",legend ="")+
  #ylab("")+
  #xlab("")+
  ylim(0,100)+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position="none")

p13fix = ggplot(data=out.comb30000fix, aes(x=time)) +
  geom_line(aes( y=C, colour = exposure.group),size = 0.8) +
  geom_line(aes( y=D, colour = exposure.group), linetype="twodash",size = 0.8) +
  scale_color_manual(values = c("#068DA9", "#7E1717"))+
  labs(title="High dose (30000)",y="",x="",legend ="")+
  #ylab("")+
  #xlab("")+
  ylim(0,100)+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position="none")


p21fix <- ggplot(outhomfix, aes(x=time))+
  geom_line(aes( y=I),colour="#7E1717",size = 0.8) +
  labs(y="Infectious population",x="Time",legend ="")+
  #ylab("Infectious population")+
  #xlab("Time")+
  ylim(0,50)+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position=c(.7,.65),
        legend.title=element_blank(),
        legend.text = element_text(size=20,face="italic"))

p22fix = ggplot(data=out.comb750fix, aes(x=time)) +
  geom_line(aes( y=I, colour = exposure.group),size = 0.8) +
  scale_color_manual(values = c("#068DA9", "#7E1717"))+
  labs(y="",x="Time",legend ="")+
  #ylab("")+
  #xlab("Time")+
  ylim(0,50)+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position="none")

p23fix = ggplot(data=out.comb30000fix, aes(x=time)) +
  geom_line(aes( y=I, colour = exposure.group),size = 0.8) +
  scale_color_manual(values = c("#068DA9", "#7E1717"))+
  labs(y="",x="Time",legend ="")+
  #ylab("")+
  #xlab("Time")+
  ylim(0,50)+
  theme_bw() +
  theme(axis.title=element_text(size=14),
        axis.text.y=element_text(size=15),
        axis.text.x=element_text(angle=45,hjust=1,size=15,face="italic"),
        panel.grid = element_blank(), 
        axis.line=element_line(),
        legend.position="none")


g2 <- arrangeGrob(p11fix, p12fix, p13fix, p21fix, p22fix, p23fix, nrow=2) #generates g2
ggsave(file="Trajectories_fixedmortality.png", g2,width=8,height=6,units="in",limitsize=FALSE) #saves g2




################################################################################
################################### PLOTTING ###################################
################################################################################

## Data for bar chart
dffixed <- data.frame(prior_exposure=rep(c("none", "low \n(750)", "high \n(30000)"),2),
                  Model = rep(c("homogeneous \n(exponential)", "heterogeneous \n(beta)"), each=3),
                  TotI2 = c(outhomfix$C[tend],outhom750fix$C[tend],outhom30000fix$C[tend],NA,out750fix$C[tend],out30000fix$C[tend]),
                  Died = c(outhomfix$D[tend],outhom750fix$D[tend],outhom30000fix$D[tend],NA,out750fix$D[tend],out30000fix$D[tend]))
dffixed




## Bar chart
bar_fixed <- ggplot(data=dffixed, aes(x=prior_exposure, y=TotI2, fill=Model,width=0.6)) +
  geom_bar(stat="identity", position=position_dodge(.7),alpha = 0.4)  + 
  scale_x_discrete(limits = c("none", "low \n(750)", "high \n(30000)")) +
  xlab("Prior exposure (dose)") +
  ylab("Total size of epidemic \n (darker color = died)") +
  geom_col(aes(x=prior_exposure,y=Died, fill = Model), position =position_dodge(.7)) + 
  ggtitle("Fixed Mortality Rate")+
  scale_fill_manual(values = c("#068DA9", "#7E1717"))+
  ylim(0,100)+
  theme_bw()+
  theme(strip.background = element_rect(fill="gray97"),
        strip.text.x = element_text (size = 15,hjust = 0.5, vjust = 0.5,face="italic"),
        axis.title=element_text(size=20),
        axis.text.y=element_text(size=15),
        plot.title = element_text(size=20),
        axis.text.x=element_text(size=15),
        #panel.grid = element_blank(), 
        panel.grid.major.y = element_line(color="lightgray", linetype="dotted"), 
        panel.grid.minor = element_blank(),
        axis.line=element_line(),
        legend.position="top",legend.title=element_blank(),
        legend.text = element_text(size=12))
bar_fixed

ggsave(file="Bar_fixedmortality.png",width=9,height=6,units="in",limitsize=FALSE)

############### Comparisons (homogeneous - heterogeneous) ###############
df_fitted$TotI2[2]-df_fitted$TotI2[5]
df_fitted$TotI2[3]-df_fitted$TotI2[6]

df_fitted$Died[2]-df_fitted$Died[5]
df_fitted$Died[3]-df_fitted$Died[6]

dffixed$TotI2[2]-dffixed$TotI2[5]
dffixed$TotI2[3]-dffixed$TotI2[6]

dffixed$Died[2]-dffixed$Died[5]
dffixed$Died[3]-dffixed$Died[6]


