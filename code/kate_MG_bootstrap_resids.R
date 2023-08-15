#bootstrapped CIs
##bootstraps for dose-response


##########run kate_gamma_integral_full_pipeline.R first!!!####
####this will feed in the functions#####
####you will want to run up to making the graphs####

store2 = read.csv("data/deviance_params_21JUL2023.csv")
z3 = read.csv("data/va_hofi_dose_response_dose_aggregated.csv")
source("code/kate_define_functions.R")

print(store2)

#models for predicting
store2 = store2 %>%
  drop_na()

#homogeneous - va0
#est.con.h
va.0.preds=pred.Exp(store2$par1.hom[store2$group=="0 va"],
         z3$dose_scale[z3$group=="0 va"])
va.0 = z3 %>%
  filter(group=="0 va")%>%
  select(frac, dose_scale, pos, tot,group)%>%
  mutate(preds = va.0.preds)%>%
  mutate(err = 
           (frac - preds)/sqrt(preds*(1-preds)/tot))


#predict va 750
#est.vax.beta.only
va.mid.preds=pred.betaPoisson(store2$par1.beta[store2$group=="750 va"],
                 store2$par2.beta[store2$group=="750 va"],
                 z3$dose_scale[z3$group=="750 va"])

va.mid = z3 %>%
  filter(group=="750 va")%>%
  select(frac, dose_scale, pos, tot,group)%>%
  mutate(preds = va.mid.preds)%>%
  mutate(err = 
         (frac - preds)/sqrt(preds*(1-preds)/tot))



#predict va 30000
va.high.preds=pred.betaPoisson(store2$par1.beta[store2$group=="30000 va"],
                                   store2$par2.beta[store2$group=="30000 va"],
                                   z3$dose_scale[z3$group=="30000 va"])

va.high = z3 %>%
  filter(group=="30000 va")%>%
  select(frac, dose_scale, pos, tot,group)%>%
  mutate(preds = va.high.preds)%>%
  #calculate chi-squared residual at each dose
  mutate(err = 
           (frac - preds)/sqrt(preds*(1-preds)/tot))

va.comb = bind_rows(va.high, va.0, va.mid)%>%
  drop_na(err)
#err.v=(pi_obs.v-pi_pred.v)/sqrt(pi_pred.v*(1-pi_pred.v)/tot.v)
#err.c=(pi_obs.c-pi_pred.c)/sqrt(pi_pred.c*(1-pi_pred.c)/tot.c)

va.comb$group.dose = paste(va.comb$bird.groups, va.comb$dose_scale, sep="_")
R=1000
pi_m.v=NA;pi_m.v.st=NA;dose.st=NA;pi_m.c=NA;pi_m.c.st=NA;

for (i in 1:length(unique(va.comb$group.dose))){
    #for each dose within each group
    for (j in 1:R) {
      #at one dose, sample all the errors for that group 1000 times 
      #and calculate a new predicted value for that dose
      # i do not understand the math on this but if you look at it,it seems to work
    pi_m.v=va.comb$preds[i]+sample(va.comb$err,1)*sqrt((va.comb$preds[i]*(1-va.comb$preds[i]))/va.comb$tot[i])
    pi_m.v[pi_m.v<=0]=0;pi_m.v[pi_m.v>=1]=1
      pi_m.v.st=c(pi_m.v.st,pi_m.v)
    }
}

length(pi_m.v.st)
#this was done a bit jenky because before I knew tidyverse so I just 
#keep everything in order and then combine it together
#obviously an issue if it get unsorted
group.dose.st = rep(va.comb$group.dose,each=1000)

x=colsplit(group.dose.st, "_",c( "bird.groups.st","dose.st"))
#separate the doses and groups for the next fitting

#need to remove the NA, the first observation
pi_m.v.st <- pi_m.v.st[!is.na(pi_m.v.st)]
length(pi_m.v.st)

#combine the vectors together to make a new dataframe for modelling
mframe=data.frame(pi_m.v.st,x, group.dose.st)

#now I need to make up the numbers of estimated positives
#based on the probabilities
#i think i originally sampled out of 20 since that was our fish sample size
#for the birds we mostly had 12 so I went with that
#i technically should have matched to the correct numbers but it shouldn't matter that much


mframe$pos.c=NA

for(i in 1:length(mframe$pi_m.v.st))  {
set1=rbinom(prob=mframe$pi_m.v.st[i],n=12,size=1)
mframe$pos.c[i]=sum(set1)
}  

mframe$tot.c=12

#this was in the original and it counts to 1000 7 times
#i think i need to do 13 because that is the number of groups
#i have no idea why i was done this way instead of grouping by the dose.groups
#basically we need to pick the same dose and the right number of values to fit to
mframe$bs=rep(1:1000,times=13)
#in order to get the groups and sequence to line up, i think i can paste together
mframe$count.me = paste(mframe$bird.groups.st, mframe$bs, sep = "_")
#this is one psuedoset for fitting

head(mframe);tail(mframe)


##fit the models to pseudo-datasets

store.pseudo.preds=data.frame(
  #par1.gamma=NA, par2.gamma=NA, lik.gamma=NA,
                 par1.beta=NA, par2.beta=NA, lik.beta=NA, 
                 par1.hom =NA, lik.hom = NA,
                 group=NA)
store.pseudo.preds2=store.pseudo.preds
unique.counts = unique(mframe$count.me)


for (i in 1:length(unique(mframe$count.me))){
#for (i in 1:nrow(mframe)){
  pos.c = mframe$pos.c[mframe$count.me==unique.counts[i]]
  tot.c = mframe$tot.c[mframe$count.me==unique.counts[i]]
  dose = mframe$dose.st[mframe$count.me==unique.counts[i]]
  group = unique.counts[i]
  #gamma.con.la<- optim(par=c(1,1), fn=nll.gamma, gr=NULL, method="Nelder-Mead", control=list(trace=10))
  #gamma.full<- optim(par=c(1,1), fn=nll.gamma.full, gr=NULL, method="Nelder-Mead", control=list(trace=10))
  beta.fit<- optim(par=c(1,3), fn=nll.bet, gr=NULL, method="Nelder-Mead", control=list(trace=10))
  hom.o<- optim(par=c(1), fn=nll.hom, gr=NULL, method="Brent",lower=0,upper=1, control=list(trace=10))
  store.pseudo.preds = data.frame(
                    #par1.gamma = gamma.full$par[1], 
                     #par2.gamma = gamma.full$par[2], 
                     #lik.gamma = gamma.full$val[1],
                     par1.beta = beta.fit$par[1], 
                     par2.beta = beta.fit$par[2], 
                     lik.beta = beta.fit$val[1],
                     par1.hom = hom.o$par[1], 
                     lik.hom = hom.o$val[1],
                     group = group)
  store.pseudo.preds2=bind_rows(store.pseudo.preds2,store.pseudo.preds)
}

p=colsplit(store.pseudo.preds2$group, "_",c( "bird.groups","counter"))
store.pseudo.preds2 = cbind(store.pseudo.preds2, p)
unique(store.pseudo.preds2$bird.groups)

ps = store.pseudo.preds2


cis = ps %>%
  drop_na()%>%
  group_by(bird.groups)%>%
  summarise(#ci_beta_par1_low = quantile(par1.beta,0.025 ),
            #ci_beta_par1_high = quantile(par1.beta,0.975 ),
            #ci_beta_par2_low = quantile(par2.beta,0.025 ),
            #ci_beta_par2_high = quantile(par2.beta,0.975 ),
            ci_hom_low = quantile(par1.hom,0.025 ),
            ci_hom_high = quantile(par1.hom,0.975 )
            )
###these so wide because alpha and beta depend on eachother so are unreliable!!!!
##that is why we use the ellipse method

library(ellipse)
library(ggplot2)
library(sp)

ps = ps %>%
  drop_na()

write.csv(ps, "figs/all_pseudosets.csv",row.names = FALSE)

#homogeneous values
hom.vals = ps %>%
  filter(bird.groups == "0 va")%>%
  select(par1.hom, bird.groups, lik.hom) %>%
  mutate(ci_hom_low = quantile(par1.hom,0.025 )) %>%
  mutate(ci_hom_high = quantile(par1.hom,0.975 ))%>%
  filter(par1.hom>=ci_hom_low&par1.hom<=ci_hom_high)

write.csv(hom.vals, "figs/va0_CIs_hom_pseudosets.csv",row.names = FALSE)


p750=ggplot(data=subset(ps, bird.groups=="750 va"),aes(x=par1.beta,y=par2.beta))+
  geom_jitter(size=1)+
  stat_ellipse(level=0.95)+
  ylab("beta")+
  xlab("alpha")+
  xlim(0,100)+
  ylim(0,100)+
  theme_bw() + 
  theme(strip.background = element_rect(fill="gray97"),strip.text.x = element_text (size = 15,hjust = 0.5, vjust = 0.5,face="italic"),axis.title=element_text(size=23),axis.text=element_text(size=15),panel.grid = element_blank(), axis.line=element_line(),legend.position="top",legend.text = element_text(size=20,face="italic"),legend.title = element_blank(),legend.background = element_blank(),legend.key=element_rect(fill="white",color="white"))
p750


# Extract components
build <- ggplot_build(p750)$data
points <- build[[1]]
ell <- build[[2]]

# Find which points are inside the ellipse, and add this to the data
dat2 <- data.frame(points[1:2], 
                  in.ell = as.logical(point.in.polygon(points$x, points$y, ell$x, ell$y)))

# Plot the result
ggplot(dat2, aes(x, y)) +
  geom_point(aes(col = in.ell)) +
  stat_ellipse()

dat2
dat3=dat2[dat2$in.ell=="TRUE",]
head(dat3)
nrow(dat3)
dat3$group = "750 va"
dat3 = dat3 %>%
  dplyr::rename(par1 = x, par2 = y)

write.csv(dat3, "figs/va750_50CIs_chisq_pseudosets.csv",row.names = FALSE)


dat3$x = dat3$par1
dat3$y= dat3$par2
###generate beta distribution
num=NA
name=NA
hx2=matrix(NA,ncol=3, nrow=1000)
hx2names=c("num","name","x")
hx2=as.data.frame(hx2)
names(hx2)=hx2names
hx=NA;hx.i.true=NA;hx.i.est=NA;hx.true=NA;hx.est=NA;hx.name.i=NA;hx.name=NA;x.i=NA
hx.x=NA;hx.alpha.st.i=NA;hx.alpha.st=NA;hx.betap.st.i=NA;hx.betap.st=NA;

out=hx2[1,]
#alpha=hom.beta$par[2],betap=hom.beta$par[3]
alpha=0.206258663;betap=0.231546872

#generate distributions
for (i in 1:nrow(dat3)) {
  x.i=seq(0.0001,0.9999,length=300)
  hx.x=c(hx.x,x.i)
  hx.i.est=dbeta(x.i,dat3$x[i],dat3$y[i])
  hx.est=c(hx.est,hx.i.est)
}

hx.est <- hx.est[!is.na(hx.est)]
hx.x <- hx.x[!is.na(hx.x)]
name2=rep("boot",times=length(hx.est))
gen.dist=data.frame(hx.x,hx.est,name2)
names(gen.dist)=c("x","hx","name")

ggplot(data=gen.dist, aes(x=hx.x,y=hx.est ))+
  geom_line()+
  coord_cartesian(ylim=c(0,100))
  
write.csv(gen.dist, "figs/va750_distributionCIs.csv",row.names = FALSE)


### va30000 ###
phigh=ggplot(data=subset(ps, bird.groups=="30000 va"),aes(x=par1.beta,y=par2.beta))+
  geom_jitter(size=1)+
  stat_ellipse(level=0.95)+
  ylab("beta")+
  xlab("alpha")+
  xlim(0,100)+
  ylim(0,100)+
  theme_bw() + 
  theme(strip.background = element_rect(fill="gray97"),strip.text.x = element_text (size = 15,hjust = 0.5, vjust = 0.5,face="italic"),axis.title=element_text(size=23),axis.text=element_text(size=15),panel.grid = element_blank(), axis.line=element_line(),legend.position="top",legend.text = element_text(size=20,face="italic"),legend.title = element_blank(),legend.background = element_blank(),legend.key=element_rect(fill="white",color="white"))
phigh


# Extract components
build <- ggplot_build(phigh)$data
points <- build[[1]]
ell <- build[[2]]

# Find which points are inside the ellipse, and add this to the data
dat2 <- data.frame(points[1:2], 
                   in.ell = as.logical(point.in.polygon(points$x, points$y, ell$x, ell$y)))

# Plot the result
ggplot(dat2, aes(x, y)) +
  geom_point(aes(col = in.ell)) +
  stat_ellipse()

dat2
dat3=dat2[dat2$in.ell=="TRUE",]
head(dat3)
dat3$group = "30000 va"
dat3 = dat3 %>%
  rename(par1 = x, par2 = y)

write.csv(dat3, "figs/va30000_CIs_chisq_pseudosets.csv",row.names = FALSE)

dat3$x = par1
dat3$y = par2
###generate beta distribution
num=NA
name=NA
hx2=matrix(NA,ncol=3, nrow=1000)
hx2names=c("num","name","x")
hx2=as.data.frame(hx2)
names(hx2)=hx2names
hx=NA;hx.i.true=NA;hx.i.est=NA;hx.true=NA;hx.est=NA;hx.name.i=NA;hx.name=NA;x.i=NA
hx.x=NA;hx.alpha.st.i=NA;hx.alpha.st=NA;hx.betap.st.i=NA;hx.betap.st=NA;

out=hx2[1,]
#alpha=hom.beta$par[2],betap=hom.beta$par[3]


#generate distributions
for (i in 1:nrow(dat3)) {
  x.i=seq(0.0001,0.9999,length=300)
  hx.x=c(hx.x,x.i)
  hx.i.est=dbeta(x.i,dat3$x[i],dat3$y[i])
  hx.est=c(hx.est,hx.i.est)
}

hx.est <- hx.est[!is.na(hx.est)]
hx.x <- hx.x[!is.na(hx.x)]
name2=rep("boot",times=length(hx.est))
gen.dist=data.frame(hx.x,hx.est,name2)
names(gen.dist)=c("x","hx","name")

ggplot(data=gen.dist, aes(x=hx.x,y=hx.est ))+
  geom_line()+
  coord_cartesian(ylim=c(0,100))

write.csv(gen.dist, "figs/va30000_distributionCIs.csv",row.names = FALSE)


################### move to other script to make figures ##############




