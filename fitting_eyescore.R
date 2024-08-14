##### Fits eyescore data to determine mortality  #####

# requires the following file:
#     eye_score_primary_dose.csv 
#     eye_score_secondary_dose.csv 

# produces the following file:
#     mortality.csv

# load important packages
library(base)
library(boot)
library(dplyr)
library(Matrix)
library(sp)
library(stats)
library(stringr)
library(tidyverse)
library(utils)

rm(list=ls()) # clear workspace
#dev.off() # close graphical device (if necessary)
cat("\014") # clear console


# load primary dose data
hofi1.bybird = read.csv("eye_score_primary_dose.csv")

hofi1.small= hofi1.bybird %>%
  group_by(population, primary_dose)%>%
  summarise(N=n(), 
            infected = sum(total.infected,na.rm=T),
            mean.eye.score = mean(mean.eye.score, na.rm=T) )%>%
  mutate(experiment="primary")%>%
  arrange(population,primary_dose)

# load secondary dose data
hofi2.bybird = read.csv("eye_score_secondary_dose.csv")

hofi2.small= hofi2.bybird %>%
  group_by(population, primary_dose,secondary_dose)%>%
  summarise(N=n(), 
            infected = sum(total.infected,na.rm=T),
            mean.eye.score = mean(mean.eye.score, na.rm=T) )%>%
  mutate(experiment="secondary")
##----hofi2.small if the secondary infection data for dose-response fittings---###
 

##some of the priming doses weren't used in the secondary challenge#
#so these birds can be treated as naive birds exposed to a large dose
hofi1.small.750 = hofi1.small %>%
  filter(primary_dose==750)%>%
  mutate(secondary_dose = primary_dose)%>%
  mutate(primary_dose = 0)
## Note: this naming may be confusing. We are using secondary dose as the name of the dose for fitting
# Thus, use secondary_dose as the name of the dose column and primary dose to 0
# For clarity, we do have a column indicating which experiment they were part of (primary)

hofi2.less.small = bind_rows(hofi1.small.750, hofi2.small)
#this is the dataset

hofi2.bybird = hofi2.bybird %>%
  filter(population=="va")
hofi.bybird<-hofi2.bybird

# Fitting to eye score data

start.vals2=list(a0=1,a1=1,a2=-1)
hofi.bybird$secondary_dose_scale = hofi.bybird$secondary_dose*0.001
hofi.bybird$primary_dose_scale = hofi.bybird$primary_dose*0.001

glm14b<-nls(mean.eye.score~6*inv.logit(a0+a1*primary_dose_scale+a2*secondary_dose_scale),
            data=hofi.bybird[hofi.bybird$total.infected==1,],start=start.vals2)

summary(glm14b)

primary_dose_scale = unique(hofi.bybird$primary_dose_scale)
secondary_dose_scale= unique(hofi.bybird$secondary_dose_scale)
newdat = expand.grid(primary_dose_scale = primary_dose_scale, secondary_dose_scale = secondary_dose_scale)
newdat$yhat = predict(glm14b, newdata = newdat)

newdat$yhat.r = round(newdat$yhat, 2)
newdat

newdat2 = expand.grid(primary_dose_scale = primary_dose_scale,
                      secondary_dose_scale = seq(min(secondary_dose_scale, na.rm = T), max(secondary_dose_scale, na.rm = T), by=.1))
newdat2$yhat = predict(glm14b, newdata = newdat2)

actual = hofi2.less.small %>%
  filter(population=="va")

actual$primary_dose_scale = actual$primary_dose*0.001
actual$secondary_dose_scale = actual$secondary_dose*0.001

dat.mod.comb = left_join(newdat, actual)
dat.mod.comb

Vmax = 0.25
mus = newdat2 %>%
  group_by(primary_dose_scale)%>%
  filter(secondary_dose_scale==7)%>%
  mutate(mu = (yhat/6)*Vmax)
mus$mu

mortality = subset(mus,select=c(primary_dose_scale,mu))
mortality$dose = primary_dose_scale*1000
mortality = subset(mortality,select=-c(primary_dose_scale))
write.csv(mortality, "mortality.csv",row.names = FALSE)
