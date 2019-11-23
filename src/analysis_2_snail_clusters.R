---
title: "Analysis 2: Number, size, and persistence of snail clusters in space and time"
author: "Skylar Hopkins and Chelsea Wood"
date: "updated 9 July 2019"
---

# We determined the number, size, and persistence of snail clusters in space and time using SatScan software, available for free at:
# https://www.satscan.org
# We then used the output of the SatScan software (site_area_clusters_numquads.csv) to test the relationship between site size and the 
# number of snail clusters. That analysis is shown below.

####################################################################################################################
######################################Load data and packages########################################################
####################################################################################################################
rm(list = ls())

# Load libraries
library(lme4)
library(gamm4)
library(ggplot2)
library(tidyverse)
library(ciTools)
library(tidyverse)
library(dplyr)
library(cowplot)

# citation(package = "MASS")
# citation(package = "stats")

#Load just the cluster and site area data
rawdata<-read.csv("data/site_area_clusters_numquads.csv")
rawdata<-rawdata[,-1]

#Take out spaces to make formatting consistent
rawdata$site<-gsub(" ", "", rawdata$sampling_site, fixed = TRUE)
table(rawdata$site, rawdata$field_mission)

#Rename some columns
names(rawdata)<-c("sampling_site","field_mission","NumQuadsSampled", "AreaSampled","num_bulinus_clusters","area", "site")

#sort by FM and site
rawdata<- rawdata %>%
  arrange(field_mission, site)

sitestoinclude<-unique(rawdata$site)

#We also want to see if total area of other/floating veg is a better predictor than just area
#We need the percent area column in this dataframe to get that
PercentOtherVegData<-read.csv("data/PercentOtherVegData.csv")
head(PercentOtherVegData)
#remove spaces in site name
PercentOtherVegData$Site<-gsub(" ", "", PercentOtherVegData$Site, fixed = TRUE)
PercentOtherVegData<- PercentOtherVegData %>%
  arrange(FM, Site)
PercentOtherVegData<-PercentOtherVegData[PercentOtherVegData$Site %in% sitestoinclude,]
dim(PercentOtherVegData)
#same order, so just merge
View(cbind(rawdata[,1:2], PercentOtherVegData[,1:2]))
rawdata<-as.data.frame(cbind(rawdata, PercentOtherVegData$PercOther))
names(rawdata)<-c(names(rawdata)[1:7], "PercentOtherVeg") #rename area_prawn to just area
rawdata$OtherVegArea<-rawdata$area*(rawdata$PercentOtherVeg/100)

write.csv(rawdata,"data/cluster_data_TESTING.csv",sep=",")

# Load just the cluster and site area data
rawdata<-read.csv("data/cluster_data_TESTING.csv")
rawdata<-rawdata[,-1]
head(rawdata)

# yes, we do need to account for area sampled
plot(rawdata$num_bulinus_clusters~rawdata$NumQuadsSampled)

####################################################################################################################
####################################################################################################################
#####################################Clusters per site per time########################################################
####################################################################################################################
####################################################################################################################

# lmer complains about predictor scale - should we transform? Yes. Log looks good.
hist(rawdata$area)
hist(log(rawdata$area)) #better
hist(rawdata$OtherVegArea)
hist(log(rawdata$OtherVegArea+1)) #better - but what constant?

# Total area - raw doesn't look good
plot<-ggplot(rawdata,aes(x=area,y=num_bulinus_clusters))+
  geom_point(position=position_dodge(0.2),size=4)+
  geom_smooth(method='lm',formula=y~x,color="black")+
  scale_color_manual(values=c("black"))+
  ylab("density of clusters detected per site over six field missions")+
  xlab("site area (m2)")+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
plot

# Other veg area
plot<-ggplot(rawdata,aes(x=OtherVegArea,y=num_bulinus_clusters))+
  geom_point(position=position_dodge(0.2),size=4)+
  geom_smooth(method='lm',formula=y~x,color="black")+
  scale_color_manual(values=c("black"))+
  ylab("density of clusters detected per site over six field missions")+
  xlab("site area (m2)")+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
plot

# log much nicer - go w/ this - but not as nice as it used to look...
# Total area
plot<-ggplot(rawdata,aes(x=log(area),y=num_bulinus_clusters))+
  geom_point(position=position_dodge(0.2),size=4)+
  geom_smooth(method='lm',formula=y~x,color="black")+
  scale_color_manual(values=c("black"))+
  ylab("density of clusters detected per site over six field missions")+
  xlab("log site area (m2)")+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
plot

# other veg area
plot<-ggplot(rawdata,aes(x=log(OtherVegArea+1),y=num_bulinus_clusters))+
  geom_point(position=position_dodge(0.2),size=4)+
  geom_smooth(method='lm',formula=y~x,color="black")+
  scale_color_manual(values=c("black"))+
  ylab("density of clusters detected per site over six field missions")+
  xlab("log site area (m2)")+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
plot

rawdata$log_area<-log(rawdata$area)
rawdata$log_veg_area<-log(rawdata$OtherVegArea+1)

# are they correlatetd?
plot(rawdata$log_veg_area~rawdata$log_area) #yes
plot(rawdata$PercentOtherVeg~rawdata$log_area) #no
plot(rawdata$log_veg_area~rawdata$NumQuadsSampled) #maybe
plot(rawdata$log_area~rawdata$NumQuadsSampled) #yes

####################################################################################################################
###############Area models - per site per field mission########################################################
####################################################################################################################

rawdata$time<-as.numeric(rawdata$field_mission)

# Gaussian model - theoretically inappropriate because we're working w/ zero truncated integers
model<-lmer(num_bulinus_clusters~log_area+(1|site), data=rawdata)
plot(residuals(model, type="pearson")~rawdata$time)
acf(residuals(model, type="pearson"))
pacf(residuals(model, type="pearson"))
summary(model)
hist(residuals(model, type="pearson"))
plot(residuals(model, type="pearson")~predict(model)) #uhhh
BIC(model) #184.03
AIC(model)
#AIC = 175.0618

# Poisson
model<-glmer(num_bulinus_clusters~log_area+(1|site), data=rawdata, family="poisson"); summary(model)
plot(residuals(model, type="pearson")~rawdata$time)
AIC(model) #166.9821
BIC(model) #162.98

# Poisson w/ offset
model<-glmer(num_bulinus_clusters~log_area+offset(log(NumQuadsSampled*0.32))+(1|site), data=rawdata, family="poisson"); summary(model)
plot(residuals(model, type="pearson")~rawdata$time)
AIC(model) #166.0014
BIC(model) #170.5008


# There are some sites with consistently low resids, but acf looks OK
c=ggplot(rawdata) +
  geom_point(aes(x=time, y=residuals(model, type="pearson"), col=site, size=6.5, alpha=0.5))+ 
  geom_line(aes(x=time, y=residuals(model, type="pearson"), col=site))+
  ylab("residual")+
  xlab("time (field mission")+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
c

acf(residuals(model, type="pearson"))
pacf(residuals(model, type="pearson"))
summary(model)
hist(residuals(model, type="pearson"))
plot(residuals(model, type="pearson")~predict(model)) 
BIC(model) #170.50 

# newdata<- expand.grid(log_area = seq(from=min(rawdata$log_area), to=max(rawdata$log_area), by=.1), site=unique(rawdata$site))
newdata<- expand.grid(log_area = seq(from=min(rawdata$log_area), to=8.7, by=.1), NumQuadsSampled=15)
newdata$phat <- predict(model, newdata, type="response", re.form=NA) #no random effects

c=ggplot(rawdata, aes(x=log_area, y=num_bulinus_clusters)) +
  geom_jitter(height = 0.03, aes(size=6.5, alpha=0.6, col=site))+ 
  #geom_jitter(height = 0.05, aes(size=6.5, alpha=0.7, col=as.numeric(as.factor(site))))+ 
  geom_line(data=newdata,aes(x=log_area,y=phat), size=1, linetype="dashed")+
  scale_color_manual(values=c("#FF0000", "#FF7F50", "#00FF00", "#008000", "#00FFFF", "#008080", "#00BFFF", "#FFFF00", "#4169E1", "#EE82EE", "#BA55D3", "#8B008B", "#FF69B4", "#2F4F4F", "#800000")) +
  scale_y_continuous(breaks=c(0,1,2)) +
  #scale_color_gradientn(colours = rainbow(15)) +
  ylab(expression(Number~of~clusters~(site^-1~time^-1)))+
  xlab(expression(Ln~site~area~(m^2)))+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
c

# (still poisson) if we take out the log(areas)<~6, does pattern persist? Yes.
rawdata$site[rawdata$log_area<5.5]
model<-glmer(num_bulinus_clusters~log_area+offset(log(NumQuadsSampled*0.32))+(1|site), data=rawdata[rawdata$log_area>5.5,], family="poisson", glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(model)
plot(residuals(model, type="pearson")~rawdata$time[rawdata$log_area>5.5])

f=ggplot(rawdata[rawdata$log_area>5.5,]) +
  geom_point(aes(x=time, y=residuals(model, type="pearson"), col=site, size=6.5, alpha=0.5))+ 
  geom_line(aes(x=time, y=residuals(model, type="pearson"), col=site))+
  ylab("residual")+
  xlab("time (field mission")+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
f

newdata<- expand.grid(log_area = seq(from=min(rawdata$log_area[rawdata$log_area>5.5]), to=max(rawdata$log_area[rawdata$log_area>5.5]), by=.1), site=unique(rawdata$site[rawdata$log_area>5.5]))
newdata<- expand.grid(log_area = seq(from=min(rawdata$log_area[rawdata$log_area>5.5]), to=max(rawdata$log_area[rawdata$log_area>5.5]), by=.1), NumQuadsSampled=15)
newdata$phat <- predict(model, newdata, type="response", re.form=NA) #no random effects

f=ggplot(rawdata[rawdata$log_area>5.5,], aes(x=log_area, y=num_bulinus_clusters)) +
  geom_point(aes(size=6.5, alpha=0.5, col=site))+ 
  geom_line(data=newdata,aes(x=log_area,y=phat),linetype="dashed")+
  ylab("density of clusters detected per site per time")+
  xlab("log site area (m2)")+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
f

# negative binomial - convergence issues
model<-glmer.nb(num_bulinus_clusters~log_area+offset(log(NumQuadsSampled))+(1|site), data=rawdata, glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))
summary(model)
plot(residuals(model, type="pearson")~rawdata$time)
BIC(model)

q=ggplot(rawdata) +
  geom_point(aes(x=time, y=residuals(model, type="pearson"), col=site, size=6.5, alpha=0.5))+ 
  geom_line(aes(x=time, y=residuals(model, type="pearson"), col=site))+
  ylab("residual")+
  xlab("time (field mission")+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
q

acf(residuals(model, type="pearson"))
pacf(residuals(model, type="pearson"))
summary(model)
hist(residuals(model, type="pearson"))
plot(residuals(model, type="pearson")~predict(model)) 
AIC(model)
#AIC = 168.3445

newdata<- expand.grid(log_area = seq(from=min(rawdata$log_area), to=max(rawdata$log_area), by=.1), site=unique(rawdata$site))
newdata<- expand.grid(log_area = seq(from=min(rawdata$log_area), to=max(rawdata$log_area), by=.1), NumQuadsSampled=15)
newdata$phat <- predict(model, newdata, type="response", re.form=NA) #no random effects

q=ggplot(rawdata, aes(x=log_area, y=num_bulinus_clusters)) +
  geom_point(aes(size=6.5, alpha=0.5, col=site))+ 
  geom_line(data=newdata,aes(x=log_area,y=phat),linetype="dashed")+
  ylab("density of clusters detected per site per time")+
  xlab("log site area (m2)")+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
q

####################################################################################################################
##################Other/Floating Veg Area models - per site per field mission########################################################
####################################################################################################################
#Poisson
model<-glmer(num_bulinus_clusters~log_veg_area+(1|site), data=rawdata, family="poisson"); summary(model)
plot(residuals(model, type="pearson")~rawdata$time)
AIC(model) #150.0294
BIC(model) #157.53

#Poisson + offset
model<-glmer(num_bulinus_clusters~log_veg_area+offset(log(NumQuadsSampled*0.32))+(1|site), data=rawdata, family="poisson")
plot(residuals(model, type="pearson")~rawdata$time)
AIC(model) #150.12
BIC(model) #150.12

#There are some sites with consistently low resids, but acf looks OK
g=ggplot(rawdata) +
  geom_point(aes(x=time, y=residuals(model, type="pearson"), col=site, size=6.5, alpha=0.5))+ 
  geom_line(aes(x=time, y=residuals(model, type="pearson"), col=site))+
  ylab("residual")+
  xlab("time (field mission")+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
g

acf(residuals(model, type="pearson"))
pacf(residuals(model, type="pearson"))
summary(model)
hist(residuals(model, type="pearson"))
plot(residuals(model, type="pearson")~predict(model)) 

#newdata<- expand.grid(log_area = seq(from=min(rawdata$log_area), to=max(rawdata$log_area), by=.1), site=unique(rawdata$site))
newdata<- expand.grid(log_veg_area = seq(from=min(rawdata$log_veg_area), to=7.5, by=.1), NumQuadsSampled=15)
newdata$phat <- predict(model, newdata, type="response", re.form=NA) #no random effects

d=ggplot(rawdata, aes(x=log_veg_area, y=num_bulinus_clusters)) +
  geom_jitter(height = 0.03, aes(size=6.5, alpha=0.6, col=site))+ 
  #geom_jitter(height = 0.05, aes(size=6.5, alpha=0.7, col=as.numeric(as.factor(site))))+ 
  geom_line(data=newdata,aes(x=log_veg_area,y=phat), size=1, linetype="dashed")+
  scale_color_manual(values=c("#FF0000", "#FF7F50", "#00FF00", "#008000", "#00FFFF", "#008080", "#00BFFF", "#FFFF00", "#4169E1", "#EE82EE", "#BA55D3", "#8B008B", "#FF69B4", "#2F4F4F", "#800000")) +
  scale_y_continuous(breaks=c(0,1,2)) +
  scale_x_continuous(limits = c(0, 8), breaks=seq(0, 8, 1)) +
  #scale_color_gradientn(colours = rainbow(15)) +
  ylab(expression(Number~of~clusters~(site^-1~time^-1)))+
  xlab(expression(Ln~vegetated~area~(m^2)))+ #it's really vegetated area+100
  theme(plot.title=element_text(size=10),axis.text.y=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
d

#stick the supplemental plots together
plot_grid(c, d, labels="auto")

test<-c(5,6,7,8)
exp(test)

####################################################################################################################
####################################################################################################################
##############################Combine by site across time points#####################################
####################################################################################################################
####################################################################################################################
#taking average site sizes seems reasonable
plot(area~as.factor(site), data=rawdata)

#taking average area sampled does NOT seem reasonable - should sum instead
ggplot(rawdata, aes(x=as.numeric(as.factor(site)), y=NumQuadsSampled)) +
  geom_jitter(size=2, alpha=0.5, width = 0.1, height=0)
table(rawdata$NumQuadsSampled)

AggData<-rawdata %>%
  group_by(site) %>%
  #summarise(meanarea = mean(area, na.rm = T))
  summarise(meanarea = mean(area), meanothervegarea=mean(OtherVegArea), totalarea=sum(area), totalotherarea=sum(OtherVegArea), totalclusters=sum(num_bulinus_clusters), totalareasampled=sum(NumQuadsSampled))

View(AggData)

#this doesn't seem as reasonable - so we sum instead of using means
plot(OtherVegArea~as.factor(site), data=rawdata)
points(meanothervegarea~as.factor(site), data=AggData, col="red")

#log area?
hist(AggData$meanarea)
hist(log(AggData$meanarea))
hist(AggData$meanothervegarea)
hist(log(AggData$meanothervegarea)) #HIP either way
hist(AggData$totalarea)
hist(log(AggData$totalarea))
hist(AggData$totalotherarea)
hist(log(AggData$totalotherarea))

#TOTAL AREA
##Poisson - raw total area + offset is log area sampled
test<-glm(totalclusters~totalarea, data=AggData, family="poisson"); summary(test)
test<-glm(totalclusters~totalarea + offset(log(log(totalareasampled*0.32))), data=AggData, family="poisson"); summary(test)
plot(residuals(test, type="pearson"), ylim=c(-3, 3))
hist(residuals(test, type="pearson"))
plot(residuals(test, type="pearson")~predict(test)) 
AIC(test)
BIC(test) #74.72
#AIC = 73.38534 to 73.30674 w/ offset

#newdata<- expand.grid(log_area = seq(from=min(rawdata$log_area), to=max(rawdata$log_area), by=.1), site=unique(rawdata$site))
newdata<- expand.grid(totalarea = seq(from=min(AggData$totalarea), to=max(AggData$totalarea), by=.1), totalareasampled=(15*6))
newdata$phat <- predict(test, newdata, type="response", re.form=NA) #no random effects

c=ggplot(AggData, aes(x=totalarea, y=totalclusters)) +
  geom_jitter(height = 0.03, aes(size=6.5, alpha=0.6, col=site))+ 
  #geom_jitter(height = 0.05, aes(size=6.5, alpha=0.7, col=as.numeric(as.factor(site))))+ 
  geom_line(data=newdata,aes(x=totalarea,y=phat), size=1, linetype="dashed")+
  #scale_color_manual(values=c("#FF0000", "#FF7F50", "#00FF00", "#008000", "#00FFFF", "#008080", "#00BFFF", "#FFFF00", "#4169E1", "#EE82EE", "#BA55D3", "#8B008B", "#FF69B4", "#2F4F4F", "#800000")) +
  #scale_y_continuous(breaks=c(0,1,2)) +
  #scale_color_gradientn(colours = rainbow(15)) +
  ylab("density of clusters detected per site per time")+
  xlab("Ln site area (m2)")+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
c

#LOG TOTAL AREA + offset is log area sampled
AggData$logtotalarea<-log(AggData$totalarea)
test<-glm(totalclusters~logtotalarea, data=AggData, family="poisson"); summary(test)
test<-glm(totalclusters~logtotalarea+offset(log(log(totalareasampled*0.32))), data=AggData, family="poisson"); summary(test)
plot(residuals(test, type="pearson"), ylim=c(-3, 3))
hist(residuals(test, type="pearson"))
plot(residuals(test, type="pearson")~predict(test)) 
AIC(test)
AICc(test)
BIC(test) #70.01
#AIC = 68.59614; AICc = 69.59614

#newdata<- expand.grid(log_area = seq(from=min(rawdata$log_area), to=max(rawdata$log_area), by=.1), site=unique(rawdata$site))
newdata<- expand.grid(logtotalarea = seq(from=6, to=11, by=.1), totalareasampled=(15*6))
newdata$phat <- predict(test, newdata, type="response", re.form=NA) #no random effects
newdata <- add_ci(newdata, test, names = c("lwr", "upr"), alpha = 0.05) %>%
  mutate(type = "parametric")
newdata$upr[newdata$upr>9]<-9 #for plotting purposes

c=ggplot(AggData, aes(x=logtotalarea)) +
  geom_ribbon(data = newdata, aes(x=logtotalarea, ymin = lwr, ymax = upr), alpha=.2) +
  geom_line(data=newdata,aes(x=logtotalarea,y=phat), size=1, col="gray50", linetype="dashed")+
  geom_jitter(data=AggData, height = 0.05, width=0.05, aes(size=6.5, y=totalclusters), col="black", alpha=0.7)+ 
  #scale_color_manual(values=c("#FF0000", "#FF7F50", "#00FF00", "#008000", "#00FFFF", "#008080", "#00BFFF", "#FFFF00", "#4169E1", "#EE82EE", "#BA55D3", "#8B008B", "#FF69B4", "#2F4F4F", "#800000")) +
  scale_x_continuous(limits=c(6, 11), breaks=seq(6, 11, 1)) +
  scale_y_continuous(limits = c(NA, 9), breaks=seq(0, 9, 1)) +
  #scale_color_gradientn(colours = rainbow(15)) +
  ylab(expression(Cluster~density~(site^-1)))+
  xlab(expression(Ln~summed~site~area~(m^2)))+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
c

#TOTAL VEG AREA - by total with offset as log area sampled
AggData$logtotalvegarea<-log(AggData$totalotherarea)
test<-glm(totalclusters~logtotalvegarea+offset(log(log(totalareasampled*0.32))), data=AggData, family="poisson"); summary(test)
plot(residuals(test, type="pearson"), ylim=c(-3, 3))
hist(residuals(test, type="pearson"))
plot(residuals(test, type="pearson")~predict(test)) 
AIC(test)
AICc(test)
BIC(test) #59.59
#AIC = 58.16935; AICc = 59.16935

#newdata<- expand.grid(log_area = seq(from=min(rawdata$log_area), to=max(rawdata$log_area), by=.1), site=unique(rawdata$site))
newdata<- expand.grid(logtotalvegarea = seq(from=min(AggData$logtotalvegarea), to=9.5, by=.1), totalareasampled=(15*16))
newdata$phat <- predict(test, newdata, type="response", re.form=NA) #no random effects
newdata <- add_ci(newdata, test, names = c("lwr", "upr"), alpha = 0.05) %>%
  mutate(type = "parametric")
newdata$upr[newdata$upr>9]<-9 #for plotting purposes

d=ggplot(AggData, aes(x=logtotalvegarea)) +
  geom_ribbon(data = newdata, aes(x=logtotalvegarea, ymin = lwr, ymax = upr), alpha=.2) +
  geom_line(data=newdata,aes(x=logtotalvegarea,y=phat), size=1, col="gray50", linetype="dashed")+
  geom_jitter(data=AggData, width=0.05, height = 0.05, aes(size=6.5, y=totalclusters), col="black", alpha=0.7)+ 
  #scale_color_manual(values=c("#FF0000", "#FF7F50", "#00FF00", "#008000", "#00FFFF", "#008080", "#00BFFF", "#FFFF00", "#4169E1", "#EE82EE", "#BA55D3", "#8B008B", "#FF69B4", "#2F4F4F", "#800000")) +
  scale_x_continuous(breaks=seq(0, 9, 1)) +
  scale_y_continuous(limits = c(NA, 9), breaks=seq(0, 9, 1)) +
  #scale_color_gradientn(colours = rainbow(15)) +
  ylab(expression(Cluster~density~(site^-1)))+
  xlab(expression(Ln~summed~vegetated~area~(m^2~site^-1)))+
  theme(plot.title=element_text(size=10),axis.text.y=element_blank(),axis.title.y=element_blank(),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
d

###stick plots together for pub
plot_grid(c, d, labels="auto")

####################
#What if we only look at sites with exactly 15 quads? Same
AggData2<-AggData[AggData$totalareasampled==(15*6),]
table(AggData2$site)
test<-glm(totalclusters~logtotalvegarea, data=AggData2, family="poisson"); summary(test)

c=ggplot(AggData2, aes(x=logtotalvegarea)) +
  geom_jitter(data=AggData2, height = 0.0, aes(size=6.5, y=totalclusters), col="black")+ 
  #scale_color_manual(values=c("#FF0000", "#FF7F50", "#00FF00", "#008000", "#00FFFF", "#008080", "#00BFFF", "#FFFF00", "#4169E1", "#EE82EE", "#BA55D3", "#8B008B", "#FF69B4", "#2F4F4F", "#800000")) +
  scale_x_continuous(breaks=seq(5, 14, 1)) +
  scale_y_continuous(limits = c(NA, 9), breaks=seq(0, 9, 1)) +
  #scale_color_gradientn(colours = rainbow(15)) +
  ylab(expression(Total~clusters~(site^-1)))+
  xlab(expression(Ln~mean~area~of~floating~vegetation~(m^2~site^-1)))+
  theme(plot.title=element_text(size=10),axis.text.y=element_text(size=10),axis.title.y=element_text(size=10),axis.text.x=element_text(size=10),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_blank(),panel.grid.minor=element_blank(),plot.margin=unit(c(0,0,0,0),"cm"))+
  theme(legend.position="none")
c
