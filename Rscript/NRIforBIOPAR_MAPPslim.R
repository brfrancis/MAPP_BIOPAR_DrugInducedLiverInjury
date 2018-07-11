library(aod)
#library(ggplot2)
library(rms)
library(Hmisc)
library(PredictABEL)
library(hsdar)
library(QuantPsyc)
library(lsr)
library(ROCR)
plotfull = function(ODlog){
  s = seq(.01,.99,length=1000)
  OUT = matrix(0,1000,4)
  for(i in 1:1000) OUT[i,]=perf(s[i],ODlog,ODlog$y)
  plot(s,OUT[,1],xlab="Cutoff",ylab="Value",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),type="l",lwd=2,axes=FALSE,col=2)
  axis(1,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
  axis(2,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
  lines(s,OUT[,2],col="darkgreen",lwd=2)
  lines(s,OUT[,3],col=4,lwd=2)
  lines(s,OUT[,4],col="darkred",lwd=2)
  box()
  legend(.1,.6,col=c(2,"darkgreen",4,"darkred"),lwd=c(2,2,2,2),c("Sensitivity","Specificity","Classification Rate","Distance"))
  ans=apply(OUT[OUT[,4]==min(OUT[,4]),],2,min)
  return(ans)
}
perf = function(cut, mod, y,nd)
{
  yhat = (predict(mod,newdata=nd,type="response")>cut)
  w = which(y==1)
  sensitivity = mean( yhat[w] == 1 ) 
  specificity = mean( yhat[-w] == 0 ) 
  c.rate = mean( y==yhat ) 
  d = cbind(sensitivity,specificity)-c(1,1)
  d = sqrt( d[1]^2 + d[2]^2 ) 
  tru=sum(y[y==1]==(yhat[y==1]>cut))
  fal=sum(y[y==0]==(yhat[y==0]>cut))
  nt=length(y[y==1])
  nf=length(y[y==0])
  out = t(as.matrix(c(cut,sensitivity, specificity, c.rate,d,tru,nt,fal,nf)))
  colnames(out) = c("probability","sensitivity", "specificity", "c.rate", "distance","correct_t","n_t","correct_f","n_f")
  return(out)
}

perfcut = function(mod,y,nd,yv,ndv){
  probs=as.matrix(seq(0.01,0.99,0.01))
  tab=apply(probs,1,function(x){perf(x,mod,y,nd)})
  mintab=min(tab[5,])
  if(sum(tab[5,]==min(tab[5,]))>1){
  res=tab[,tab[5,]==min(tab[5,])][,1]
  prob=min(probs[tab[5,]==min(tab[5,])])}
  if(sum(tab[5,]==min(tab[5,]))==1){
    res=tab[,tab[5,]==min(tab[5,])]
    prob=probs[tab[5,]==min(tab[5,])]}
  res=perf(prob,mod,y,nd)
  resv=perf(prob,mod,yv,ndv)
  return(rbind(res,resv))
  }

compCalc = function(ODlogitN,ODlogitO){
  new=predict(ODlogitN,type="response")
  old=predict(ODlogitO,type="response")
  comp=improveProb(old,new,ODlogitN$y)
  return(comp)  
}

#Sample size calculation

library(pwr)
#sensitivity of ALT: 0.52
#lowest sensitivity of novels: 0.56 
#highest sens of novels: 0.84

h1 <- ES.h(0.63, 0.52)
pwr.p.test(h = h1, n = NULL, sig.level = 0.05, power = 0.9, alt = "greater")

h1 <- ES.h(0.84, 0.52)
pwr.p.test(h = h1, n = NULL, sig.level = 0.05, power = 0.9, alt = "greater")

library(pROC)

probs_dist=sapply(probs,function(x){sqrt(((sum(new[data3$LateALI==1]>x)/length(new[data3$LateALI==1])-1)^2)+((sum(new[data3$LateALI==0]<x)/length(new[data3$LateALI==0])-1)^2))})
cut=median(probs[probs_dist==min(sapply(probs,function(x){sqrt(((sum(new[data3$LateALI==1]>x)/length(new[data3$LateALI==1])-1)^2)+((sum(new[data3$LateALI==0]<x)/length(new[data3$LateALI==0])-1)^2))}))])
sum(new[data3$LateALI==1]>cut)
sum(new[data3$LateALI==0]<cut)

expit = function(x) exp(predict(x))/(1+exp(predict(x)))
dataOD_Combo<-read.csv("/home/ben/Dropbox/CDSS/Dan/COMBO/COMBO.csv", header=TRUE, sep=",", na.strings=c("","NA"))
dataOD_Combo = read.csv("C://Users/brfra/Dropbox/CDSS/Dan/COMBO/COMBO.csv", header=TRUE, sep=",", na.strings=c("","NA"))

dataOD<-dataOD_Combo[substr(dataOD_Combo$ID,1,1)=="M"& !is.na(dataOD_Combo$LateALI),]
dataOD$SexN<-ifelse(dataOD$Sex=="F",1,0)
dataOD$StaggeredN<-ifelse(dataOD$Staggered=="Y",1,0)
colSums(is.na(dataOD))

#Set up data
data3<-na.omit(cbind(dataOD$LateALI,dataOD$ALT,dataOD$HMGB1,dataOD$necrosis_K18,dataOD$apoptosis_K18,dataOD$miR122_Let_7d,dataOD$GLDH))
data3=as.data.frame(data3)
colnames(data3)<-c("LateALI","ALT","HMGB1","necrosis_K18","apoptosis_K18","miR122_Let_7d","GLDH")

#Fit multiple regression model
ODlogit1 <- glm(data3$LateALI ~ scale(data3$ALT,center=TRUE) + scale(data3$HMGB1,center=TRUE) + scale(data3$necrosis_K18,center=TRUE) + scale(data3$apoptosis_K18,center=TRUE) + scale(data3$miR122_Let_7d,center=TRUE) + scale(data3$GLDH,center=TRUE), family = "binomial")
summary(ODlogit1)
ODlogitNothing <- glm(data3$LateALI ~  scale(data3$ALT,center=TRUE), data = data3, family = "binomial")
summary(ODlogitNothing)
both = step(ODlogit1,direction="both",scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)))
back = step(ODlogit1,direction="backward",scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)))
forward = step(ODlogitNothing,scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)),direction="forward")
formula(both)
formula(back)
formula(forward)

#NRI analysis in MAPP
data3<-na.omit(cbind(dataOD$LateALI,scale(dataOD$ALT,scale=TRUE,center = TRUE),scale(dataOD$HMGB1,scale=TRUE,center = TRUE),scale(dataOD$necrosis_K18,scale=TRUE,center = TRUE),scale(dataOD$apoptosis_K18,scale=TRUE,center = TRUE),scale(dataOD$miR122_Let_7d,scale=TRUE,center = TRUE),scale(dataOD$GLDH,scale=TRUE,center = TRUE)))
colnames(data3)<-c("LateALI","ALT","HMGB1","necrosis_K18","apoptosis_K18","miR122_Let_7d","GLDH")
data3=as.data.frame(data3)
ODlogit1 <- glm(LateALI ~ ALT+HMGB1+necrosis_K18+apoptosis_K18+miR122_Let_7d,data=data3, family = "binomial")
summary(ODlogit1)
ODlogit2 <- glm(LateALI ~ ALT+miR122_Let_7d,data=data3, family = "binomial")
summary(ODlogit2)
ODlogit3 <- glm(LateALI ~ ALT+apoptosis_K18+miR122_Let_7d,data=data3, family = "binomial")
summary(ODlogit3)
ODlogit4 <- glm(LateALI ~ ALT+necrosis_K18+apoptosis_K18+miR122_Let_7d,data=data3, family = "binomial")
summary(ODlogit4)
ODlogitNothing <- glm(LateALI ~ ALT,data=data3, family = "binomial")
summary(ODlogitNothing)


compCalc(ODlogitN = ODlogit2,ODlogitO = ODlogitNothing)#mir 
compCalc(ODlogitN = ODlogit3,ODlogitO = ODlogitNothing)#+necrosis
compCalc(ODlogitN = ODlogit4,ODlogitO = ODlogitNothing)#+apoptosis
compCalc(ODlogitN = ODlogit1,ODlogitO = ODlogitNothing)#+HMGB1

##### For new study sample size #####

roc1<-roc(data3$LateALI, data3$miR122_Let_7d,ci=TRUE)

#BIOPAR data
dataOD_B<-dataOD_Combo[substr(dataOD_Combo$ID,1,1)=="B"& !is.na(dataOD_Combo$LateALI),]
colSums(is.na(dataOD_B))
data4<-na.omit(cbind(dataOD_B$LateALI,scale(dataOD_B$ALT,scale=TRUE,center = TRUE),scale(dataOD_B$HMGB1,scale=TRUE,center = TRUE),scale(dataOD_B$necrosis_K18,scale=TRUE,center = TRUE),scale(dataOD_B$apoptosis_K18,scale=TRUE,center = TRUE),scale(dataOD_B$miR122_Let_7d,scale=TRUE,center = TRUE),scale(dataOD_B$GLDH,scale=TRUE,center = TRUE)))
data4=as.data.frame(data4)
colnames(data4)<-c("LateALI","ALT","HMGB1","necrosis_K18","apoptosis_K18","miR122_Let_7d","GLDH")

#NRI analysis for BIOPAR
valid3=predict.glm(ODlogit2,newdata=data4,type="response")
valid4=predict.glm(ODlogit3,newdata=data4,type="response")
valid5=predict.glm(ODlogit4,newdata=data4,type="response")
valid1=predict.glm(ODlogit1,newdata=data4,type="response")
valid2=predict.glm(ODlogitNothing,newdata=data4,type="response")
valid3_comp=improveProb(valid2,valid3,data4$LateALI)
valid4_comp=improveProb(valid2,valid4,data4$LateALI)
valid5_comp=improveProb(valid2,valid5,data4$LateALI)
valid1_comp=improveProb(valid2,valid1,data4$LateALI)
valid3_comp
valid4_comp
valid5_comp
valid1_comp

#Table of Correctly Identified

perfcut(ODlogitNothing,data3$LateALI,data3,data4$LateALI,data4)
perfcut(ODlogit2,data3$LateALI,data3,data4$LateALI,data4)
perfcut(ODlogit3,data3$LateALI,data3,data4$LateALI,data4)
perfcut(ODlogit4,data3$LateALI,data3,data4$LateALI,data4)
perfcut(ODlogit1,data3$LateALI,data3,data4$LateALI,data4)

#Plot
png("C://Users/Ben/Dropbox/CDSS/Dan/LateALI_MultivariablePlots.png", height=900,width=1200,res=120)
layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
par(mai=rep(0.5, 4))
new=predict(ODlogit1,type="response")
old=predict(ODlogitNothing,type="response")
probs=as.matrix(seq(0,1,0.01))
plot(probs,sapply(probs,function(x){sum(new[data3$LateALI==1]>x)/length(new[data3$LateALI==1])}),type="l",xlab="Calculated Risk",ylab="Sensitivity (Black) 1-Specificity (Red)",main="MAPP dataset",ylim=c(0,1),xlim=c(0,1),lwd=2)
lines(probs,sapply(probs,function(x){sum(new[data3$LateALI==0]>x)/length(new[data3$LateALI==0])}),col="red",lwd=2)
new=predict(ODlogit2,type="response")
lines(probs,sapply(probs,function(x){sum(new[data3$LateALI==1]>x)/length(new[data3$LateALI==1])}),lty=2,type="l",xlab="Calculated Risk",ylab="Sensitivity (Black) 1-Specificity (Red)",main="MAPP dataset",ylim=c(0,1),xlim=c(0,1),lwd=2)
lines(probs,sapply(probs,function(x){sum(new[data3$LateALI==0]>x)/length(new[data3$LateALI==0])}),lty=2,col="red",lwd=2)
lines(probs,sapply(probs,function(x){sum(old[data3$LateALI==1]>x)/length(old[data3$LateALI==1])}),lty=3,type="l",xlab="Calculated Risk",ylab="Sensitivity (Black) 1-Specificity (Red)",main="MAPP dataset",ylim=c(0,1),xlim=c(0,1),lwd=2)
lines(probs,sapply(probs,function(x){sum(old[data3$LateALI==0]>x)/length(old[data3$LateALI==0])}),lty=3,col="red",lwd=2)

new=predict(ODlogit1,newdata=data4,type="response")
old=predict(ODlogitNothing,newdata=data4,type="response")
probs=as.matrix(seq(0,1,0.01))
plot(probs,sapply(probs,function(x){sum(new[data4$LateALI==1]>x)/length(new[data4$LateALI==1])}),type="l",xlab="Calculated Risk",ylab="Sensitivity (Black) 1-Specificity (Red)",main="BIOPAR dataset",ylim=c(0,1),xlim=c(0,1),lwd=2)
lines(probs,sapply(probs,function(x){sum(new[data4$LateALI==0]>x)/length(new[data4$LateALI==0])}),col="red",lwd=2)
new=predict(ODlogit2,newdata=data4,type="response")
lines(probs,sapply(probs,function(x){sum(new[data4$LateALI==1]>x)/length(new[data4$LateALI==1])}),lty=2,type="l",xlab="Calculated Risk",ylab="Sensitivity (Black) 1-Specificity (Red)",main="MAPP dataset",ylim=c(0,1),xlim=c(0,1),lwd=2)
lines(probs,sapply(probs,function(x){sum(new[data4$LateALI==0]>x)/length(new[data4$LateALI==0])}),lty=2,col="red",lwd=2)
lines(probs,sapply(probs,function(x){sum(old[data4$LateALI==1]>x)/length(old[data4$LateALI==1])}),lty=3,type="l",xlab="Calculated Risk",ylab="Sensitivity (Black) 1-Specificity (Red)",main="MAPP dataset",ylim=c(0,1),xlim=c(0,1),lwd=2)
lines(probs,sapply(probs,function(x){sum(old[data4$LateALI==0]>x)/length(old[data4$LateALI==0])}),lty=3,col="red",lwd=2)

par(mai=c(0,0,0,0))
plot.new()
legend("center",legend=c("Full","miR-122 Let 7d","ALT only"),lty=c(1,2,3),lwd=c(2,2,2))
dev.off()



#SevereALI work
dataOD<-dataOD_Combo[substr(dataOD_Combo$ID,1,1)=="M"& !is.na(dataOD_Combo$LateALI),]
dataOD$SexN<-ifelse(dataOD$Sex=="F",1,0)
dataOD$StaggeredN<-ifelse(dataOD$Staggered=="Y",1,0)
colSums(is.na(dataOD))

#Set up data
data3<-na.omit(cbind(dataOD$SLALI,scale(dataOD$ALT,scale=TRUE,center = TRUE),scale(dataOD$HMGB1,scale=TRUE,center = TRUE),scale(dataOD$necrosis_K18,scale=TRUE,center = TRUE),scale(dataOD$apoptosis_K18,scale=TRUE,center = TRUE),scale(dataOD$miR122_Let_7d,scale=TRUE,center = TRUE),scale(dataOD$GLDH,scale=TRUE,center = TRUE)))
colnames(data3)<-c("SLALI","ALT","HMGB1","necrosis_K18","apoptosis_K18","miR122_Let_7d","GLDH")
data3=as.data.frame(data3)

#Fit multiple regression model
ODlogit1 <- glm(data3$SLALI ~ data3$ALT+data3$HMGB1+data3$necrosis_K18+data3$apoptosis_K18+data3$miR122_Let_7d+data3$GLDH, family = "binomial")
summary(ODlogit1)
ODlogitNothing <- glm(data3$SLALI ~  1, data = data3, family = "binomial")
summary(ODlogitNothing)
both = step(ODlogit1,direction="both",scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)))
back = step(ODlogit1,direction="backward",scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)))
forward = step(ODlogitNothing,scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)),direction="forward")
formula(both)
formula(back)
formula(forward)
data3<-na.omit(cbind(dataOD$SLALI,scale(dataOD$HMGB1,scale=TRUE,center = TRUE)))
data3=as.data.frame(data3)
colnames(data3)<-c("SLALI","HMGB1")
ODlogit1 <- glm(SLALI ~ HMGB1, data=data3, family = "binomial")
summary(ODlogit1)
compCalc(ODlogitN = ODlogit1,ODlogitO = ODlogitNothing)

##### For new study sample size #####

roc2<-roc(data3$SLALI, data3$HMGB1,ci=TRUE)

#Validation
dataOD_B<-dataOD_Combo[substr(dataOD_Combo$ID,1,1)=="B"& !is.na(dataOD_Combo$LateALI),]
colSums(is.na(dataOD_B))
data4<-na.omit(cbind(dataOD_B$SLALI,scale(dataOD_B$HMGB1,scale=TRUE,center = TRUE)))
data4=as.data.frame(data4)
colnames(data4)<-c("SLALI","HMGB1")

#NRI analysis for BIOPAR
valid1=predict.glm(ODlogit1,newdata=data4,type="response")
valid2=predict.glm(ODlogitNothing,newdata=data4,type="response")
valid1_comp=improveProb(valid2,valid1,data4$SLALI)

#Proceed with HMGB1 as indicator of Severe Late ALI
perfcut(ODlogit1,data3$SLALI,data3,data4$SLALI,data4)

#Plot
png("C://Users/Ben/Dropbox/CDSS/Dan/SLALI_MultivariablePlots.png", height=900,width=1200,res=120)
layout(matrix(c(1,2,3,3), ncol=2, byrow=TRUE), heights=c(4, 1))
par(mai=rep(0.5, 4))
new=predict(ODlogit1,type="response")
old=predict(ODlogitNothing,type="response")
probs=as.matrix(seq(0,1,0.01))
plot(probs,sapply(probs,function(x){sum(new[data3$SLALI==1]>x)/length(new[data3$SLALI==1])}),type="l",xlab="Calculated Risk",ylab="Sensitivity (Black) 1-Specificity (Red)",main="MAPP dataset",ylim=c(0,1),xlim=c(0,1),lwd=2)
lines(probs,sapply(probs,function(x){sum(new[data3$SLALI==0]>x)/length(new[data3$SLALI==0])}),col="red",lwd=2)

new=predict(ODlogit1,newdata=data4,type="response")
old=predict(ODlogitNothing,newdata=data4,type="response")
probs=as.matrix(seq(0,1,0.01))
plot(probs,sapply(probs,function(x){sum(new[data4$SLALI==1]>x)/length(new[data4$SLALI==1])}),type="l",xlab="Calculated Risk",ylab="Sensitivity (Black) 1-Specificity (Red)",main="BIOPAR dataset",ylim=c(0,1),xlim=c(0,1),lwd=2)
lines(probs,sapply(probs,function(x){sum(new[data4$SLALI==0]>x)/length(new[data4$SLALI==0])}),col="red",lwd=2)

par(mai=c(0,0,0,0))
plot.new()
legend("center",legend="HMGB1",lty=1,lwd=c(2))
dev.off()







