library(aod)
library(ggplot2)
library(rms)
library(Hmisc)
library(PredictABEL)
library(hsdar)
library(QuantPsyc)
library(lsr)
perf = function(cut, mod, y)
{
  yhat = (mod$fit>cut)
  w = which(y==1)
  sensitivity = mean( yhat[w] == 1 ) 
  specificity = mean( yhat[-w] == 0 ) 
  c.rate = mean( y==yhat ) 
  d = cbind(sensitivity,specificity)-c(1,1)
  d = sqrt( d[1]^2 + d[2]^2 ) 
  out = t(as.matrix(c(sensitivity, specificity, c.rate,d)))
  colnames(out) = c("sensitivity", "specificity", "c.rate", "distance")
  return(out)
}

expit = function(x) exp(predict(x))/(1+exp(predict(x)))
#dataOD<-read.csv("/home/ben/Dropbox/CDSS/ParaOverdose.csv",sep=",",header=TRUE)
dataOD_Combo = read.csv("C://Users/Ben/Dropbox/CDSS/Dan/COMBO/COMBO.csv", header=TRUE, sep=",", na.strings=c("","NA"))
#dataOD_Combo = read.csv("/home/ben/Dropbox/CDSS/Dan/COMBO/COMBO.csv", header=TRUE, sep=",", na.strings=c("","NA"))

plot(dataOD_Combo$Paracetamol,dataOD_Combo$ALT)
points(dataOD_Combo$Paracetamol[dataOD_Combo$LateALI==0],dataOD_Combo$ALT[dataOD_Combo$LateALI==0],col="blue")
points(dataOD_Combo$Paracetamol[dataOD_Combo$LateALI==1],dataOD_Combo$ALT[dataOD_Combo$LateALI==1],col="red")

plot(dataOD$Paracetamol,dataOD$ALT)
points(dataOD$Paracetamol[dataOD$LateALI==0],dataOD$ALT[dataOD$LateALI==0],col="blue")
points(dataOD$Paracetamol[dataOD$LateALI==1],dataOD$ALT[dataOD$LateALI==1],col="red")

plot(dataOD_Combo$Paracetamol,dataOD_Combo$HMGB1)
points(dataOD_Combo$Paracetamol[dataOD_Combo$LateALI==0],dataOD_Combo$HMGB1[dataOD_Combo$LateALI==0],col="blue")
points(dataOD_Combo$Paracetamol[dataOD_Combo$LateALI==1],dataOD_Combo$HMGB1[dataOD_Combo$LateALI==1],col="red")

plot(dataOD$Paracetamol,dataOD$HMGB1)
points(dataOD$Paracetamol[dataOD$LateALI==0],dataOD$HMGB1[dataOD$LateALI==0],col="blue")
points(dataOD$Paracetamol[dataOD$LateALI==1],dataOD$HMGB1[dataOD$LateALI==1],col="red")

plot(dataOD_Combo$Paracetamol,dataOD_Combo$necrosis_K18)
points(dataOD_Combo$Paracetamol[dataOD_Combo$LateALI==0],dataOD_Combo$necrosis_K18[dataOD_Combo$LateALI==0],col="blue")
points(dataOD_Combo$Paracetamol[dataOD_Combo$LateALI==1],dataOD_Combo$necrosis_K18[dataOD_Combo$LateALI==1],col="red")

plot(dataOD$Paracetamol,dataOD$necrosis_K18)
points(dataOD$Paracetamol[dataOD$LateALI==0],dataOD$necrosis_K18[dataOD$LateALI==0],col="blue")
points(dataOD$Paracetamol[dataOD$LateALI==1],dataOD$necrosis_K18[dataOD$LateALI==1],col="red")

plot(dataOD_Combo$Paracetamol,dataOD_Combo$apoptosis_K18)
points(dataOD_Combo$Paracetamol[dataOD_Combo$LateALI==0],dataOD_Combo$apoptosis_K18[dataOD_Combo$LateALI==0],col="blue")
points(dataOD_Combo$Paracetamol[dataOD_Combo$LateALI==1],dataOD_Combo$apoptosis_K18[dataOD_Combo$LateALI==1],col="red")

plot(dataOD$Paracetamol,dataOD$apoptosis_K18)
points(dataOD$Paracetamol[dataOD$LateALI==0],dataOD$apoptosis_K18[dataOD$LateALI==0],col="blue")
points(dataOD$Paracetamol[dataOD$LateALI==1],dataOD$apoptosis_K18[dataOD$LateALI==1],col="red")

plot(dataOD_Combo$Paracetamol,dataOD_Combo$miR122_Let_7d)
points(dataOD_Combo$Paracetamol[dataOD_Combo$LateALI==0],dataOD_Combo$miR122_Let_7d[dataOD_Combo$LateALI==0],col="blue")
points(dataOD_Combo$Paracetamol[dataOD_Combo$LateALI==1],dataOD_Combo$miR122_Let_7d[dataOD_Combo$LateALI==1],col="red")

plot(dataOD$miR122_Let_7d)
points(dataOD$Paracetamol[dataOD$LateALI==0],dataOD$miR122_Let_7d[dataOD$LateALI==0],col="blue")
points(dataOD$Paracetamol[dataOD$LateALI==1],dataOD$miR122_Let_7d[dataOD$LateALI==1],col="red")

dataOD<-dataOD_Combo[substr(dataOD_Combo$ID,1,1)=="M"& !is.na(dataOD_Combo$LateALI),]
dataOD$SexN<-ifelse(dataOD$Sex=="F",1,0)
dataOD$StaggeredN<-ifelse(dataOD$Staggered=="Y",1,0)
colSums(is.na(dataOD))

ODlogitAge <- glm(LateALI ~ scale(Age,center=TRUE), data = dataOD, family = "binomial")
ODlogitGender <- glm(LateALI ~ scale(Sex,center=TRUE), data = dataOD, family = "binomial")
ODlogitS <- glm(LateALI ~ scale(StaggeredN,center=TRUE), data = dataOD, family = "binomial")
ODlogitConc <- glm(LateALI ~ scale(Paracetamol,center=TRUE), data = dataOD, family = "binomial")

ODlogitAI <- glm(LateALI ~ scale(AmountIng,center=TRUE), data = dataOD, family = "binomial")
ODlogitPOD <- glm(LateALI ~ scale(ParaOD,center=TRUE), data = dataOD, family = "binomial")
ODlogitS <- glm(LateALI ~ scale(StaggeredN,center=TRUE), data = dataOD, family = "binomial")
ODlogitConc <- glm(LateALI ~ scale(Paracetamol,center=TRUE), data = dataOD, family = "binomial")
ODlogitCr <- glm(LateALI ~ scale(Creatinine,center=TRUE), data = dataOD, family = "binomial")
ODlogitBil <- glm(LateALI ~ scale(Bilirubin,center=TRUE), data = dataOD, family = "binomial")
ODlogitALT <- glm(LateALI ~ scale(ALT,center=TRUE), data = dataOD, family = "binomial")
ODlogitAlk <- glm(LateALI ~ scale(Alk.Phos,center=TRUE), data = dataOD, family = "binomial")
ODlogitINR <- glm(LateALI ~ scale(INR,center=TRUE), data = dataOD, family = "binomial")
ODlogitHMG <- glm(LateALI ~ scale(HMGB1,center=TRUE), data = dataOD, family = "binomial")
ODlogitNec <- glm(LateALI ~ scale(necrosis_K18,center=TRUE), data = dataOD, family = "binomial")
ODlogitApop <- glm(LateALI ~ scale(apoptosis_K18,center=TRUE), data = dataOD, family = "binomial")
ODlogitmiR <- glm(LateALI ~ scale(miR122_Let_7d,center=TRUE), data = dataOD, family = "binomial")
ODlogitGLDH <- glm(LateALI ~ scale(GLDH,center=TRUE), data = dataOD, family = "binomial")

summary(ODlogitAI)#(LateALI ~  AmountIng, data = dataOD, family = "binomial")
#summary(ODlogitPOD)#(LateALI ~  ParaOD, data = dataOD, family = "binomial")
#summary(ODlogitS)#(LateALI ~  Staggered, data = dataOD, family = "binomial")
summary(ODlogitConc)#(LateALI ~  Paracetamol, data = dataOD, family = "binomial")
#summary(ODlogitCr)#(LateALI ~  Creatinine, data = dataOD, family = "binomial")
summary(ODlogitBil)#(LateALI ~  Bilirubin, data = dataOD, family = "binomial")
summary(ODlogitALT)#(LateALI ~  ALT, data = dataOD, family = "binomial")
summary(ODlogitAlk)#(LateALI ~  Alk.Phos, data = dataOD, family = "binomial")
summary(ODlogitINR)#(LateALI ~  INR, data = dataOD, family = "binomial")
summary(ODlogitHMG)#(LateALI ~  HMGB1, data = dataOD, family = "binomial")
summary(ODlogitNec)#(LateALI ~  necrosis_K18, data = dataOD, family = "binomial")
summary(ODlogitApop)#(LateALI ~  apoptosis_K18, data = dataOD, family = "binomial")
summary(ODlogitmiR)#(LateALI ~  miR122_Let_7d, data = dataOD, family = "binomial")
summary(ODlogitGLDH)#(LateALI ~  GLDH, data = dataOD, family = "binomial")

#Make dataset for stepwise selection with no "NA"
data1<-na.omit(cbind(dataOD$LateALI,dataOD$Paracetamol,dataOD$Bilirubin,dataOD$ALT,dataOD$Alk.Phos,dataOD$INR,dataOD$HMGB1,dataOD$necrosis_K18,dataOD$apoptosis_K18,dataOD$miR122_Let_7d,dataOD$GLDH))
data1=as.data.frame(data1)
colnames(data1)<-c("LateALI","Paracetamol","Bilirubin","ALT","Alk.Phos","INR","HMGB1","necrosis_K18","apoptosis_K18","miR122_Let_7d","GLDH")

ODlogit1 <- glm(data1$LateALI ~ scale(data1$Paracetamol,center=TRUE) + scale(data1$Bilirubin,center=TRUE) + scale(data1$ALT,center=TRUE) + scale(data1$Alk.Phos,center=TRUE) + scale(data1$INR,center=TRUE) + scale(data1$HMGB1,center=TRUE) + scale(data1$necrosis_K18,center=TRUE) + scale(data1$apoptosis_K18,center=TRUE)+ scale(data1$miR122_Let_7d,center=TRUE) + scale(data1$GLDH,center=TRUE), family = "binomial")
summary(ODlogit1)
ODlogitNothing <- glm(data1$LateALI ~  scale(data1$ALT,center=TRUE), data = data1, family = "binomial")
summary(ODlogitNothing)
both = step(ODlogit1,direction="both",scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)))
back = step(ODlogit1,direction="backward",scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)))
forward = step(ODlogitNothing,scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)),direction="forward")
formula(both)
formula(back)
formula(forward)

#Make different dataset with less "NA" omits - remove AmountIng
data2<-na.omit(cbind(dataOD$LateALI,dataOD$Paracetamol,dataOD$Bilirubin,dataOD$ALT,dataOD$Alk.Phos,dataOD$INR,dataOD$HMGB1,dataOD$necrosis_K18,dataOD$apoptosis_K18,dataOD$miR122_Let_7d,dataOD$GLDH))
data2=as.data.frame(data2)
colnames(data2)<-c("LateALI","Paracetamol","Bilirubin","ALT","Alk.Phos","INR","HMGB1","necrosis_K18","apoptosis_K18","miR122_Let_7d","GLDH")

ODlogit1 <- glm(data2$LateALI ~ scale(data2$Paracetamol,center=TRUE) + scale(data2$Bilirubin,center=TRUE) + scale(data2$ALT,center=TRUE) + scale(data2$Alk.Phos,center=TRUE) + scale(data2$INR,center=TRUE) + scale(data2$HMGB1,center=TRUE) + scale(data2$necrosis_K18,center=TRUE) + scale(data2$apoptosis_K18,center=TRUE) + scale(data2$miR122_Let_7d,center=TRUE) + scale(data2$GLDH,center=TRUE), family = "binomial")
summary(ODlogit1)
ODlogitNothing <- glm(data2$LateALI ~  scale(data2$ALT,center=TRUE), data = data2, family = "binomial")
summary(ODlogitNothing)
both = step(ODlogit1,direction="both",scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)))
back = step(ODlogit1,direction="backward",scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)))
forward = step(ODlogitNothing,scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)),direction="forward")
formula(both)
formula(back)
formula(forward)

#Make dataset of just significant multiple model data 
data3<-na.omit(cbind(dataOD$LateALI,dataOD$ALT,dataOD$HMGB1,dataOD$necrosis_K18,dataOD$apoptosis_K18,dataOD$miR122_Let_7d,dataOD$GLDH))
data3=as.data.frame(data3)
colnames(data3)<-c("LateALI","ALT","HMGB1","necrosis_K18","apoptosis_K18","miR122_Let_7d","GLDH")

ODlogit1 <- glm(data3$LateALI ~ scale(data3$ALT,center=TRUE) + scale(data3$HMGB1,center=TRUE) + scale(data3$necrosis_K18,center=TRUE) + scale(data3$apoptosis_K18,center=TRUE) + scale(data3$miR122_Let_7d,center=TRUE) + scale(data3$GLDH,center=TRUE), family = "binomial")
summary(ODlogit1)
ODlogitNothing <- glm(data3$LateALI ~  scale(data3$ALT,center=TRUE), data = data2, family = "binomial")
summary(ODlogitNothing)
both = step(ODlogit1,direction="both",scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)))
back = step(ODlogit1,direction="backward",scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)))
forward = step(ODlogitNothing,scope=list(lower=formula(ODlogitNothing),upper=formula(ODlogit1)),direction="forward")
formula(both)
formula(back)
formula(forward)

#For word table

ODlogit1 <- glm(data3$LateALI ~ scale(data3$ALT,center=TRUE,scale=TRUE) + scale(data3$HMGB1,center=TRUE,scale=TRUE) + scale(data3$necrosis_K18,center=TRUE,scale=TRUE) + scale(data3$apoptosis_K18,center=TRUE,scale=TRUE) + scale(data3$miR122_Let_7d,center=TRUE,scale=TRUE), family = "binomial")
summary(ODlogit1)

new=expit(ODlogit1)
new[is.na(new)]<-1
old=expit(ODlogitNothing)

comp=improveProb(old[!is.nan(new)],new[!is.nan(new)],ODlogit1$y[!is.nan(new)])
comp
plot(new[ODlogit1$y==0],old[ODlogit1$y==0],xlab="Augmented Probability for Late ALI",ylab="ALT-based Probability for Late ALI",main ="MAPP Patients\n who did not experience Late ALI",xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)
plot(new[ODlogit1$y==1],old[ODlogit1$y==1],xlab="Augmented Probability for Late ALI",ylab="ALT-based Probability for Late ALI",main="MAPP Patients\n who did experience Late ALI",xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)

ODlogit2 <- glm(data3$LateALI ~ scale(data3$ALT,center=TRUE,scale=TRUE) + scale(data3$miR122_Let_7d,center=TRUE,scale=TRUE), family = "binomial")
summary(ODlogit2)

new=expit(ODlogit2)
new[is.na(new)]<-1
old=expit(ODlogitNothing)

comp=improveProb(old[!is.nan(new)],new[!is.nan(new)],ODlogit2$y[!is.nan(new)])
comp

ODlogit3 <- glm(data3$LateALI ~ scale(data3$ALT,center=TRUE,scale=TRUE) + scale(data3$necrosis_K18,center=TRUE,scale=TRUE) + scale(data3$miR122_Let_7d,center=TRUE,scale=TRUE), family = "binomial")
summary(ODlogit3)

new=expit(ODlogit3)
new[is.na(new)]<-1
old=expit(ODlogitNothing)

comp=improveProb(old[!is.nan(new)],new[!is.nan(new)],ODlogit3$y[!is.nan(new)])
comp

ODlogit4 <- glm(data3$LateALI ~ scale(data3$ALT,center=TRUE,scale=TRUE) + scale(data3$necrosis_K18,center=TRUE,scale=TRUE) +scale(data3$apoptosis_K18,center=TRUE,scale=TRUE) + scale(data3$miR122_Let_7d,center=TRUE,scale=TRUE), family = "binomial")
summary(ODlogit4)

new=expit(ODlogit4)
new[is.na(new)]<-1
old=expit(ODlogitNothing)

comp=improveProb(old[!is.nan(new)],new[!is.nan(new)],ODlogit4$y[!is.nan(new)])
comp



s = seq(.01,.99,length=1000)
OUT = matrix(0,1000,4)
for(i in 1:1000) OUT[i,]=perf(s[i],ODlogit1,ODlogit1$y)
plot(s,OUT[,1],xlab="Cutoff",ylab="Value",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),type="l",lwd=2,axes=FALSE,col=2)
axis(1,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
axis(2,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
lines(s,OUT[,2],col="darkgreen",lwd=2)
lines(s,OUT[,3],col=4,lwd=2)
lines(s,OUT[,4],col="darkred",lwd=2)
box()
legend(.1,.6,col=c(2,"darkgreen",4,"darkred"),lwd=c(2,2,2,2),c("Sensitivity","Specificity","Classification Rate","Distance"))

s = seq(.01,.99,length=1000)
OUTold = matrix(0,1000,4)
for(i in 1:1000) OUTold[i,]=perf(s[i],ODlogitNothing,ODlogitNothing$y)
plot(s,OUTold[,1],xlab="Cutoff",ylab="Value",cex.lab=1.5,cex.axis=1.5,ylim=c(0,1),type="l",lwd=2,axes=FALSE,col=2)
axis(1,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
axis(2,seq(0,1,length=5),seq(0,1,length=5),cex.lab=1.5)
lines(s,OUTold[,2],col="darkgreen",lwd=2)
lines(s,OUTold[,3],col=4,lwd=2)
lines(s,OUTold[,4],col="darkred",lwd=2)
box()
legend(.1,.6,col=c(2,"darkgreen",4,"darkred"),lwd=c(2,2,2,2),c("Sensitivity","Specificity","Classification Rate","Distance"))


probs=as.matrix(seq(0.01,0.99,0.01))
new[is.na(new)]<-1
probs_dist=sapply(probs,function(x){sqrt(((sum(new[data3$LateALI==1]>x)/length(new[data3$LateALI==1])-1)^2)+((sum(new[data3$LateALI==0]<x)/length(new[data3$LateALI==0])-1)^2))})
probs[probs_dist==min(sapply(probs,function(x){sqrt(((sum(new[data3$LateALI==1]>x)/length(new[data3$LateALI==1])-1)^2)+((sum(new[data3$LateALI==0]<x)/length(new[data3$LateALI==0])-1)^2))}))]
sum(new[data3$LateALI==1]>0.44)
sum(new[data3$LateALI==0]<0.44)


new=expit(ODlogit1)
new[is.na(new)]<-1
probs_dist=sapply(probs,function(x){sqrt(((sum(new[data3$LateALI==1]>x)/length(new[data3$LateALI==1])-1)^2)+((sum(new[data3$LateALI==0]<x)/length(new[data3$LateALI==0])-1)^2))})
cut=median(probs[probs_dist==min(sapply(probs,function(x){sqrt(((sum(new[data3$LateALI==1]>x)/length(new[data3$LateALI==1])-1)^2)+((sum(new[data3$LateALI==0]<x)/length(new[data3$LateALI==0])-1)^2))}))])
sum(new[data3$LateALI==1]>cut)
sum(new[data3$LateALI==0]<cut)

probs_dist=sapply(probs,function(x){sqrt(((sum(old[data3$LateALI==1]>x)/length(old[data3$LateALI==1])-1)^2)+((sum(old[data3$LateALI==0]<x)/length(old[data3$LateALI==0])-1)^2))})
probs[probs_dist==min(sapply(probs,function(x){sqrt(((sum(old[data3$LateALI==1]>x)/length(old[data3$LateALI==1])-1)^2)+((sum(old[data3$LateALI==0]<x)/length(old[data3$LateALI==0])-1)^2))}))]
sum(old[data3$LateALI==1]>0.05)
sum(old[data3$LateALI==0]<0.05)

png("C:\\Users/Ben/Dropbox/CDSS/Dan/PlotforJo.png",width = 900,height=900,res=120)
plot(probs,sapply(probs,function(x){sum(new[data3$LateALI==1]>x)/length(new[data3$LateALI==1])}),type="l",xlab="Calculated Risk",ylab="Sensitivity (Black) 1-Specificity (Red)",main="MAPP dataset",ylim=c(0,1),xlim=c(0,1),lwd=2)
lines(probs,sapply(probs,function(x){sum(new[data3$LateALI==0]>x)/length(new[data3$LateALI==0])}),col="red",lwd=2)
lines(probs,sapply(probs,function(x){sum(old[data3$LateALI==1]>x)/length(old[data3$LateALI==1])}),type="l",lty=2)
lines(probs,sapply(probs,function(x){sum(old[data3$LateALI==0]>x)/length(old[data3$LateALI==0])}),col="red",lty=2)

ODlogitmiR <- glm(data3$LateALI ~ scale(data3$miR122_Let_7d,center=TRUE,scale=TRUE), family = "binomial")
summary(ODlogitmiR)
new=expit(ODlogitmiR)
new[is.na(new)]<-1

lines(probs,sapply(probs,function(x){sum(new[data3$LateALI==1]>x)/length(new[data3$LateALI==1])}),type="l",lty=4,lwd=1.75,xlab="Calculated Risk",ylab="Sensitivity (Black) 1-Specificity (Red)",main="MAPP dataset",ylim=c(0,1),xlim=c(0,1))
lines(probs,sapply(probs,function(x){sum(new[data3$LateALI==0]>x)/length(new[data3$LateALI==0])}),col="red",lty=4,lwd=1.75)
legend(.35,.6,legend=c("Full","miR-122 let 7d","ALT"),lty=c(1,4,2),lwd =c(2,1.75,1))
dev.off()

new=expit(ODlogit1)
newmiR=expit(ODlogitmiR)
old=expit(ODlogitNothing)


############BIOPARR dataset

dataOD_B<-dataOD_Combo[substr(dataOD_Combo$ID,1,1)=="B"& !is.na(dataOD_Combo$LateALI),]
colSums(is.na(dataOD_B))
data4<-na.omit(cbind(dataOD_B$LateALI,scale(dataOD_B$ALT,scale=TRUE,center = TRUE),scale(dataOD_B$HMGB1,scale=TRUE,center = TRUE),scale(dataOD_B$necrosis_K18,scale=TRUE,center = TRUE),scale(dataOD_B$apoptosis_K18,scale=TRUE,center = TRUE),scale(dataOD_B$miR122_Let_7d,scale=TRUE,center = TRUE)))
data4=as.data.frame(data4)
colnames(data4)<-c("LateALI","ALT","HMGB1","necrosis_K18","apoptosis_K18","miR122_Let_7d")

coef=ODlogit1$coefficients
pred1=coef[1]+coef[2]*data4$ALT+coef[3]*data4$HMGB1+coef[4]*data4$necrosis_K18+coef[5]*data4$apoptosis_K18+coef[6]*data4$miR122_Let_7d
valid1=(exp(pred1)/(1+exp(pred1)))
valid1[is.na(valid1)]<-1
coef2=ODlogitNothing$coefficients
pred2=coef2[1]+coef2[2]*data4$ALT
valid2=(exp(pred2)/(1+exp(pred2)))
valid_comp=improveProb(valid2,valid1,data4$LateALI)
valid_comp
new=valid1
probs_dist=sapply(probs,function(x){sqrt(((sum(new[data4$LateALI==1]>x)/length(new[data4$LateALI==1])-1)^2)+((sum(new[data4$LateALI==0]<x)/length(new[data4$LateALI==0])-1)^2))})
probs[probs_dist==min(sapply(probs,function(x){sqrt(((sum(new[data3$LateALI==1]>x)/length(new[data3$LateALI==1])-1)^2)+((sum(new[data3$LateALI==0]<x)/length(new[data3$LateALI==0])-1)^2))}))]
sum(new[data3$LateALI==1]>0.44)
sum(new[data3$LateALI==0]<0.44)

coef=ODlogit2$coefficients
pred1=coef[1]+coef[2]*data4$ALT+coef[3]*data4$miR122_Let_7d
valid1=(exp(pred1)/(1+exp(pred1)))
valid1[is.na(valid1)]<-1
coef2=ODlogitNothing$coefficients
pred2=coef2[1]+coef2[2]*data4$ALT
valid2=(exp(pred2)/(1+exp(pred2)))
valid_comp=improveProb(valid2,valid1,data4$LateALI)
valid_comp

coef=ODlogit3$coefficients
pred1=coef[1]+coef[2]*data4$ALT+coef[3]*data4$necrosis_K18+coef[4]*data4$miR122_Let_7d
valid1=(exp(pred1)/(1+exp(pred1)))
valid1[is.na(valid1)]<-1
coef2=ODlogitNothing$coefficients
pred2=coef2[1]+coef2[2]*data4$ALT
valid2=(exp(pred2)/(1+exp(pred2)))
valid_comp=improveProb(valid2,valid1,data4$LateALI)
valid_comp

coef=ODlogit4$coefficients
pred1=coef[1]+coef[2]*data4$ALT+coef[3]*data4$necrosis_K18+coef[4]*data4$apoptosis_K18+coef[5]*data4$miR122_Let_7d
valid1=(exp(pred1)/(1+exp(pred1)))
valid1[is.na(valid1)]<-1
coef2=ODlogitNothing$coefficients
pred2=coef2[1]+coef2[2]*data4$ALT
valid2=(exp(pred2)/(1+exp(pred2)))
valid_comp=improveProb(valid2,valid1,data4$LateALI)
valid_comp

plot(valid1[data4$LateALI==0],valid2[data4$LateALI==0],xlab="Augmented Probability for Late ALI",ylab="ALT-based Probability for Late ALI",main ="BIOPAR Patients\n who did not experience Late ALI",xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)
plot(valid1[data4$LateALI==1],valid2[data4$LateALI==1],xlab="Augmented Probability for Late ALI",ylab="ALT-based Probability for Late ALI",main="BIOPAR Patients\n who did experience Late ALI",xlim=c(0,1),ylim=c(0,1))
abline(a=0,b=1)

probs=as.matrix(seq(0,1,0.01))
plot(probs,sapply(probs,function(x){sum(valid1[data4$LateALI==1]>x)/length(valid1[data4$LateALI==1])}),type="l",xlab="Calculated Risk",ylab="Sensitivity (Black) 1-Specificity (Red)",main="BIOPAR dataset")
lines(probs,sapply(probs,function(x){sum(valid1[data4$LateALI==0]>x)/length(valid1[data4$LateALI==0])}),col="red")
lines(probs,sapply(probs,function(x){sum(valid2[data4$LateALI==1]>x)/length(valid2[data4$LateALI==1])}),type="l",lty=2)
lines(probs,sapply(probs,function(x){sum(valid2[data4$LateALI==0]>x)/length(valid2[data4$LateALI==0])}),col="red",lty=2)


 
# # specify column number of outcome variable
# cOutcome <- 1
# # specify column numbers of non-genetic predictors
# cNonGenPred <- c(2:7)
# # specify column numbers of non-genetic predictors that are categorical
# cNonGenPredCat <- c(0)
# # specify column numbers of genetic predictors
# cGenPred <- c(0)
# # specify column numbers of genetic predictors that are categorical
# cGenPredCat <- c(0)
# # fit logistic regression model
# riskmodel <- fitLogRegModel(data=data4, cOutcome=cOutcome,cNonGenPreds=cNonGenPred, cNonGenPredsCat=cNonGenPredCat,cGenPreds=cGenPred, cGenPredsCat=cGenPredCat)
# # show summary details for the fitted risk model
# summary(riskmodel)

