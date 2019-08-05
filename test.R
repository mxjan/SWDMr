library(SWDMr)

# Data
# ****** Prepare Sleep-Wake data ****** 
files <- list.files(path="../data/Sleep_States/",pattern = paste("^","BL6",sep=""),full.names = T)
SWdf<-Read_SW(files,concattimesec = 360) # 300 = 5min, 180 = 3 min, 360 = 6min
SWdf<-SWdf_AddLD(SWdf) # Add Light and Dark
SWdf<-SWdf_DayMerging(SWdf,Daysformat=list(c(1,2),c(1,2),3,4,c(1,2)),concattimesec=360)
SWdf<-SWdf_AddSD(SWdf,c(48,54))

# ***** Prepare Gene expression data ******
# Gene expression
load("../data/CHor2018_TimeCourse_NormalizedData_Genes.RData")
# Meta data
load("../data/CHor2018_TimeCourse_MetaData.RData")
# Remove data sampled 7days after SD
rna_expr<-rna_expr[,rna_meta$time<100]
rna_meta<-rna_meta[rna_meta$time<100,]
# Row = samples, column = gene
rna_expr<-t(rna_expr)
#rna_expr<-2^rna_expr
rna_expr<-as.data.frame(rna_expr)
# Add time factor
rna_expr$Time<-rna_meta$time+24

swdmr <- SWDMr(SWdist=SWdf, Gexp=rna_expr)
model<-initDDHOmodel(swdmr,VarExp = "Arntl")
summary(swdmr)
summary(model)
model
FixNaturalFrequency(model,2*pi/24)
model@omega
#model<-FixNaturalFrequency(model,2*pi/24)
#model<-FixDamping(model,log(0.015))
model@omega
model<-AddForce(model,"Wake",1)
model<-AddForce(model,"Wake")
model<-AddEffect(model,"SD",1)
model<-SetYinitMode(model,mode = "Free")
model<-SetYinitMode(model,mode = "Fixed",c(1,2))
model<-SetYinitMode(model,mode = "Intercept_0")
model<-SetYinitMode(model,mode = "CircadianFit",c(0,48))
model@initpos
model@initspeed

model<-ReplicateDrivingForce(model,c(0.1,24.0),10)
plot(model@SWdist$Time,model@SWdist$NREM,type="l")

model<-AddSinF(model,FixPer = 24)
model@SinForce
model@PerSin

model<-SetFittingValue(model,value = "RSS")

model<-PenalizeUnstableFit(model,value = T,PredictedValueInterval = c(0,24), StabilityDayCheck = 10)

model


params=c(intercept=1,omega=1,loggamma=1,Wake=1,AmpSin=1,PhiSin=1)

library(microbenchmark)
microbenchmark(SWDMrFit(model,params),times=1000)

out<-SWDMrFit(model,params)
plot(out$time,out$y1,type="l")

model<-SetYinitMode(model,mode = "Intercept_0")
out<-SWDMrFit(model,params)
plot(out$time,out$y1,type="l")

SWDMrStats(model,out,detailed = T)

stat<-SWDMrStats(model,out,FittingValue = "RSS")
stat$val
RSSpen<-AddUnstableFitPenalization(object = model,fitted = out,FittingValue = "RSS",val = stat$val, var = stat$var)
RSSpen

stat<-SWDMrStats(model,out,FittingValue = "NLL")
stat$val
NLLpen<-AddUnstableFitPenalization(object = model,fitted = out,FittingValue = "NLL",val = stat$val, var = stat$var)
NLLpen


library(optimx)
objfun<-SWDMrGetEvalFun(model)
params["intercept"]<-5.5
params["omega"]<-2*pi/23.75
params["loggamma"]<-log(0.01)
params["Wake"]<-0
params["AmpSin"]<-0
optimxres<-optimx(params,objfun,method="nlminb")
out<-SWDMrFit(model,params = optimxres[1,])
plot(out$time,out$y1,type="l")


################################################################################
Gene<-"Arntl"
MeanPerTime<-aggregate(rna_expr[rna_meta$SD_NSD == "NSD",Gene],list(rna_meta$time[rna_meta$SD_NSD == "NSD"]),mean)
MeanGeneExprInBaseline<-(max(MeanPerTime$x)+min(MeanPerTime$x))/2

############# Build model ###############
model<-initDDHOmodel(swdmr,VarExp = Gene)

# Mean expression in baseline
MeanPerTime<-aggregate(rna_expr[rna_meta$SD_NSD == "NSD",Gene],list(rna_meta$time[rna_meta$SD_NSD == "NSD"]),mean)
MeanGeneExprInBaseline<-(max(MeanPerTime$x)+min(MeanPerTime$x))/2
# Fix the intercepts
model<-FixIntercept(model,MeanGeneExprInBaseline)
# Add sleep-wake force
model<-AddForce(model,"Wake")
model<-AddForce(model,"Sleep")
# Start is set at intercept with speed of 0
model<-SetYinitMode(model,mode = "Intercept_0",values = c(0,48))
# We replicate baseline for 20 day
model<-ReplicateDrivingForce(model,c(0.1,24.0),500)
# A sin-wave force is applied with a period of 24h
model<-AddSinF(model,FixPer = 24)
# Compute the fit using RSS
model<-SetFittingValue(model,value = "RSS")
# Penalize the fitting for unstable value for 10 replicated days
model<-PenalizeUnstableFit(model,value = T,PredictedValueInterval = c(0,48), StabilityDayCheck = 400)

# Get objective function
objfun<-SWDMrGetEvalFun(model)
# Limits of the model
params<-c(Wake=0,Sleep=0,loggamma=log(1e-1),omega=2*pi/23.8,AmpSin=1,PhiSin=5.2)
#params<-c(Wake=0.056834519219219,Sleep=0.00814929150840472,loggamma=-5.22,omega=0.296,AmpSin=0.00204,PhiSin=2.48179822476184)
lower<-c(Wake=-Inf,Sleep=-Inf,loggamma=-Inf,omega=2*pi/27,AmpSin=0,PhiSin=0)
upper<-c(Wake=Inf,Sleep=Inf,loggamma=Inf,omega=2*pi/20,AmpSin=Inf,PhiSin=2*pi)
# Fit
optimxres<-optimx(params,objfun,method=c("nlminb"),lower = lower,upper=upper)
optimxres

out<-SWDMrFit(model,params = optimxres[1,])
#out<-SWDMrFit(model,params = params)
plot(out$time,out$y1,type="l",ylim=c(min(c(out$y1,model@Gexp[,Gene])) , max(c(out$y1,model@Gexp[,Gene]))))
plot(out$time,out$y1,type="l",xlim=c(0,120),ylim=c(min(c(out$y1,model@Gexp[,Gene])) , max(c(out$y1,model@Gexp[,Gene]))))
points(model@Gexp$Time,model@Gexp[,Gene])
SWDMrStats(model,out,detailed = T)


SWDMr:::StandardFittingPlot(model,optimxres[1,])

par(mfrow=c(2,1))
Fa<-ForceApplied(model,optimxres[1,])
mF<-PctAbsForceApplied(model,optimxres[1,])
max<-apply(mF[,c(2,3,4)],1,function(x) as.character(which(x == max(x))))
max[max == "1"]<-"blue"
max[max == "2"]<-"red"
max[max == "3"]<-"green"

plot(Fa$time,Fa$force,xlim=c(0,120),xaxt="n",col=max,pch=19,cex=1,xlab="Time",ylab="Force [expr / h / h]",main="Total Force applied to the model")
abline(h=0)
lines(Fa$time,Fa$force)
axis(1,seq(0,120,by=12),seq(0,120,by=12))

plot(-1,-1,xlim=c(0,120),ylim=c(0,1),xaxt="n",xlab="Time",ylab="Percentage of absolute force",main="Contributors to force applied")
axis(1,seq(0,120,by=12),seq(0,120,by=12))
polygon(c(seq(min(mF$Time),120,by=.1),seq(120,min(mF$Time),by=-.1)), c(rep(0,length(mF[,"Wake"])),rev(mF[,"Wake"])),col = "blue")
polygon(c(seq(min(mF$Time),120,by=.1),seq(120,min(mF$Time),by=-.1)), c(mF[,"Wake"],rev(mF[,"Wake"]+mF[,"Sleep"])),col = "red")
polygon(c(seq(min(mF$Time),120,by=.1),seq(120,min(mF$Time),by=-.1)), c(mF[,"Wake"]+mF[,"Sleep"],rev(mF[,"Wake"]+mF[,"Sleep"]+mF[,"SinF"])),col = "green")


