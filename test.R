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
model<-SetYinitMode(model,mode = "Intercept_0",c(1,2))
model<-SetYinitMode(model,mode = "CircadianFit",c(0,48))
model@initpos
model@initspeed

model<-ReplicateDrivingForce(model,c(0.1,24.0),10)
plot(model@SWdist$Time,model@SWdist$NREM,type="l")

model<-AddSinF(model,FixPer = 24)
model@SinForce
model@PerSin

model<-SetFittingValue(model,value = "RSS")

object<-model

model<-PenalizeUnstableFit(model,value = T,PredictedValueInterval = c(0,24), StabilityDayCheck = 10)

model

ControlParams(model,params=c(Wake=1,AmpSin=1,PhiSin=1)) # miss core param on purpose
ControlParams(model,params=c(intercept=1,omega=1,loggamma=1,Wake=1,AmpSin=1,PhiSin=1))

microbenchmark(GetFit(model,params),times=1000)