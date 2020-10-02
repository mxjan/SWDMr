############# SET CLASS ################

setClass(
  # Set the name for the class
  Class = "SWDMr_ProcS",
  
  contains="SWDMr",
  
  # Define the slots
  representation = representation(
    
    # Name of variable explained
    VarExp = "character",
    
    # Core parameters of process S
    AsympWake = "numeric", # Asymptote for wake
    AsympSleep = "numeric", # Asymptote for sleep
    TauWake = "numeric", # time constant wake
    TauSleep = "numeric", # time constant sleep
    
    # How the initial value of process S is calculated (fixed or free)
    initmod = "character",
    initpos = "numeric",

    # Objective function is optimized for FitttingValue, can be RSS or Negative Log Likelihood (NLL)
    FittingValue = "character",
    
    verbose = "numeric"
  ),
  
  prototype=list(
    initmod = "Free",
    initpos = 0,

    FittingValue = "RSS",
    
    verbose = 1
  ),
  
  validity=SWDMr_ProcScheck
  
)



############# Summaries ################

setMethod("GetFreeFixedParams","SWDMr_ProcS",function(object,...) {
  
  freeparams<-matrix(ncol=3,nrow=0)
  fixedparams<-matrix(ncol=3,nrow=0)

  # Core parameters
  CoreParams<-c("AsympWake","AsympSleep","TauWake","TauSleep")
  
  for (param in CoreParams){
    if ( length(slot(object,param)) == 0){
      freeparams<-rbind(freeparams,c("Core parameter",param,NA))
    }else{
      fixedparams<-rbind(fixedparams,c("Core parameter",param,slot(object,param)))
    }
  }
  
  # Yinit mode
  if (object@initmod == "Free"){
    freeparams<-rbind(freeparams,c("Yinit","start_pos",NA))
  }
  if (object@initmod == "Fixed"){
    fixedparams<-rbind(fixedparams,c("Yinit","start_pos",object@initpos))
  }
  
  freeparams<-data.frame(freeparams,stringsAsFactors = F)
  colnames(freeparams)<-c("parameter","subparameter","value")
  fixedparams<-data.frame(fixedparams,stringsAsFactors = F)
  colnames(fixedparams)<-c("parameter","subparameter","value")
  
  return(list(FreeParams = freeparams,FixedParams = fixedparams))
  
})


setMethod("ShowFreeParams","SWDMr_ProcS",function(object) {
  
  cat("~~~~~~~~~~~ Current parameter setting ~~~~~~~~~~ \n\n")
  
  parsedparams<-GetFreeFixedParams(object)
  
  freemsg<-"* [free parameters] "
  
  if (nrow(parsedparams$FreeParams) > 0){
    for (i in 1:nrow(parsedparams$FreeParams)){
      cat(freemsg,parsedparams$FreeParams[i,2]," (",parsedparams$FreeParams[i,1],") ","\n",sep = "") 
    }
  }
  
  cat("\n\n")
  
  freemsg<-"* [fixed parameters] "
  
  if (nrow(parsedparams$FixedParams) > 0){
    for (i in 1:nrow(parsedparams$FixedParams)){
      cat(freemsg,parsedparams$FixedParams[i,2]," (",parsedparams$FixedParams[i,1],") : ",parsedparams$FixedParams[i,3],"\n",sep = "") 
    }
  }
  
  cat("\n")
  
})

setMethod("summary", "SWDMr_ProcS", function(object) {
  cat("~~~~~~~~ This is a S4 SWDMr_ProcS object ~~~~~~~~ \n\n")
  cat("Display the current setting for your fitting\n\n")
  ShowFreeParams(object)
})



############# Set model free/fixed parameters ################
setMethod("FixAsymptWake",signature="SWDMr_ProcS", function(object,value){
  object@AsympWake <- value 
  return(object)
})

setMethod("FixAsymptSleep",signature="SWDMr_ProcS", function(object,value){
  object@AsympSleep <- value 
  return(object)
})
setMethod("FixTimeConstWake",signature="SWDMr_ProcS", function(object,value){
  object@TauWake <- value 
  return(object)
})
setMethod("FixTimeConstSleep",signature="SWDMr_ProcS", function(object,value){
  object@TauSleep <- value 
  return(object)
})


setMethod("SetYinitMode",signature="SWDMr_ProcS", function(object,mode="Free",values=numeric(0)){
  
  if (mode == "Free"){
    object@initmod <- "Free"
    object@initpos <- numeric(0)
  } else if (mode == "Fixed"){
    object@initmod <- "Fixed"
    object@initpos <- values[1]
  } else {
    stop("Mode not found")
  }
  return(object)
})

setMethod("SetFittingValue",signature="SWDMr_ProcS", function(object,value = "RSS"){
  
  if (! value %in% c("RSS","NLL","LL")){
    stop("The returned value should be RSS, NLL or LL")
  }
  object@FittingValue <- value
  return(object)
})


setMethod("SWDMrFit",signature="SWDMr_ProcS", function(object,params){
  
  # Control parameters given
  allparams<-GetAllParams(object,params)
  
  # Get Core params
  AsympWake_idx<-which(allparams$subparameter == "AsympWake")
  AsympWake<-allparams$value[AsympWake_idx]
  
  AsympSleep_idx<-which(allparams$subparameter == "AsympSleep")
  AsympSleep<-allparams$value[AsympSleep_idx]
  
  TauWake_idx<-which(allparams$subparameter == "TauWake")
  TauWake<-allparams$value[TauWake_idx]
  
  TauSleep_idx<-which(allparams$subparameter == "TauSleep")
  TauSleep<-allparams$value[TauSleep_idx]
  
  
  # Initial position
  if (object@initmod == "Fixed"){y.init = object@initpos}
  if (object@initmod == "Free"){y.init = as.numeric(params["start_pos"])}
  
  # time
  time <- object@SWdist$Time
  
  # Wake
  Wake <- object@SWdist$Wake
  
  # Sleep
  Sleep <- object@SWdist$Sleep
  
  # Fit
  out<-SWDMr:::CptSprocess(Sleep=Sleep,Wake=Wake,time=time,U=AsympWake,L=AsympSleep,t_w=TauWake,t_s=TauSleep,init=y.init)
  
  outl<-list(time=time,y1=out[-1])
 
  return(outl)
  
})



setMethod("SWDMrStats",signature="SWDMr_ProcS", function(object,fitted,FittingValue="RSS",detailed=F){
  
  predv<-approxfun(fitted$time,fitted$y1)
  predval<-predv(object@Gexp$Time)
  
  idx<- ! is.na(predval) & ! is.na(object@Gexp[,object@VarExp])
  
  GeneExp<-object@Gexp[idx,object@VarExp]
  predval<-predval[idx]
  
  n <- length(GeneExp)
  
  # Return only RSS | NLL | LL 
  if (detailed == F){
    if (length(predval) == 0){
      return(Inf)
    }else{
      RSS<-sum((GeneExp-predval)^2)
    }
    
    # Return RSS 
    if (FittingValue == "RSS"){
      return(list(val=RSS,var=RSS/n))
      
      # Return Negative Log likelihood or Log Likelihood  
    }else{
      
      # biased estimator of sigma^2
      var<-RSS/n
      
      # Compute Negative log likelihood
      NLL <- (n/2)*(log(2*pi)+log(var)+1)
      
      if (FittingValue == "NLL"){
        return(list(val=NLL,var=var))
      }else if (FittingValue == "LL") {
        return(list(val=-1*NLL,var=var))
      }
    }
    
    
  }else{
    
    # Number of parameters
    k <- nrow(GetFreeFixedParams(object)$FreeParams)
    # residuals
    residualsV<-(GeneExp-predval)
    # Residuals sum of square
    RSS<-sum(residualsV^2)
    # biased estimator of sigma^2
    var<-RSS/n
    # Negative Log Likelihood
    NLL<- (n/2)*(log(2*pi)+log(var)+1)
    # Bayesian Information Criterion
    BIC <- -2*(-NLL)+(k)*log(n)
    # Akaike information criterion
    AIC <- 2*(k) - 2*(-NLL)
    
    # FLAT MODEL
    lm_flat<-lm(GeneExp~1)
    RSS_flat<-sum(resid(lm_flat)^2)
    var_flat<-RSS_flat/n
    NLL_flat<- (n/2)*(log(2*pi)+log(var_flat)+1)
    BIC_flat <- -2*(-NLL_flat)+(1)*log(n)
    
    # Bayes Factor
    BF_DDHOvFlat<-exp((BIC_flat-BIC)/2)
    
    # Kendal Tau
    KendalTau<-cor(GeneExp,predval,method="kendall")
    
    return(list(stats = as.data.frame(list(Variable=object@VarExp,RSS=RSS,NLL=NLL,
                                           BIC=BIC,BIC_flat=BIC_flat,BayesFactor=BF_DDHOvFlat,
                                           AIC=AIC,n=n,k=k,ErrorVariance=var, KendalTau = KendalTau)), residuals = residualsV, fitted = predval ))
    
  }
  
})


setMethod("SWDMrGetEvalFun",signature="SWDMr_ProcS", function(object){
  
  # Return an objective function
  objfun<-function(params){
    
    # Step 1. Fit function given params
    out<-SWDMrFit(object,params)
    
    # Step 2. Get fitting value
    stat<-SWDMrStats(object,out,FittingValue = object@FittingValue)
    
    
    return(stat$val)
    
  }
  
  return(objfun)
  
})


setMethod("GetAllParams",signature="SWDMr_ProcS",function(object,params){
  
  paramofmodel<-GetFreeFixedParams(object)
  expectedparams<-paramofmodel$FreeParams[,"subparameter"]
  if (any(! expectedparams %in% names(params))){
    stop("Missing values in parameter given")
  }
  paramofmodel$FreeParams$value<-as.numeric(params[paramofmodel$FreeParams$subparameter])
  allparams<-rbind(paramofmodel$FreeParams,paramofmodel$FixedParams)
  allparams$value<-as.numeric(allparams$value)
  return(allparams)
})

setMethod("StatsPerTimePoint",signature="SWDMr_ProcS", function(object){
  meansd<-matrix(ncol=3,nrow=length(unique(object@Gexp$Time)))
  colnames(meansd)<-c("mean","sd","se")
  rownames(meansd)<-unique(object@Gexp$Time)
  for (i in unique(object@Gexp$Time)){
    vals<-object@Gexp[object@Gexp$Time == i,object@VarExp]
    meansd[as.character(i),1]<-mean(vals)
    meansd[as.character(i),2]<-sd(vals)
    meansd[as.character(i),3]<-sd(vals)/sqrt(length(vals))
  }
  return(meansd)
})