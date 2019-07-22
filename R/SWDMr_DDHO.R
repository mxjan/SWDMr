############# SET CLASS ################

setClass(
  # Set the name for the class
  Class = "SWDMr_DDHO",
  
  contains="SWDMr",
  
  # Define the slots
  representation = representation(
    
    VarExp = "character",
    
    intercept = "numeric",
    omega = "numeric",
    loggamma = "numeric",
    
    Forces = "list",
    
    AddEffects = "list",
    
    SinForce = "logical",
    AmpSin = "numeric",
    PhiSin = "numeric",
    PerSin = "numeric",
    
    initmod = "character",
    initpos = "numeric",
    initspeed = "numeric",
    
    PenalizeUnstableFit = "logical",
    PredictedValueInterval = "numeric",
    StabilityDayCheck = "numeric",
    
    
    FittingValue = "character",
    
    verbose = "numeric"
  ),
  
  prototype=list(
    initmod = "Free",
    initpos = 0,
    initspeed = 0,
    
    SinForce = F,
    AmpSin = 0,
    PhiSin = 0,
    PerSin = 0,
    
    PenalizeUnstableFit = F,
    
    FittingValue = "RSS",
    
    verbose = 1
  ),
  
  validity=SWDMr_DDHOcheck
  
)

############# Summaries ################

setGeneric("GetFreeFixedParams", function(object)
  standardGeneric("GetFreeFixedParams") )

setMethod("GetFreeFixedParams","SWDMr_DDHO",function(object) {
  
  freeparams<-matrix(ncol=3,nrow=0)
  fixedparams<-matrix(ncol=3,nrow=0)
  #freeparams<-data.frame(parameter = character(0),subparameter=character(0),value=numeric(0),stringsAsFactors = F)
  #fixedparams<-data.frame(parameter = character(0),subparameter=character(0),value=numeric(0),stringsAsFactors = F)
  
  # Core params
  for (param in c("intercept","omega","loggamma")){
    if ( length(slot(object,param)) == 0){
      #freeparams<-rbind(freeparams,data.frame(parameter="Core parameter",subparameter=param,value=NA,stringsAsFactors = F))
      freeparams<-rbind(freeparams,c("Core parameter",param,NA))
    }else{
      #fixedparams<-rbind(fixedparams,data.frame(parameter="Core parameter",subparameter=param,value=slot(object,param),stringsAsFactors = F))
      fixedparams<-rbind(fixedparams,c("Core parameter",param,slot(object,param)))
    }
    #if ( length(slot(object,param)) == 0){freeparams[param] = T}else{fixedparams[param]<-slot(object,param)}
  }
  
  # Forces
  if (length(slot(object,"Forces")) > 0){
    for (ForceName in names(slot(object,"Forces"))){
      if (length(object@Forces[[ForceName]]) == 0){
        freeparams<-rbind(freeparams,c("Forces",ForceName,NA))
        #freeparams<-rbind(freeparams,data.frame(parameter="Forces",subparameter=ForceName,value=NA,stringsAsFactors = F))
      }else{
        #fixedparams<-rbind(fixedparams,data.frame(parameter="Forces",subparameter=ForceName,value=object@Forces[[ForceName]],stringsAsFactors = F))
        fixedparams<-rbind(fixedparams,c("Forces",ForceName,object@Forces[[ForceName]]))
      }
    }
  }
  
  # Add Effect
  if (length(slot(object,"AddEffects")) > 0){
    for (AddEffect in names(slot(object,"AddEffects"))){
      if (length(object@AddEffects[[AddEffect]]) == 0){
        freeparams<-rbind(freeparams,c("AddEffects",AddEffectAddEffect,NA))
        #freeparams<-rbind(freeparams,data.frame(parameter="AddEffects",subparameter=AddEffect,value=NA,stringsAsFactors = F))
      }else{
        #fixedparams<-rbind(fixedparams,data.frame(parameter="AddEffects",subparameter=AddEffect,value=object@AddEffects[[AddEffect]],stringsAsFactors = F))
        fixedparams<-rbind(fixedparams,c("AddEffects",AddEffect,object@AddEffects[[AddEffect]]))
      }
    }
  }
  
  # Sin force
  if (slot(object,"SinForce") == T){
    for (SinFparam in c("AmpSin","PhiSin","PerSin")){
      if (length(slot(object,SinFparam)) == 0){
        freeparams<-rbind(freeparams,c("SinForce",SinFparam,NA))
        #freeparams<-rbind(freeparams,data.frame(parameter="SinForce",subparameter=SinFparam,value=NA,stringsAsFactors = F))
      }else{
        fixedparams<-rbind(fixedparams,c("SinForce",SinFparam,slot(object,SinFparam)))
        #fixedparams<-rbind(fixedparams,data.frame(parameter="SinForce",subparameter=SinFparam,value=slot(object,SinFparam),stringsAsFactors = F))
      }
    }
  }
  
  # Yinit mode
  if (object@initmod == "Free"){
    freeparams<-rbind(freeparams,c("Yinit","start_pos",NA))
    freeparams<-rbind(freeparams,c("Yinit","start_speed",NA))
  }
  if (object@initmod == "Fixed"){
    fixedparams<-rbind(fixedparams,c("Yinit","start_pos",object@initpos))
    fixedparams<-rbind(fixedparams,c("Yinit","start_speed",object@initspeed))
  }
  
  freeparams<-data.frame(freeparams,stringsAsFactors = F)
  colnames(freeparams)<-c("parameter","subparameter","value")
  fixedparams<-data.frame(fixedparams,stringsAsFactors = F)
  colnames(fixedparams)<-c("parameter","subparameter","value")
  
  return(list(FreeParams = freeparams,FixedParams = fixedparams))
  
})


setGeneric("ShowFreeParams", function(object)
  standardGeneric("ShowFreeParams") )
setMethod("ShowFreeParams","SWDMr_DDHO",function(object) {
  
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


setMethod("summary", "SWDMr_DDHO", function(object) {
  cat("~~~~~~~~ This is a S4 SWDMr_DDHO object ~~~~~~~~ \n\n")
  cat("Display the current setting for your fitting\n\n")
  ShowFreeParams(object)
})


############# Set model free/fixed parameters ################

setMethod("FixNaturalFrequency",signature="SWDMr_DDHO", function(object,value){
 object@omega <- value 
 return(object)
})

setMethod("FixDamping",signature="SWDMr_DDHO", function(object,value){
  object@loggamma <- value 
  return(object)
})


setMethod("FixIntercept",signature="SWDMr_DDHO", function(object,value){
  object@intercept <- value 
  return(object)
})


setMethod("AddForce",signature="SWDMr_DDHO", function(object,ForceName,value=numeric(0)){
  
  if (! ForceName %in% colnames(object@SWdist)){
   stop("Force not in SWdist data.frame !") 
  }
  object@Forces[[ForceName]] <- value 
  return(object)
})

setMethod("AddEffect",signature="SWDMr_DDHO", function(object,AddName,value=numeric(0)){
  
  if (! AddName %in% colnames(object@SWdist)){
    stop("Force not in SWdist data.frame !") 
  }
  object@AddEffects[[AddName]] <- value 
  return(object)
})

setMethod("SetYinitMode",signature="SWDMr_DDHO", function(object,mode="Free",values=numeric(0)){
  
  if (mode == "Free"){
    object@initmod <- "Free"
    object@initpos <- numeric(0)
    object@initspeed <- numeric(0)
  } else if (mode == "Fixed"){
    object@initmod <- "Fixed"
    object@initpos <- values[1]
    object@initspeed <- values[2]
  } else if (mode == "CircadianFit"){
    object@initmod <- "CircadianFit"
    
    TimeCircaidx<-object@Gexp$Time >= values[1] & object@Gexp$Time <= values[2]
    TimeCirca<-object@Gexp$Time[TimeCircaidx] # Time in given interval
    GexpCirca<-object@Gexp[,object@VarExp][TimeCircaidx] # Variable in given interval
    
    omega<-(2*pi)/24
    
    reslm <- lm(GexpCirca ~ sin(omega*TimeCirca)+cos(omega*TimeCirca)) # Fit
    t1<-object@SWdist[1,"Time"]
    t2<-object@SWdist[2,"Time"]
    timestep<-t2-t1
    t0<-t1-timestep
    
    post0<-as.numeric(reslm$coefficients[1]+reslm$coefficients[2]*sin(omega*t0)+ # position at t0
                        reslm$coefficients[3]*cos(omega*t0))
    # postm1<-as.numeric(reslm$coefficients[1]+reslm$coefficients[2]*sin(((2*pi)/24)*(t0-timestep))+ # position at t0-t1
    #                      reslm$coefficients[3]*cos(((2*pi)/24)*(t0-timestep)))
    # speedt0<- (post0 - postm1) / timestep
    # 
    
    if (any(c(is.na(reslm$coefficients[2]),is.na(reslm$coefficients[3])))){stop("Fit failed in the given interval")}
    
    # Initial position @ T0
    object@initpos <- post0
    # First derivative @ T0
    object@initspeed <- as.numeric(reslm$coefficients[2] * (cos(omega*t0)) * omega + reslm$coefficients[3] * (-sin(omega*t0)) * omega)
    
  } else if (mode == "Intercept_0"){
    object@initmod <- "Intercept_0"
    object@initpos <- 0
    object@initspeed <- 0
  } else {
    
    stop("Mode not found")
    
  }
  
  return(object)
})


setMethod("AddSinF",signature="SWDMr_DDHO", function(object,FixAmp=numeric(0),FixPhi=numeric(0),FixPer=numeric(0)){
  
  object@SinForce<-T
  object@AmpSin<-FixAmp
  object@PhiSin<-FixPhi
  object@PerSin<-FixPer

  return(object)
})


setMethod("SetFittingValue",signature="SWDMr_DDHO", function(object,value = "RSS"){
  
  if (! value %in% c("RSS","NLL","LL")){
    stop("The returned value should be RSS, NLL or LL")
  }
  
  object@FittingValue <- value
  
  return(object)
})


setMethod("PenalizeUnstableFit",signature="SWDMr_DDHO", function(object,value = T,PredictedValueInterval,StabilityDayCheck){
  
  # Set penalization to True
  object@PenalizeUnstableFit<-T
  
  if (length(PredictedValueInterval) != 2){stop("The PredictedValueInterval provided must be 2 values (e.g. c(0,24))")}
  object@PredictedValueInterval<-PredictedValueInterval
  
  object@StabilityDayCheck<-StabilityDayCheck
  
  return(object)
})

setGeneric("SumForces",  function(object,params,allparams=NULL)
  standardGeneric("SumForces") )
setMethod("SumForces",signature="SWDMr_DDHO", function(object,params,allparams=NULL){
  
  if (is.null(allparams)){
    paramofmodel$FreeParams$value<-params[paramofmodel$FreeParams$subparameter]
    allparams<-rbind(paramofmodel$FreeParams,paramofmodel$FixedParams)
    allparams$value<-as.numeric(allparams$value)
  }
  
  idxForces<-which(allparams$parameter == "Forces")
  if (length(idxForces) > 0){
    Forces<-allparams[idxForces,"value"]
    SWf<-object@SWdist[,allparams[idxForces,"subparameter"],drop=F]
    SWf<-as.matrix(mapply(`*`,SWf,Forces))
    if (length(Forces) > 1){
      Force<-rowSums(SWf)
    }else{
      Force<-SWf[,1]
    }
  }else{
    Force<-rep(0,nrow(GEMS$SWdistr))
  }
  
  return(Force)
  
})

setMethod("GetFit",signature="SWDMr_DDHO", function(object,params){
  
  # Control parameters given
  paramofmodel<-GetFreeFixedParams(object)
  expectedparams<-paramofmodel$FreeParams[,"subparameter"]
  if (any(! expectedparams %in% names(params))){
    stop("Missing values in parameter given")
  }
  paramofmodel$FreeParams$value<-params[paramofmodel$FreeParams$subparameter]
  allparams<-rbind(paramofmodel$FreeParams,paramofmodel$FixedParams)
  allparams$value<-as.numeric(allparams$value)
  
  # Get Core params
  intercept_idx<-which(allparams$subparameter == "intercept")
  intercept<-allparams$value[intercept_idx]
  
  omega_idx<-which(allparams$subparameter == "omega")
  omega<-allparams$value[omega_idx]
  
  gamma_idx<-which(allparams$subparameter == "loggamma")
  gamma<-exp(allparams$value[gamma_idx])
  
  # Initial position
  if (object@initmod == "Fixed" || object@initmod == "CircadianFit" || object@initmod == "Free"){y.init = c(object@initpos-intercept, object@initspeed)}
  if (object@initmod == "Intercept_0"){y.init = c(0, 0)}
  
  # Sum forces
  force<-SumForces(object,params,allparams)
  
  # time
  time <- object@SWdist$Time
  
  # string constant
  k<-omega*omega
  
  # Sin F
  if (object@SinForce == T){
    AmpSin<-allparams$value[which(allparams$subparameter == "AmpSin")]
    PhiSin<-allparams$value[which(allparams$subparameter == "PhiSin")]
    PerSin<-allparams$value[which(allparams$subparameter == "PerSin")]
  }else{
    AmpSin<-0
    PhiSin<-0
    PerSin<-24
  }
  
  # Fit
  out<-SDDHO_SinF_RungeKutta(y=y.init,time = time,force = force,gamma = gamma, k=k, AmpSin = AmpSin, PhiSin = PhiSin, PerSin = PerSin) # Solve ODE, see src/sddho_sinforce_RK4.cpp
  
  # Add intercept to data
  out$y1 <- out$y1 + intercept
  
  # Add additive parameters
  Addparams_idx<-which(allparams$parameter == "AddEffects")
  if (length(Addparams_idx) > 0){
    Addeffs<-allparams[Addparams_idx,"value"]
    SWf<-object@SWdist[,allparams[Addparams_idx,"subparameter"],drop=F]
    SWf<-as.matrix(mapply(`*`,SWf,Addeffs))
    if (length(Addeffs) > 1){
      Addeff<-rowSums(SWf)
    }else{
      Addeff<-SWf[,1]
    }
    out$y1[2:length(out$y1)]<-out$y1[2:length(out$y1)]+Addeff
  }
  
  return(out)
  
})

