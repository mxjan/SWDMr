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
  
  freeparams<-data.frame(parameter = character(0),subparameter=character(0),value=numeric(0),stringsAsFactors = F)
  fixedparams<-data.frame(parameter = character(0),subparameter=character(0),value=numeric(0),stringsAsFactors = F)
  
  # Core params
  for (param in c("intercept","omega","loggamma")){
    if ( length(slot(object,param)) == 0){
      freeparams<-rbind(freeparams,data.frame(parameter="Core parameter",subparameter=param,value=NA,stringsAsFactors = F))
    }else{
      fixedparams<-rbind(fixedparams,data.frame(parameter="Core parameter",subparameter=param,value=slot(object,param),stringsAsFactors = F))
    }
    #if ( length(slot(object,param)) == 0){freeparams[param] = T}else{fixedparams[param]<-slot(object,param)}
  }
  
  # Forces
  if (length(slot(object,"Forces")) > 0){
    for (ForceName in names(slot(object,"Forces"))){
      if (length(object@Forces[[ForceName]]) == 0){
        freeparams<-rbind(freeparams,data.frame(parameter="Forces",subparameter=ForceName,value=NA,stringsAsFactors = F))
      }else{
        fixedparams<-rbind(fixedparams,data.frame(parameter="Forces",subparameter=ForceName,value=object@Forces[[ForceName]],stringsAsFactors = F))
      }
    }
  }
  
  # Add Effect
  if (length(slot(object,"AddEffects")) > 0){
    for (AddEffect in names(slot(object,"AddEffects"))){
      if (length(object@AddEffects[[AddEffect]]) == 0){
        freeparams<-rbind(freeparams,data.frame(parameter="AddEffects",subparameter=AddEffect,value=NA,stringsAsFactors = F))
      }else{
        fixedparams<-rbind(fixedparams,data.frame(parameter="AddEffects",subparameter=AddEffect,value=object@AddEffects[[AddEffect]],stringsAsFactors = F))
      }
    }
  }
  
  # Sin force
  if (slot(object,"SinForce") == T){
    for (SinFparam in c("AmpSin","PhiSin","PerSin")){
      if (length(slot(object,SinFparam)) == 0){
        freeparams<-rbind(freeparams,data.frame(parameter="SinForce",subparameter=SinFparam,value=NA,stringsAsFactors = F))
      }else{
        fixedparams<-rbind(fixedparams,data.frame(parameter="SinForce",subparameter=SinFparam,value=slot(object,SinFparam),stringsAsFactors = F))
      }
    }
  }
  
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

#' Fix natural frequency of the model
#'
#' @param object An SWDMr_DDHO object
#' @param value natural frequency 
#' @export
#' @docType methods
#' @examples
#' model<-FixNaturalFrequency(model,2*pi/24)
setGeneric("FixNaturalFrequency", function(object,value)
  standardGeneric("FixNaturalFrequency") )

setMethod("FixNaturalFrequency",signature="SWDMr_DDHO", function(object,value){
 object@omega <- value 
 return(object)
})

#' Fix damping of the model
#'
#' @param object An SWDMr_DDHO object
#' @param value log gamma 
#' @export
#' @docType methods
#' @examples
#' model<-FixNaturalFrequency(model,-2)
setGeneric("FixDamping", function(object,value)
  standardGeneric("FixDamping") )

setMethod("FixDamping",signature="SWDMr_DDHO", function(object,value){
  object@loggamma <- value 
  return(object)
})


#' Fix intercept of the model
#'
#' @param object An SWDMr_DDHO object
#' @param value intercept
#' @export
#' @docType methods
#' @examples
#' model<-FixIntercept(model,0)
setGeneric("FixIntercept", function(object,value)
  standardGeneric("FixIntercept") )

setMethod("FixIntercept",signature="SWDMr_DDHO", function(object,value){
  object@intercept <- value 
  return(object)
})


#' Add Force to the model
#'
#' @param object An SWDMr_DDHO object
#' @param ForceName character present in SWdist data.frame
#' @param value if this is a free parameter, let it be empty.
#' @export
#' @docType methods
#' @examples
#' model<-AddForce(model,"Wake")
setGeneric("AddForce", function(object,ForceName,value=numeric(0))
  standardGeneric("AddForce") )

setMethod("AddForce",signature="SWDMr_DDHO", function(object,ForceName,value=numeric(0)){
  
  if (! ForceName %in% colnames(object@SWdist)){
   stop("Force not in SWdist data.frame !") 
  }
  object@Forces[[ForceName]] <- value 
  return(object)
})


#' Add Additive factor to the model
#'
#' @param object An SWDMr_DDHO object
#' @param AddName character present in SWdist data.frame
#' @param value if this is a free parameter, let it be empty.
#' @export
#' @docType methods
#' @examples
#' model<-AddEffect(model,"SD")
setGeneric("AddEffect", function(object,AddName,value=numeric(0))
  standardGeneric("AddEffect") )

setMethod("AddEffect",signature="SWDMr_DDHO", function(object,AddName,value=numeric(0)){
  
  if (! AddName %in% colnames(object@SWdist)){
    stop("Force not in SWdist data.frame !") 
  }
  object@AddEffects[[AddName]] <- value 
  return(object)
})


#' Set the start position and speed of the oscillator
#'
#' @param object An SWDMr_DDHO object
#' @param mode character can be Free (default), Fixed, Intercept_0 or CircadianFit
#' @param values if the mode is "Fixed", the values is start and speed of the oscillator (e.g. c(0,0) for start at 0 and speed of 0. If the mode is CircadianFit, values are the time interval for the circadian fit (e.g. c(0,48) for time 0 to time 48)
#' @export
#' @docType methods
#' @examples
#' model<-SetYinitMode(model,"Intercept_0")
setGeneric("SetYinitMode", function(object,mode="Free",values=numeric(0))
  standardGeneric("SetYinitMode") )

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
    object@initspeed <- reslm$coefficients[2] * (cos(omega*t0)) * omega + reslm$coefficients[3] * (-sin(omega*t0)) * omega
    
  } else if (mode == "Intercept_0"){
    object@initmod <- "Intercept_0"
    object@initpos <- 0
    object@initspeed <- 0
  } else {
    
    stop("Mode not found")
    
  }
  
  return(object)
})


#' Add Sin Force to the model
#'
#' @param object An SWDMr_DDHO object
#' @param FixAmp Fix Sin Force Amplitude to the given value
#' @param FixPhi Fix the Sin Force phase to the given value
#' @param FixPer Fix the period (in hour) of the Sin Force
#' @export
#' @docType methods
#' @examples
#' model<-AddSinF(model,FixPer = 24)
setGeneric("AddSinF", function(object,FixAmp=numeric(0),FixPhi=numeric(0),FixPer=numeric(0))
  standardGeneric("AddSinF") )

setMethod("AddSinF",signature="SWDMr_DDHO", function(object,FixAmp=numeric(0),FixPhi=numeric(0),FixPer=numeric(0)){
  
  object@SinForce<-T
  object@AmpSin<-FixAmp
  object@PhiSin<-FixPhi
  object@PerSin<-FixPer

  return(object)
})


#' Fitting value
#'
#' @param object An SWDMr_DDHO object
#' @param value Set the fitting value return, can be Residual Sum of Square (RSS), Log likelihood (LL) or Negative Log Likelihood (NLL)
#' @export
#' @docType methods
#' @examples
#' model<-SetFittingValue(model,value = "RSS")
setGeneric("SetFittingValue", function(object,value = "RSS")
  standardGeneric("SetFittingValue") )

setMethod("SetFittingValue",signature="SWDMr_DDHO", function(object,value = "RSS"){
  
  if (! value %in% c("RSS","NLL","LL")){
    stop("The returned value should be RSS, NLL or LL")
  }
  
  object@FittingValue <- value
  
  return(object)
})

#' Unstable Fit penalization
#' 
#' Penalize the fit that are unstable through replicated points
#'
#' @param object An SWDMr_DDHO object
#' @param value Set the penalization method (TRUE or FALSE)
#' @param PredictedValueInterval A two value vector for min and max value observed
#' @param StabilityDayCheck The number of day control for similiar min and max value compared to PredictedValueInterval
#' @export
#' @docType methods
#' @examples
#' # Example of 20 day of baseline replicated
#' model<-ReplicateDrivingForce(model,c(0,24),20)
#' # Then we control that the last 10 replicated day have similar min and max value compared to true baseline
#' model<-PenalizeUnstableFit(model,value = T,PredictedValueInterval = c(0,24), StabilityDayCheck = 10)
setGeneric("PenalizeUnstableFit", function(object,value = T,PredictedValueInterval,StabilityDayCheck)
  standardGeneric("PenalizeUnstableFit") )

setMethod("PenalizeUnstableFit",signature="SWDMr_DDHO", function(object,value = T,PredictedValueInterval,StabilityDayCheck){
  
  # Set penalization to True
  object@PenalizeUnstableFit<-T
  
  if (length(PredictedValueInterval) != 2){stop("The PredictedValueInterval provided must be 2 values (e.g. c(0,24))")}
  object@PredictedValueInterval<-PredictedValueInterval
  
  object@StabilityDayCheck<-StabilityDayCheck
  
  return(object)
})

