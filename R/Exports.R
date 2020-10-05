#' Initiate DDHO model
#'
#' Initiate a driven damped harmonic oscillator model
#'
#' @param object A SWDMr object 
#' @param VarExp The name of the variable to explain by the model
#' @param UseDampingRatio The core parameters of the model change for intercept, damping ratio and natural frequency instead of intercept, damping constant and natural frequency
#' 
#' @return SWDMr model object
#'
#' @examples
#' model <- initDDHOmodel(swdmr,VarExp = "Arntl")
#'
#'@export
initDDHOmodel<-function(object,VarExp, UseDampingRatio = F){
  new(Class = "SWDMr_DDHO", SWdist = object@SWdist, 
      Gexp = object@Gexp[c(VarExp,"Time")],
      VarExp = VarExp,
      UseDampingRatio = UseDampingRatio)
}


#' Create a SWDMr object
#'
#' Create a SWDMr object using a sleep-wake data.frame and a gene expression data.frame
#'
#' @param SWdist sleep-wake dataframe, must contain the column: Wake, Sleep and Time
#' @param Gexp Gene expression dataframe which must contain genes as column and a Time column corresponding to the sleep - wake dataframe
#'
#' @return SWDMr object
#'
#' @examples
#' swdmr <- SWDMr(SWdist=SWdf, Gexp=Gexpdf)
#'
#'@export
SWDMr <- function(SWdist,Gexp,verbose=1){
  new (Class="SWDMr",SWdist=SWdist,Gexp=Gexp,verbose=verbose)
}

#' Unstable Fit penalization
#' 
#' Penalize the fit that are unstable through replicated points
#'
#' @param object An SWDMr_DDHO object
#' @param value Set the penalization method (TRUE or FALSE)
#' @param PredictedValueInterval A two value vector for min and max value observed
#' @param StabilityDayCheck The number of day control for similiar min and max value compared to PredictedValueInterval
#' @param weight Weight applied on unstable points
#' @export
#' @docType methods
#' @examples
#' # Example of 20 day of baseline replicated
#' model<-ReplicateDrivingForce(model,c(0,24),20)
#' # Then we control that the last 10 replicated day have similar min and max value compared to true baseline
#' model<-PenalizeUnstableFit(model,value = T,PredictedValueInterval = c(0,24), StabilityDayCheck = 10)
setGeneric("PenalizeUnstableFit", function(object,value = T,PredictedValueInterval,StabilityDayCheck,weight=1)
  standardGeneric("PenalizeUnstableFit") )


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


#' Replicate baseline driving force
#' 
#' Replicate the baseline (given time interval) and happend it at the begining of your Force dataset
#'
#' @param object An SWDMr_DDHO object
#' @param interval A two value time interval that is a 24h interval (time given is included, i.e. <= or >=)
#' @param Nrep Number of repliated time interval
#' @export
#' @docType methods
#' @examples
#' swdmr<-ReplicateDrivingForce(swdmr,c(0,24),10)
setGeneric("ReplicateDrivingForce", function(object,interval,Nrep)
  standardGeneric("ReplicateDrivingForce") )

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

#' Fit data given model
#' @param object An SWDMr model object
#' @param params a named vector with required free parameter value
#' @export
#' @docType methods
#' @examples 
#' fitted<-SWDMrFit(model,params)
setGeneric("SWDMrFit", function(object,params)
  standardGeneric("SWDMrFit") )


#' Return statistics of the fitting
#' @param object An SWDMr model object
#' @param fitted fitted value
#' @param FittingValue What value to return for the fitting (RSS or NLL or LL)
#' @param detailed If set to True, return all value (RSS,NLL,LL,BIC), N, nparams
#' @export
#' @docType methods
#' @examples 
#' fitted<-SWDMrFit(model,params)
#' stats<-SWDMrStats(model,fitted)
setGeneric("SWDMrStats", function(object,fitted,FittingValue="RSS",detailed=F)
  standardGeneric("SWDMrStats") )

#' Return function to evaluate
#' @param object An SWDMr model object
#' @export
#' @docType methods
#' @examples 
#' objfun<-SWDMrGetEvalFun(model)
#' optimx(objfunc,params)
setGeneric("SWDMrGetEvalFun", function(object)
  standardGeneric("SWDMrGetEvalFun") )



##########################################################################
######################## PROCESS S #######################################

#' Initiate Process S model
#'
#' Initiate a Process S model of sleep
#'
#' @param object A SWDMr object 
#' @param VarExp The name of the variable to explain by the model
#' 
#' @return SWDMr model object
#'
#' @examples
#' model <- initProcessSmodel(swdmr,VarExp = "Arntl")
#'
#'@export
initProcessSmodel<-function(object,VarExp){
  new(Class = "SWDMr_ProcS", SWdist = object@SWdist, 
      Gexp = object@Gexp[c(VarExp,"Time")],
      VarExp = VarExp)
}


#' Fix Asymptote of wake in the model
#'
#' @param object An SWDMr_ProcS object
#' @param value top value reached by wake
#' @export
#' @docType methods
#' @examples
#' model<-FixAsymptWake(model,300)
setGeneric("FixAsymptWake", function(object,value)
  standardGeneric("FixAsymptWake") )


#' Fix Asymptote of sleep in the model
#'
#' @param object An SWDMr_ProcS object
#' @param value top value reached by wake
#' @export
#' @docType methods
#' @examples
#' model<-FixAsymptSleep(model,80)
setGeneric("FixAsymptSleep", function(object,value)
  standardGeneric("FixAsymptSleep") )

#' Fix time constant of wake in the model
#'
#' @param object An SWDMr_ProcS object
#' @param value top value reached by wake
#' @export
#' @docType methods
#' @examples
#' model<-FixTimeConstWake(model,16)
setGeneric("FixTimeConstWake", function(object,value)
  standardGeneric("FixTimeConstWake") )


#' Fix Time constant of sleep in the model
#'
#' @param object An SWDMr_ProcS object
#' @param value top value reached by wake
#' @export
#' @docType methods
#' @examples
#' model<-FixTimeConstSleep(model,8)
setGeneric("FixTimeConstSleep", function(object,value)
  standardGeneric("FixTimeConstSleep") )





######################### GENERIC
setGeneric(name="ShowFreeParams", def = function(object,...)
  standardGeneric("ShowFreeParams"))
setGeneric("GetFreeFixedParams", function(object,...)
  standardGeneric("GetFreeFixedParams") )
setGeneric("GetAllParams", function(object,params)
  standardGeneric("GetAllParams") )
setGeneric("SumForces",  function(object,params,allparams=NULL)
  standardGeneric("SumForces") )
setGeneric("AddUnstableFitPenalization",  function(object,fitted,FittingValue,val,var,weight=1)
  standardGeneric("AddUnstableFitPenalization") )
setGeneric("StatsPerTimePoint", function(object) 
  standardGeneric("StatsPerTimePoint") )
setGeneric("PctAbsForceApplied", function(object,params,pct=T) 
  standardGeneric("PctAbsForceApplied") )
setGeneric("AllForceApplied", function(object,params) 
  standardGeneric("AllForceApplied") )
setGeneric("MatchPoints", function(object)
  standardGeneric("MatchPoints") )

