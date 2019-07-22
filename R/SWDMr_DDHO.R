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
    AmpSin = "numeric",
    PhiSin = "numeric",
    PerSin = "numeric",
    Forces = "list",
    AddEffects = "list",
    
    initmod = "character",
    initpos = "numeric",
    initspeed = "numeric",
    
    verbose = "numeric"
  ),
  
  prototype=list(
    initmod = "Free",
    initpos = 0,
    initspeed = 0,
    verbose = 1
  ),
  
  validity=SWDMr_DDHOcheck
  
)

############# Summaries ################

setMethod("summary", "SWDMr_DDHO", function(object) {
  cat("This is a SWDMr_DDHO object\n")
  cat("This are you default free parameters\n")
})

setMethod("show", "SWDMr_DDHO", function(object) {
  summary(object)
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
