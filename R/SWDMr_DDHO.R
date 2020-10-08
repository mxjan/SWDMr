############# SET CLASS ################

setClass(
  # Set the name for the class
  Class = "SWDMr_DDHO",
  
  contains="SWDMr",
  
  # Define the slots
  representation = representation(
    
    # Name of variable explained
    VarExp = "character",
    
    # Core parameters of damped oscillator
    intercept = "numeric", # intercept of the model (0 position of the oscillator)
    omega = "numeric", # natural frequency 
    loggamma = "numeric", # damping constant (in log as it is easier to optimize with derivative free method)
    ## Use damping ration for optimization
    UseDampingRatio = "logical", # If True, then use only omega and get loggamma = log(dampratio/2*omega)
    dampratio = "numeric", # damping ration (gamma/2*omega), > 1 = overdamping, 1 = critical damping,  < 1 = underdamnping
    
    # List of forces applied to the model
    Forces = "list", 
    
    # List of additive effects applied to the model
    AddEffects = "list",
    
    # A second force added (True or False) in the form of a sin-wave 
    SinForce = "logical", # T || F
    AmpSin = "numeric", # Amplitude
    PhiSin = "numeric", # Phase
    PerSin = "numeric", # Period
    
    # How the initial value of position and speed are calculated
    initmod = "character",
    initpos = "numeric", # position
    initspeed = "numeric", # speed
    
    # The fit is penalized for unstable fit
    PenalizeUnstableFit = "logical",
    PredictedValueInterval = "numeric",
    StabilityDayCheck = "numeric",
    PenalizeUnstableFitWeight = "numeric",
    
    # Objective function is optimized for FitttingValue, can be RSS or Negative Log Likelihood (NLL)
    FittingValue = "character",
    
    verbose = "numeric"
  ),
  
)

setMethod("initialize","SWDMr_DDHO",function(.Object,initmod = "Free",initpos = 0,initspeed = 0,
                                             SinForce = F, AmpSin = 0,PhiSin = 0,PerSin = 0,
                                             PenalizeUnstableFit = F, 
                                             FittingValue = "RSS", verbose=1,...){
  
  .Object <- callNextMethod(.Object, ...)
  
  .Object@initmod <- initmod
  .Object@initpos <- initpos
  .Object@initpos <- initspeed
  
  .Object@SinForce <- SinForce
  .Object@AmpSin <- AmpSin
  .Object@PhiSin <- PhiSin  
  .Object@PerSin <- PerSin  
  
  .Object@PenalizeUnstableFit <- PenalizeUnstableFit
  
  .Object@FittingValue <- FittingValue
  
  .Object@verbose <- verbose
  
  .Object
})

######################################
############# METHODS ################
######################################


####Return a matrix of free and fixed parameters
setMethod("GetFreeFixedParams","SWDMr_DDHO",function(object,...) {

  # Matrix of free, fixed parameters  
  freeparams<-matrix(ncol=3,nrow=0)
  fixedparams<-matrix(ncol=3,nrow=0)

  # Core params as intercept [equilibrium], omega [nat. freq.], dampratio | loggamma
  if (object@UseDampingRatio == T){
    CoreParams<-c("intercept","omega","dampratio")
  }else{
    CoreParams<-c("intercept","omega","loggamma")
  }
  for (param in CoreParams){
    if ( length(slot(object,param)) == 0){
      freeparams<-rbind(freeparams,c("Core parameter",param,NA))
    }else{
      fixedparams<-rbind(fixedparams,c("Core parameter",param,slot(object,param)))
    }
  }
  
  # add Forces of the model
  if (length(slot(object,"Forces")) > 0){
    for (ForceName in names(slot(object,"Forces"))){
      if (length(object@Forces[[ForceName]]) == 0){
        freeparams<-rbind(freeparams,c("Forces",ForceName,NA))
      }else{
        fixedparams<-rbind(fixedparams,c("Forces",ForceName,object@Forces[[ForceName]]))
      }
    }
  }
  
  # Simple Additive Effects
  if (length(slot(object,"AddEffects")) > 0){
    for (AddEffect in names(slot(object,"AddEffects"))){
      if (length(object@AddEffects[[AddEffect]]) == 0){
        freeparams<-rbind(freeparams,c("AddEffects",AddEffectAddEffect,NA))
      }else{
        fixedparams<-rbind(fixedparams,c("AddEffects",AddEffect,object@AddEffects[[AddEffect]]))
      }
    }
  }
  
  # Sinewave force 
  if (slot(object,"SinForce") == T){
    for (SinFparam in c("AmpSin","PhiSin","PerSin")){
      if (length(slot(object,SinFparam)) == 0){
        freeparams<-rbind(freeparams,c("SinForce",SinFparam,NA))
      }else{
        fixedparams<-rbind(fixedparams,c("SinForce",SinFparam,slot(object,SinFparam)))
      }
    }
  }
  
  # Yinit mode (initial position and speed)
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

# Display the parameter matrix
setMethod("ShowFreeParams","SWDMr_DDHO",function(object,...) {
  
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

# Summary function
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
    
    # Fit a cosine to phenotype, will be used to estimate start and speed of the oscillator
    reslm <- lm(GexpCirca ~ sin(omega*TimeCirca)+cos(omega*TimeCirca))
    t1<-object@SWdist[1,"Time"]
    t2<-object@SWdist[2,"Time"]
    timestep<-t2-t1
    t0<-t1-timestep
    
    post0<-as.numeric(reslm$coefficients[1]+reslm$coefficients[2]*sin(omega*t0)+ # position at t0
                        reslm$coefficients[3]*cos(omega*t0))
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

# Add a force in the form of a sinewave
setMethod("AddSinF",signature="SWDMr_DDHO", function(object,FixAmp=numeric(0),FixPhi=numeric(0),FixPer=numeric(0)){
  object@SinForce<-T
  object@AmpSin<-FixAmp
  object@PhiSin<-FixPhi
  object@PerSin<-FixPer
  return(object)
})

# Value returned by the objective function as Residual Sum of Square (RSS), Negative log-likelihood or Likelihood
setMethod("SetFittingValue",signature="SWDMr_DDHO", function(object,value = "RSS"){
  if (! value %in% c("RSS","NLL","LL")){
    stop("The returned value should be RSS, NLL or LL")
  }
  object@FittingValue <- value
  return(object)
})

# Set a penalization of the objective function for unstable oscillation
setMethod("PenalizeUnstableFit",signature="SWDMr_DDHO", function(object,value = T,PredictedValueInterval,StabilityDayCheck,weight=1){
  
  # Set penalization to True
  object@PenalizeUnstableFit<-T
  
  if (length(PredictedValueInterval) != 2){stop("The PredictedValueInterval provided must be 2 values (e.g. c(0,24))")}
  object@PredictedValueInterval<-PredictedValueInterval
  
  object@StabilityDayCheck<-StabilityDayCheck
  
  object@PenalizeUnstableFitWeight<-weight
  
  return(object)
})


# Add forces given in the model
setMethod("SumForces",signature="SWDMr_DDHO", function(object,params,allparams=NULL){
  
  if (is.null(allparams)){
    allparams<-GetAllParams(object,params)
  }
  idxForces<-which(allparams$parameter == "Forces")
  if (length(idxForces) > 0){
    Forces<-allparams[idxForces,"value"]
    SWf<-object@SWdist[,allparams[idxForces,"subparameter"],drop=F]
    Force <- (as.matrix(SWf) %*% Forces)[,1]
  }else{
    Force<-rep(0,nrow(object@SWdist))
  }
  return(Force)
  
})

# function to return the fitted values given parameters
setMethod("SWDMrFit",signature="SWDMr_DDHO", function(object,params){
  
  # Control parameters given
  allparams<-GetAllParams(object,params)
  
  # Get Core params
  intercept_idx<-which(allparams$subparameter == "intercept")
  intercept<-allparams$value[intercept_idx]
  
  omega_idx<-which(allparams$subparameter == "omega")
  omega<-allparams$value[omega_idx]
  
  if (object@UseDampingRatio == F){
    gamma_idx<-which(allparams$subparameter == "loggamma")
    gamma<-exp(allparams$value[gamma_idx])
  }else{
    dampratio_idx<-which(allparams$subparameter == "dampratio")
    gamma <- allparams$value[dampratio_idx] / (2*allparams$value[omega_idx])
  }

  # Initial position
  if (object@initmod == "Fixed" || object@initmod == "CircadianFit"){y.init = c(object@initpos-intercept, object@initspeed)}
  if (object@initmod == "Intercept_0"){y.init = c(0, 0)}
  if (object@initmod == "Free"){y.init = c(as.numeric(params["start_pos"]),as.numeric(params["start_speed"]))}
  
  
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
  
  # Fit# Solve ODE, see src/sddho_sinforce_RK4.cpp
  out<-SWDMr:::SDDHO_SinF_RungeKutta(y=y.init,time = time,force = force,gamma = gamma, k=k, AmpSin = AmpSin, PhiSin = PhiSin, PerSin = PerSin) 
  
  # Add intercept to data
  out$y1 <- out$y1 + intercept
  
  # Add potential additive effects
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

# Apply unstable oscillation penalization of the fit
setMethod("AddUnstableFitPenalization",signature="SWDMr_DDHO",function(object,fitted,FittingValue,val,var,weight=1){
  
  # Weight
  if (object@PenalizeUnstableFitWeight != 1){
    weight<-object@PenalizeUnstableFitWeight
  }
  
  # Get min and max in intervals
  interval_eval<-fitted$time >= object@PredictedValueInterval[1] & fitted$time <= object@PredictedValueInterval[2]
  sodeinterval<-fitted$y1[interval_eval]
  idxminv<-which(sodeinterval == min(sodeinterval))[1]
  idxmaxv<-which(sodeinterval == max(sodeinterval))[1]
  
  modultimeminv<-fitted$time[interval_eval][idxminv] %% 24
  modultimemaxv<-fitted$time[interval_eval][idxmaxv] %% 24
  
  evaldaytimemin<-cumsum(c(modultimeminv-24,rep(-24,object@StabilityDayCheck)))
  evaldaytimemax<-cumsum(c(modultimemaxv-24,rep(-24,object@StabilityDayCheck)))

  
  obsmin<-fitted$y1[sapply(evaldaytimemin,function(x) which.min(abs(fitted$time-x)))]
  obsmax<-fitted$y1[sapply(evaldaytimemax,function(x) which.min(abs(fitted$time-x)))]
  obs<-c(obsmin,obsmax )
  pred<-c(rep(min(sodeinterval),length(obsmin)),rep(max(sodeinterval),length(obsmax)))
  
  penalization<-sum((pred-obs)^2)*weight
  
  if (FittingValue == "RSS"){
    return(val+penalization)
  }
  else if (FittingValue == "NLL"){
    val <- val+penalization
    return(val)
  }
  else if (FittingValue == "LL"){
    val <- val-penalization
    return(val)
  }
})



setMethod("SWDMrStats",signature="SWDMr_DDHO", function(object,fitted,FittingValue="RSS",detailed=F,match="exact"){
  
  if (match == "exact"){
    idx<-na.omit(object@Match_Tgexp_Tswd)
    GeneExp<-object@Gexp[which(! is.na(object@Match_Tgexp_Tswd)),object@VarExp]
    predval<-fitted$y1[idx]
  }else if (match == "approx"){
    predv<-approxfun(fitted$time,fitted$y1)
    predval<-predv(object@Gexp$Time)
    idx<- ! is.na(predval) & ! is.na(object@Gexp[,object@VarExp])
    GeneExp<-object@Gexp[idx,object@VarExp]
    predval<-predval[idx]
  }

  
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


setMethod("SWDMrGetEvalFun",signature="SWDMr_DDHO", function(object,match="exact"){
  
  # Return an objective function
  objfun<-function(params){
    
    # Step 1. Fit function given params
    out<-SWDMrFit(object,params)
    
    # Step 2. Get fitting value
    stat<-SWDMrStats(object,out,FittingValue = object@FittingValue,match=match)
    
    # Step 3. If Penalization is set to True, compute 
    if (object@PenalizeUnstableFit == T){
      retval <- AddUnstableFitPenalization(object = object,
                                           fitted = out,
                                           FittingValue = object@FittingValue,
                                           val = stat$val,
                                           var = stat$var)
    }else{
      retval<-stat$val
    }
    
    return(retval)
    
  }
  
  return(objfun)
  
})



setMethod("StatsPerTimePoint",signature="SWDMr_DDHO", function(object){
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


setMethod("GetAllParams",signature="SWDMr_DDHO",function(object,params){
  
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

setMethod("SWdrivenValue",signature="SWDMr_DDHO",function(object,params,interval){
  
  params<-SWDMr:::GetAllParams(model,params)
  
  # Time interval given
  idxInterval<- object@SWdist$Time >= (interval[1] - object@tol) & object@SWdist$Time <= (interval[2] + object@tol)
  
  # Get variance of the force applied by sleep & wake
  coefS<-params[params$subparameter == "Sleep","value"]
  coefW<-params[params$subparameter == "Wake","value"]
  Fsw<-object@SWdist$Sleep[idxInterval]*coefS+object@SWdist$Wake[idxInterval]*coefW
  # Get variance of the force applied by the circadian system
  AmpSin<-params[params$subparameter == "AmpSin","value"]
  PhiSin<-params[params$subparameter == "PhiSin","value"]  
  PerSin<-params[params$subparameter == "PerSin","value"]
  Fsin<-AmpSin*sin(2*pi/PerSin*object@SWdist$Time[idxInterval]+PhiSin)
  
  return(var(Fsw)/(var(Fsw)+var(Fsin)))
  
})



setGeneric("ForceApplied", function(object,params) 
  standardGeneric("ForceApplied") )
setMethod("ForceApplied",signature="SWDMr_DDHO",function(object,params){
  
  allparams<-GetAllParams(object,params)
  
  force<-rep(0,length(object@SWdist$Time))
  time<-(object@SWdist$Time)  
  
  if (any(allparams$parameter %in% c("Forces","SinForce"))){
    
    if ("Forces" %in% allparams$parameter){
      for (i in allparams$subparameter[allparams$parameter == "Forces"]){
        force <- force + (object@SWdist[,i]*allparams$value[allparams$subparameter == i])
      }
    }
    
    if ("SinForce" %in% allparams$parameter){
      amp<-allparams$value[allparams$subparameter == "AmpSin"]
      phi<-allparams$value[allparams$subparameter == "PhiSin"]
      per<-allparams$value[allparams$subparameter == "PerSin"]
      force <- force + amp * sin(2*pi/per * time + phi) 
    }
    
    return(list(force=force,time=time))
    
  }else{
    stop("No force applied to the model")
  }
  
})


setMethod("PctAbsForceApplied",signature="SWDMr_DDHO",function(object,params,pct=T){
  
  allparams<-GetAllParams(object,params)
  
  time<-(object@SWdist$Time)  
  
  if (any(allparams$parameter %in% c("Forces","SinForce"))){
    
    mF<-matrix(time,ncol=1)
    colnames(mF)<-"Time"
    
    if ("Forces" %in% allparams$parameter){
      for (i in allparams$subparameter[allparams$parameter == "Forces"]){
        #force <- force + (object@SWdist[,i]*allparams$value[allparams$subparameter == i])
        mF<-cbind(mF,(object@SWdist[,i]*allparams$value[allparams$subparameter == i]))
        colnames(mF)[ncol(mF)]<-i
      }
    }
    
    if ("SinForce" %in% allparams$parameter){
      amp<-allparams$value[allparams$subparameter == "AmpSin"]
      phi<-allparams$value[allparams$subparameter == "PhiSin"]
      per<-allparams$value[allparams$subparameter == "PerSin"]
      #force <- force + amp * sin(2*pi/per * time + phi) 
      mF<-cbind(mF,amp * sin(2*pi/per * time + phi) )
      colnames(mF)[ncol(mF)]<-"SinF"
    }
    
    totabsF<-apply(abs(mF[,2:ncol(mF)]),1,sum)
    
    if (pct==T){
      for (i in 2:ncol(mF)){
        mF[,i]<-abs(mF[,i]/totabsF)
      }
    }

    return(as.data.frame(mF))
    
  }else{
    stop("No force applied to the model")
  }
  
})

setMethod("AllForceApplied",signature="SWDMr_DDHO",function(object,params){
  
  allparams<-SWDMr:::GetAllParams(object,params)
  
  time<-(object@SWdist$Time)  
  
  mF<-matrix(time,ncol=1)
  colnames(mF)<-"Time"
  
  # External forces
  if ("Forces" %in% allparams$parameter){
    for (i in allparams$subparameter[allparams$parameter == "Forces"]){
      #force <- force + (object@SWdist[,i]*allparams$value[allparams$subparameter == i])
      mF<-cbind(mF,(object@SWdist[,i]*allparams$value[allparams$subparameter == i]))
      colnames(mF)[ncol(mF)]<-i
    }
  }
  
  # Zeitgeber force
  if ("SinForce" %in% allparams$parameter){
    amp<-allparams$value[allparams$subparameter == "AmpSin"]
    phi<-allparams$value[allparams$subparameter == "PhiSin"]
    per<-allparams$value[allparams$subparameter == "PerSin"]
    #force <- force + amp * sin(2*pi/per * time + phi) 
    mF<-cbind(mF,amp * sin(2*pi/per * time + phi) )
    colnames(mF)[ncol(mF)]<-"SinF"
  }
  
  # Inner force by string constant
  out <- SWDMr::SWDMrFit(object,params)
  k<-allparams[allparams$subparameter == "omega","value"]^2
  pos<-out$y1[-1]-allparams[allparams$subparameter == "intercept","value"]
  mF<-cbind(mF,pos*k)
  colnames(mF)[ncol(mF)]<-"FstringConst"
  
  # Inner force by damping constant
  if ("loggamma" %in% allparams$subparameter){
    mF<-cbind(mF,out$y2[-1]*exp(allparams[allparams$subparameter == "loggamma","value"]))
    colnames(mF)[ncol(mF)]<-"Fdamp"
  }else if ("dampratio" %in% allparams$subparameter){
    omegav<-allparams[allparams$subparameter == "omega","value"]
    damprv<-allparams[allparams$subparameter == "dampratio","value"]
    gamma<-damprv*2*omegav
    mF<-cbind(mF,out$y2[-1]*gamma)
    colnames(mF)[ncol(mF)]<-"Fdamp"
  }

  
  return(as.data.frame(mF))
  
})
