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
  if (object@UseDampingRatio == T){
    CoreParams<-c("intercept","omega","dampratio")
  }else{
    CoreParams<-c("intercept","omega","loggamma")
  }
  for (param in CoreParams){
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


setMethod("PenalizeUnstableFit",signature="SWDMr_DDHO", function(object,value = T,PredictedValueInterval,StabilityDayCheck,weight=1){
  
  # Set penalization to True
  object@PenalizeUnstableFit<-T
  
  if (length(PredictedValueInterval) != 2){stop("The PredictedValueInterval provided must be 2 values (e.g. c(0,24))")}
  object@PredictedValueInterval<-PredictedValueInterval
  
  object@StabilityDayCheck<-StabilityDayCheck
  
  object@PenalizeUnstableFitWeight<-weight
  
  return(object)
})

setGeneric("SumForces",  function(object,params,allparams=NULL)
  standardGeneric("SumForces") )
setMethod("SumForces",signature="SWDMr_DDHO", function(object,params,allparams=NULL){
  
  if (is.null(allparams)){
    allparams<-GetAllParams(object,params)
  }
  
  idxForces<-which(allparams$parameter == "Forces")
  if (length(idxForces) > 0){
    Forces<-allparams[idxForces,"value"]
    SWf<-object@SWdist[,allparams[idxForces,"subparameter"],drop=F]
    Force <- (as.matrix(SWf) %*% Forces)[,1]
    # SWf<-as.matrix(mapply(`*`,SWf,Forces))
    # if (length(Forces) > 1){
    #   Force<-rowSums(SWf)
    # }else{
    #   Force<-SWf[,1]
    # }
  }else{
    Force<-rep(0,nrow(object@SWdist))
  }
  
  return(Force)
  
})

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
  
  # Fit
  out<-SWDMr:::SDDHO_SinF_RungeKutta(y=y.init,time = time,force = force,gamma = gamma, k=k, AmpSin = AmpSin, PhiSin = PhiSin, PerSin = PerSin) # Solve ODE, see src/sddho_sinforce_RK4.cpp
  
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

setGeneric("AddUnstableFitPenalization",  function(object,fitted,FittingValue,val,var,weight=1)
  standardGeneric("AddUnstableFitPenalization") )
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
  
  if (FittingValue == "RSS"){
    return(val+sum((pred-obs)^2)*weight)
  }
  if (FittingValue == "NLL"){
    # transform in LL
    val <- -val
    # Add likelihood
    val <- val+sum(dnorm(obs,mean=pred,sd=sqrt(var),log=T))*weight
    # Return NLL
    return(-val)
  }
  if (FittingValue == "LL"){
    val <- val+sum(dnorm(obs,mean=pred,sd=sqrt(var),log=T))*weight
    return(val)
  }
})



setMethod("SWDMrStats",signature="SWDMr_DDHO", function(object,fitted,FittingValue="RSS",detailed=F){
  
  predv<-approxfun(fitted$time,fitted$y1)
  predval<-predv(object@Gexp$Time)
  
  idx<- ! is.na(predval) & ! is.na(object@Gexp[,object@VarExp])
  
  GeneExp<-object@Gexp[idx,object@VarExp]
  predval<-predval[idx]
  
  n <- length(GeneExp)
  
  if (detailed == F){
    if (length(predval) == 0){
      return(1e30)
    }else{
      RSS<-sum((GeneExp-predval)^2)
    }
    
    if (FittingValue == "RSS"){
      return(list(val=RSS,var=RSS/n))
    }else{
      var<-RSS/n
      LL <- sum(dnorm(GeneExp,mean=predval,sd=sqrt(var),log=T))
      if (FittingValue == "NLL"){
        return(list(val=-LL,var=var))
      }else if (FittingValue == "LL") {
        return(list(val=LL,var=var))
      }
    }
  }else{
    
    # N parameters
    k <- nrow(GetFreeFixedParams(object)$FreeParams)
    # residuals
    residualsV<-(GeneExp-predval)
    # Residuals sum of square
    RSS<-sum(residualsV^2)
    # Variance
    var<-RSS/n
    # Negative Log Likelihood
    NLL<- -1* sum(dnorm(GeneExp,mean=predval,sd=sqrt(var),log=T))
    # Bayesian Information Criterion
    BIC <- -2*(-NLL)+(k)*log(n)
    # Akaike information criterion
    AIC <- 2*(k) + n*log(RSS/n)
    
    # FLAT MODEL
    lm_flat<-lm(GeneExp~1)
    RSS_flat<-sum(resid(lm_flat)^2)
    var_flat<-RSS_flat/n
    NLL_flat<- -1* ( ((-n/2)*log(2*pi)) - ((n/2)*log(var_flat)) - (RSS_flat/(2*var_flat)) )
    BIC_flat <- -2*(-NLL_flat)+(1)*log(n)
    
    # Bayes Factor
    BF_DDHOvFlat<-exp((BIC_flat-BIC)/2)
    
    # Kendal Tau
    KendalTau<-cor(GeneExp,predval,method="kendall")
    
    return(list(stats = as.data.frame(list(Variable=object@VarExp,RSS=RSS,NLL=NLL,
                              BIC=BIC,BIC_flat=BIC_flat,BayesFactor=BF_DDHOvFlat,
                              AIC=AIC,n=n,k=k,ErrorVariance=var, KendalTau = KendalTau)), residuals = residualsV, fitted = predval ))
    
    
    # k <- nrow(GetFreeFixedParams(object)$FreeParams)
    # 
    # # R2 NOT VALID FOR NON LINEAR MODELS ! DOES NOT ADD UP TO 100 !!
    # # https://statisticsbyjim.com/regression/r-squared-invalid-nonlinear-regression/
    # # Use Standard Error of regression https://statisticsbyjim.com/regression/standard-error-regression-vs-r-squared/
    # # SEE SOME STATS HERE: http://sia.webpopix.org/nonlinearRegression.html
    # # See Note on the R2 measure of goodness of fit for nonlinear models (TARALD O. KVALSETH, 1983)
    # RSS<-sum((GeneExp-predval)^2)
    # MSS<-sum((predval - mean(predval))^2)
    # R2<-MSS/(MSS+RSS)
    # 
    # rdf<-n-k
    # AdjR2<- 1 - (1-R2) * ((n - 1)/rdf)
    # 
    # resvar <- RSS/rdf
    # Fstat<-(MSS/(k - 1))/resvar
    # pvalF<-1-pf(Fstat,k-1,rdf)
    # 
    # var<-RSS/n
    # NLL<- -1* sum(dnorm(GeneExp,mean=predval,sd=sqrt(var),log=T))
    # # Faster ?
    # #((-n/2)*log(2*pi)) - ((n/2)*log(sigma2)) - (RSS/(2*var))
    # 
    # BIC <- -2*(-NLL)+(k+1)*log(n) # add variance as parameter
    # AIC <- 2*(k+1) + n*log(RSS/n)
    # 
    # # Flat model
    # lm_flat<-lm(GeneExp~1)
    # RSS_flat<-sum(resid(lm_flat)^2)
    # var_flat<-RSS_flat/n
    # NLL_flat<- -1* ( ((-n/2)*log(2*pi)) - ((n/2)*log(var_flat)) - (RSS_flat/(2*var_flat)) )
    # BIC_flat <- -2*(-NLL_flat)+(1+1)*log(n)
    # BF_DDHOvFlat<-exp((BIC_flat-BIC)/2)
    # 
    # 
    # return(as.data.frame(list(Variable=object@VarExp,RSS=RSS,NLL=NLL,
    #                           BIC=BIC,BIC_flat=BIC_flat,BayesFactor=BF_DDHOvFlat,
    #                           AIC=AIC,R2=R2,AdjR2=AdjR2,Fstat=Fstat,pvalF=pvalF,numdf=k-1,rdf=rdf,n=n,k=k,ErrorVariance=var)))
  
  }


})


setMethod("SWDMrGetEvalFun",signature="SWDMr_DDHO", function(object){
  
  # Return an objective function
  objfun<-function(params){
    
    # Step 1. Fit function given params
    out<-SWDMrFit(object,params)
    
    # Step 2. Get fitting value
    stat<-SWDMrStats(object,out,FittingValue = object@FittingValue)
    
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


setGeneric("StatsPerTimePoint", function(object) 
  standardGeneric("StatsPerTimePoint") )
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


setGeneric("GetAllParams", function(object,params)
  standardGeneric("GetAllParams") )
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

setGeneric("PctAbsForceApplied", function(object,params,pct=T) 
  standardGeneric("PctAbsForceApplied") )
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


setGeneric("AllForceApplied", function(object,params) 
  standardGeneric("AllForceApplied") )
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
  mF<-cbind(mF,out$y2[-1]*exp(allparams[allparams$subparameter == "loggamma","value"]))
  colnames(mF)[ncol(mF)]<-"Fdamp"
  
  return(as.data.frame(mF))
  
})
