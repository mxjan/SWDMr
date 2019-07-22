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

model