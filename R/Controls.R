SWDMr_DDHOcheck<-function(object){
  
  errors<-character()
  
  if (length(object@VarExp) == 0){
    msg<-"No VarExp given"
    return(msg)
  }
  
  if (! object@VarExp %in% colnames(object@Gexp)){
    msg<-"The VarExp is not in data.frame"
    errors<-c(errors,msg)
  }
  
  if (length(errors) == 0) TRUE else errors
}

SWDMrcheck<-function(object){
  
  errors<-character()
  
  if (! "Time" %in% colnames(object@SWdist)){
    msg<-"\"Time\" is missing in SWdist data frame"
    errors<-c(errors,msg)
  }
  if (! "Time" %in% colnames(object@Gexp)){
    msg<-"\"Time\" is missing in Gexp data frame"
    errors<-c(errors,msg)
  }
  
  if (length(errors) == 0) TRUE else errors
  
}


SWDMr_ProcScheck<-function(object){
  
  errors<-character()
  
  if (length(object@VarExp) == 0){
    msg<-"No VarExp given"
  }
  
  if (! "Sleep" %in% colnames(object@SWdist)){
    msg<-"No Sleep columns in data"
    errors<-c(errors,msg)
  }
  
  if (! "Wake" %in% colnames(object@SWdist)){
    msg<-"No Wake columns in data"
    errors<-c(errors,msg)
  }
  
  if (! object@VarExp %in% colnames(object@Gexp)){
    msg<-"The VarExp is not in data.frame"
    errors<-c(errors,msg)
  }
  
  if (length(errors) == 0) TRUE else errors
}