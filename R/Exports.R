#' Initiate DDHO model
#'
#' Initiate a driven damped harmonic oscillator model
#'
#' @param object A SWDMr object 
#' @param VarExp The name of the variable to explain by the model
#' 
#' @return SWDMr model object
#'
#' @examples
#' model <- initDDHOmodel(swdmr,VarExp = "Arntl")
#'
#'@export
initDDHOmodel<-function(object,VarExp){
  new(Class = "SWDMr_DDHO", SWdist = object@SWdist, Gexp = object@Gexp[c(VarExp,"Time")],VarExp = VarExp)
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