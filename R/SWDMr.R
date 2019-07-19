setClass(
  # Set the name for the class
  Class = "SWDMr",
  
  # Define the slots
  representation = representation(
    SWdist = "data.frame",
    Gexp = "data.frame",
    verbose = "numeric"
  ),
  
  prototype=list(
    verbose = 1
  ),
  
  validity=SWDMrcheck
  
)


setMethod("summary", "SWDMr", function(object) {
  
  cat("This is a SWDMr object\n")
  
  nG<-ncol(object@Gexp)-1
  cat("This object contains: ",nG, "Genes",sep = " ")
  nTp<-nrow(object@Gexp)
  cat(" Over",nTp, "Time points",sep=" ")
  
  SWdf_nrow<-nrow(object@SWdist)
  SWdf_colnames<-colnames(object@SWdist)
  
  cat(paste("\nYour force data frame contain",SWdf_nrow, "values with the following possible forces:"))
  cat(paste(SWdf_colnames,collapse=";"))
  
})

setMethod("show", "SWDMr", function(object) {
  summary(object)
})