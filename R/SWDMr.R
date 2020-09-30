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

setMethod("ReplicateDrivingForce",signature="SWDMr", function(object,interval,Nrep){
    
    # Some Controls
    if (length(interval) != 2){stop("The interval provided must be 2 values (e.g. c(0,24))")}
    
    # Control that the interval is modulo 24h
    t1<-object@SWdist[1,"Time"]
    t2<-object@SWdist[2,"Time"]
    timestep<-t2-t1
    idxInterval<- object@SWdist$Time >= interval[1] & object@SWdist$Time <= interval[2]
    nvals <- nrow(object@SWdist[idxInterval, ])
    if (round(nvals*timestep,5) %% 24 != 0){stop("The interval must be a multiple of 24h, here it is: ",nvals*timestep, "(modulo 24 is:",nvals*timestep %% 24,")")}
    
    disttmp<-do.call("rbind", replicate(Nrep, object@SWdist[idxInterval, ] , simplify = FALSE))
    disttmp$Time<- t1-rev(cumsum(rep(timestep,nrow(disttmp))))
    object@SWdist<-rbind(disttmp,object@SWdist)
    
    return(object)
})


setGeneric("SDSimulation",  function(object,StartTime,Duration)
    standardGeneric("SDSimulation") )
setMethod("SDSimulation",signature="SWDMr",function(object,StartTime,Duration){
    
    timestep<-object@SWdist$Time[2]-object@SWdist$Time[1]
    idxSimu<-which(object@SWdist$Time >= StartTime & object@SWdist$Time <= StartTime+Duration-timestep)
    
    # Keep only time points 
    newSWdf<-object@SWdist[seq(1,min(idxSimu)-1),]
    SndSDSWdf<-object@SWdist[idxSimu,]
    
    SndSDSWdf$NREM<-rep(0,nrow(SndSDSWdf))
    SndSDSWdf$REM<-rep(0,nrow(SndSDSWdf))
    SndSDSWdf$Sleep<-rep(0,nrow(SndSDSWdf))
    SndSDSWdf$Wake<-rep(max(SWdf$Wake),nrow(SndSDSWdf))
    SndSDSWdf$SD<-rep(1,nrow(SndSDSWdf))
    
    object@SWdist<-rbind(newSWdf,SndSDSWdf)
    
    return(object)
    
})