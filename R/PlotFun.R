# Plot Functions
StandardFittingPlot<-function(object,params){
    
    FittedData<-SWDMrFit(object = object,params = params)
    
    Gene<-object@VarExp
    
    dfd<-as.data.frame(cbind(object@Gexp$Time,object@Gexp[,object@VarExp]))
    colnames(dfd)<-c("Time","Gene")
    
    p_Gene <- ggplot(dfd,aes(x=Time,y=Gene))
    # LD rect
    for (i in seq(12,108,by=24)){
        p_Gene <- p_Gene + annotate("rect", xmin = i, xmax = i+12, ymin = -Inf, ymax = Inf, fill="black",alpha = 0.2)
    }
    
    
    # Boxplot + Data points
    #p_Gene <- p_Gene + geom_point(color = "forestgreen")
    #p_Gene <- p_Gene + geom_boxplot(fill="forestgreen",aes(group = cut_width(Time, 2))) #+ geom_jitter(shape=16, position=position_jitter(0.2))
    # Mean point and sd per time points
    meansd<-StatsPerTimePoint(object)
    
    p_Gene <- p_Gene + ylab("log2 CPM")+xlab("Time") + ggtitle(Gene)
    p_Gene <- p_Gene + scale_x_continuous(breaks=seq(0,120,12))
    p_Gene <- p_Gene + theme(panel.background = element_blank(),panel.border = element_rect(fill=NA)) 
    
    # SD rect
    p_Gene <- p_Gene + annotate("line",x=FittedData$time[FittedData$time>=0],y=FittedData$y1[FittedData$time>=0],color="red",size=2,alpha=.8)
    p_Gene <- p_Gene + annotate(x=unique(object@Gexp$Time),"errorbar",ymin=meansd[,"mean"]-meansd[,"sd"],ymax=meansd[,"mean"]+meansd[,"sd"],colour="black", width=1,size=1)
    rangev<-ggplot_build(p_Gene)$layout$panel_params[[1]]$y.range
    p_Gene <- p_Gene + annotate("rect", xmin = 48, xmax = 48+6, ymin = -Inf, ymax = rangev[[1]]+(rangev[[2]]-rangev[[1]])*0.1, fill="red",alpha = 0.3)
    p_Gene <- p_Gene + annotate("point",x=unique(object@Gexp$Time),y=meansd[,"mean"],size=3,col="black",shape=21,fill="forestgreen")
    
    return(p_Gene)
}