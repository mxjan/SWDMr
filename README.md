Package
=======

load packages

``` r
library(SWDMr)
library(optimx)
library(ggplot2)
library(doSNOW)
```

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Loading required package: snow

Prepare data
============

Sleep-wake data per hour

``` r
# Data
# ****** Prepare Sleep-Wake data ****** 
files <- list.files(path="../data/Sleep_States/",pattern = paste("^","BL6",sep=""),full.names = T)
SWdf<-Read_SW(files,concattimesec = 360) # 300 = 5min, 180 = 3 min, 360 = 6min
SWdf<-SWdf_AddLD(SWdf) # Add Light and Dark
SWdf<-SWdf_DayMerging(SWdf,Daysformat=list(c(1,2),c(1,2),3,4,c(1,2)),concattimesec=360)
SWdf<-SWdf_AddSD(SWdf,c(48,54))
head(SWdf)
```

    ##         NREM         REM       Wake      Sleep LenW LenS Day Time Light
    ## 1 0.03273148 0.001296296 0.06597222 0.03402778 59.5 30.5   1  0.1     1
    ## 2 0.03648148 0.002824074 0.06069444 0.03930556 55.0 35.0   1  0.2     1
    ## 3 0.04069444 0.004351852 0.05495370 0.04504630 49.5 40.5   1  0.3     1
    ## 4 0.04986111 0.002824074 0.04731481 0.05268519 42.5 47.5   1  0.4     1
    ## 5 0.06240741 0.003611111 0.03398148 0.06601852 30.5 59.5   1  0.5     1
    ## 6 0.06763889 0.009907407 0.02245370 0.07754630 20.5 69.5   1  0.6     1
    ##   Dark SD
    ## 1    0  0
    ## 2    0  0
    ## 3    0  0
    ## 4    0  0
    ## 5    0  0
    ## 6    0  0

Explained variable is gene expression

``` r
# ***** Prepare Gene expression data ******
# Gene expression
load("../data/CHor2018_TimeCourse_NormalizedData_Genes.RData")
# Meta data
load("../data/CHor2018_TimeCourse_MetaData.RData")
# Remove data sampled 7days after SD
rna_expr<-rna_expr[,rna_meta$time<100]
rna_meta<-rna_meta[rna_meta$time<100,]
# Row = samples, column = gene
rna_expr<-t(rna_expr)
#rna_expr<-2^rna_expr
rna_expr<-as.data.frame(rna_expr)
# Add time factor
rna_expr$Time<-rna_meta$time+24
head(rna_expr[,c("Arntl","Dbp","Time")])
```

    ##          Arntl      Dbp Time
    ## ZT0A  5.944392 4.898669   24
    ## ZT0B  5.977196 4.750183   24
    ## ZT0H  5.849927 5.281751   24
    ## J70N1 5.860580 4.980111   24
    ## T3N2  5.835763 5.282515   27
    ## T3N3  5.900950 5.159076   27

``` r
#rna_expr<-rna_expr[-which(rna_expr$Time %in% c(51,54)),]
```

Run fitting
===========

Build swdmr object

``` r
swdmr <- SWDMr(SWdist=SWdf, Gexp=rna_expr)
swdmr
```

    ## This is a SWDMr object
    ## This object contains:  17185 Genes Over 56 Time points
    ## Your force data frame contain 1200 values with the following possible forces:NREM;REM;Wake;Sleep;LenW;LenS;Day;Time;Light;Dark;SD

Initiate a Driven Damped Harmonic Oscillator \[DDHO\] model for a gene

``` r
Gene<-"Nr1d1"
model<-initDDHOmodel(swdmr,VarExp = Gene)
```

Set some parameter of our model

``` r
# Mean expression in baseline between highest and lowest value
MeanPerTime<-aggregate(rna_expr[rna_meta$SD_NSD == "NSD",Gene],list(rna_meta$time[rna_meta$SD_NSD == "NSD"]),mean)
MeanGeneExprInBaseline<-(max(MeanPerTime$x)+min(MeanPerTime$x))/2

# Fix the intercepts
model<-FixIntercept(model,MeanGeneExprInBaseline)
# Add sleep-wake force
model<-AddForce(model,"Wake")
model<-AddForce(model,"Sleep")


# Start is set at intercept with speed of 0
model<-SetYinitMode(model,mode = "Intercept_0",values = c(0,48))
# We replicate baseline for 20 day
model<-ReplicateDrivingForce(model,c(0.1,24.0),20)
# A sin-wave force is applied with a period of 24h
model<-AddSinF(model,FixPer = 24)
# Compute the fit using RSS
model<-SetFittingValue(model,value = "RSS")
# Penalize the fitting for unstable value for 10 replicated days
model<-PenalizeUnstableFit(model,value = T,PredictedValueInterval = c(0,48), StabilityDayCheck = 10)
```

summary of the model

``` r
model
```

    ## ~~~~~~~~ This is a S4 SWDMr_DDHO object ~~~~~~~~ 
    ## 
    ## Display the current setting for your fitting
    ## 
    ## ~~~~~~~~~~~ Current parameter setting ~~~~~~~~~~ 
    ## 
    ## * [free parameters] omega (Core parameter) 
    ## * [free parameters] loggamma (Core parameter) 
    ## * [free parameters] Wake (Forces) 
    ## * [free parameters] Sleep (Forces) 
    ## * [free parameters] AmpSin (SinForce) 
    ## * [free parameters] PhiSin (SinForce) 
    ## 
    ## 
    ## * [fixed parameters] intercept (Core parameter) : 6.84537089043882
    ## * [fixed parameters] PerSin (SinForce) : 24

Fit data with optimx

``` r
# Get objective function
objfun<-SWDMrGetEvalFun(model)
# Limits of the model
params<-c(Wake=0,Sleep=0,loggamma=log(1e-1),omega=2*pi/23.8,AmpSin=1e-2,PhiSin=2.5)
lower<-c(Wake=-Inf,Sleep=-Inf,loggamma=-Inf,omega=2*pi/30,AmpSin=0,PhiSin=0)
upper<-c(Wake=Inf,Sleep=Inf,loggamma=Inf,omega=2*pi/18,AmpSin=Inf,PhiSin=2*pi)

# Fit
optimxres<-optimx(params,objfun,method=c("nlminb"),lower = lower,upper=upper)
optimxres
```

    ##              Wake     Sleep loggamma     omega      AmpSin   PhiSin
    ## nlminb -0.1338645 0.2065329 -3.63269 0.2437751 0.008298693 2.665784
    ##           value fevals gevals niter convcode kkt1 kkt2 xtime
    ## nlminb 1.545931     63    276    41        0 TRUE TRUE  0.61

Get fit

``` r
out<-SWDMrFit(model,params = optimxres[1,])
par(mfrow=c(2,1))
plot(out$time,out$y1,type="l",ylab="Position",xlab="Time")
plot(out$time,out$y2,type="l",ylab="Speed",xlab="Time")
```

![](README_files/figure-markdown_github/unnamed-chunk-9-1.png)

Get some statistics for the fit

``` r
SWDMrStats(model,out,detailed = T)
```

    ##   Variable     RSS       NLL      BIC BIC_flat  BayesFactor       AIC
    ## 1    Nr1d1 1.54593 -21.05178 -13.9261 46.52304 7.475457e-14 -187.0247
    ##          R2     AdjR2    Fstat       pvalF numdf rdf  n k ErrorVariance
    ## 1 0.7628008 0.7390809 32.15866 1.64313e-14     5  50 56 6    0.02760588

Plot fit and gene expression using ggplot2

``` r
SWDMr:::StandardFittingPlot(model,optimxres[1,])
```

![](README_files/figure-markdown_github/unnamed-chunk-11-1.png)

Compute confidence interval using empirical bootstrap

``` r
# See http://sia.webpopix.org/nonlinearRegression.html
EmpBoot<-function(model, optimxres, params, upper, lower, nboot = 100, NCORES = 7){
  
  out<-SWDMrFit(model,params = optimxres[1,])
  
  F <- matrix(nrow=nboot,ncol=length(out$y1))
  Y <- matrix(nrow=nboot,ncol=length(out$y1))
  
  pred<-SWDMrFit(model,params = optimxres[1,])
  predv<-pred$y1
  names(predv)<-pred$time
  predval<-predv[as.character(model@Gexp[,"Time"])]
  
  stats<-SWDMrStats(model,out,detailed = T)
  model2<-model
  
  paramsboot<-as.numeric(optimxres[1,names(params)])
  names(paramsboot)<-names(params)
  
  cl <- makeCluster(NCORES)
  #clusterExport(cl,c("lower","upper","model2","stats"))
  registerDoSNOW(cl)
  
  B <- foreach(i = 1:nboot,.packages=c("SWDMr","optimx"),.combine="rbind") %dopar% {
    n<-length(model2@Gexp[,1])
    model2@Gexp[,1]<- predval + sqrt(stats[["ErrorVariance"]])*rnorm(n)
    objfunboot<-SWDMrGetEvalFun(model2)
    optimxresboot<-optimx(paramsboot,objfunboot,method=c("nlminb"),lower = lower,upper=upper)
    
    return(as.numeric(optimxresboot[1,names(params)]))
  }
  stopCluster(cl)
  
  for (l in 1:nrow(B)){
    bootpar<-B[l,]
    names(bootpar)<-names(params)
    outboot<-SWDMrFit(model,params = bootpar)
    F[l,] <- outboot$y1
    Y[l,] <- outboot$y1 + rnorm(1,0,sqrt(stats[["ErrorVariance"]]))
  }
  
  level <- 0.95
  alpha <- 1 - level
  df.mc<-data.frame(Time=out$time)
  df.mc[c("lwr.conf","upr.conf")] <- t(apply(F,MARGIN=2,function(x) quantile(x,c(alpha/2,1-alpha/2))))
  df.mc[c("lwr.pred","upr.pred")] <- t(apply(Y,MARGIN=2,function(x) quantile(x,c(alpha/2,1-alpha/2))))
  df.mc<-df.mc[df.mc$Time>0,]
  V.mc <- cov(B)
  
  se.mc <- sqrt(diag(V.mc))
  
  t.stat <- paramsboot/se.mc
  p.value <- 2*(1 - pt(abs(t.stat),51))
  
  Est.pvals<-cbind(Estimate=paramsboot,SE=se.mc,Tstat=t.stat,p.value)
  
  
  b <- c(apply(B,MARGIN=2,function(x) quantile(x,alpha/2)),
         apply(B,MARGIN=2,function(x) quantile(x,1-alpha/2)))
  ci.mc=matrix(b,ncol=2)
  row.names(ci.mc) <- names(paramsboot)
  colnames(ci.mc) <- c(paste0((1-level)/2*100,"%"),paste0((1+level)/2*100,"%"))

  return(list(B=B,Y=Y,F=F,CI=ci.mc,Est.pvals=Est.pvals,df.mc=df.mc))
}

EmpBootres<-EmpBoot(model,optimxres,params,upper,lower,nboot=100)
EmpBootres$CI
```

    ##                 2.5%       97.5%
    ## Wake     -0.24219532 -0.09455252
    ## Sleep     0.14951609  0.32472725
    ## loggamma -4.22821840 -2.74486673
    ## omega     0.23503532  0.25340632
    ## AmpSin    0.00630117  0.01035181
    ## PhiSin    2.32597939  2.97287549

``` r
EmpBootres$Est.pvals
```

    ##              Estimate          SE     Tstat      p.value
    ## Wake     -0.133864464 0.038961462 -3.435817 1.182974e-03
    ## Sleep     0.206532924 0.046494522  4.442092 4.807501e-05
    ## loggamma -3.632690235 0.423803002 -8.571648 1.877543e-11
    ## omega     0.243775139 0.004983053 48.920844 0.000000e+00
    ## AmpSin    0.008298693 0.001135900  7.305833 1.790384e-09
    ## PhiSin    2.665783946 0.167902374 15.876988 0.000000e+00

``` r
p<-SWDMr:::StandardFittingPlot(model,optimxres[1,])
p <- p + annotate("ribbon",x=EmpBootres$df.mc$Time, ymin=EmpBootres$df.mc$lwr.pred, ymax=EmpBootres$df.mc$upr.pred, alpha=0.2, fill="blue")
p <- p + annotate("ribbon",x=EmpBootres$df.mc$Time, ymin=EmpBootres$df.mc$lwr.conf, ymax=EmpBootres$df.mc$upr.conf, alpha=0.4, fill="#339900")
p
```

![](README_files/figure-markdown_github/unnamed-chunk-13-1.png)

``` r
#Mean Force applied

Fmod<-SWDMr:::AllForceApplied(model,optimxres[1,])
pctFmod<-Fmod
sumF<-apply(abs(pctFmod[,seq(2,5)]),1,sum)
pctFmod[,seq(2,5)]<-abs(pctFmod[,seq(2,5)])/sumF
head(pctFmod)
```

    ##     Time      Wake     Sleep      SinF FstringConst        Fdamp
    ## 1 -479.9 0.4536661 0.3610213 0.1852851 2.750075e-05 4.762287e-06
    ## 2 -479.8 0.4133653 0.4130123 0.1734895 1.329944e-04 1.373346e-05
    ## 3 -479.7 0.3700768 0.4680347 0.1615250 3.635297e-04 2.730707e-05
    ## 4 -479.6 0.3129342 0.5376112 0.1486841 7.705054e-04 4.716335e-05
    ## 5 -479.5 0.2164140 0.6486834 0.1334894 1.413160e-03 7.838720e-05
    ## 6 -479.4 0.1386800 0.7389431 0.1199814 2.395454e-03 1.193020e-04

``` r
apply(pctFmod[pctFmod$Time>0,],2,mean)
```

    ##         Time         Wake        Sleep         SinF FstringConst 
    ## 6.005000e+01 2.336547e-01 2.591884e-01 1.444320e-01 3.627249e-01 
    ##        Fdamp 
    ## 6.060075e-05

``` r
sessionInfo()
```

    ## R version 3.5.1 (2018-07-02)
    ## Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## Running under: Windows 10 x64 (build 17763)
    ## 
    ## Matrix products: default
    ## 
    ## locale:
    ## [1] LC_COLLATE=French_Switzerland.1252  LC_CTYPE=French_Switzerland.1252   
    ## [3] LC_MONETARY=French_Switzerland.1252 LC_NUMERIC=C                       
    ## [5] LC_TIME=French_Switzerland.1252    
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ## [1] doSNOW_1.0.16    snow_0.4-3       iterators_1.0.10 foreach_1.4.4   
    ## [5] ggplot2_3.1.0    optimx_2018-7.10 SWDMr_1.0       
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] Rcpp_1.0.0        pillar_1.3.1      compiler_3.5.1   
    ##  [4] plyr_1.8.4        tools_3.5.1       digest_0.6.18    
    ##  [7] evaluate_0.13     tibble_2.0.1      gtable_0.2.0     
    ## [10] pkgconfig_2.0.2   rlang_0.4.0       parallel_3.5.1   
    ## [13] yaml_2.2.0        xfun_0.6          withr_2.1.2      
    ## [16] stringr_1.4.0     dplyr_0.8.0.1     knitr_1.22       
    ## [19] grid_3.5.1        tidyselect_0.2.5  glue_1.3.1       
    ## [22] R6_2.4.0          rmarkdown_1.12    purrr_0.3.1      
    ## [25] magrittr_1.5      scales_1.0.0      codetools_0.2-15 
    ## [28] htmltools_0.3.6   assertthat_0.2.1  colorspace_1.4-0 
    ## [31] numDeriv_2016.8-1 labeling_0.3      stringi_1.4.3    
    ## [34] lazyeval_0.2.1    munsell_0.5.0     crayon_1.3.4
