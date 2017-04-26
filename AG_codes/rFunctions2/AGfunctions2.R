## This file needs to be sourced in R
## functions to fit ancestral graphs from fmri data
## tailored to Sara data
## first you need to decide which conditions to test en which to aggregate

# check whether package "ggm" is available. If not, install and obtain
if(!("ggm" %in% .packages(all=TRUE))) install.packages(ggm)
require(ggm)

# make object "brain"
as.brain <- function(yData,rois=NA,cond=NA,...)
{
  if(!is.na(rois[1])) yData <- yData[,rois]
  if(!is.na(cond[1])) yData <-yData[cond,]
  attr(yData,"n") <- length(cond)
  class(yData) <- c("brain","data.frame")
  return(yData)
}


aic.dev <- function(fit,S,n)
{
  p <- length(fit$Shat[1,])
  q <- p*(p-1)/2 - fit$df
  npars <- p+q
  fitS <- -(n/2)*log(2*pi) -(n/2)*log(det(S)) -(n*p)/2
  aic <- fit$dev + fitS + 2*npars
  return(list(aic=aic,npars=npars))
}





# function to fit all subjects using an AG model
# result is aic value for each subject on the model
# or if <detail="LR"> then a likelihood ratio test is performed
# or if <detail="both"> then both results are given
fitAG <- function(yData,amat,rois=NA,cond=NA,detail="AIC",...)
{
  N <- length(yData)
  n <- length(cond)
  p <- length(rois)
  results <- numeric(0)
  for(i in 1:N)
  {
    Data <- as.brain(yData[[i]],rois,cond)
    S <- makeCov(Data)
    fitS <- fitAncestralGraph(amat,S,attr(S,"n"))
    if(detail=="AIC") results <- c(results,aic=aic.dev(fitS,S,n)$aic)
    if(detail=="LR")  results <- c(results, list(chi2=T2(fitS$dev,n),p=pchisq(T2(fitS$dev,n),fitS$df,lower.tail=FALSE)))
    if(detail=="both") results <- c(results,list(fit=list(chi2=T2(fitS$dev,n),p=pchisq(T2(fitS$dev,n),fitS$df,lower.tail=FALSE))),aic=aic.dev(fitS,S,n)$aic,npars=aic.dev(fitS,S,n)$npars)
  }
  return(results)
}

makeCov <- function(yData)
{
    S <- cov(yData)
    attr(S,"n") <- attr(yData,"n")
    return(S)
}


T2 <- function(fit,n)
{
    fit/(1+(fit/n))
}


