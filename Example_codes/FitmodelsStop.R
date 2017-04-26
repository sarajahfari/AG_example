
# prepare matrix to save AIC values for individual subject
AIC_direct = round(matrix(,length(yDat),length(models)),digits=2)
colnames(AIC_direct) =c(paste('aic',1:length(models),sep=''))
rownames(AIC_direct)= names(yDat)
	
# prepare matrix to save chi (yuan bentler chi test) for individuals	
Fit_direct = round(matrix(,length(yDat),length(models)),digits=2)
colnames(Fit_direct) = c(paste('chi',1:length(models),sep=''))
rownames(Fit_direct)= names(yDat)
	
# save P-val chi test, if p<0.05 model is no fit for subject compared to saturated model
P_direct = round(matrix(,length(yDat),length(models)),digits=2)
colnames(P_direct) = c(paste('chi-p',1:length(models),sep=''))
rownames(P_direct)= names(yDat)

# make matrix to save log-likelihood per models en per ppn
logL_direct = round(matrix(,length(yDat),length(models)),digits=2)
colnames(logL_direct) = c(paste('logL',1:length(models),sep=''))
rownames(logL_direct)= names(yDat)
	
# save number of observations for each subject to compute BIC
nObs = round(matrix(,length(yDat),length(models)),digits=2) 
colnames(nObs) = c(paste('nobs',1:length(models),sep=''))
rownames(nObs)= names(yDat)
	
# save individual AIC values
AMT=round(matrix(,length(cond),length(models)*2),digits=2)
colnames(AMT)=c(names(models),paste('nfit-',names(models),sep=''))
rownames(AMT)=names(cond)
	
# save individual BIC values
BMT=round(matrix(,length(cond),length(models)*2),digits=2)
colnames(BMT)=c(names(models),paste('nfit-',names(models),sep=''))
rownames(BMT)=names(cond)
	
for (k in 1:length(cond)) # loop over the conditions in cond
	{
	
	  for (i in 1:length(yDat)) # do AG for each subject in yDat
	
	    {
		    Aic =c()
		    Fit =c()
		    P =c()
        nPars = c() # vecotr met nr par voor ieder models
		    nobs=c() # vector with  number of observations for each model
		
		    # var models contains the models
		    # yDat contains the data for all subjects in a list
		    # C.Lab indicates the Region of input for the models
		    # cond contains the index numbers for each Condition, each subject in yDat
		
		    for (j in 1:length(models)) # loop over the different models in var models
		
			  {
				Aic[j]=fitAG(yDat[i], models[[j]],C.Lab, cond[[k]][[i]],detail='both')$aic
            nPars[j]=fitAG(yDat[i], models[[j]],C.Lab, cond[[k]][[i]],detail='both')$npars # bewaar per models nr par
				      Fit[j]=fitAG(yDat[i], models[[j]],C.Lab, cond[[k]][[i]],detail='both')$fit$chi2
				        P[j]=fitAG(yDat[i], models[[j]],C.Lab, cond[[k]][[i]],detail='both')$fit$p
				          nobs[j] = length(cond[[k]][[i]])
				  }# close j loop over models
		
		AIC_direct[i,] = round(Aic,digits=2)
		Fit_direct[i,] = round(Fit,digits=2)
		P_direct[i,] = round(P,digits=2)
		logL_direct[i,] = Aic - 2*nPars # de log-likelihood is dan de aic - de factor van de aic
    nObs[i,]=round(nobs[j],digits=2)
                
		}	# close i loop over subjects


		A=cbind(AIC_direct,P_direct)
		
		AIC.direct <- apply(logL_direct,2,sum)+2*nPars # radom effecs AIC with pooled LLH over subjects
		BIC.totaal <- apply(logL_direct,2,sum)+(nPars*log(apply(nObs,2,sum))) # radom effecs BIC with pooled LLH over subjects
    nfit.rfx = apply(ifelse(A[,(length(models)+1):(2*length(models))]>0.049,1,0),2,sum) # get the number of subjects where the model fits
		
    # save the random effects AIC and BIC but also individual subjest (from chi test) where the model did not fit (key!)
		AMT[k,]=c(AIC.direct,nfit.rfx) 
		BMT[k,]=c(BIC.totaal,nfit.rfx)
	
} # close cond loop

