# r script to see how stopping is affected by visual information

rm(list=ls())
# set datadir Base
dir_base='/Users/sarajahfari/Github/AG_stopexample'

#dir AG_codes, example_codes,
dir_agcodes=paste(dir_base,'/AG_codes',sep='')
dir_excodes=paste(dir_base,'/Example_codes',sep='')
dir_data=paste(dir_base,'/Example_data',sep='')

# source checked packages
source("http://bioconductor.org/biocLite.R")
biocLite("graph")
source("http://bioconductor.org/biocLite.R")
biocLite("RBGL")
install.packages(paste(dir_agcodes,'/Rpackage_checked/ggm-1.995-3.tar.gz',sep=''), repos = NULL, type="source")

# source AG_codes
file.sources = list.files(paste(dir_agcodes,'/rFunctions2',sep=''),pattern="*.R",full.names = TRUE,recursive=T)
sapply(file.sources,source,.GlobalEnv)

# Read yDat list with single subject beta's for each trial/roi, and make index list for the conditions of interest
source(paste(dir_excodes,'/index_yDatstop.R',sep=''), chdir = F)

#########################################################################

# these are the conditions that are evalauted seperatly in FitmodelsStop.R	
cond=list(ST,SR,Go,GoL,GoR)
names(cond)=c('ST', 'SR', 'Go','GoL','GoR')

# if responses are made with the left and right hand it might be though to fit the model with the trials collapsed
# but alos note that this wil reduce your number of trials by half
#cond = list(STL, STR, SRL, SRR,GoL,GoR)
#names(cond)=c('STL', 'STR', 'SRL', 'SRR','GoL','GoR')

# select the graph nodes for connectivity network
M0roi=c("CaudateR40exc" ,"PreSMARsmall","IFGR",
		"maxSTNR25exc","maxGPeR30exc","maxGPiR30exc",
		"ThalamusR40exc")

C.Lab=M0roi

########################################################################

# source the models to use, these are the defined PFC-BG models to evaluate for connectivity/fits 
source(paste(dir_excodes,'/StopModel.R',sep=''), chdir = F) 
# Now do the connectivity evaluation for each subject, each model/condition, to compute aic/bic and n-fits
source(paste(dir_excodes,'/FitmodelsStop.R',sep=''), chdir = F)

# print results based on random effects aic, + number of fits based on yuan and bentler chi (n= number of subject where model fits)
AMT
# print results based on random effects bic, + number of fits based on yuan and bentler chi (n= number of subject where model fits)
BMT

# the hyperdirect-Indirect (hindi) model fits all the 43 ppn for succesful stop (ST), and failed stop (SR) trials and  has the lowes AIC and BIC


#------------------------------------------------------------------------------#
#           Compute (indvidual) connectivity strengths winning model           #
#------------------------------------------------------------------------------#

# ag0 was the hindi model in StopModel.R
hindi=ag0

#conditions to get beta for [these are the conditions where the model fitted all ppn]
cond2=list(ST,SR)
names(cond2)=c('ST', 'SR')

# Make succesfull stop and failed stop lists using yDat
ST=list()
SR=list()

for (i in 1:length(yDat)){
  ST[[i]] =yDat[[i]][cond2[[1]][[i]],C.Lab]
  SR[[i]] =yDat[[i]][cond2[[2]][[i]],C.Lab]
}	


# compute covarience matrix			
covList.ST =lapply(ST,cov)
covList.SR=lapply(SR,cov)

fitST = list()
fitSR = list()

for (i in 1:length(yDat)){
  fitST[[i]] <- fitAncestralGraph(hindi, covList.ST[[i]],dim(ST[[i]])[1])
  fitSR[[i]] <- fitAncestralGraph(hindi, covList.SR[[i]],dim(SR[[i]])[1])
  }

# compute the variance of beta per subject
var.g.ST <-ag.var.group(fitST, covList.ST, ST)
var.g.SR <-ag.var.group(fitSR, covList.SR, SR)

# save matrix for var beta and beta itself (zijn nog leeg)
varbeta.ST = matrix(,dim(var.g.ST)[3],dim(var.g.ST)[1])
varbeta.SR = matrix(,dim(var.g.SR)[3],dim(var.g.SR)[1])
beta.ST = matrix(,length(yDat), length(ag.theta(fitST[[1]])))
beta.SR = matrix(,length(yDat),length(ag.theta(fitSR[[1]])))

for (i in 1:length(yDat)){	
  varbeta.ST[i,] = diag(var.g.ST[,,i])
  varbeta.SR[i,] = diag(var.g.SR[,,i])
  beta.ST[i,] = ag.theta(fitST[[i]])
  beta.SR[i,] = ag.theta(fitSR[[i]])
  }

# these are the standarized connections (to use for group comparisons)
sbeta.ST= round(beta.ST/varbeta.ST,dig=8)
sbeta.SR= round(beta.SR/varbeta.SR,dig=8)

# note here beta and sbeta contain the estimated for all the defined connections in the model

#------------------------------------------------------------------------------#
#           give names to the rows and columns of the beta matrix              #                               
#------------------------------------------------------------------------------#			

# select the defined directed and undirected connections
beta.ST=beta.ST[,c(1:8,10)]
beta.SR=beta.SR[,c(1:8,10)]
sbeta.ST=sbeta.ST[,c(1:8,10)]
sbeta.SR=sbeta.SR[,c(1:8,10)]


colnames(beta.ST) =paste('ST-',c('caud->gpe','presma->caud','presma->stn','ifg->caud',
                                 'ifg->stn','stn->gpi','gpe->gpi','gpi->thalamus',
                                 'presma-ifg'),sep='')

colnames(beta.SR) =paste('SR-',c('caud->gpe','presma->caud','presma->stn','ifg->caud',
                                 'ifg->stn','stn->gpi','gpe->gpi','gpi->thalamus',
                                 'presma-ifg'),sep='')

colnames(sbeta.ST) =colnames(beta.ST)
colnames(sbeta.SR) =colnames(beta.SR)

Beta_output=list(Connectivity_ST=beta.ST, Connectivity_SR=beta.SR,standarized_con_ST=sbeta.ST,standarized_con_SR=sbeta.SR)

for (rn in 1:length(Beta_output)) {rownames(Beta_output[[rn]])=names(yDat)}


