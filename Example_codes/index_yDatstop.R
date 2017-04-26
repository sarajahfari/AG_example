# define conditions for AG model fitting davinci

	# load library
	load(paste(dir_data,"/yDat_stopexample.Rdat",sep='')) # f-test, and bg no-ftest combined

	yDat=yDat.stopexample
	names(yDat)=names(yDat.stopexample)
	rm(yDat.stopexample)
#################################################
  
	GoL= list();GoR= list()
	STL= list();STR= list()
	SRL= list();SRR= list()
	
	error= list();miss= list()
	
	# collapsed over hand
	Go= list()
	ST= list()
	SR= list()
######################################################
	
	
	#yDat[[1]]$type
	
	for (c in 1:length(yDat))
	
		{
			GoL[[c]] =grep('GoL', yDat[[c]]$type)
			GoR[[c]] =grep('GoR', yDat[[c]]$type)

			STL[[c]] =grep('StL', yDat[[c]]$type)
			STR[[c]] =grep('StR', yDat[[c]]$type)

			SRL[[c]] =grep('SrL', yDat[[c]]$type)
			SRR[[c]] =grep('SrR', yDat[[c]]$type)
	
			error[[c]] =grep('error', yDat[[c]]$type)
			miss[[c]] =grep('miss', yDat[[c]]$type)
			
	# collapsed over hand
	
			Go[[c]] =c(GoL[[c]],GoR[[c]])
			ST[[c]] =c(STL[[c]],STR[[c]])
			SR[[c]] =c(SRL[[c]],SRR[[c]])

		}
		
	
	# give names of subjects to each list	
	names(GoL)= names(yDat)
	names(GoR)= names(yDat)
	names(STL)= names(yDat)
	names(STR)= names(yDat)
	names(SRL)= names(yDat)
	names(SRR)= names(yDat)
	
	names(error)= names(yDat)
	names(miss)= names(yDat)
	names(Go)= names(yDat)
	names(ST)= names(yDat)
	names(SR)= names(yDat)
	
	
		