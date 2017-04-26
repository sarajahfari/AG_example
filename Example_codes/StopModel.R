
	# classic hyper-indirect model (hindi)
	ag0=makeMG(dg=DAG(
					CaudateR40exc~PreSMARsmall+IFGR,
					maxSTNR25exc ~PreSMARsmall+IFGR,
					maxGPeR30exc~CaudateR40exc,
					maxGPiR30exc~maxGPeR30exc+ maxSTNR25exc,
					ThalamusR40exc~maxGPiR30exc			
					),ug=UG(~PreSMARsmall*IFGR))

	# hyperdirect (hyp)
	ag1=makeMG(dg=DAG(
					#CaudateR40exc~PreSMARsmall+IFGR,
					maxSTNR25exc ~PreSMARsmall+IFGR,
					#maxGPeR30exc~CaudateR40exc,
					maxGPiR30exc~maxSTNR25exc,
					ThalamusR40exc~maxGPiR30exc			
					),ug=UG(~PreSMARsmall*IFGR+CaudateR40exc+maxGPeR30exc))

	# indirect (indir)
	ag2=makeMG(dg=DAG(
					CaudateR40exc~PreSMARsmall+IFGR,
					#maxSTNR25exc ~PreSMARsmall+IFGR,
					maxGPeR30exc~CaudateR40exc,
					maxGPiR30exc~maxGPeR30exc,
					ThalamusR40exc~maxGPiR30exc			
					),ug=UG(~PreSMARsmall*IFGR+maxSTNR25exc))

	# direct (dir)
	ag3=makeMG(dg=DAG(
					CaudateR40exc~PreSMARsmall+IFGR,
					#maxSTNR25exc ~PreSMARsmall+IFGR,
					#maxGPeR30exc~CaudateR40exc,
					maxGPiR30exc~CaudateR40exc,
					ThalamusR40exc~maxGPiR30exc			
					),ug=UG(~PreSMARsmall*IFGR+maxSTNR25exc+maxGPeR30exc))


	models=list(ag0,ag1,ag2,ag3)
   
  
	# 
	names(models)=c('hindi','hyp','indir','dir')
	