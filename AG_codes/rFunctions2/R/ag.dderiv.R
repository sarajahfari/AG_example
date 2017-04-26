ag.dderiv <-
## second-order derivatives of ancestral graphs
function(ag.fit,S,n)
{
  require(MASS)
	Sigma <- ag.fit$Shat
	Bhat <- ag.fit$Bhat
	Lhat <- ag.fit$Lhat
	Ohat <- ag.fit$Ohat
        v <- nrow(S)
	if(identical(Ohat,Lhat)) Lhat <- 0
        beta <- Bhat[Bhat!=0 & Bhat!=1]
        Llow.diag <- Lhat[lower.tri(Lhat,diag=TRUE)]
        lambda <- Llow.diag[Llow.diag!=0]
        Olow.diag <- Ohat[lower.tri(Ohat,diag=TRUE)]
        omega <- Olow.diag[Olow.diag!=0]
        pB <- length(beta)
        pL <- length(lambda)
        pO <- length(omega)
	    theta <- c(beta,lambda,omega)
	    p <- length(theta)
        SigmaI <- ginv(Sigma)
        BhatI <- ginv(Bhat)
        if(pL!=0) LhatI <- ginv(Lhat)
        OhatI <- ginv(Ohat)
        dderiv <- array(0,dim=c(p,p))

	for(i in 1:p)
	{
	  ## first-order derivatives of Sigma
	  if(i<=pB)
            {
              DBi <- Bhat
              DBi[theta[i]!=DBi] <- 0
              DBi[DBi!=0] <- 1
              Di <- -BhatI%*%DBi%*%Sigma - Sigma%*%t(DBi)%*%t(BhatI)
            }
          if(i>pB & i<=(pB+pL))
            {
              DLi <- Lhat
              DLi[theta[i]!=DLi] <- 0
              DLi[DLi!=0] <- 1
              Di <- -BhatI%*%DLi%*%t(BhatI)
            }
          if(i>(pB+pL))
            {
              DOi <- Ohat
              DOi[theta[i]!=DOi] <- 0
              DOi[DOi!=0] <- 1
              Di <- BhatI%*%DOi%*%t(BhatI)
            }
          for(j in 1:i)
          {
          if(j<=pB)
            {
              DBj <- Bhat
              DBj[theta[j]!=DBj] <- 0
              DBj[DBj!=0] <- 1
              Dj <- -BhatI%*%DBj%*%Sigma - Sigma%*%t(DBj)%*%t(BhatI)
            }
          if(j>pB & j<=(pB+pL))
            {
              DLj <- Lhat
              DLj[theta[j]!=DLj] <- 0
              DLj[DLj!=0] <- 1
              Dj <- -BhatI%*%DLj%*%t(BhatI)
            }
          if(j>(pB+pL))
            {
              DOj <- Ohat
              DOj[theta[j]!=DOj] <- 0
              DOj[DOj!=0] <- 1
              Dj <- BhatI%*%DOj%*%t(BhatI)
            }
          ## second-rder derivatives of Sigma
          if(i<=pB & j<=pB) # beta,beta
            {
              D2 <- 2*( BhatI%*%Di%*%BhatI%*%Dj%*%Sigma 
                      + BhatI%*%Di%*%Sigma%*%t(Dj)%*%t(BhatI)
                      + Sigma%*%Di%*%t(BhatI)%*%t(Dj)%*%t(BhatI) )
            }
          if(i<=pB & j>pB & j<=(pB+pL)) # beta,lambda
            {
                    D2 <- 2*( BhatI%*%Di%*%BhatI%*%LhatI%*%Dj%*%LhatI%*%t(BhatI)
                      + BhatI%*%LhatI%*%Dj%*%LhatI%*%t(BhatI%*%Di%*%BhatI) )
            }
          if(j<=pB & i>pB & i<=(pB+pL)) # lambda,beta  
            {
                    D2 <- 2*( BhatI%*%Dj%*%BhatI%*%LhatI%*%Di%*%LhatI%*%t(BhatI)
                      + BhatI%*%LhatI%*%Di%*%LhatI%*%t(BhatI)%*%Dj%*%t(BhatI) )
            }
          if(i<=pB & j>(pB+pL)) # beta,omega
            {
                    D2 <- -2*( BhatI%*%Di%*%BhatI%*%Dj%*%t(BhatI)
                      + BhatI%*%Dj%*%t(BhatI%*%Di%*%BhatI) )
            }
          if(j<=pB & i>(pB+pL)) # omega,beta
            {
                    D2 <- -2*( BhatI%*%Dj%*%BhatI%*%Di%*%t(BhatI)
                      + BhatI%*%Di%*%t(BhatI%*%Dj%*%BhatI) )
            }
          if(i>pB &i<=(pB+pL) & j>(pB+pL) || j>pB &j<=(pB+pL) & i>(pB+pL)) # lambda,omega
            {
                    D2 <- array(0,dim=c(v,v))
            }
          if(i>pB & i<=(pB+pL) & j>pB & j<=(pB+pL)) # lambda,lambda
            {
                    D2 <- 2*BhatI%*%LhatI%*%Di%*%LhatI%*%Dj%*%LhatI%*%t(BhatI)
            }
          if(i>(pB+pL) & j>(pB+pL)) # omega,omega
            {   
                    D2 <- array(0,dim=c(v,v))
            }
        dderiv[i,j] <- -0.5*n*sum(diag( SigmaI%*%(2*S - Sigma)%*%SigmaI%*%Di%*%SigmaI%*%Dj )) #- SigmaI%*%(S - Sigma)%*%SigmaI%*%D2 )) 
	   }
    }
    dderiv <- dderiv + t(dderiv) - diag(diag(dderiv))
    return(dderiv)
}
