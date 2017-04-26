ag.deriv <- 
## first-order derivatives of ancestral graphs
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
  deriv <- numeric(0)
        
    for(i in 1:p)
	{
          if(i<=pB)
            {
              DB <- Bhat
              DB[theta[i]!=DB] <- 0
              DB[DB!=0] <- 1
              D <- -BhatI%*%DB%*%Sigma - Sigma%*%t(DB)%*%t(BhatI)
            }
          if(i>pB & i<=(pB+pL))
            {
              DL <- Lhat
              DL[theta[i]!=DL] <- 0
              DL[DL!=0] <- 1
              D <- -BhatI%*%LhatI%*%DL%*%LhatI%*%t(BhatI)
            }
          if(i>(pB+pL))
            {
              DO <- Ohat
              DO[theta[i]!=DO] <- 0
              DO[DO!=0] <- 1
              D <- BhatI%*%OhatI%*%DO%*%OhatI%*%t(BhatI)
            }
          deriv <- c(deriv,0.5*n*sum(diag( SigmaI%*%(S - Sigma)%*%SigmaI%*%D )))
	}
        return(deriv)
}
