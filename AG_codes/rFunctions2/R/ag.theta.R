ag.theta <- 
## test on ancestral graph parameters
function(ag.fit)
{
  require(MASS)
  Sigma <- ag.fit$Shat
  Bhat <- ag.fit$Bhat
  Lhat <- ag.fit$Lhat
  Ohat <- ag.fit$Ohat
  #v <- nrow(S)
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
  
  return(theta)
}
