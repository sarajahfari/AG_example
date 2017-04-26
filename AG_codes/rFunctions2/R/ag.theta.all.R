ag.theta.all <- 
## test on ancestral graph parameters
function(fit)
{
  N <- length(fit)
  q <- length(ag.theta(fit[[1]]))
  theta <- array(NA,dim=c(N,q))
  for(i in 1:N)
  {
    theta[i,] <- ag.theta(fit[[i]])
  }
  return(theta)
}
