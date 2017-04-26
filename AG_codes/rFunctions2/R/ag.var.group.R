ag.var.group <- 
## variance at group level determined by Xg (Liang & Zeger, 1986)
function(fit,S,Y)
{
  N <- length(fit)
  n <- length(Y[[1]][,1])
  q <- length(ag.deriv(fit[[1]],S[[1]],n))
  sandwich.g <- array(NA,dim=c(q,q,N))
  for(i in 1:N)
  {
    n <- length(Y[[i]][,1])
    sandwich.g[,,i] <- ag.sandwich(fit[[i]],S[[i]],Y[[i]],n)
  }
  return(sandwich.g)
}
