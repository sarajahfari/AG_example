ag.sandwich <- 
## variance at individual level
function(fit,S,Y,n)
{
  p <- length(ag.deriv(fit,S,n))
  H <- ag.dderiv(fit,S,n)
  uu <- array(0,dim=c(p,p))
  for(j in 1:n)
  {    	     
    yy <- t(as.matrix(Y[j,]))%*%as.matrix(Y[j,])
    score <- ag.deriv(fit,yy,n)
    uu <- uu + score%*%t(score)/(n^2)
  }
  HI <- solve(H)
  V <- HI %*% uu %*% HI
  attr(V,"hessian") <- H
  attr(V,"fisher") <- uu
  return(V)
}
