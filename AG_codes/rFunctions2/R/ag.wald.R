ag.wald <- 
## test on ancestral graph parameters
function(theta,psi,contr,u=0,N,n)
{
  if(!is.na(n[2]))
    {
      if(sum(n)/N==n[1]) n <- n[1]
      else
        {
          if(!length(n)==N) stop("Number of subjects is different than number of subjects entered in trial vecotr")
          N <- 1
          n <- sum(n)
        }
    }
    k <- nrow(as.matrix(t(contr)))
    wald <- ((N*n - k)/(N*n*k))*t(contr%*%theta - u)%*%solve(contr%*%psi%*%contr)%*%(contr%*%theta - u)
    df1 <- k
    df2 <- N*n - k
    pw <- pf(wald,df1,df2,lower.tail=FALSE)
    return(list(w=wald,p=pw,df1=df1,df2=df2))
}

# ag.wald(theta1,var.a,c(1,-1,0,0,0,0,0,0,0,0),N=10,n=c(5,6,7,8,9,10,20,20,20,4))
