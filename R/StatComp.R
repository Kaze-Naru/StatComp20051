#' @title The mean vector for multinormal distribution.
#' @description Generate the mean vector for multinormal distribution.
#' @param p Dimension of the data.
#' @param ter length of nonzero components.
#' @param co  co*p is the mean of the nonzero components.
#' @return a vector
#' @export
ftheta <- function(p,ter,co)
{
  ttheta <- rep(0,p)
  for(ii in 1:ter)
  {
    ttheta[ii] <- log(co*p)
  }
  return(ttheta)  
}
#' @title Correlation matrix.
#' @description Generate correlation matrix for multinormal distribution.
#' @param rho cor(xi,xj) = rho^abs(i-j).
#' @param p Dimension of the data.
#' @return A matrix
#' @export
fSIGMA <- function(rho,p)
{
  repzero <- rep(0,p*p)
  SIGMA <- matrix(data = repzero, nrow = p, ncol = p)
  for(ii in 1:p)
  {
    for(jj in 1:p)
      SIGMA[ii,jj] <- rho^(abs(ii-jj))
  }
  return(SIGMA)  
}
#' @title Log of compositional dataset
#' @description Generate compositional dataset, and apply log() on each component
#' @param m the number of samples.
#' @param p Dimension of the data.
#' @param ter length of nonzero components.
#' @param co the mean of the nonzero components.
#' @param rho cor(xi,xj) = rho^abs(i-j).
#' @return A matrix
#' @importFrom mvtnorm rmvnorm
#' @export
fZ <- function(m,p,ter,co,rho)
{
  SIGMAtemp <- fSIGMA(rho,p)
  thetatemp <- ftheta(p,ter,co)
  wtemp <- rmvnorm(m, mean = thetatemp, sigma = SIGMAtemp)  #m rows and p columns
  rrepzero <- rep(0,m*p)
  xtemp <- matrix(data = rrepzero, nrow = m, ncol = p)
  for(ii in 1:m)
  {
    xexptemp <- exp(wtemp[ii,])
    sumexptemp <- sum(xexptemp)
    for(jj in 1:p)
      xtemp[ii,jj] <- xexptemp[jj]/sumexptemp
  }
  ztemp <- log(xtemp)
  return(ztemp)
}
#' @title Response vector
#' @description Generate response vector.
#' @param zk Log of componsional dataset.
#' @param beta Vector of coefficient. Sum(beta)=0.
#' @param mean mean vector when generating composional dataset.
#' @param sd standard error for error term in linear regression model.
#' @return A vector
#' @importFrom stats rnorm
#' @export
fY <- function(zk,beta,mean,sd)
{
  q <- length(zk[,1])
  epsilontemp <- rnorm(q, mean , sd)
  y <- zk%*%beta+epsilontemp
  return(y)  
}
#' @title Soft thresholding operator
#' @description Function of soft thresholding operator.
#' @param t A number of input.
#' @param lambda Threshold
#' @return A number
#' @export
fsto <- function(t,lambda)
{
  if(abs(t)<=lambda)
    return(0)
  else
  {
    if(t>0)
      return(t-lambda)
    else
      return(t+lambda)
  }  
}
#' @title Estimation of coefficients
#' @description Find the estimation of coefficents for high dimensional compositional data, L1 and L2 penalized functions are used
#' @param yy a vector of response variable y
#' @param Z a matrix of independent variables. It should be log of the compositional data matrix, in which each row should sum to 1, and each component should be between 0 and 1.
#' @param betastart betastart is starting vector for the estimation of coefficient beta.
#' @param alphastart alphastart is starting vector for the estimation of lagrange multiplier for sum of the coefficients.
#' @param ac1 ac1 is accuracy when getting an estimation of parameters. High accuracy means more accuracy results.
#' @param ac2 ac2 is accuracy when getting an estimation of parameters. High accuracy means more accuracy results.
#' @param lambda1 coefficient for L1 penlized function
#' @param lambda2 coefficient for L2 penlized function
#' @param mu mu is multiplier of the square of the sum of the coefficients in Lagrange method.
#' @return A matrix
#' @import utils 
#' @export
#' @useDynLib StatComp20051
festimate <- function(yy, Z, betastart, alphastart, ac1, ac2, lambda1, lambda2, mu)
{
  pp <- length(Z[1, ]) 
  nn <- length(Z[, 1])  
  beta1 <- rep(betastart, pp)  # betastart is starting vector for the estimation of coefficient beta.
  beta2 <- beta1+100
  alpha1 <- alphastart  # alphastart is starting vector for the estimation of lagrange multiplier for sum of the coefficients.
  alpha2 <- alphastart+100
  vv <- rep(0, pp)
  for(kk in 1:pp)
    vv[kk] <- sum(Z[, kk]*Z[, kk])/nn ##squared sum of each column
  
  while(sqrt(sum((beta1-beta2)*(beta1-beta2))+(alpha1-alpha2)^2)>ac2)
  {
    beta2 <- beta1+100   
    alpha2 <- alpha1
    while(sqrt(sum((beta1-beta2)*(beta1-beta2)))>ac1)
    {
      beta2 <- beta1
      for(jj in 1:pp)
      {
        tempj <- 0
        for(i in 1:nn)
        {
          tempj <- tempj+1/nn*(yy[i]-sum(Z[i,]*beta1)+Z[i,jj]*beta1[jj])*Z[i,jj]
        }
        tj <- tempj-mu*(sum(beta1)-beta1[jj]+alpha2)
        beta1[jj] <- 1/(vv[jj]+2*lambda2+mu)*fsto(tj,lambda1)  # mu is multiplier of the square of the sum of the coefficients in Lagrange method.
      }    
    }
    alpha2 <- alpha1
    alpha1 <- alpha2+sum(beta1)      
  }
  
  vec <- c(beta1, betastart, alphastart, ac1, ac2, lambda1, lambda2, mu)
  write.csv(vec,"temptempest.csv")
  file.append("tempest.csv", "temptempest.csv") 
  return(beta1)  # return Estimation of coefficients
}