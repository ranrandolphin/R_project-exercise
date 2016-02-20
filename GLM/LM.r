##### 1. Write an R function that takes as parameters a matrix X and a vector Y and returns the vector of least squares estimates, the fitted values and the residuals.
r.function <- function(matrix.X, vector.Y) {
  
  # Least square estimator
  LSE <- solve(t(matrix.X) %*% matrix.X) %*% t(matrix.X) %*% vector.Y
 
  # Fitted values
  fitted <- matrix.X %*% LSE
  
  # Residuals
  residuals <- vector.Y - fitted
  
  final <- list(LSE, fitted, residuals)
  return(final)
  
}

#### Test
X <- matrix(c(1,2,4,5,6,7,8,9,12,34,56,26,12,11,9,10), nrow = 4, ncol = 4, byrow = TRUE)
Y <- matrix(c(12,32,54,67), nrow = 4, ncol = 1, byrow = TRUE)
r.function(X,Y)

##### 2. Write an R function that takes as parameters a vector of means $\mu$ a covariane matrix $S$ and a number $m$ and generates a sample of size m from the multivariate normal distribution with mean $\mu$ and covariance $S$. Plot the observations for 2 different bivariate normals (show both plots) with your choice of means $\mu$ and covariances $S$.

library(matrixcalc)

r.function <- function(mean.mu, cov.S, size.m){
  if(is.positive.definite(cov.S) == TRUE){
    p <- qr(cov.S)$rank 
    mu <- t(t(mean.mu))
    # generates a size.m by p random normal matrix
    randn.x <- matrix(rnorm(p * size.m), size.m, p)
    
    # find cholesky decomposition of cov.S, chol(covS)
    eval <- eigen(cov.S)$values
    evec <- eigen(cov.S, symmetric = TRUE)$vector
    diag <- diag(sqrt(pmax(eval, 0)), p)
    new.sigma <- evec %*% diag %*% t(randn.x)
    
    # set a null matrix to multivariate normal distribution
    mnd <- matrix(NA,nrow = p, ncol = size.m) 
    for( i in 1 : size.m){
      mnd[,i] <- mu + new.sigma[,i]
    }
    return(mnd)
    }
  else {
    stop("Covariance matrix is not positive definite!")
    }
}

Size <- 10000
Mean1 <- rep(1,2)
Cov1 <- matrix(c(10,3,3,4),2,2)
multivariate1 <- r.function(Mean1, Cov1, Size)
p1 <- plot(multivariate1[1,], multivariate1[2,], main ="First Plot of Bivariate Normals")

##### 3.1  Plot the power of the F-test as a function of $c/\sigma$ for 4 differents values of $n$
# First Case with noncentrality parameter is (4/3) * n* c^2/sigma^2
power1 <- function(size, c_sigma1){
  noncentrality <- (4/3) * (c_sigma1^2) * size
  probability <- 1 - pf(qf(0.95, 2, size - 2), 2, size - 4, ncp = noncentrality)
  final <- list(c_sigma1, probability)
  return(final)
}

c_sigma <- 5 / seq(1, 150, by = 0.1)
power1_5 <- power1(5, c_sigma) 
power1_10 <- power1(10, c_sigma)
power1_15 <- power1(15, c_sigma)
power1_100 <- power1(100, c_sigma)

plot(power1_5[[1]],power1_5[[2]], ylim = c(0, 1), type = "l", xlab = "c/sigma1", 
     ylab = "Power of the F-test", 
     main = "The Power of F-test associated with c/sigma1", col = "red", lty = 1, lwd=2)
lines(power1_10[[1]], power1_10[[2]], col = "blue", lty = 2, lwd=2)
lines(power1_15[[1]], power1_15[[2]], col = "orangered", lty = 3, lwd=3)
lines(power1_100[[1]], power1_100[[2]], col = "green", lty = 4, lwd=2)

legend("bottomright", c(expression(paste(F[list(2,1)]("c_sigma1"))),
 expression(paste(F[list(2,6)]("c_sigma1"))),
 expression(paste(F[list(2,11)]("c_sigma1"))),
 expression(paste(F[list(2,96)]("c_sigma1")))),
 lty = 1 : 4,
 col = c("red", "blue", "orangered", "green"), lwd=3)
 
 ##### 3.2 For one of the values of $n$ above verify the finding by a simulation
 
 # chose n = 100 as the value to verify and noncentrality parameter in both case1 and case2
c_sigma <- (4/3) * 100 * (1/2) ^ 2
U1 <- rchisq(100, df = 2, ncp = c_sigma)
U2 <- rchisq(100, df = 98)
F <- (U1 / 2) / (U2 / 98) 
hist(F, main = "Noncentrality F distribution by rchisq function", xlab = " ")

#Althernatively,
simf <- rf(100, 2, 98, c_sigma)
hist(simf, main = "Noncentrality F distribution by rf function", xlab = " ")

#Calculate power of F-test by a simulation

sigmaval <- seq(1, 100, by = 1)
cvalue <- rep(5, length(sigmaval))
powerval <- NULL

csigma <- function(c, sigma)
{
  samplesize <- c(100, 100, 100) #sample size for three groups, each n = 100
  
  mu <- c(c, c, 3*c) # true mean for each group, miu1 = miu2 = 1, miu3 = 3
  sigma <- c(sigma, sigma, sigma) #true standard deviations for three groups
  
  mus <- rep(mu, times = samplesize)
  sigmas <- rep(sigma, times = samplesize)
  
  group <- factor(rep(1:3, times = samplesize)) #factor three groups which will be used in ANOVA
  nsims <- 100
  
  simulation <- function(){
    Y.data <- rnorm(sum(samplesize), mus, sigmas)
    model <- aov(lm(Y.data ~ group))
    p_value <- summary(model)[[1]][["Pr(>F)"]][1]
  }
  pVals <- replicate(nsims, simulation())
  power <- mean(pVals < 0.05)
  return(power)
}


for ( i in 1:length(sigmaval)){
 
 powerval[i] <- csigma(cvalue[i], sigmaval[i])

}

cs <- cvalue/sigmaval
plot(cs, powerval, main = "Noncentrality F distribution", type = "l", col = "red")
