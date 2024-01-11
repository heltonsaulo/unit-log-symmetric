#####################################################
Unit-log-symmetric models
#####################################################

Example with the unit-lognormal model
eta1 <- 1

eta2 <- 1

sigma1 <- 0.5

sigma2 <- 0.5

rho <- 0.5

xi <- -0.9

Sigma <- matrix(c(sigma1^2,rho * sigma1 * sigma2,rho * sigma1 * sigma2,sigma2^2), 2,2, byrow = TRUE)

data = rmlogsym.Choleski(n=1000, etas=c(eta1,eta2), sigmas=c(sigma1,sigma2), rho=rho, xi=xi, family="Powerexp")

est <- mle.bivlogsym(data=data,family="Powerexp", xi=c(-0.9), print = T)
