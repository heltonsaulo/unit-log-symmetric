
####################################################################################
################################### Packages #######################################
####################################################################################

require(ssym)
require(resample)
require(normalp) 
require(maxLik)
require(numDeriv)
require(extraDistr)

####################################################################################
################################### PDF and CDF ####################################
####################################################################################

dunitlogsym1 <- function(w, 
                         sigma,
                         eta,
                         xi, 
                         family,
                         log = FALSE){
  
  
  if(family=="Normal"){
    xi  <- 0
     pdf <- 1/(w*(1-w)*sigma) * dnorm(log((w/(eta*(1-w)))^(1/sigma)))
  }
 
  if(family=="Student"){
    if(xi[1]<=0) stop("nu must be positive!!",call.=FALSE)
    nu <- xi[1]
    pdf <- 1/(w*(1-w)*sigma) * dt(log((w/(eta*(1-w)))^(1/sigma)),df=nu)
  }
  
  if(family=="Laplace"){
    xi  <- 0
    pdf <- 1/(w*(1-w)*sigma) * dlaplace(log((w/(eta*(1-w)))^(1/sigma)))
  }
  

  
  if (log==TRUE){pdf <-log(pdf)}
  
  return(pdf) 
}
 

punitlogsym1 <- function(w, 
                         sigma, 
                         eta,
                         xi, 
                         family){
  
  
  if(family=="Normal"){
    xi  <- 0
    cdf <- pnorm(log((w/(eta*(1-w)))^(1/sigma)))
  }
  
  if(family=="Student"){
    if(xi[1]<=0) stop("nu must be positive!!",call.=FALSE)
    nu <- xi[1]
    cdf <-  pt(log((w/(eta*(1-w)))^(1/sigma)),df=nu)
  }
  
  
  if(family=="Laplace"){
    xi  <- 0
    cdf <- plaplace(log((w/(eta*(1-w)))^(1/sigma)))
  }
  
  
  return(cdf)
}




####################################################################################
############################### Random Numbers #####################################
####################################################################################


runitlogsym.inverse.sigma <- function(n, sigma,  eta, xi, family) {
  

  
  if(family=="Normal"){ 
    
    inte <- qnorm(runif(n)) *sigma
    
    W <- (eta*exp(inte)) / (1+(eta*exp(inte)))
    
  }
  
  if(family=="Student"){ 
    
    inte <- qt(runif(n), df = xi[1] ) *sigma 
    
    W <- (eta*exp(inte)) / (1+(eta*exp(inte)))
    
  }
  
  if(family=="Laplace"){ 
    
    inte <- qlaplace(runif(n)) *sigma
    
    W <- (eta*exp(inte)) / (1+(eta*exp(inte)))
  }
    
  return(W)
}



####################################################################################
################################ ML estimation #####################################
####################################################################################


mle.unitlogsym <- function(dat,family,print=TRUE) {
  
  if(family!="Normal" & family!="Student"  & family!="Laplace")
    stop("family of distributions specified by the user is not supported!!",call.=FALSE)
   

    ## log-likelihood
    loglik <- function(theta,dat,xi,family) { 
      
      sigma <- theta[1]
      eta   <- theta[2] 
            
      loglikvalue <- numeric()
      
      for(i in 1:length(dat)){
        
        loglikvalue[i] <- ( log( dunitlogsym1(w = dat[i], 
                                              sigma = sigma,
                                              eta = eta,
                                              xi = xi, 
                                              family = family,
                                              log = FALSE) ) )
      }
      
      return(sum(loglikvalue))  
    }
    
    n <-  length(dat)
    sigma_guess <- sqrt( 1/n * sum( (log(dat/(1-dat)))^2 ) ) 
    eta_guess   <- median(dat)
      
    startvalues <- c(sigma_guess,eta_guess)
    
    if(family == "Normal"){
      xi = 0
      opt = maxBFGS(fn = loglik, start = startvalues, dat=dat,xi=xi,family="Normal",control = list(iterlim = 5000))
      log.lik.est = opt$maximum #-opt$value
      estimates   = opt$estimate #opt$par
      fisher = -opt$hessian
      se     = sqrt(diag(solve(fisher)))
    }
    
    
    if(family == "Laplace"){
      xi = 0
      opt = maxBFGS(fn = loglik, start = startvalues, dat=dat,xi=xi,family="Laplace",control = list(iterlim = 5000))
      log.lik.est = opt$maximum #-opt$value
      estimates   = opt$estimate #opt$par
      fisher = -opt$hessian
      se     = sqrt(diag(solve(fisher)))
      
    }
    
    
    if(family == "Student"){
    nus           <- seq(1,10,by=1)
    logLikelihood <- rep(NA,length(nus))
    
    for(i in 1:length(nus)){
      opt = maxBFGS(fn = loglik, start = startvalues, dat=dat,xi=nus[i],family="Student",control = list(iterlim = 5000))
      logLikelihood[i] = opt$maximum #-opt$value
    }
    
    maximum <- max(logLikelihood)
    best    <- which(logLikelihood == maximum)
    xi      <- nus[best]
    
    
    opt = maxBFGS(fn = loglik, start = startvalues, dat=dat,xi=xi,family="Student",control = list(iterlim = 5000))
    log.lik.est = opt$maximum #-opt$value
    estimates   = opt$estimate #opt$par
    fisher = -opt$hessian
    se     = sqrt(diag(solve(fisher)))
    
    
    }
    
   
      
    ## Information criteria
    kp     = length(startvalues)
    AIC   = - 2 * log.lik.est + 2 * (kp)
    AICc  = AIC + (2 * (kp) * ((kp) + 1)) / (n - (kp) - 1)
    BIC   = - 2 * log.lik.est + log(n) * (kp)
    HIC   = - 2 * log.lik.est + log(log(n)) * (kp)
  
    
    
    confidence.int <- cbind(estimates - qnorm(0.975)*se, estimates + qnorm(0.975)*se)
    colnames(confidence.int) <- c("2.5%","97.5%")
    
    
    zstatpars  = estimates / se
    pvalorpars = 2 * pnorm(abs(zstatpars), lower.tail = F)
    

    tb <- miscTools::coefTable(estimates,se, df=(n-(kp)))
    
    ## residuals
    zhat        <- punitlogsym1(w=dat, 
                                sigma=estimates[1],
                                eta=estimates[2],
                                xi=xi, 
                                family=family)
      
    GCSresidual <- -log(1-zhat)
    RQresidual  <- qnorm(zhat)
    
    
    if(print == TRUE & family=="Normal" | family=="Laplace"){
      cat("\n")
      cat("--------------------------------------------------------------\n")
      cat("                 Unit log-symmetric distribution              \n")
      cat("--------------------------------------------------------------\n")
      cat("--------------------------------------------------------------\n")
      cat("Maximum Likelihood estimation \n")
      cat("Log-Likelihood:", log.lik.est, "\n")
      cat("AIC:", AIC, "AICc:", AICc, "BIC:", BIC, "\n")
      cat("Number of observations:", n, "\n")
      cat("Family:", family, "\n")
      cat("--------------------------------------------------------------\n")
      cat("Coefficients:\n")
      printCoefmat(tb, signif.stars = TRUE, signif.legend = TRUE, digits = 4)
      cat("--------------------------------------------------------------\n")
    }
    
    if(print == TRUE & family=="Student"){
      cat("\n")
      cat("--------------------------------------------------------------\n")
      cat("                 Unit log-symmetric distribution              \n")
      cat("--------------------------------------------------------------\n")
      cat("--------------------------------------------------------------\n")
      cat("Maximum Likelihood estimation \n")
      cat("Log-Likelihood:", log.lik.est, "\n")
      cat("AIC:", AIC, "AICc:", AICc, "BIC:", BIC, "\n")
      cat("Number of observations:", n, "\n")
      cat("Family:", family, "\n")
      cat("Xi (best):", xi, "\n")
      cat("--------------------------------------------------------------\n")
      cat("Coefficients:\n")
      printCoefmat(tb, signif.stars = TRUE, signif.legend = TRUE, digits = 4)
      cat("--------------------------------------------------------------\n")
    }
    

    
    rval = list(coefficients = estimates, 
                nu  = xi, 
                se = se, 
                conf.int = confidence.int,
                pvalor =   pvalorpars,
                information.criterions = list(aic = AIC,bic = BIC,aicc=AICc,hic=HIC), 
                loglik = log.lik.est,
                n = n, 
                family = family,
                GCSresidual = GCSresidual,
                RQresidual = RQresidual)
    
    
    return(rval)
    
    
}

  

