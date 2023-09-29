#####R-Project03 Niccolo Ducci 7036789#####

#Import Data

xobs = c(3560,3739,2784,2571,2729,
         3952,993, 1908,948,1172,
         1047,3138,5485,5554,2943,
         4969,4828)

x=c(NA,xobs)

yobs=c( 3,4,1,1,3,1,2,0,2,0,
        1,3,5,4,6,2,5,4)

cbind(x,yobs)

###ES1-Develop an EM algorithm

log.lik.incomplete<-function(beta,tau1,taui,yobs,xobs){
  sumFact=0
  for(i in xobs){
    for(j in 1:i){
      sumFact <- sumFact + log(j)
    }
  }
  tau <- c(tau1,taui)
  ( sum( log( (exp(-beta*tau)*(beta*tau)^(yobs))/factorial(yobs) ) ) + 
      sum( -taui+xobs*log(taui) ) - sumFact  )
}

EM.poisson<-function(xobs, yobs, tau1=NULL, beta=NULL, tol = 1e-06, maxit =1000){
  
  ##starting values
  if(is.null(tau1) == TRUE){
    tau1<- mean(xobs)
  }
  
  if(is.null(beta) == TRUE){
    beta<- mean(yobs)/tau1
  }
  
  
  taui=xobs
  Phi=c(beta, tau1, taui ) 
  lik.obs <- log.lik.incomplete(beta,tau1,taui,yobs,xobs)
  n.it <- 0
  converged <- FALSE
  
  while(converged==FALSE){
    tau1S<- tau1
    betaS<- beta
    tauiS <- taui
    
    ##E-step 
    #we do not need this step as we already have 
    #the solutions(anyway x1=tau1)
    
    ##M-step
    beta <- sum(yobs)/( tau1+sum(xobs) )
    tau1 <- (tau1 + yobs[1]) / (beta+1)
    for(i in 2:length(yobs)){
      taui[i-1] <- (xobs[i-1]+yobs[i]) / (beta+1)
    }
    
    Phi <- rbind(Phi,c(beta, tau1, taui))
    lik.obs <-c(lik.obs, log.lik.incomplete(beta,tau1,taui,yobs,xobs) )
    
    ##Convergence check

    if(abs(tau1S-tau1)< tol & abs(betaS-beta)< tol &
       all(abs(tauiS-taui)<tol) & n.it <= maxit){
      converged <- TRUE
    }
    
    if(converged == FALSE & n.it < maxit){
      n.it <- n.it +1
    }
    
    if(converged == FALSE & n.it == maxit){
      warning("The maximum number of iterations has been reached without reaching convergence")
      break
    }
    
    
  }#end loop over while 

  out <- list(phi=c(beta=beta, tau1=tau1, taui=taui), Phi=Phi,
              converged = converged, n.it = n.it, lik.obs=lik.obs)
  out    
  
}

EM.result <- EM.poisson(xobs, yobs, maxit=20000)
Phi.em <- EM.result$phi
options(scipen=5)
Phi.em

#plottin the likelihood
plot(1:length(ll),ll,type="l",
     col="navy", xlab="Number of Iteration", ylab="log-likelihood")

###ES2-Compare the direct solution with the EM solution

beta.hat <- sum(yobs)/sum(xobs)

tau1.hat <- yobs[1]/beta.hat

taui.hat <- (xobs+yobs[2:length(yobs)]) / (beta.hat+1)

Phi.hat <- c(beta.hat, tau1.hat, taui.hat)

cbind(Phi.em, Phi.hat, diff= abs(Phi.em-Phi.hat))

#log-likelihood comparison
log.lik.incomplete(beta.hat,tau1.hat,taui.hat,yobs,xobs)
-112.2879

log.lik.incomplete(Phi.em[1],Phi.em[2],Phi.em[3:19],yobs,xobs)
-112.1901

#EM solutions slightly better

#####

