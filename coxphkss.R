## Median bias reduction for Cox regression ##
library(survival)
library(coxphf)
library(flexsurv)


###################################################
###           score function                   ####
###################################################

# beta : parameter to be estimated
# data : dataset(time,cens,covariates)
score <- function(beta,data){
  p <- length(beta)
  dim <- ncol(data) 
  index <- order(data[,1])
  data <- data[index,]
  status <- data[,2]
  X <- as.matrix(data[,3:dim])
  n <- nrow(data)
  ofInterest <- which(status==1)
  m <- length(ofInterest)
  score_k <- function(k){
    res <-0
    for (i in 1:m) {
      R_i <- ofInterest[i]:n
      weights <- exp(X[R_i,]%*%beta)/sum(exp(X[R_i,]%*%beta))
      res <- res+ (X[ofInterest[i],k]-sum(X[R_i,k]*weights))
    }
    res
  }
  sapply(1:p,score_k)
}

####################################################
####        observed    information             ####
####################################################

obsInfo <- function(beta,data){
  p <- length(beta)
  dim <- ncol(data) 
  index <- order(data[,1])
  data <- data[index,]
  status <- data[,2]
  X <- as.matrix(data[,3:dim])
  n <- nrow(data)
  ofInterest <- which(status==1)
  m <- length(ofInterest)
  result <- matrix(NA,p,p)
  obsinfo_kl <- function(k,l){
    res <- 0
    for (i in 1:m) {
      R_i <- ofInterest[i]:n
      weights <- exp(X[R_i,]%*%beta)/sum(exp(X[R_i,]%*%beta))
      res <- (res+ (sum(X[R_i,k]*X[R_i,l]*weights)
                    -sum(X[R_i,k]*weights)*sum(X[R_i,l]*weights)) )
    }
    res
  }
  for (k in 1:p) {
    for (l in k:p) {
      result[k,l] <-result[l,k]<- obsinfo_kl(k,l)
    }
  }
  return(result)
}


##########################################
###------ third order cumulant   --##
##########################################

obsThOrdcum <- function(beta,data){
  p <- length(beta)
  dim <- ncol(data) 
  index <- order(data[,1])
  data <- data[index,]
  status <- data[,2]
  X <- as.matrix(data[,3:dim])
  n <- nrow(data)
  ofInterest <- which(status==1)
  m <- length(ofInterest)
  result <- array(NA,c(p,p,p))
  obscum_rst <- function(r,s,t){
    res <- 0
    for (i in 1:m) {
      R_i <- ofInterest[i]:n
      weights <- exp(X[R_i,]%*%beta)/sum(exp(X[R_i,]%*%beta))
      res <- (res+ (sum(X[R_i,r]*X[R_i,s]*X[R_i,t]*weights)-sum(X[R_i,t]*weights)*sum(X[R_i,r]*X[R_i,s]*weights)
                    -sum(X[R_i,s]*weights)*sum(X[R_i,r]*X[R_i,t]*weights)-sum(X[R_i,r]*weights)*sum(X[R_i,s]*X[R_i,t]*weights)
                    +2*sum(X[R_i,r]*weights)*sum(X[R_i,s]*weights)*sum(X[R_i,t]*weights)) )
    }
    res
  }
  for (r in 1:p) {
    for (s in 1:p) {
      for (t in s:p) {
        result[r,s,t] <- result[r,t,s]<- obscum_rst(r,s,t)
      }
    }
  }
  return(result)
}
####################################################
####         Mean adjusted score               ####
####################################################

score_meanBR <- function(beta,data){
  invInfo  <- solve(obsInfo(beta,data))
  cum <- obsThOrdcum(beta,data)
  score(beta,data)+0.5*sapply(1:length(beta),function(r) sum(diag(invInfo%*%cum[r,,])))
}

####################################################
####       Median  adjusted score               ####
####################################################

score_medianBR <- function(beta,data){
  
  adjustment <- function(beta,data,nu_r_s_t )
  {
    p <- length(beta)
    b_vector <- numeric(p)
    LKsum <- function(t){
      return ( sum(diag(invInfo%*%(nu_r_s_t[t,,]))) )
    }
    LKsum2 <- function(t){
      return ( sum(diag(K_r%*%( nu_r_s_t[t,,]/3))) )
    }
    for (r in 1:p) 
    {
      inverse_info_r <- invInfo[r, ]
      K_r <- tcrossprod(inverse_info_r) / inverse_info_r[r]
      b_vector[r] <- sum(inverse_info_r*sapply(1:p,LKsum2)) 
    }
    return( 0.5*sapply(1:p,LKsum) - info%*%b_vector )
  }
  
  info <- obsInfo(beta,data)
  invInfo  <- solve(info)
  cum <- obsThOrdcum(beta,data)
  score(beta,data)+adjustment(beta,data,cum)
}

##################################################
## function to get median BR estimates         ###
##################################################

# beta : parameter to be estimated
# data : dataset(time,cens,covariates)
# type : either "AS_median" for medianBR estimator or "AS_mean" for meanBR estimator
# maxit: maximum number of iterations
# tol  : stopping criterion
coxphkss<- function(beta,data,type="AS_median",maxit=50,tol=1e-08) {
  inverse_info <- try(solve(obsInfo(beta, data)),silent = TRUE)
  if (failed <- inherits(inverse_info, "try-error")) {
    warning("failed to invert the information matrix: iteration stopped prematurely")
    break
  }
  adjScore <-  score_medianBR(beta,data)
  if(type=="AS_mean") adjScore <-  score_meanBR(beta,data)
  step <-  inverse_info%*%adjScore
  
  if(maxit <=0 ) 
    warning("meanBR and medianBR  cannot be performed with maxit <= 0")
  iter <- 0
  if( maxit > 0 ) {
    for (iter in seq.int(maxit)) {
      stepPrev <- step
      stepFactor <- 0
      testhalf <- TRUE
      while (testhalf & stepFactor < 11) {
        beta <- beta + 2^(-stepFactor)*step 
        inverse_info <- try(solve(obsInfo(beta, data)),silent = TRUE)
        if (failed <- inherits(inverse_info, "try-error")) {
          warning("failed to invert the information matrix: iteration stopped prematurely")
          break
        }
        adjScore <-  score_medianBR(beta,data)
        if(type=="AS_mean") adjScore <-  score_meanBR(beta,data)
        step <-  inverse_info%*%adjScore
        stepFactor <- stepFactor + 1
        testhalf <- drop(crossprod(stepPrev) < crossprod(step))  
      }
      if (failed | (all(abs(step) < tol))) {
        break
      }
    }
  }
  #print(iter)
  beta <- drop(beta)
  converged <- maxit > iter
  exp.beta <- exp(beta)
  # print(beta)
  jacobian <- diag(exp.beta,length(beta))
  # print(jacobian)
  # print(inverse_info)
  var.exp.coef <- crossprod(jacobian,inverse_info)%*%jacobian
  list(coefficients=beta,se=sqrt(diag(inverse_info)),
       exp.coefficients=exp.beta,exp.se=sqrt(diag(var.exp.coef)),
       nIter=iter,adjustedScore=adjScore,convergence=converged
  )
}

