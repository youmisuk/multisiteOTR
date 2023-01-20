###########################################################################
# Designing Optimal, Data-Driven Policies from Multisite Randomized Trials
# by Youmi Suk and Chan Park

# :: create a function for creating data from a multisite randomized trial.
multilevel_DGP <- function(J=50, n_j=50, n_j_test = 50, beta1= 0, overall=0.5, 
                           beta2 = 2.5, nonlinear.effect= 0, R.sd=0.9) {
  
  id <- rep(1:J, each=n_j) # cluster ids
  n.obs  <- J * n_j # number of total observations in training data
  n.obs_test  <- J * n_j_test # number of total observations in test data
  n.vars <- 5 # number of covariates
  x <- matrix(rnorm(n.obs * n.vars, sd = 1), n.obs, n.vars) # generate covariates in training data
  x_test <- matrix(rnorm(n.obs_test * n.vars, sd = 1), n.obs_test, n.vars)  # generate covariates in test data
  colnames(x) <- colnames(x_test) <- paste0("X", 1:n.vars) # give names for covariates
  
  R <- rnorm(J, 0, sd=R.sd) # generate a cluster-level covariate that affects treatment assignment
  U1 <- rnorm(J, 0, 1)  # generate a cluster-level covariate that is assumed to be unmeasured and affects the decision rule
  
  # simulate a randomized treatment but randomization may be affected by clusters in multisite randomized trials
  xbetat   <-  log(overall/(1-overall)) + 0.8*R[id]
  trt.prob <- exp(xbetat) / (1 + exp(xbetat))
  trt      <- rbinom(n.obs, 1, prob = trt.prob) # treatment status in training data
  trt_test      <- rbinom(n.obs_test, 1, prob = trt.prob)   # treatment status test data
  
  # simulate delta (benefit scores)
  delta <- (0.2 + 0.5*x[,3] - 0.2 * x[,4] - 0.2 * x[,5] + beta1*U1[id]) # delta in training data
  delta_test <- (0.2 + 0.5*x_test[,3] - 0.2 * x_test[,4] - 0.2 * x_test[,5]  + beta1*U1[id]) # delta in test data

  # simulate potential outcomes and the observed outcome
  e <- rnorm(n.obs) # random error in training data
  xbeta <- 0.7 + x[,1] + x[,2] -2* x[,3] + x[,4] + 0.5* x[,5] + 
    nonlinear.effect*(exp(-0.6*x[,1]*x[,2]+0.3*x[,1]))  + beta2*U1[id] # the main effects of covariates in the outcome model in training data
  
  y0 <- xbeta  + e  # potential control outcome in training data
  y1 <- xbeta + delta  + e  # potential treatment outcome in training data  
  y <- xbeta + delta * trt + e # observed outcome in training data
  
  e_test <- rnorm(n.obs_test) # random error in test data
  xbeta_test <- 0.7 + x_test[,1] + x_test[,2] -2* x_test[,3] + x_test[,4] + 0.5* x_test[,5] + 
    nonlinear.effect*(exp(-0.6*x_test[,1]*x_test[,2]+0.3*x_test[,1])) + beta2*U1[id]  # the main effects of covariates in the outcome model in test data
  
  y0_test <- xbeta_test  + e_test  # potential control outcome in test data
  y1_test <- xbeta_test + delta_test  + e_test  # potential treatment outcome in test data     
  y_test <- xbeta_test + delta_test * trt_test + e_test # observed outcome in test data
  
  # save data
  dat <- data.frame(id=as.factor(id), Y=y, Y0=y0, Y1=y1, trt=trt, trt.prob, x, R=R[id], U1=U1[id], delta) # training data
  
  dat_test <- data.frame(id=as.factor(id), Y=y_test, Y0=y0_test, Y1=y1_test, trt=trt_test, trt.prob, 
                         x_test, R=R[id], U1=U1[id], delta=delta_test) # test data
  
  
  return(list(dat=dat, dat_test=dat_test))
}
