###########################################################################
# Designing Optimal, Data-Driven Policies from Multisite Randomized Trials
# by Youmi Suk and Chan Park

# :: load packages
library(psych)
library(personalized)
library(lme4)
source("Qlearn.R")
source("DataGeneratingModels.R")

# :: function for R package ''personalized''
propensity.const <- function(x, trt)
{
  rep(prop.table(table(trt))[2], length(trt))
}

propensity.clustermean <- function(x, trt)
{
  prop.table(table(dat$id, trt), 1)[,2][dat$id]
}

aug.func_l <- function(x, y) {
  aug <- lm(y ~ covmat_aug)
  return(predict(aug))
}

aug.func_s <- function(x, y) {
  aug <- lm(y ~ covmat_aug + dat$id)
  return(predict(aug))
}

# :: conduct simulations
reps <- 2 # num. of replications

rslt_1 <- rslt_2 <- rslt_3 <- rslt_4 <- list() # lists for saving benefit scores under sample size conditions 1, 2, 3, and 4 for (25, 25), (25, 150), (150, 25), and (150, 150), respectively

for (smpl.size in list(c(25, 25), c(25, 150), c(150, 25), c(150, 150))) {
  print(smpl.size)  
  
  for (i in 1:reps) {
    
    print(i)  
    
    # :: create a simulated data
    repeat{
      
      # :: choose a design that you like to replicate.
      df <- multilevel_DGP(J=smpl.size[1], n_j=smpl.size[2], n_j_test = smpl.size[2], beta1= 0, R.sd=0.9) # for Design 1
      #df <- multilevel_DGP(J=smpl.size[1], n_j=smpl.size[2], n_j_test = smpl.size[2], beta1= 0.2, R.sd=0.9) # for Design 2
      #df <- multilevel_DGP(J=smpl.size[1], n_j=smpl.size[2], n_j_test = smpl.size[2], beta1= 0, R.sd=0.9, nonlinear.effect=1) # for Design 3
      #df <- multilevel_DGP(J=smpl.size[1], n_j=smpl.size[2], n_j_test = smpl.size[2], beta1= 0, R.sd=1.5) # for Design 4
      
      dat <- df$dat
      dat_test <- df$dat_test
      
      temp.p1 <- prop.table(table(dat$id, dat$trt), 1)[,2]
      
      if ((sum(range(temp.p1) == 0) + sum(range(temp.p1) == 1)) == 0){
        break
      }
    }
    
    # ::::::::::::::::::::::::
    # :: Q-learning methods ::
    # ::::::::::::::::::::::::
    
    # :: implment Q-learning methods
    out1 <- Qlearn(Yname="Y", Trtname="trt", Xname=c("X1", "X2", "X3", "X4", "X5"), interTrt=formula(~ X3 + X4 + X5), IDname="id", 
                   traindata=dat, testdata=dat_test, method = c("Q-lm")) # baseline

    out2 <- Qlearn(Yname="Y", Trtname="trt", Xname=c("X1", "X2", "X3", "X4", "X5"), interTrt=formula(~ X3 + X4 + X5), IDname="id", 
                   traindata=dat, testdata=dat_test, method = c("Q-fe"))

    out3 <- Qlearn(Yname="Y", Trtname="trt", Xname=c("X1", "X2", "X3", "X4", "X5"), interTrt=formula(~ X3 + X4 + X5), IDname="id", 
                   traindata=dat, testdata=dat_test, method = c("Q-ri"))
    
    out4 <- Qlearn(Yname="Y", Trtname="trt", Xname=c("X1", "X2", "X3", "X4", "X5"), interTrt=formula(~ X3 + X4 + X5), IDname="id", 
                   traindata=dat, testdata=dat_test, method = c("Q-rs"))
    
    out5 <- Qlearn(Yname="Y", Trtname="trt", Xname=c("X1", "X2", "X3", "X4", "X5"), interTrt=formula(~ X3 + X4 + X5), IDname="id", 
                   traindata=dat, testdata=dat_test, method = c("Q-dm"))
    
    out6 <- Qlearn(Yname="Y", Trtname="trt", Xname=c("X1", "X2", "X3", "X4", "X5"), interTrt=formula(~ X3 + X4 + X5), IDname="id", 
                   traindata=dat, testdata=dat_test, method = c("Q-dm-ri"))
    
    out7 <- Qlearn(Yname="Y", Trtname="trt", Xname=c("X1", "X2", "X3", "X4", "X5"), interTrt=formula(~ X3 + X4 + X5), IDname="id", 
                   traindata=dat, testdata=dat_test, method = c("Q-dm-rs"))
    
    # :: get the benefit scores from Q-learning methods
    out1_score <- out1$benefit.scores
    out2_score <- out2$benefit.scores
    out3_score <- out3$benefit.scores
    out4_score <- out4$benefit.scores
    out5_score <- out5$benefit.scores
    out6_score <- out6$benefit.scores
    out7_score <- out7$benefit.scores
    
    # :::::::::::::::::::::::
    # :: weighting methods ::
    # :::::::::::::::::::::::
    
    # :: create covariate matrices
    covmat <- as.matrix(dat[, c(paste0("X", c(3, 4, 5)))])
    covschmat <- model.matrix(~ X3 + X4 + X5 + id + 0, dat)
    covschmat <- covschmat[, !colnames(covschmat) %in% "id1"]
    covmat_demean <- apply(covmat, 2, function(x) x - ave(x, dat$id))
    dat$Y_dm <- dat$Y - ave(dat$Y, dat$id)
    covmat_aug <- model.matrix(~ X1 + X2 + X3 + X4 + X5, dat)[,-1]
    
    schpen_only <- c(0, rep(0, ncol(covmat)), rep(1, length(levels(dat$id))-1))
    
    # :: implment weighting methods
    m0 <- fit.subgroup(x = covmat, y = dat$Y,
                       trt = dat$trt, propensity.func = propensity.const,
                       loss   = "sq_loss_lasso", 
                       nfolds = 5) # the baseline weighting method that uses the constant propensity score
    
    m1 <- fit.subgroup(x = covmat, y = dat$Y,
                       trt = dat$trt, propensity.func = propensity.clustermean,
                       loss   = "sq_loss_lasso", 
                       nfolds = 5) # the modification that uses cluster-specific propensity scores
    
    m1_aug_l <- fit.subgroup(x = covmat, y = dat$Y,
                             trt = dat$trt, propensity.func = propensity.clustermean,
                             loss   = "sq_loss_lasso", 
                             nfolds = 5,
                             augment.func = aug.func_l) # the modification that uses cluster-specific propensity scores and augmentation terms with covariates only
    
    m1_aug_s <- fit.subgroup(x = covmat, y = dat$Y,
                             trt = dat$trt, propensity.func = propensity.clustermean,
                             loss   = "sq_loss_lasso", 
                             nfolds = 5,
                             augment.func = aug.func_s) # the modification that uses cluster-specific propensity scores and augmentation terms with covariates and cluster dummies

    m2 <- fit.subgroup(x = covschmat, y = dat$Y,
                       trt = dat$trt, propensity.func = propensity.clustermean,
                       loss   = "sq_loss_lasso", penalty.factor= schpen_only,
                       nfolds = 5) # the modification that uses cluster-specific propensity scores and adds cluster dummies as moderator variables
    
    m2_aug_l <- fit.subgroup(x = covschmat, y = dat$Y,
                             trt = dat$trt, propensity.func = propensity.clustermean,
                             loss   = "sq_loss_lasso",  penalty.factor= schpen_only,
                             augment.func = aug.func_l, 
                             nfolds = 5) # the modification that uses cluster-specific propensity scores, adds cluster dummies as moderator variables, and adds augmentation terms with covariates only
    
    m2_aug_s <- fit.subgroup(x = covschmat, y = dat$Y,
                             trt = dat$trt, propensity.func = propensity.clustermean,
                             loss   = "sq_loss_lasso", penalty.factor= schpen_only,
                             augment.func = aug.func_s, 
                             nfolds = 5) # the modification that uses cluster-specific propensity scores, adds cluster dummies as moderator variables, and adds augmentation terms with covariates  and cluster dummies
    
    # :: get the benefit scores from weighting methods
    covmat_test <- as.matrix(dat_test[, c(paste0("X", c(3, 4, 5)))])
    covschmat_test <- model.matrix(~ X3 + X4 + X5 + id + 0, dat_test)
    covschmat_test <- covschmat_test[, !colnames(covschmat_test) %in% "id1"]
    covmat_demean_test <- apply(covmat_test, 2, function(x) x - ave(x, dat_test$id))
    
    m0_score <- predict(m0, covmat_test) # get the benefit scores
    
    m1_score <- predict(m1, covmat_test)
    m1_aug_l_score <- predict(m1_aug_l, covmat_test)
    m1_aug_s_score <- predict(m1_aug_s, covmat_test)
    
    m2_score <- predict(m2, covschmat_test)
    m2_aug_l_score <- predict(m2_aug_l, covschmat_test)
    m2_aug_s_score <- predict(m2_aug_s, covschmat_test)
    
    # :: save benefit scores in the simulated data
    temp_score <- cbind(delta=dat_test$delta, # true benefit scores
                        out1_score, out2_score, out3_score, out4_score, 
                        out5_score, out6_score, out7_score,
                        m0_score, 
                        m1_score, m1_aug_l_score, m1_aug_s_score,
                        m2_score, m2_aug_l_score, m2_aug_s_score)
    
    # :: save the results for the respective sample size conditions.
    if (sum(smpl.size == c(25, 25)) ==2) {
      rslt_1[[i]] <- temp_score
    } else if (sum(smpl.size == c(25, 150)) ==2) {
      rslt_2[[i]] <- temp_score
     } else if (sum(smpl.size == c(150, 25)) ==2) {
      rslt_3[[i]] <- temp_score
     } else if (sum(smpl.size == c(150, 150)) ==2) {
      rslt_4[[i]] <- temp_score
    }
    
  }
}

