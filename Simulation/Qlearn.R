###########################################################################
# Designing Optimal, Data-Driven Policies from Multisite Randomized Trials
# by Youmi Suk and Chan Park

# The Qlearn function is a tool for implementing our proposed modifications for Q-learning as well as the baseline Q-learning method. 
Qlearn <- function(Yname, Trtname, Xname, interTrt=formula(~1), IDname, traindata, testdata=NULL,
                   method = c("Q-lm", "Q-fe", "Q-ri", "Q-rs", "Q-dm", "Q-dm-ri", "Q-dm-rs")) {

  # The function takes several parameters:
  # Yname: The name of the outcome variable.
  # Trtname: The name of the treatment variable.
  # Xname: The names of the covariates.
  # interTrt: A formula for the moderators.
  # IDname: The name of the cluster identifiers.
  # traindata: Training data.
  # testdata: Test data. If not provided, the training data will be used as test data.
  # method: The Q-learning method you want to use. You can choose from "Q-lm", "Q-fe", "Q-ri", "Q-rs", "Q-dm", "Q-dm-ri", "Q-dm-rs". 
  #         Note that "Q-lm" is the baseline Q-learning estimator.
  
  if(is.null(testdata)) {
    testdata <- traindata
  }
  
  # variables
  IDfactor = as.factor(traindata[, IDname])
  IDnum = as.numeric(IDfactor)

  Xtraindata <- traindata[, Xname]
  Xtestdata <- testdata[, Xname]
  
  lvl2var <- colSums(apply(as.matrix(Xtraindata), 2, function(x) tapply(x, IDnum, function(x) length(unique(x)))) == 1) == length(levels(IDfactor))
  X_lvl1_train <- as.matrix(Xtraindata)[, !lvl2var]; X_lvl1_test <- as.matrix(Xtestdata)[, !lvl2var]
  W_lvl2_train <- as.matrix(Xtraindata)[, lvl2var]; W_lvl2_test <- as.matrix(Xtestdata)[, lvl2var]

  qdat = data.frame(Y=traindata[, Yname], trt=traindata[, Trtname], X_lvl1_train, IDfactor=IDfactor)
  qdat_test = data.frame(Y=testdata[, Yname], trt=testdata[, Trtname], X_lvl1_test, IDfactor=as.factor(testdata[, IDname]))
  
  interTrt.mat <- model.matrix(interTrt, data=testdata)
  interTrt.name <- paste0(colnames(interTrt.mat)[!colnames(interTrt.mat) %in% "(Intercept)"], ":trt")
  
  if (method == "Q-lm") {
    
    lmFrml <- formula(paste("Y ~", paste(c(colnames(X_lvl1_train), "trt", interTrt.name), collapse="+")))
    out <- lm(lmFrml, data=qdat) 
    out_score <- as.numeric(interTrt.mat%*%matrix(coef(out)[c("trt", interTrt.name)]))
    
  }
  
  if (method == "Q-fe") {
    
    IDFrml <- formula(paste("Y ~", paste(c(colnames(X_lvl1_train), "trt", interTrt.name, "IDfactor"), collapse="+")))
    out <- lm(IDFrml, data=qdat) 
    out_score <- as.numeric(interTrt.mat%*%matrix(coef(out)[c("trt", interTrt.name)]))
    
  }
  
  if (method == "Q-ri") {
    
    rlFrml <- formula(paste("Y ~", paste(c(colnames(X_lvl1_train), "trt", interTrt.name), collapse="+") , " + (1 | IDfactor)"))
    out <- lmer(rlFrml, data=qdat) 
    out_score <- as.numeric(interTrt.mat %*% t(as.matrix(coef(out)$IDfactor[1, c("trt", interTrt.name)])))
    
  }
  
  if (method == "Q-rs") {
    
    rsFrml <- formula(paste("Y ~", paste(c(colnames(X_lvl1_train), "trt", interTrt.name), collapse="+") , " + (trt | IDfactor)"))
    out <- lmer(rsFrml, data=qdat)
    out_score <- as.numeric(matrix(coef(out)$IDfactor[, "trt"][qdat_test$IDfactor]) + interTrt.mat[, -1]%*%t(as.matrix(coef(out)$IDfactor[1, interTrt.name])))
    
  }

  if (method %in% c("Q-dm", "Q-dm-ri", "Q-dm-rs")) {
    
    dmvarFrml <- formula(paste("~", paste(c("Y", colnames(X_lvl1_train), "trt", interTrt.name, "0"), collapse="+")))
    m <- model.matrix(dmvarFrml, data=qdat)
    m_test <- model.matrix(dmvarFrml, data=qdat_test)
    
    dmvar_mat <- data.frame(apply(m, 2, function(x) x - ave(x, qdat$IDfactor)))
    dmvar_mat_test <- data.frame(apply(m_test, 2, function(x) x - ave(x, qdat_test$IDfactor)))
    dmvar_mat$IDfactor <- qdat$id; dmvar_mat_test$IDfactor <- qdat_test$id
    
    interTrt.name.dot <- paste0(colnames(interTrt.mat)[!colnames(interTrt.mat) %in% "(Intercept)"], ".trt")
    
    if (method == "Q-dm") {
      
      dmFrml <- formula(paste("Y ~", paste(c(colnames(X_lvl1_train), "trt", interTrt.name.dot), collapse="+")))
      out <- lm(dmFrml, data=dmvar_mat)
      out_score <- as.numeric(interTrt.mat%*%matrix(coef(out)[c("trt", interTrt.name.dot)]))
      
    } else if (method == "Q-dm-ri") {
      
      dm_rlFrml <- formula(paste("qdat$Y ~", paste(c(colnames(X_lvl1_train), "trt", interTrt.name.dot), collapse="+") , " + (1 | IDfactor)"))
      out <- lmer(dm_rlFrml, data=dmvar_mat)
      out_score <- as.numeric(interTrt.mat %*% t(as.matrix(coef(out)$IDfactor[1, c("trt", interTrt.name.dot)])))
      
    } else {
      
      dm_rsFrml <- formula(paste("qdat$Y ~", paste(c(colnames(X_lvl1_train), "trt", interTrt.name.dot), collapse="+") , " + (trt | IDfactor)"))
      out <- lmer(dm_rsFrml, data=dmvar_mat)
      out_score <- as.numeric(matrix(coef(out)$IDfactor[, "trt"][qdat_test$IDfactor]) + interTrt.mat[, -1]%*%t(as.matrix(coef(out)$IDfactor[1, interTrt.name.dot])))
    }

  }
  
  return(list(summary=summary(out), coef=summary(out)$coef, benefit.scores=out_score,  recommended.trtrts=as.numeric(I(out_score>0))))
  
}
