###########################################################################
# Designing Optimal, Data-Driven Policies from Multisite Randomized Trials
# by Youmi Suk and Chan Park

library(readstata13)
library(data.table)
library(lme4)
library(tidyverse)
library(personalized)
source("https://github.com/youmisuk/multisiteOTR/blob/main/Simulation/Qlearn.R")

# :: download data from https://doi.org/10.3886/E113783V1.

# :: load data
dat <- read.dta13("[...changeDirectoryName...]/CCTproject_shortterm/AEJApp_2010-0132_Data/Public_Data_AEJApp_2010-0132.dta", generate.factors=TRUE, nonint.factors = TRUE)
dat$s_sexo <- ifelse(dat$s_sexo == "0", 0, 1) # change the sex variable into dummy variable.

# :: subset data in San Cristobal
dat_San <- dat %>% filter(suba == 0 & (control == 1 | treatment == 1)) 
dat_San$trt <- dat_San$treatment

# :: evaluate whether there are any duplicated cases or not. 
# If duplicated cases, delete those cases because they are thought to come from the same household. 
# Otherwise, cases are considered as if coming from a separate household.
dat_San$caseID <- 10000:(9999+nrow(dat_San)) # give individual identifiers.
dat_San_NA_dup <- data.table(dat_San[is.na(dat_San$fu_nim_hogar),  c("s_durables", "s_edadhead", "s_estcivil", "s_estrato", "s_infraest_hh", "s_ingtotal", "s_num18", "s_puntaje", "s_single", "s_tpersona", "s_teneviv", "s_utilities", "s_yrshead")], key= "s_durables,s_edadhead,s_estcivil,s_estrato,s_infraest_hh,s_ingtotal,s_num18,s_puntaje,s_single,s_tpersona,s_teneviv,s_utilities,s_yrshead")
dupRowsIndex <- dat_San_NA_dup[unique(dat_San_NA_dup[duplicated(dat_San_NA_dup)]), which=T]
dat_San_NA <- dat_San[is.na(dat_San$fu_nim_hogar),]
dupCaseID <- dat_San_NA[dupRowsIndex, "caseID"]

MoreOneS_houseID <- as.numeric(names(table(dat_San$fu_nim_hogar))[table(dat_San$fu_nim_hogar) > 1])

# :: create data where (i) case ID is not duplicated ID, 
#                      (ii) households do not have two children who registered,
#                      (iii) schools have more than one student who registered,
#                      (iv) schools have survey data about attendance.
df <- dat_San[!dat_San$caseID %in% dupCaseID &  !(dat_San$fu_nim_hogar %in% MoreOneS_houseID) & 
                dat_San$school_code %in% names(table(dat_San$school_code))[table(dat_San$school_code) > 1] &
                dat_San$survey_selected==1, ]

# :: estimate the treatment proportions: (i) overall and (ii) cluster-specific
trt_prop_overall <- with(df, prop.table(table(trt)))
trt_prop_sch <- with(df, prop.table(table(school_code, trt), 1))
trt_prop_sch
hist(trt_prop_sch[,2])

# :: re-define factor variables
df$school_code <- factor(df$school_code)
df$inc_380 <- as.numeric(df$s_ingtotal <= 380)
df$s_estcivil <- ifelse(df$s_estcivil =="Single", 1, 0)
df$s_estrato <- factor(df$s_estrato)
df$grade <- factor(df$grade)


# :: create functions for R pakcage "personalized"

propensity.const <- function(x, trt)
{
  rep(prop.table(table(trt))[2], length(trt))
}

propensity.clustermean <- function(x, trt)
{
  trt_prop_sch[,2][df$school_code]
}


# create a covariate matrix for augmentation
covmat_aug <- model.matrix(~  s_teneviv + s_utilities + s_durables + s_infraest_hh + s_age_sorteo  + s_years_back + s_sexo + 
                             s_estcivil + s_single + s_edadhead + s_yrs + s_yrshead + s_tpersona + s_num18 + s_estrato + inc_380 +
                             s_puntaje + s_ingtotal + grade +  s_over_age, df)[, -1]

aug.func_l <- function(x, y) {
  aug <- lm(y ~ covmat_aug)
  return(predict(aug))
}

aug.func_s <- function(x, y) {
  aug <- lm(y ~ covmat_aug + df$school_code)
  return(predict(aug))
}

# ::::::::::::::::::::
# :::: Q-Learning :::: ####
# ::::::::::::::::::::

# :: create a matrix that contains covariates working as moderators 
covmat <- model.matrix(~  s_teneviv + s_years_back + grade + s_estrato + inc_380, df)[, -1]

# :: create matrices that contains school dummies
sch <- model.matrix( ~ school_code, df)[, -1]

covmat_sch <- cbind(covmat, sch) # a matrix that contains covariates and school dummies that potentially work as moderators 
covmat_sch_aug <- model.matrix(~  s_teneviv + s_utilities + s_durables + s_infraest_hh + s_age_sorteo  + s_years_back + s_sexo + 
                                 s_estcivil + s_single + s_edadhead + s_yrs + s_yrshead + s_tpersona + s_num18 + s_estrato + inc_380 +
                                 s_puntaje + s_ingtotal + grade +  s_over_age + school_code, df)[, -1] # a matrix that contains both all the covariates and school dummies

# :: replace " " with "." in column names
colnames(covmat) <- gsub(" ", ".", colnames(covmat))
colnames(covmat_sch) <- gsub(" ", ".", colnames(covmat_sch))
colnames(covmat_aug) <- gsub(" ", ".", colnames(covmat_aug))
colnames(covmat_sch_aug) <- gsub(" ", ".", colnames(covmat_sch_aug))

# :: create data for Q-learning
out_dat <- data.frame(at_msamean=df$at_msamean, trt = df$trt, covmat_sch_aug, school_code=df$school_code)

# :: implement Q-learning methods
out1 <- Qlearn(Yname="at_msamean", Trtname="trt", Xname=colnames(covmat_aug), interTrt=formula(paste("~", paste(colnames(covmat), collapse=" + "))), IDname="school_code", 
               traindata=out_dat, testdata=out_dat, method = c("Q-lm"))

out2 <- Qlearn(Yname="at_msamean", Trtname="trt", Xname=colnames(covmat_sch_aug), interTrt=formula(paste("~", paste(colnames(covmat), collapse=" + "))), IDname="school_code", 
               traindata=out_dat, testdata=out_dat, method = c("Q-fe"))

out3 <- Qlearn(Yname="at_msamean", Trtname="trt", Xname=colnames(covmat_aug), interTrt=formula(paste("~", paste(colnames(covmat), collapse=" + "))), IDname="school_code", 
               traindata=out_dat, testdata=out_dat, method = c("Q-ri"))

out4 <- Qlearn(Yname="at_msamean", Trtname="trt", Xname=colnames(covmat_aug), interTrt=formula(paste("~", paste(colnames(covmat), collapse=" + "))), IDname="school_code", 
               traindata=out_dat, testdata=out_dat, method = c("Q-rs"))

out5 <- Qlearn(Yname="at_msamean", Trtname="trt", Xname=colnames(covmat_aug), interTrt=formula(paste("~", paste(colnames(covmat), collapse=" + "))), IDname="school_code", 
               traindata=out_dat, testdata=out_dat, method = c("Q-dm"))

out6 <- Qlearn(Yname="at_msamean", Trtname="trt", Xname=colnames(covmat_aug), interTrt=formula(paste("~", paste(colnames(covmat), collapse=" + "))), IDname="school_code", 
               traindata=out_dat, testdata=out_dat, method = c("Q-dm-ri"))

out7 <- Qlearn(Yname="at_msamean", Trtname="trt", Xname=colnames(covmat_aug), interTrt=formula(paste("~", paste(colnames(covmat), collapse=" + "))), IDname="school_code", 
               traindata=out_dat, testdata=out_dat, method = c("Q-dm-rs"))

# save benefit scores from Q-learning methods
out1_score <- out1$benefit.scores
out2_score <- out2$benefit.scores
out3_score <- out3$benefit.scores
out4_score <- out4$benefit.scores
out5_score <- out5$benefit.scores
out6_score <- out6$benefit.scores
out7_score <- out7$benefit.scores


# :::::::::::::::::::::::::::
# :::: Weighting methods :::: ####
# :::::::::::::::::::::::::::

# create a covariate matrix that is demeaned
covmat_demean <- apply(covmat, 2, function(x) x - ave(x, df$school_code))

set.seed(123)
m0 <- fit.subgroup(x = covmat, y = df$at_msamean,
                   trt = df$trt, propensity.func = propensity.const,
                   loss   = "sq_loss_lasso", penalty.factor= c(0, rep(1, 3), rep(0, 9)),
                   nfolds = 5) # the baseline weighting method that uses the constant propensity score

set.seed(123)
m1 <- fit.subgroup(x = covmat, y = df$at_msamean,
                   trt = df$trt, propensity.func = propensity.clustermean,
                   loss   = "sq_loss_lasso", penalty.factor= c(0, rep(1, 3), rep(0, 9)),
                   nfolds = 5) # the modification that uses cluster-specific propensity scores

set.seed(123)
m1_aug_l <- fit.subgroup(x = covmat, y = df$at_msamean,
                         trt = df$trt, propensity.func = propensity.clustermean,
                         loss   = "sq_loss_lasso", penalty.factor= c(0, rep(1, 3), rep(0, 9)),
                         nfolds = 5,
                         augment.func = aug.func_l) # the modification that uses cluster-specific propensity scores and augmentation terms with covariates only

set.seed(123)
m1_aug_s <- fit.subgroup(x = covmat, y = df$at_msamean,
                         trt = df$trt, propensity.func = propensity.clustermean,
                         loss   = "sq_loss_lasso", penalty.factor= c(0, rep(1, 3), rep(0, 9)),
                         nfolds = 5,
                         augment.func = aug.func_s) # the modification that uses cluster-specific propensity scores and augmentation terms with covariates and cluster dummies

set.seed(123)
m2 <- fit.subgroup(x = covmat_sch, y = df$at_msamean,
                   trt = df$trt, propensity.func = propensity.clustermean,
                   loss   = "sq_loss_lasso", penalty.factor= c(0, rep(1, 3), rep(0, 9), rep(1, ncol(sch))),
                   nfolds = 5) # the modification that uses cluster-specific propensity scores and adds cluster dummies as moderator variables

set.seed(123)
m2_aug_l <- fit.subgroup(x = covmat_sch, y = df$at_msamean,
                         trt = df$trt, propensity.func = propensity.clustermean,
                         loss   = "sq_loss_lasso", penalty.factor= c(0, rep(1, 3), rep(0, 9), rep(1, ncol(sch))),
                         augment.func = aug.func_l,
                         nfolds = 5) # the modification that uses cluster-specific propensity scores, adds cluster dummies as moderator variables, and adds augmentation terms with covariates only

set.seed(123)
m2_aug_s <- fit.subgroup(x = covmat_sch, y = df$at_msamean,
                         trt = df$trt, propensity.func = propensity.clustermean,
                         loss   = "sq_loss_lasso", penalty.factor= c(0, rep(1, 3), rep(0, 9), rep(1, ncol(sch))),
                         augment.func = aug.func_s,
                         nfolds = 5) # the modification that uses cluster-specific propensity scores, adds cluster dummies as moderator variables, and adds augmentation terms with covariates  and cluster dummies

# :: save benefit scores from weighting methods
m0_score <- predict(m0, covmat)

m1_score <- predict(m1, covmat)
m1_aug_l_score <- predict(m1_aug_l, covmat)
m1_aug_s_score <- predict(m1_aug_s, covmat)

m2_score <- predict(m2, covmat_sch)
m2_aug_l_score <- predict(m2_aug_l, covmat_sch)
m2_aug_s_score <- predict(m2_aug_s, covmat_sch)

# :::::::::::::::::
# :::: Summary :::: ####
# :::::::::::::::::

# save all the beneift scores
att_score <- data.frame(out1_score, out2_score, out3_score, out4_score, out5_score, out6_score, out7_score,
                        m0_score, m1_score, m1_aug_l_score, m1_aug_s_score, m2_score, m2_aug_l_score, m2_aug_s_score)

colnames(att_score) <-index <- c("Q-base", "Q-fe", "Q-ri", "Q-rs", "Q-dm", "Q-dm-ri", "Q-dm-rs",
                                 "W-base", "W-noaug", "W-aug", "W-augID", "W[S]-noaug", "W[S]-aug", "W[S]-augID")

att_score2 <- att_score[, c("Q-base", "Q-fe", "Q-ri", "Q-rs", "Q-dm", "Q-dm-ri", "Q-dm-rs",
                            "W-base", "W-noaug", "W[S]-noaug", "W-aug", "W[S]-aug", "W-augID", "W[S]-augID")]

att_score2_opt <- att_score2 > 0

apply(att_score2_opt, 2, mean) # the proportion of the recommended group in each estimator
