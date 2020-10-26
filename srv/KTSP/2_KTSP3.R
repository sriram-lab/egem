library(dplyr)
library(tidyr)
library(caret)
library(pracma)
require(switchBox)

indices <- sample(10,nrow(HM_Data), replace = TRUE, prob = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))
X <- HM_Data %>% select(3:227)
Y <- HM_Data[2]
acc <- rep(10,1)
prec <- rep(10,1)
rec <- rep(10,1)
F1 <- rep(10,1)
MCC <- rep(10,1)

for (k in 1:10)
{
  test <- indices == k 
  train <- indices != k
  Xtrain <- X[train,]
  Ytrain <- Y[train,]
  Xtest <- X[test,]
  Ytest <- Y[test,]
  
  # TAKING TRANSPOSE TO BRING DATA IN RIGHT FORMAT
  Xtrain <- t(Xtrain)
  Xtest <- t(Xtest)
  
  aa <- data.frame(Ytrain)
  aa$Target <- as.factor(aa$Target)
  Ytrain_fac <- aa$Target
  
  aa_test <- data.frame(Ytest)
  aa_test$Target <- as.factor(aa_test$Target)
  Ytest_fac <- aa_test$Target
  
  # Train a classifier using default filtering function based on the Wilcoxon test
  #classifier <- SWAP.Train.KTSP(xtrain, ytrain_fac, FilterFunc=NULL, featureNo = 30)
  classifier <- SWAP.Train.KTSP(Xtrain, Ytrain_fac, krange=c(3:15))
  
  testPrediction <- SWAP.KTSP.Classify(Xtest, classifier)
  # Show
  #table(testPrediction)
  
  # Resubstitution performance in the TEST set
  tab <- table(testPrediction, Ytest_fac)
  acc[k] <- sum(Ytest_fac == testPrediction)/length(Ytest_fac)
  prec[k] <- sum(Ytest_fac == 1 & testPrediction == 1)/sum(testPrediction == 1)
  rec[k] <- sum(Ytest_fac == 1 & testPrediction == 1)/sum(Ytest_fac == 1)
  F1[k] <- 2*rec[k]*prec[k]/(rec[k]+prec[k])
  
  # MCC Calculation
  TP <- tab[1]
  TN <- tab[4]
  FP <- tab[3]
  FN <- tab[2]
  N <- TN + TP + FN + FP
  S <- (TP + FN)/N
  P <- (TP + FP)/N
  MCC[k] <- ((TP/N) - (S*P))/(sqrt(P*S*(1-S)*(1-P)))
}

# Random guess accuracy
guess_acc <- rep(0,100)
for (j in 1:100)
{
  guess <- Ytest_fac[randperm(length(Ytest_fac))]
  guess_acc[j] =sum(Ytest_fac == guess)/length(Ytest_fac)
}

mean_acc <- mean(acc)
mean_prec <- mean(prec)
mean_rec <- mean(rec)
mean_F1 <- mean(F1)
mean_mcc <- mean(MCC)
mean_guess_acc <- mean(guess_acc)
p_val <- t.test(acc, guess_acc)$p.value
