"""
Regressor

This script trains regression-based statistical learning models on histone post-translational modifications to
predict metabolomics data.

This module trains several regression-based models. First, hyperparameter tuning is performed with k-fold cross
validation. Both the best k-fold and hyperparameters are to be selected for getting the best model with the highest
adjusted coefficient of determination.

TODO:
  * Incorporate wild-type metabolic fluxes as the target variable
  * Incorporate different datasets other than the CCLE
  * Remove histone markers with a lot of NaNs and remove specific cell lines with few NaNs in a given histone marker
  * Iterate through each histone marker individually to compute the coefficients
  * Use a subset of the cancer cell lines (ie TCGA-110, NCI-60)

# EXPECTATIONS
  * From EDA, there are highly correlated histone markers. Thus, I would expect linear models to not perform as well.

@author: Scott Campit
"""

# Essential data object manipulations
import pandas as pd
import xlsxwriter
import numpy as np

# Regression Models
from sklearn.model_selection import RandomizedSearchCV
from sklearn.utils import shuffle
from sklearn import linear_model
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import AdaBoostRegressor
from sklearn.neural_network import MLPRegressor
#import torch
#from torch.autogrid import Variable
#import torch.nn as nn
#import torch.nn.functional as F

# Model evaluation
from sklearn.model_selection import train_test_split
from sklearn.model_selection import KFold

# Custom scripts
import preprocess
import evaluateMLmodel

# Load CCLE data
print("LOAD DATA")
db = '/home/scampit/Data/Expression/Combined/mapped_ccle_data.xlsx'
metabolomics = pd.read_excel(db,
                             sheet_name='CCLE_Metabolomics',
                             index_col='CCL')
gcp = pd.read_excel(db,
                    sheet_name='CCLE_GCP',
                    index_col='CCL')

# Create validation set
Xtrain, Xtest, Ytrain, Ytest = train_test_split(gcp, metabolomics, test_size=0.3)

# Set up variables to store results from 10-fold Cross Validation
kcv = KFold(n_splits=10, shuffle=True)

regressors = [linear_model.LinearRegression(),
              linear_model.Ridge(),
              linear_model.Lasso(),
              linear_model.ElasticNet(),
              DecisionTreeRegressor(),
              RandomForestRegressor(oob_score=True),
              MLPRegressor(early_stopping=True, max_iter=1000)
              ]
final_names = ["LR", "Ridge", "LASSO", "ElasticNet",
               "DT", "RF", "MLP"]
# Hyperparameters to sample
reg_param = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1, 1, 10]
learning_rate = [1E-6, 1E-5, 1E-4, 1E-3, 1E-2, 1E-1]
l1rat = [0.1, 0.25, 0.50, 0.75, 0.90, 1]
lr_schedule = ['constant', 'optimal', 'invscaling', 'adaptive']
max_feats = ['auto', 'sqrt', 'log2']
cp = [1E-3, 0.005, 1E-2, 0.05, 1E-1, 0.5, 1]
lossFnc = ['mse', 'mae']
actFnc = ['identity', 'logistic', 'tanh', 'relu']
solve = ['lbfgs', 'sgd', 'adam']
hidden_layer_size = list(range(5, 1000))
param_list = [dict(),
              dict(alpha=reg_param),
              dict(alpha=reg_param),
              dict(alpha=reg_param, l1_ratio=l1rat),
              dict(criterion=lossFnc, max_features=max_feats, ccp_alpha=cp),
              dict(criterion=lossFnc, max_features=max_feats, ccp_alpha=cp),
              dict(hidden_layer_sizes=(hidden_layer_size), activation=actFnc,
                   alpha=reg_param, learning_rate=lr_schedule,
                   learning_rate_init=learning_rate,
                   solver=solve)]

# Perform 10-fold cross validation on all regressors to get R, R2, MSE, and MAE
allMdls = []
print("Random Grid Search")
p = 0
for mdl in regressors:
    bestMdl = RandomizedSearchCV(estimator=mdl,
                                 param_distributions=param_list[p],
                                 n_iter=50)
    p += 1
    allMdls.append(bestMdl)

print("CROSS VALIDATION")
bestMdls = []
mdlPerformance = []
for mdl in allMdls:
    tmpMdl = []
    tmpPerformance = []
    for train_idx, test_idx in kcv.split(Xtrain):
        X, Y = shuffle(Xtrain, Ytrain)

        # Create cross validation indices and data
        Xtrain2, Xtest2 = X.iloc[train_idx, :], X.iloc[test_idx, :]
        Ytrain2, Ytest2 = Y.iloc[train_idx, :], Y.iloc[test_idx, :]

        # Scale the data using min-max scaling
        Xtrain2, xmax, xmin = preprocess.scale(Xtrain2)
        Xtest2 = (Xtest2 - xmin) / (xmax - xmin)
        Ytrain2, ymax, ymin = preprocess.scale(Ytrain2)
        Ytest2 = (Ytest2 - ymin) / (ymax - ymin)

        trainedMdl, mdlScore = evaluateMLmodel.fitMLModel(bestMdl,
                                                          Xtrain2, Ytrain2,
                                                          Xtest2, Ytest2)
        tmpMdl.append(trainedMdl)
        tmpPerformance.append(mdlScore)
    mdl, mdlPerf = evaluateMLmodel.get_best_model_metrics(tmpMdl, tmpPerformance, Xtest2, Ytest2)
    bestMdls.append(mdl)
    mdlPerformance.append(mdlPerf)
print("FINISHED CROSS VALIDATION")

# Scale the test data
Xtest, _, _ = preprocess.scale(Xtest)
Ytest, _, _ = preprocess.scale(Ytest)

best_mdl_metrics = []
writer = pd.ExcelWriter('predictMetabolismFromGCP_paramOpt.xlsx', engine='xlsxwriter')
for i in range(0, len(bestMdls)):
    _, mdlPerf = evaluateMLmodel.get_best_model_metrics(tmpMdl[i], tmpPerformance[i], Xtest, Ytest)
    best_mdl_metrics.append(mdlPerf)
    best_mdl_metrics[i].to_excel(writer, sheet_name=final_names[i])
writer.save()