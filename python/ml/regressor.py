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
import numpy as np
from numba import jit

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

ols_mdls = [];     ols_performance = []
ridge_mdls = [];   ridge_performance = []
lasso_mdls = [];   lasso_performance = []
elastic_mdls = []; elastic_performance = []
dt_mdls = [];      dt_performance = []
rf_mdls = [];      rf_performance = []
snn_mdls = [];     snn_performance = []

# Perform 10-fold cross validation on all regressors to get R, R2, MSE, and MAE
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

    # 1. Ordinary least squares
    ols = linear_model.LinearRegression()
    ols.fit(Xtrain2, Ytrain2)
    ols_pred = pd.DataFrame(ols.predict(Xtest2),
                            index=Ytest2.index,
                            columns=Ytest2.columns)
    ols_scores = evaluateMLmodel.regression_eval_metrics(ols_pred, Ytest2)
    ols_mdls.append(ols); ols_performance.append(ols_scores.mean(axis=0))

    # 2. Ridge Regression
    ridger = linear_model.Ridge(alpha=0.01)
    ridger.fit(Xtrain2, Ytrain2)
    ridger_pred = pd.DataFrame(ridger.predict(Xtest2),
                               index=Ytest2.index,
                               columns=Ytest2.columns)
    ridger_scores = evaluateMLmodel.regression_eval_metrics(ridger_pred, Ytest2)
    ridge_mdls.append(ridger); ridge_performance.append(ridger_scores.mean(axis=0))

    # 3. LASSO
    lassor = linear_model.Lasso(alpha=0.01)
    lassor.fit(Xtrain2, Ytrain2)
    lassor_pred = pd.DataFrame(lassor.predict(Xtest2),
                               index=Ytest2.index,
                               columns=Ytest2.columns)
    lassor_scores = evaluateMLmodel.regression_eval_metrics(lassor_pred, Ytest2)
    lasso_mdls.append(lassor); lasso_performance.append(lassor_scores.mean(axis=0))

    # 4. Elastic Net
    elasticr = linear_model.ElasticNet(alpha=0.50)
    elasticr.fit(Xtrain2, Ytrain2)
    elasticr_pred = pd.DataFrame(elasticr.predict(Xtest2),
                                 index=Ytest2.index,
                                 columns=Ytest2.columns)
    elasticr_scores = evaluateMLmodel.regression_eval_metrics(elasticr_pred, Ytest2)
    elastic_mdls.append(elasticr); elastic_performance.append(elasticr_scores.mean(axis=0))

    # 5. Decision Tree
    dtr = DecisionTreeRegressor(criterion='mse')
    dtr.fit(Xtrain2, Ytrain2)
    dtr_pred = pd.DataFrame(dtr.predict(Xtest2),
                            index=Ytest2.index,
                            columns=Ytest2.columns)
    dtr_scores = evaluateMLmodel.regression_eval_metrics(dtr_pred, Ytest2)
    dt_mdls.append(dtr); dt_performance.append(dtr_scores.mean(axis=0))

    # 6. Random Forest
    rfr = RandomForestRegressor(criterion='mse',
                                max_features='sqrt',
                                oob_score=True)
    rfr.fit(Xtrain2, Ytrain2)
    rfr_pred = pd.DataFrame(rfr.predict(Xtest2),
                            index=Ytest2.index,
                            columns=Ytest2.columns)
    rfr_scores = evaluateMLmodel.regression_eval_metrics(rfr_pred, Ytest2)
    rf_mdls.append(rfr); rf_performance.append(rfr_scores.mean(axis=0))

    # 7. AdaBoosting
    #adar = AdaBoostRegressor(loss='exponential')
    #adar.fit(Xtrain, Ytrain)
    #adar_pred = adar.predict(Xtest)
    #print(adar_pred)

    # 8. Shallow neural networks
    snnr = MLPRegressor(hidden_layer_sizes=(100),
                         activation='relu',
                         solver='adam',
                         alpha=0.001,
                         learning_rate='adaptive')
    snnr.fit(Xtrain2, Ytrain2)
    snnr_pred = pd.DataFrame(snnr.predict(Xtest2),
                             index=Ytest2.index,
                             columns=Ytest2.columns)
    snnr_scores = evaluateMLmodel.regression_eval_metrics(snnr_pred, Ytest2)
    snn_mdls.append(snnr); snn_performance.append(snnr_scores.mean(axis=0))

# Get best model and get final metrics from the validation set
all_mdls = [ols_mdls, ridge_mdls, lasso_mdls, elastic_mdls, dt_mdls, rf_mdls, snn_mdls]
all_perf = [ols_performance, ridge_performance, lasso_performance, elastic_performance,
            dt_performance, rf_performance, snn_performance]
final_names = ['OLS', 'Ridge', 'LASSO', 'ElasticNet', 'DecisionTree', 'RandomForest', 'SNN']

# Scale the test data
Xtest, _, _ = preprocess.scale(Xtest)
Ytest, _, _ = preprocess.scale(Ytest)

best_mdl_metrics = []
writer = pd.ExcelWriter('predictMetabolismFromGCP.xlsx', engine='xlsxwriter')
for i in range(0, len(all_mdls)):
    best_mdl_metrics.append(evaluateMLmodel.get_best_model_metrics(all_mdls[i], all_perf[i], Xtest, Ytest))
    best_mdl_metrics[i].to_excel(writer, sheet_name=final_names[i])
writer.save()