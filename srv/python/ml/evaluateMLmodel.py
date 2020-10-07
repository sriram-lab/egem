"""
evaluateMLmodel

This script evaluates machine learning models that were trained using the sklearn framework.

TODO:
  * Set up script to also evaluate classifiers as well as regressors.

  
"""

import numpy as np
import pandas as pd
import math
from numba import jit

from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error

def regression_eval_metrics(ypred, validation_set):
    """

    :return:
    """
    r = []
    r2 = []
    mse = []
    mae = []

    m = 0
    for i in validation_set:
        r.append(validation_set[i].corr(ypred[i], method="pearson"))
        r2.append(r[m]**2)
        mse.append(mean_squared_error(validation_set[i], ypred[i], squared=True))
        mae.append(mean_absolute_error(validation_set[i], ypred[i]))
        m += 1

    r = pd.DataFrame(r, index=validation_set.columns, columns=['R'])
    r2 = pd.DataFrame(r2, index=validation_set.columns, columns=['R2'])
    mse = pd.DataFrame(mse, index=validation_set.columns, columns=['MSE'])
    mae = pd.DataFrame(mae, index=validation_set.columns, columns=['MAE'])
    return pd.concat([r, r2, mse, mae], axis=1)

def get_best_model_metrics(mdlStrct, performance_df, Xtest, Ytest):
    """

    :param mdlStrct:
    :param performance_df:
    :param Xtest:
    :param Ytest:
    :return:
    """
    df = pd.concat(performance_df, axis=1)
    idx = df.iloc[1, :].idxmax(axis=1)
    final_model = mdlStrct[idx]
    ypred = pd.DataFrame(final_model.predict(Xtest),
                             index=Ytest.index,
                             columns=Ytest.columns)
    return final_model, regression_eval_metrics(ypred, Ytest)

def fitMLModel(mdlObj, Xtrain, Ytrain, Xtest, Ytest):
    """

    :param mdlObj:
    :return:
    """
    mdlObj.fit(Xtrain, Ytrain)
    mdl_pred = pd.DataFrame(mdlObj.predict(Xtest),
                            index=Ytest.index,
                            columns=Ytest.columns)
    mdlScore = regression_eval_metrics(mdl_pred, Ytest)
    return mdlObj, mdlScore.mean(axis=0)