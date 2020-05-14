"""

"""

import numpy as np
import pandas as pd
import math

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
    return regression_eval_metrics(ypred, Ytest)
