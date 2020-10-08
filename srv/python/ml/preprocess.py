"""

"""

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import MinMaxScaler

def data_imputation(df):
    """

    :param df:
    :return:
    """
    imp = SimpleImputer(missing_values=np.nan, strategy='median')
    #print(imp.fit_transform(df))
    return pd.DataFrame(imp.fit_transform(df.values),
                        index=df.index.values,
                        columns=df.columns.values)

def scale(df):
    """

    :param df:
    :return:
    """
    scaler = MinMaxScaler()
    return pd.DataFrame(scaler.fit_transform(df),
                        index=df.index.values,
                        columns=df.columns.values), \
           df.values.max(), df.values.min()