"""
"""
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import train_test_split

import pandas as pdb
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

file = ""
df = pd.DataFrame(file, index='Reactions')
FeatureNames = df.columns.tolist()
RxnNames = df.index.tolist()

data = df.values[:, [1, 2]]
predictors = df.values[:, 3]


metTrain, metTest, histTrain, histTest = train_test_split(
    data, predictors, test_size=0.2)

lm = LinearRegression()
lm.fit(metTrain, histTrain)

predictions = lm.predict(metTest)
plt.scatter(histTest, predictions)

holdOutAcc = lm.score(metTest, histTest)
