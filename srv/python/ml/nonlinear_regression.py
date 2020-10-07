# -*- coding: utf-8 -*-
"""Nonlinear_Regression.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/19uJBMazJfLyPw9zKaBy6sOsh5tRXve07

# Training tree-based algorithms to predict epigenomic-metabolic interactions
**Author**: Scott Campit

## Mount Google Drive to Colab
This bit of code mounts Drive to the Colab notebook, and writes in an accessory function that is needed to read in Google Sheets as Pandas dataframes.
"""

# Import necessary data science libraries
import pandas as pd
import numpy as np

"""# Load datasets
Now let's load the GCP datasets we'll be computing ratios for. Right now, we'll compute the following ratios:
  * Cancer Cell Line Encyclopedia
  * LeRoy et al., 2012
"""

<<<<<<< HEAD
gcp_path = '~/Data/Proteomics/CCLE/CCLE Global Chromatin Profiles.xlsx'
met_path = '~/Data/Metabolomics/CCLE/CCLE_ALL_Ratios.csv'
=======
gcp_path = '/nfs/turbo/umms-csriram/scampit/Data/Proteomics/CCLE/CCLE Global Chromatin Profiles.xlsx'
met_path = '/nfs/turbo/umms-csriram/scampit/Data/Metabolomics/CCLE/CCLE_ALL_Ratios.csv'

>>>>>>> 6460e3458cffd5e531b83fa268d6c2a5083b1b59
GCP = pd.read_excel(gcp_path, 'All Ratios')
MET = pd.read_csv(met_path)

"""To preprocess the data, we'll do a couple of things, including:
  * Match by cell lines
  * Sort by index
  * Remove unncessary columns
  * Z-score the metabolomics data
"""

idx = list(set(GCP['index']) & set(MET['index']))
GCP = GCP[GCP['index'].isin(idx)]
MET = MET[MET['index'].isin(idx)]
GCP = GCP.drop_duplicates(subset='index', keep='first')

GCP = GCP.sort_values(by=['index'])
MET = MET.sort_values(by=['index'])

cell_lines = GCP['index'].values
gcpcol_to_drop = ['index'] 
metcol_to_drop = ['index', 'Unnamed: 0', 'Cell Lines']

GCP = GCP.drop(labels=gcpcol_to_drop, axis=1)
MET = MET.drop(labels=metcol_to_drop, axis=1)

"""Save the column names, which will be used later when constructing dataframes for evaluating model performance.

metabolites = list(MET.columns)
gcps = list(GCP.columns)
"""

from sklearn.preprocessing import quantile_transform
from sklearn.preprocessing import robust_scale
from scipy.stats import zscore
GCP_norm = GCP
MET_norm = zscore(MET, axis=1)

"""# Cancer cell line encyclopedia GCP -> Metabolism models
First, let's split the data into training and test sets.
"""

from sklearn.model_selection import train_test_split

"""Convert to Numpy array."""

GCP_norm = GCP_norm.to_numpy()

"""Split the data into validation and training data."""

# Split the CCLE data into a validation set
Xtrain, Xval, Ytrain, Yval = train_test_split(
    GCP_norm, MET_norm, test_size=0.3, random_state=0
)


"""## 3-fold cross validation for non-linear regressor selection
Now let's train a bunch of non-linear regressors and evaluate their performance.

We'll train the following ML models:
  * Random Forests
  * Gradient boosting
  * Gaussian Process Regression
  * XGBoost
"""

# ML models
from sklearn.ensemble import RandomForestRegressor
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.ensemble import ExtraTreesRegressor
import xgboost as xgb

# Accessory functions
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error
from skopt import BayesSearchCV
from skopt.space import Real, Categorical, Integer

# Suppress annoying warnings
import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
from joblib import Parallel, delayed
from skopt import dump, load
from skopt.utils import use_named_args

"""### Load some hyperparameters to sample

Let's define the search parameters we'll use for hyperparameter optimization
"""

# Gradient boosting
gb_params = {
    'max_depth': Integer(1, 5),
    'max_features': Categorical(['auto', 'sqrt', 'log2']),
    'min_samples_split': Integer(2, 100),
    'min_samples_leaf': Integer(1, 100),
    'learning_rate': Real(10**-5, 10**0, "log-uniform")
}

# Random Forests
rf_params = {
    'max_depth': Integer(1, 5),
    'max_features': Integer(1, Xtrain.shape[1]-1),
    'min_samples_split': Integer(2, 100),
    'min_samples_leaf': Integer(1, 100)
}

# Extra Trees
et_params = {
    'max_depth': Integer(1, 5),
    'max_features': Integer(1, Xtrain.shape[1]-1),
    'min_samples_split': Integer(2, 100),
    'min_samples_leaf': Integer(1, 100),
}

# XGBoost
xgb_params ={
    'gamma': Integer(1, 10),
    'learning_rate': Real(10**-5, 0.99, prior="log-uniform"),
    'max_depth': Integer(3, 10),
    'reg_alpha': Real(10**-5, 1, prior="log-uniform"),
    'reg_lambda':Real(10**-5, 1, prior="log-uniform"),
    'max_delta_step': Integer(0, 10),
}

"""### Random forest regression

The following line below performs Bayesian hyperparameter optimization using a random forest regressor. By default, 3-fold cross validation is used
"""

from sklearn.model_selection import KFold
<<<<<<< HEAD
from sklearn.model_selection import KFold
from time import sleep
import progressbar
=======
from time import sleep
import progressbar

bar = progressbar.ProgressBar(maxval=Ytrain.shape[1], \
    widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
>>>>>>> da8992e4440baf2110e8a01ef6ff283617765fe4

# Set the kfold operator to split 3 times with shuffle
kfold = KFold(n_splits=3, 
              shuffle=True, 
              random_state=0)

# Set the bayesian hyperparameter tuning to also include kfold cross validation
<<<<<<< HEAD
opt1 = BayesSearchCV(
    estimator=RandomForestRegressor(),
    search_spaces=rf_params,
    cv=kfold,
    n_iter=30,
    n_jobs=-1,
    random_state=0
)

bar.start()
=======
#opt1 = BayesSearchCV(
#    estimator=RandomForestRegressor(),
#    search_spaces=rf_params,
#    cv=kfold,
#    n_iter=30,
#    n_jobs=-1,
#    random_state=0
#)

#bar.start()
>>>>>>> da8992e4440baf2110e8a01ef6ff283617765fe4
# Construct univariate random forest models and append to mdls list
<<<<<<< HEAD
#mdls = []
#model_path='/home/scampit/Data/Models/GCP2Met/rf.pkl'
#print("Starting hyperparameter optimization and cross validation")
#for i in range(Ytrain.shape[1]):
#    _ = opt1.fit(Xtrain, Ytrain[:, i])
#    mdls.append(opt1)
#    dump(res=mdls, filename=model_path)
#    
#    bar.update(i+1)
#    sleep(0.1)
#print("Finished hyperparameter optimization and cross validation")
#bar.finish()
=======
mdls = []
model_path='/nfs/turbo/umms-csriram/scampit/Data/Models/GCP2Met/rf.pkl'
print("Starting random forest hyperparameter optimization and cross validation")
for i in range(Ytrain.shape[1]):
    _ = opt1.fit(Xtrain, Ytrain[:, i])
    mdls.append(opt1)

<<<<<<< HEAD
    model_path='/home/scampit/Data/Models/GCP2Met/rf.pkl'
    dump(res=mdls, filename=model_path)
    
    bar.update(i+1)
    sleep(0.1)
print("Finished random forest hyperparameter optimization and cross validation")
bar.finish()
=======
"""### Save RF GCP -> MET model 
Ensure model persistence by saving the serialized version of the model
"""

model_path='/nfs/turbo/umms-csriram/scampit/Data/Models/GCP2Met/rf.pkl'
dump(mdls, model_path)
>>>>>>> 6460e3458cffd5e531b83fa268d6c2a5083b1b59
>>>>>>> da8992e4440baf2110e8a01ef6ff283617765fe4

"""### Gradient boosting
Now let's train the remaining GCP -> MET models and save them as well.
"""

# Gradient Boosting
<<<<<<< HEAD
#opt2 = BayesSearchCV(
#    estimator=GradientBoostingRegressor(),
#    search_spaces=gb_params,
#    cv=kfold,
#    n_iter=30,
#    n_jobs=-1,
#    random_state=0
#)

#bar.start()
# Construct univariate gradient boosting models and append to mdls list
#mdls = []
#model_path='/home/scampit/Data/Models/GCP2Met/gb.pkl'
#print("Starting hyperparameter optimization and cross validation")
#for i in range(Ytrain.shape[1]):
#    _ = opt2.fit(Xtrain, Ytrain[:, i])
#    mdls.append(opt2)
#    dump(res=mdls, filename=model_path)
#    
#    bar.update(i+1)
#    sleep(0.1)
    
#print("Finished hyperparameter optimization and cross validation")
#bar.finish()
=======
opt2 = BayesSearchCV(
    estimator=GradientBoostingRegressor(),
    search_spaces=gb_params,
    cv=kfold,
    n_iter=30,
    n_jobs=-1,
    random_state=0
)

mdls = []
for i in range(Ytrain.shape[1]):
  _ = opt2.fit(Xtrain, Ytrain)
  mdls.append(opt2)
model_path='/nfs/turbo/umms-csriram/scampit/Data/Models/GCP2Met/gb.pkl'
dump(mdls, model_path)
>>>>>>> 6460e3458cffd5e531b83fa268d6c2a5083b1b59

"""### Extra Trees
Same for the extra trees -> train them and save them.
"""

# Extra Trees
<<<<<<< HEAD
#opt3 = BayesSearchCV(
#    estimator=ExtraTreesRegressor(),
#    search_spaces=et_params,
#    cv=kfold,
#    n_iter=30,
#    n_jobs=-1,
#    random_state=0
#)

#bar.start()
#mdls = []
#model_path='/home/scampit/Data/Models/GCP2Met/et.pkl'
#for i in range(Ytrain.shape[1]):
#    _ = opt3.fit(Xtrain, Ytrain[:, i])
#    mdls.append(opt3)
#    dump(res=mdls, filename=model_path)
#    
#    bar.update(i+1)
#    sleep(0.1)
    
#print("Finished hyperparameter optimization and cross validation")
#bar.finish()
=======
opt3 = BayesSearchCV(
    estimator=ExtraTreesRegressor(),
    search_spaces=et_params,
    cv=kfold,
    n_iter=30,
    n_jobs=-1,
    random_state=0
)

mdls = []
for i in range(Ytrain.shape[1]):
  _ = opt3.fit(Xtrain, Ytrain)
  mdls.append(opt3)
model_path='/nfs/turbo/umms-csriram/scampit/Data/Models/GCP2Met/et.pkl'
dump(mdls, model_path)
>>>>>>> 6460e3458cffd5e531b83fa268d6c2a5083b1b59

"""### XGBoost
Same for gradient boosting -> train them and save them.
"""

# Create object that will perform Bayesian hyperparameter tuning with 30 different iterations
<<<<<<< HEAD
opt4 = BayesSearchCV(
    estimator=xgb.XGBRegressor(),
    search_spaces=xgb_params,
    cv=kfold,
    n_iter=30,
    n_jobs=-1,
    random_state=0
)
=======
#opt4 = BayesSearchCV(
#    estimator=xgb.XGBRegressor(),
#    search_spaces=xgb_params,
#    cv=kfold,
#    n_iter=30,
#    n_jobs=-1,
#    random_state=0
#)
>>>>>>> da8992e4440baf2110e8a01ef6ff283617765fe4

# Create an object that will store all models
<<<<<<< HEAD
#bar.start()
#mdls = []
#model_path='/home/scampit/Data/Models/GCP2Met/xgb.pkl'
#for i in range(Ytrain.shape[1]):
#    _ = opt4.fit(Xtrain, Ytrain[:, i])
#    mdls.append(opt4)
#    dump(res=mdls, filename=model_path)
#    
#    bar.update(i+1)
#    sleep(0.1)
    
#print("Finished hyperparameter optimization and cross validation")
#bar.finish()
=======
mdls = []
for i in range(Ytrain.shape[1]):
  _ = opt4.fit(Xtrain, Ytrain[:, i])
  mdls.append(opt4)
model_path='/nfs/turbo/umms-csriram/scampit/Data/Models/GCP2Met/xgb.pkl'
dump(mdls, model_path)
>>>>>>> 6460e3458cffd5e531b83fa268d6c2a5083b1b59

"""# Cancer cell line encyclopedia Metabolism -> GCP models
Now we'll try to learn models that do the reverse problem: predicting chromatin profiles using metabolic data.

Now that I have some sense on how to write up K-fold cross validation with hyperparameter tuning, I'll create a function that does all of the bits above.
"""

# Split the CCLE data into a validation set
Xtrain, Xval, Ytrain, Yval = train_test_split(
    MET_norm, GCP_norm, test_size=0.3, random_state=0
)

models = [
          RandomForestRegressor(),
          GradientBoostingRegressor(),
          ExtraTreesRegressor(),
          xgb.XGBRegressor()
]

params = [
          gb_params,
          rf_params,
          et_params,
          xgb_params
]
names = [
<<<<<<< HEAD
         '/home/scampit/Data/Models/Met2GCP/rf.pkl',
         '/home/scampit/Data/Models/Met2GCP/gb.pkl',
         '/home/scampit/Data/Models/Met2GCP/et.pkl',
         '/home/scampit/Data/Models/Met2GCP/xgb.pkl'
=======
         '/nfs/turbo/umms-csriram/scampit/Data/Models/Met2GCP/rf.pkl',
         '/nfs/turbo/umms-csriram/scampit/Data/Models/Met2GCP/gb.pkl',
         '/nfs/turbo/umms-csriram/scampit/Data/Models/Met2GCP/et.pkl',
         '/nfs/turbo/umms-csriram/scampit/Data/Models/Met2GCP/xgb.pkl'
>>>>>>> 6460e3458cffd5e531b83fa268d6c2a5083b1b59
]

"""Here is the function that I have defined in order to do the following steps for each model:
  1. Define the BayesOpt object.
  2. Compute a model for each feature.
  3. Save all models based on a designated path.
"""

def train_models(models, params, Xtrain, Ytrain, kfold, filename):
  """
  train_models performs kfold bayesian hyperparameter tuning for different 
  models, and saves the output for model persistence.

  :param models: A single sklearn model object or list of sklearn model objects.
  :param params: A dictionary or list of dictionaries containing hyperparameters 
                 to tune.
  :param Xtrain: A numpy array or pandas dataframe containing the training data.
  :param Ytrain: A numpy array or pandas dataframe containing the output data.
  :param kfold:  An integer or sklearn object determining the kfold operation 
                 performed.
  :param filename: A string or list of paths to save the models (pickle).

  """
  print("Starting hyperparameter optimization and cross validation")
  for i in range(len(models)):
    model_path = filename[i]
    opt = BayesSearchCV(
                          estimator=models[i],
                          search_spaces=params[i],
                          cv=kfold,
		          n_iter=30,
		          n_jobs=-1,
		          random_state=0
    )
    bar.start()
    mdls =[]
    for j in range(Ytrain.shape[1]):
      _ = opt.fit(Xtrain, Ytrain[:, j])
      mdls.append(opt)
<<<<<<< HEAD
      dump(mdls, model_path)
=======
      dump(res=mdls, filename=model_path)
      
      bar.update(j+1)
      sleep(0.1)
    print("Finished hyperparameter optimization and cross validation")
    bar.finish()
>>>>>>> da8992e4440baf2110e8a01ef6ff283617765fe4

"""Finally, let's train the models."""

train_models(models, params, Xtrain, Ytrain, kfold, names)

"""# Evaluating the models

## GCP to MET models
Now that I have some models trained up, it's time to create some data structures that will have the metrics I want. First, let's grab the validatoin set again from the `train_test_split()` function. Because the seed is set, it should get me the same entries.
"""

# Split the CCLE data into a validation set
Xtrain, Xval, Ytrain, Yval = train_test_split(
    GCP_norm, MET_norm, test_size=0.3, random_state=0
)

"""Next, we'll load some libraries we'll be using to evaluate the predicted value against the true value."""

from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error

"""Let's now load the models we have trained to predict metabolite values from chromatin profiles."""

<<<<<<< HEAD
mdls = [load('/home/scampit/Data/Models/GCP2Met/rf.pkl'),
        load('/home/scampit/Data/Models/GCP2Met/gb.pkl'),
        load('/home/scampit/Data/Models/GCP2Met/et.pkl'),
        load('/home/scampit/Data/Models/GCP2Met/xgb.pkl')
=======
mdls = [load('/nfs/turbo/umms-csriram/scampit/Data/Models/GCP2Met/rf.pkl'),
        load('/nfs/turbo/umms-csriram/scampit/Data/Models/GCP2Met/gb.pkl'),
        load('/nfs/turbo/umms-csriram/scampit/Data/Models/GCP2Met/et.pkl'),
        load('/nfs/turbo/umms-csriram/scampit/Data/Models/GCP2Met/xgb.pkl')
>>>>>>> 6460e3458cffd5e531b83fa268d6c2a5083b1b59
]

"""The `evaluate_models()` function will compute evaluation metrics and spit out the final metrics of interest."""

def evaluate_models(models, Xval, Yval):
  """
  evaluate_models returns metrics from the model predictions, include the pearson
  correlation coefficient, coefficient of determination, MSE, and MAE.

  :param models:         A scikit-learn model object or list of model objects.
  :param Xval:           A numpy array or pandas dataframe containing 
                         validation set input data.
  :param Yval:           A numpy array or pandas dataframe containing 
                         validation set output data.
  :return final_metrics: A pandas dataframe or list of dfs containing the final 
                         evaluation metrics
  """

  final_metrics = []
  for j in range(len(models)):
    # Iterate through model objects
    m = models[j]

    r_values = list()
    p_values = list()
    mse_values = list()
    mae_values = list()

    # Iterate through features
    for i in range(len(m)):
      mdl = m[i]
      ypred = mdl.predict(Xval)
      r, pvalue = pearsonr(ypred, Yval[:, i])
      mse = mean_squared_error(ypred, Yval[:, i])
      mae = mean_absolute_error(ypred, Yval[:, i])

      r_values.append(r)
      p_values.append(pvalue)
      mae.append(mae)
      mse.append(mse)

    # Save the metrics in a dataframe
    pre_df = {
              "Pearson": r_values, 
              "Pvalue":  p_values,
              "MSE":     mse,
              "MAE":     mae
              }
    df = pd.DataFrame(pre_df)
    final_metrics.append(df)

    return final_metrics

"""Let's run the function on the list of regressors I have trained. Then we'll perform the following operations:

  1. Concatenate the results into a single dataframe
  2. Append the metabolite names to the list
  3. Sort the values in ascending alphabetical order by metabolite name
  4. Save the final results to the Google Sheet.
"""

final_metrics = evaluate_models(mdls, Xval, Yval)

# Flatten the array so that 
final_metrics = pd.concat(final_metrics, axis=1)
final_metrics["Metabolites"] = metabolites
final_metrics = final_metrics.sort_values(by=["Metabolites"], 
                                          axis=1, 
                                          ascending=True)


<<<<<<< HEAD
path = '/home/scampit/Data/Models/GCP2Met/gcp2met_metrics.csv'
=======
path = '/nfs/turbo/umms-csriram/scampit/Data/Models/GCP2Met/gcp2met_metrics.csv'
>>>>>>> 6460e3458cffd5e531b83fa268d6c2a5083b1b59
final_metrics.to_csv(final_metrics)

"""## MET to GCP models
Now let's do the reverse using the same operations described above.
"""

# Split the CCLE data into a validation set
Xtrain, Xval, Ytrain, Yval = train_test_split(
    MET_norm, GCP_norm, test_size=0.3, random_state=0
)

<<<<<<< HEAD
mdls = [load('/home/scampit/Data/Models/Met2GCP/rf.pkl'),
        load('/home/scampit/Data/Models/Met2GCP/gb.pkl'),
        load('/home/scampit/Data/Models/Met2GCP/et.pkl'),
        load('/home/scampit/Data/Models/Met2GCP/xgb.pkl')
=======
mdls = [load('/nfs/turbo/umms-csriram/scampit/Data/Models/Met2GCP/rf.pkl'),
        load('/nfs/turbo/umms-csriram/scampit/Data/Models/Met2GCP/gb.pkl'),
        load('/nfs/turbo/umms-csriram/scampit/Data/Models/Met2GCP/et.pkl'),
        load('/nfs/turbo/umms-csriram/scampit/Data/Models/Met2GCP/xgb.pkl')
>>>>>>> 6460e3458cffd5e531b83fa268d6c2a5083b1b59
]

final_metrics = evaluate_models(mdls, Xval, Yval)

# Flatten the array so that 
final_metrics = pd.concat(final_metrics, axis=1)
final_metrics["GCP"] = gcps
final_metrics = final_metrics.sort_values(by=["GCP"], 
                                          axis=1, 
                                          ascending=True)

<<<<<<< HEAD
path = '/home/scampit/Data/Models/Met2GCP/met2gcp_metrics.csv'
=======
path = '/nfs/turbo/umms-csriram/scampit/Data/Models/Met2GCP/met2gcp_metrics.csv'
>>>>>>> 6460e3458cffd5e531b83fa268d6c2a5083b1b59
final_metrics.to_csv(final_metrics)
