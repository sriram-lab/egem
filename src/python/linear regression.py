#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
from scipy.stats import pearsonr
from sklearn.metrics import mean_squared_error
from sklearn.metrics import mean_absolute_error
from sklearn.model_selection import KFold
from scipy.stats import median_absolute_deviation as mad
from scipy.stats import norm
### Fit train sets in Lasso with best alpha


# ### Merge datasets based on cancer cell line name;
# ### concatenate the two dataframes and match based on cancer cell lines; 
# ### save them into two new dataframes- MET and GCP; 
# ### split GCP and MET into train and test set respectively.

# In[27]:


md = pd.read_excel(io='CCLE metabolomics dataset.xlsx',sheet_name="All")
mt = md.drop(['Tissue', 'Medium','Culture'], axis=1)
hm = pd.read_csv('GCP_proteomics_remapped.csv')
merge_tb = mt.merge(hm,how='inner',left_on='CCL', right_on='Cell Line')
MET = merge_tb.iloc[:,1:226]
GCP = merge_tb.iloc[:,227:269]


# ### Plot the H3K4me1 and K3K9ac1 z-score from originial histone marker data set.

# In[31]:


ax1 = hm.plot.scatter(x='H3K4me1',
                      y='Cell Line',
                      c='DarkBlue')
# ax1.figure.savefig('H3K4me1.pdf')
ax2 = hm.plot.scatter(x='H3K9ac1K14ac0',
                      y='Cell Line',
                      c='Red')


# ### Eliminate NaN
# ### Split GCP and MET into train and test set respectively.

# In[32]:


GCP = np.nan_to_num(GCP, nan=0)
MET = np.nan_to_num(MET, nan=0)
Xtrain, Xtest, Ytrain, Ytest =train_test_split(GCP, MET, test_size=0.3, random_state=0)

print(Xtest.shape)

print(Ytrain.shape)


# ### Fit train sets with Linear regression to generate models

# In[4]:


GCP2MET_models = []
for i in range(Ytrain.shape[1]):
    mdl_G2M = LinearRegression().fit(Xtrain, Ytrain[:, i])
    GCP2MET_models.append(mdl_G2M)


MET2GCP_models = []

for j in range(Xtrain.shape[1]):
    mdl_M2G = LinearRegression().fit(Ytrain, Xtrain[:,j])
    MET2GCP_models.append(mdl_M2G)


# In[5]:


# def robust_zscore(df):
#     """
#     robust_zscore computes a modified version of the Z-score centered around the media and scaled by the median absolute deviation.
#     INPUT:
#       :param df: A data frame corresponding to numerical values where rows correspond to cells and columns correspond to genes.
#     OUTPUTS:
#       :output Z: A data frame of the robust Z-scores.
#       :output pval: A data frame of the p-values corresponding to the robust Z-scores.
#     # Compute model parameters
#     """
#     med = df.median(axis=0)
#     n = df.shape[1]
#     med_abs_dev = mad(df)
#     scale_factor = 1.4826
#     # Compute robust Z-score
#     Z = (df - med) / (med_abs_dev * scale_factor)
#     Z[~np.isfinite(Z)] = 0
#     # Compute P-value based on Z-score
#     pval = norm.sf(abs(Z))*2 
#     pval = pd.DataFrame(pval,
#                     index=df.index,
#                     columns=df.columns)
#     return(Z, pval)


# ### To evaluate models and get correlation coefficient, pvalue, and mse value.

# In[42]:


def evaluate_models(models, Xtest, Ytest):
    """
    evaluate_models returns results from the model predictions, including the pearson
    correlation coefficient, p-Values, and MSE.

    :param models:         A list of scikit-learn model objects.
    :param Xtest:          A numpy array or pandas dataframe containing validation set input data.
    :param Ytest:          A numpy array or pandas dataframe containing validation set output data.
    :return pred_resul:    A dictionary containing the final MSE, pValue, or rValue.
    
    """
    
    predictions = []
    rValue = list()
    pValue = list()
    MSE = list()
    for i in range(len(models)):
        mdl = models[i]
        Ypred = mdl.predict(Xtest)
        predictions.append(Ypred)   
    
        r, p_value = pearsonr(Ypred, Ytest[:, i])
        # Compute P-value based on r-value
        pval = norm.sf(abs(r))*2 
        rValue.append(r)
        pValue.append(pval)
        
        mse = mean_squared_error(Ypred, Ytest[:, i])
        MSE.append(mse)
        
    df_MSE = pd.DataFrame(MSE)
    df_pValue = pd.DataFrame(pValue)
    df_rValue = pd.DataFrame(rValue)
    
    return df_MSE, df_pValue, df_rValue


# ## GCP to MET
# ### perform 3-fold cross validation to find model that generate the smallest mse value.

# In[43]:


mse1,p1,pearson1 = evaluate_models(GCP2MET_models, Xtest, Ytest)
mse1.columns = ["GCP2MET_models"]

print(mse1)


# In[8]:



kf = KFold(n_splits=3)

KFold(n_splits=3, random_state=None, shuffle=True)
mse_result = pd.DataFrame()
pvalue_result = pd.DataFrame()
rvalue_result = pd.DataFrame()
for train_index, test_index in kf.split(Xtrain, Ytrain):

    X_train, X_test = Xtrain[train_index], Xtrain[test_index]
    y_train, y_test = Ytrain[train_index], Ytrain[test_index]
    v_GCP2MET_models = []
    final_metrics = []
    for i in range(y_train.shape[1]):
        mdl_G2M = LinearRegression().fit(X_train, y_train[:, i])
        v_GCP2MET_models.append(mdl_G2M)
        df_MSE, df_pValue, df_rValue = evaluate_models(v_GCP2MET_models, X_test, y_test)
    mse_result = pd.concat([mse_result,df_MSE],axis = 1)
    pvalue_result = pd.concat([pvalue_result,df_pValue],axis = 1)
    rvalue_result = pd.concat([rvalue_result,df_rValue],axis = 1)

    
mse_result.columns = ["mse1","mse2","mse3"]
mse_result = mse_result.T
print(mse_result)
       


# In[9]:


minMSE = mse_result.idxmin()
print(minMSE)


# In[10]:


mse_min = pd.DataFrame()
MSE_min = list()
for index in mse_result:
    v = mse_result[index][minMSE[index]]
    MSE_min.append(v)
#     print(v)
mse_min = pd.DataFrame(MSE_min)
mse_min.columns = ["v_GCP2MET_models"]
# print(mse_min)
raw = mse1.join(mse_min)
print(raw)


# ### Generate new list of model

# In[11]:


new_GCP2MET = list()
for index, row in raw.iterrows():
    if row['GCP2MET_models'] < row['v_GCP2MET_models']:
        
        new_GCP2MET.append(GCP2MET_models[index])
    else:
        new_GCP2MET.append(v_GCP2MET_models[index])
        
    
print(new_GCP2MET)


# ### Evaluate if the correlation coefficient generated with the list of best model is significant or not

# In[12]:


mse_G2M,p_G2M,pearson_G2M = evaluate_models(new_GCP2MET, Xtest, Ytest)
p_G2M.columns = ["GCP2MET_p"]
pearson_G2M.columns = ["GCP2MET_r"]
raw_G2M = p_G2M.join(pearson_G2M)
pearson_G2M.columns = ["GCP2MET_r_0.05"]
new_G2M = raw_G2M.join(pearson_G2M)

for index, row in new_G2M.iterrows():
    if row['GCP2MET_p'] >= 0.05:
        row['GCP2MET_r_0.05'] = 'NaN'

print(new_G2M)        
new_G2M.to_csv('Linear_Regression_G2M.csv')


# ## MET TO GCP (repeat the same process above)

# In[14]:


mse2,p2,pearson2 = evaluate_models(MET2GCP_models, Ytest, Xtest)
mse2.columns = ["MET2GCP_models"]

kf = KFold(n_splits=3)
KFold(n_splits=3, random_state=None, shuffle=True)
mse2_result = pd.DataFrame()
# pvalue_result = pd.DataFrame()
# rvalue_result = pd.DataFrame()
for train_index, test_index in kf.split(Xtrain, Ytrain):
    X_train, X_test = Xtrain[train_index], Xtrain[test_index]
    y_train, y_test = Ytrain[train_index], Ytrain[test_index]
    v_MET2GCP_models = []
    metrics = []
    for i in range(X_train.shape[1]):
        mdl_M2G = LinearRegression().fit(y_train, X_train[:, i])
        v_MET2GCP_models.append(mdl_M2G)
        df_MSE, df_pValue, df_rValue = evaluate_models(v_MET2GCP_models, y_test, X_test)
    mse2_result = pd.concat([mse2_result,df_MSE],axis = 1)

mse2_result.columns = ["mse1","mse2","mse3"]
mse2_result = mse2_result.T
print(mse2_result) 


# In[15]:


minMSE2 = mse2_result.idxmin()
# print(minMSE)

mse2_min = pd.DataFrame()
MSE2_min = list()
for index in mse2_result:
    v = mse2_result[index][minMSE2[index]]
    MSE2_min.append(v)
#     print(v)
mse2_min = pd.DataFrame(MSE2_min)
mse2_min.columns = ["v_MET2GCP_models"]
# print(mse2_min)
raw2 = mse2.join(mse2_min)

new_MET2GCP = list()
for index, row in raw2.iterrows():
    if row['MET2GCP_models'] < row['v_MET2GCP_models']:
        
        new_MET2GCP.append(MET2GCP_models[index])
    else:
        new_MET2GCP.append(v_MET2GCP_models[index])
        
    
print(new_MET2GCP)


# In[16]:


mse_M2G,p_M2G,pearson_M2G = evaluate_models(new_MET2GCP, Ytest, Xtest)
p_M2G.columns = ["MET2GCP_p"]
pearson_M2G.columns = ["MET2GCP_r"]
raw_M2G = p_M2G.join(pearson_M2G)
pearson_M2G.columns = ["MET2GCP_r_0.05"]
new_M2G = raw_M2G.join(pearson_M2G)

for index, row in new_M2G.iterrows():
    if row['MET2GCP_p'] >= 0.05:
        row['MET2GCP_r_0.05'] = 'NaN'

print(new_M2G)        
new_M2G.to_csv('Linear_Regression_M2G.csv')


# In[ ]:




