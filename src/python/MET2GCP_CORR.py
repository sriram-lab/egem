#!/usr/bin/env python
# coding: utf-8

# # Chromatin Profile and Metabolomics Correlation
# **Author**: Rachael Jin

# ## Summary
# This notebook computes the Pearson correlation coefficient and p-value between each histone marker and metabolite pair for 800+ cancer cell lines.

# ## Import data
# First, we'll read in the metabolomics and chromatin profile data, which are saved as Excel workbooks.

# In[1]:


import pandas as pd
import numpy as np
import scipy.stats
from scipy.stats import kendalltau, pearsonr, spearmanr
from scipy import stats
md = pd.read_excel(io='CCLE metabolomics dataset.xlsx',sheet_name="All")
md.head()


# In[2]:


mt = md.drop(['Tissue', 'Medium','Culture'], axis=1)
print('\n\nmd after deleting column\n--------------')
print(mt)


# In[3]:


hm = pd.read_csv('GCP_proteomics_remapped.csv')
hm.head()


# ## Merge datasets based on unique cancer cell line name
# Next, we'll concatenate the two dataframes and match based on cancer cell lines.

# In[4]:


merge_tb = mt.merge(hm,how='inner',left_on='CCL', right_on='Cell Line')
merge_tb.head()


# In[5]:


mt.info()
print('\n')
hm.info()
print('\n')
merge_tb.info()


# ## Separate dataframes 
# Now that the data is matched by cell lines, we can separate the dataframes again.

# In[6]:


metabolites = merge_tb.iloc[:,1:226]
metabolites.head()


# In[7]:


histone_markers = merge_tb.iloc[:,227:269]
histone_markers.head()


# ## Compute pearson correlation coefficient and pvalue
# Finally, we'll compute the correlation coefficients between metabolites and histone markers and the p-value of correlation. Note that we're also computing metabolite-metabolite and histone-histone correlations. While those are interesting as well, we'll ignore those for downstream analyses.
# 

# In[8]:


correlation = merge_tb.corr(method ='pearson')
correlation.to_csv('correlation.csv')


# In[9]:


corr = pd.concat([metabolites, histone_markers], axis=1, keys=['metabolites', 'histone_markers']).corr().loc['metabolites', 'histone_markers']
corr


# In[10]:


corr.to_csv('corr.csv')


# In[11]:


def calculate_pvalues(df,mt,hm):
    """
    :param df: A pandas dataframe containing the merged table of "CCLE metabolomics dataset" and "GCP_proteomics_remapped". Rows correspond to X, Columns correspond to Y.
    :param mt: A pandas dataframe containing just the metabolites portion of the merged table.
    :param hm: A pandas dataframe containing just the histone markers portion of the mergerd table.
    :return newpvalues: A pandas dataframe containg correlations's p-values that are less or equal to 0.05
    
    """
    df = df.dropna()._get_numeric_data()
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    df_mt = mt.dropna()._get_numeric_data()
    df_hm = hm.dropna()._get_numeric_data()
    df_mtcols = pd.DataFrame(columns=df_mt.columns)
    df_hmcols = pd.DataFrame(columns=df_hm.columns)
    newpvalues = df_mtcols.transpose().join(df_hmcols, how='outer')
    for r in (df.columns):
        for c in (df.columns):
            pvalues[r][c] = round(pearsonr(df[r], df[c])[1],4)
    for r in (df_mt.columns):
        for c in (df_hm.columns):
            if pvalues[r][c] <= 0.05:
                newpvalues[c][r] = pvalues[r][c]
    newpvalues.to_csv('pvalues_0.05.csv')
    return newpvalues

calculate_pvalues(merge_tb,metabolites,histone_markers)


# In[14]:


def generate_corrtable(df,mt,hm,corr):
    """
    :param df: A pandas dataframe containing the merged table of "CCLE metabolomics dataset" and "GCP_proteomics_remapped". Rows correspond to X, Columns correspond to Y.
    :param mt: A pandas dataframe containing just the metabolites portion of the merged table.
    :param hm: A pandas dataframe containing just the histone markers portion of the mergerd table.
    :param corr: A pandas dataframe containing previous generated correlation value between metabolites and histone markers.
    :return new_corr: A pandas dataframe containg correlation that has p-values <= 0.05
    
    """
    df = df.dropna()._get_numeric_data()
    dfcols = pd.DataFrame(columns=df.columns)
    pvalues = dfcols.transpose().join(dfcols, how='outer')
    df_mt = mt.dropna()._get_numeric_data()
    df_hm = hm.dropna()._get_numeric_data()
    df_mtcols = pd.DataFrame(columns=df_mt.columns)
    df_hmcols = pd.DataFrame(columns=df_hm.columns)
    newpvalues = df_mtcols.transpose().join(df_hmcols, how='outer')
    new_corr = df_mtcols.transpose().join(df_hmcols, how='outer')
    for r in (df.columns):
        for c in (df.columns):
            pvalues[r][c] = round(pearsonr(df[r], df[c])[1],4)
    for r in (df_mt.columns):
        for c in (df_hm.columns):
            if pvalues[r][c] <= 0.05:
                newpvalues[c][r] = pvalues[r][c] 
    for r in (df_mt.columns):
        for c in (df_hm.columns):
            if np.isnan(newpvalues[c][r]) == False:
                new_corr[c][r] = corr[c][r]
                if (new_corr[c][r]) > 0:
                    new_corr[c][r] = 1
                elif (new_corr[c][r] < 0):
                    new_corr[c][r] = -1
    new_corr.to_csv('new_corr.csv')
    return new_corr
    
generate_corrtable(merge_tb,metabolites, histone_markers, corr)


# In[ ]:




