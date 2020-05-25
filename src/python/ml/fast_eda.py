"""
Fast exploratory data analysis

This module outputs several simple plots that are *not* publication worthy, but are meant to serve diagnostic functions.

@author: Scott Campit
"""

# Essential
import pandas as pd
import numpy as np

# Data visualization
import matplotlib.pyplot as plt
import seaborn as sns

# Custom scripts
import preprocess

# Load CCLE data
db = '/home/scampit/Data/Expression/Combined/mapped_ccle_data.xlsx'
metabolomics = pd.read_excel(db,
                             sheet_name='CCLE_Metabolomics',
                             index_col='CCL')
gcp = pd.read_excel(db,
                    sheet_name='CCLE_GCP',
                    index_col='CCL')

# 1. In-sample correlation plots
gcp_corr = gcp.corr(method='pearson')
#plt.figure()
#ax1 = sns.heatmap(gcp_corr, annot=True,
#                  xticklabels=1, yticklabels=1,
#                  cmap='RdBu_r')
#plt.show()

#metabolomics_corr =  metabolomics.corr(method='pearson')
#plt.figure()
#ax2 = sns.heatmap(metabolomics_corr, annot=True,
#                  xticklabels=1, yticklabels=1)
#plt.show()

# 2. Out-of-sample correlation
#epimet = gcp.corrwith(metabolomics, axis=1)
#print(epimet)

# 3. Dendrogram
#gcp_dist = 1 - gcp_corr
#sns.clustermap(gcp_dist, method='average',
#               metric='euclidean', annot=True,
#               xticklabels=1, yticklabels=1,
#               cmap='RdBu_r')
#plt.show()