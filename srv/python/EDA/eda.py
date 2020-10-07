"""
eda.py is a Python script to make fast static plots that quickly visualizes data sets.
@author: Scott Campit
"""

import pandas as pd
import scipy
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

# File paths
gcpFile = '~/Data/Expression/Proteomics/GCP/MappedGCP_subset.csv'
metabolomicsFile = '~/Data/Expression/Metabolomics/CCLE/CCLE_metabolomics_subset.csv'
transcriptomicsFile = '~/Data/Expression/Microarray/CCLE/CCLE_Microarray.xlsx'
corrmatFile = '~/Data/Expression/Combined/GCP_metabolomics_corrmat.csv'

# Visualize the entire dataset

## Metabolomics
#df = pd.read_csv(metabolomicsFile)
#df = df.drop(['CCLE_ID'], axis=1)
#tmp = scipy.stats.zscore(df)
#df =  pd.DataFrame(tmp, columns=df.columns)
#print(df)

#sns.set(style="darkgrid")
#sns.stripplot(data=df,
#              jitter=0.25,
#              alpha=0.25)
#ax = sns.violinplot(data=df,
#                    orient='v')
#ax.set_xticklabels(df.columns,
#                   rotation=90)
#ax.set_title('Metabolomics distribution for subset of metabolites in CCLE dataset', loc='left', fontsize=30)
#ax.set_xlabel('Metabolites (subset)', fontsize=30)
#ax.set_ylabel('Z-score distribution', fontsize=30)
#ax.tick_params(axis='both', which='major', labelsize=18)
#ax.figure.subplots_adjust(bottom = 0.3)
#plt.show()

# Heatmap
#ax = sns.heatmap(df,
#            cmap='RdYlBu_r',
#            )
#ax.set_yticks([])
#ax.set_xticks([])
#plt.show()

## GCP
#df = pd.read_csv(gcpFile)
#df = df.drop(['CellLine'], axis=1)

# Violin plots
#sns.set(style="darkgrid")
#sns.stripplot(data=df,
#              jitter=0.25,
#              alpha=0.25)
#ax = sns.violinplot(data=df,
#                    orient='v')
#ax.set_xticklabels(df.columns,
#                   rotation=90)
#ax.set_title('GCP distribution for subset of histone markers in CCLE dataset', loc='left', fontsize=30)
#ax.set_xlabel('GCP (subset)', fontsize=30)
#ax.set_ylabel('Z-score distribution', fontsize=30)
#ax.tick_params(axis='both', which='major', labelsize=18)
#ax.figure.subplots_adjust(bottom = 0.3)
#plt.show()

# Heatmap
#ax = sns.heatmap(df,
#            cmap='RdYlBu_r',
#            )
#ax.set_yticks([])
#ax.set_xticks([])
#plt.show()

## Transcriptomics
#df = pd.read_excel(transcriptomicsFile)
#df = df.drop(['Gene'], axis=1)

# Heatmap
#ax = sns.heatmap(df,
#            cmap='RdYlBu_r',
#            mask=df.isnull()
#            )
#ax.set_yticks([])
#ax.set_xticks([])
#plt.show()

## Correlation map
#sns.set(style="darkgrid")
#df = pd.read_csv(corrmatFile)
#df = df.drop('GCP', axis=1)
#ax = sns.heatmap(df, cmap='RdBu_r')
#ax.set_yticks([])
#ax.set_xticks([])
#plt.show()

#corrmatFile2 = '~/Data/Expression/Combined/GCP_metabolomics_corrmat_subset.csv'
#sns.set(style="darkgrid")
#df = pd.read_csv(corrmatFile2)

# Subset analysis
#df = df[df['GCP'].isin(['H3K4me1', 'H3K4me2', 'H3K9ac1', 'H3K36me1', 'H3K36me2'])]
#GCP = df['GCP']
#df = df.drop('GCP', axis=1)

#f, ax = plt.subplots(1, 1)
#sns.heatmap(df, cmap='RdBu_r', yticklabels=GCP)
#ax.set_title('Correlation matrix for metabolite-histone marker pairs', loc='left', fontsize=45)
#ax.set_xlabel('Metabolites', fontsize=30)
#ax.set_ylabel('Global Chromatin Profiles', fontsize=30)
#ax.tick_params(axis='both', which='major', labelsize=45)
#ax.figure.subplots_adjust(bottom = 0.4)
#cbar = ax.collections[0].colorbar
#cbar.ax.tick_params(labelsize=45)
#plt.show()

## LASSO
lassoFile = '~/Downloads/PredictCCLEMetabolomicsFromGCP.xlsx'

sns.set(style="darkgrid")
#df = pd.read_excel(lassoFile, sheet_name='LASSO Coef All No P')
df = pd.read_excel(lassoFile, sheet_name='LASSO Metab Subset 2')

#ax = sns.heatmap(df, cmap='RdBu_r')
#ax.set_yticks([])
#ax.set_xticks([])
#plt.show()
df = df[df.columns.drop(list(df.filter(regex='Unnamed')))]
#print(df)

# All LASSO Heatmap
metabolites = df['Metabolites']
df = df.drop('Metabolites', axis=1)
ax = sns.heatmap(df, yticklabels=metabolites, xticklabels=df.columns,  cmap='RdBu_r')
ax.set_title('LASSO coefficients',
             loc='left',
             fontsize=45)
#ax.set_yticks([])
#ax.set_xticks([])
ax.set_xlabel('Global Chromatin Profiles', fontsize=30)
ax.set_ylabel('Metabolites', fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=45)
ax.figure.subplots_adjust(bottom = 0.35, left=0.2)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=45)
plt.show()

# Box Plots of accuracy
#shistacc = '~/Data/Models/eGEM/predict_metabolomics/supplementation_example_subset.csv'
#df = pd.read_csv(shistacc)
#df = pd.melt(df)
#print(df)

#sns.stripplot(x='variable', y='value', data=df,
#             jitter=0.25,
#              alpha=0.75,
#              s=20)
#ax = sns.boxplot(x='variable', y='value', data=df)
#ax.set_title('Accuracy distribution for nutrient supplementation',
#             loc='left',
#             fontsize=45)
#ax.set_xticklabels(df['variable'].unique(),
#                   rotation=90)
#ax.set_xlabel('Global Chromatin Profiles', fontsize=30)
#ax.set_ylabel('Accuracy distribution', fontsize=30)
#ax.axhline(y=0.33)
#ax.tick_params(axis='both', which='major', labelsize=45)
#ax.figure.subplots_adjust(bottom = 0.3)
#plt.show()
