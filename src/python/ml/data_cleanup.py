"""
Data cleanup

This script specifically formats the data to have matching entries between different datasets. This is to make loading
data for the statistical learning models less convoluted.

TO DO:
  * Perform this data cleaning step with the CCLE microarray data
  * Perform this data cleaning step with the CCLE RNASeq data
  * Perform this data cleaning step with the LeRoy histone marker data
@author: Scott Campit
"""

import pandas as pd

# Load CCLE metabolomics data
metabolomics = pd.read_csv('~/Data/Expression/Metabolomics/CCLE/CCLE_metabolomics.csv')
del metabolomics['DepMap_ID']
metabolomics["CCL"], metabolomics["Tissue"] = metabolomics["CCLE_ID"].str.split('_', 1).str

# Load CCLE GCP data
gcp = pd.read_csv('~/Data/Expression/Proteomics/GCP/CCLE_GCP.csv')
gcp["CCL"], gcp["Tissue"] = gcp["CellLineName"].str.split('_', 1).str

# Get intersection between metabolomics and GCP datasets
combine_df = pd.merge(metabolomics, gcp, how='inner', left_on='CCL', right_on='CCL')
matched_metabolomics = combine_df.loc[:, list(metabolomics.columns)]
matched_metabolomics = matched_metabolomics.set_index('CCL')
matched_metabolomics = matched_metabolomics.drop(['CCLE_ID', 'Tissue'], axis=1)

matched_gcp = combine_df.loc[:, list(gcp.columns)]
matched_gcp = matched_gcp.set_index('CCL')
matched_gcp = matched_gcp.drop(['CellLineName', 'BroadID', 'Tissue'], axis=1)

# Remove specific columns with mostly NaN in GCP dataset
matched_gcp = matched_gcp.drop(['H3K18ac0K23ub1', 'H3K56me1'], axis=1)
matched_metabolomics = matched_metabolomics[~matched_gcp.isnull().any(axis=1)]
matched_gcp = matched_gcp[~matched_gcp.isnull().any(axis=1)]


# Save as a single file for easy mapping
with pd.ExcelWriter('/home/scampit/Data/Expression/Combined/mapped_ccle_data.xlsx', engine='xlsxwriter') as writer:
    matched_metabolomics.to_excel(writer, sheet_name='CCLE_Metabolomics')
    matched_gcp.to_excel(writer, sheet_name='CCLE_GCP')
writer.save()

