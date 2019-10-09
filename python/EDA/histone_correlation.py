#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: Marc Di Meo & Scott Campit
"""

import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.impute import SimpleImputer
from scipy.stats.stats import pearsonr

import pearson_matrix

#geneExp = pd.read_excel('./../../data/CCLE/Microarray/CCLE_Microarray.xlsx')
geneExp = pd.read_csv(
    './../../data/CCLE/RNASeq/CCLE_934_TPM_EBI.tsv', delimiter='\t')
gcpExp = pd.read_csv('./../../data/CCLE/GCP/GCP_proteomics_remapped.csv')

# remove extra stuff before comma
tmp = []
geneCellLineList = list(geneExp.columns)
for CCL in geneCellLineList:
    cleanedCCL = CCL.split(',')[-1]
    cleanedCCL = cleanedCCL.replace('-', '')
    cleanedCCL = cleanedCCL.replace(' ', '')
    cleanedCCL = cleanedCCL.replace('.', '')
    cleanedCCL = cleanedCCL.upper()
    tmp.append(cleanedCCL)

geneCellLineList = tmp
geneExp.columns = geneCellLineList

gcpCellLineList = list(gcpExp['Cell Line'])
commonCCL = set(geneCellLineList).intersection(set(gcpCellLineList))
commonCCL.add("GENES")

GeneExpImputer = SimpleImputer(missing_values=np.nan, strategy='median')
gcpExpImputer = SimpleImputer(missing_values=np.nan, strategy='median')

# Get common cancer cell lines
geneExp = geneExp.loc[:, list(commonCCL)]
#geneExp = geneExp.groupby('Genes', as_index=False).median()

# From this set, get histone modifying genes only.
ptmModifiers = pd.read_csv('./../../data/EpiFactors/knownEpiFactors.csv')
ptmModifierGenes = list(ptmModifiers['HGNC'])
ptmModifierGenes = [x.strip(' ') for x in ptmModifierGenes]

# Gene expression dataframe
geneExpGenes = list(geneExp['GENES'])
commonGenes = set(geneExpGenes).intersection(set(ptmModifierGenes))
geneExp = geneExp.reset_index()
geneExp = geneExp.drop('index', axis=1)
geneExp = geneExp[geneExp['GENES'].isin(list(commonGenes))]

geneExp = geneExp.rename(columns={'GENES': 'Cell line'})
geneExp = geneExp.set_index('Cell line')
geneExp = geneExp.transpose()

geneExp[geneExp == np.inf] = np.nan
GeneExpImputer.fit(geneExp)
geneExpTransformed = GeneExpImputer.transform(geneExp)
geneMatrix = pd.DataFrame(
    geneExpTransformed, columns=geneExp.columns, index=geneExp.index)
geneMatrix = geneMatrix.reset_index()
#meltedGeneMatrix = pd.melt(
#    geneMatrix, id_vars=['Genes'], var_name='Cell line', value_name='TPM')
#meltedGeneMatrix = meltedGeneMatrix.set_index('Cell line')
geneMatrix = geneMatrix.set_index('index')

# GCP Proteomics Dataframe
gcpExp = gcpExp[gcpExp['Cell Line'].isin(list(commonCCL))]
gcpExp = gcpExp.transpose()
gcpExp.columns = gcpExp.iloc[0]
gcpExp.drop(gcpExp.index[0], inplace=True)
gcpExp[gcpExp == np.inf] = np.nan
gcpExpImputer.fit(gcpExp)
gcpExpTransformed = gcpExpImputer.transform(gcpExp)
gcpMatrix = pd.DataFrame(
    gcpExpTransformed, index=gcpExp.index, columns=gcpExp.columns)
gcpMatrix = gcpMatrix.transpose()

#meltedGCPMatrix = pd.melt(
#    gcpMatrix, id_vars=['index'], var_name='Cell line', value_name='[Marker]')
#meltedGCPMatrix = gcpMatrix.set_index('Cell line')
#meltedGCPMatrix = meltedGCPMatrix.rename(columns={'index': 'Histone  Marker'})
# Merged dataframe
merged = pd.merge(geneMatrix, gcpMatrix,
                  how='inner', left_index=True, right_index=True)
merged.to_csv('expressionDFNoMelt.csv', index=True)

"""
#H3K4me1 --> MLL3 and SETD7
h3k4me1 = (2, len(common_celllines))
h3k4me1  = np.zeros(h3k4me1)

h3k4me1[0] = ccle_common_cellline_expression[no_repeat_genes.index('MLL3')]
h3k4me1[1] = ccle_common_cellline_expression[no_repeat_genes.index('SETD7')]

h3k4me1_averaged = np.mean(h3k4me1 ,axis = 0)


#H3K9me1K14ac0 --> EHMT1, EHMT2, and SETDB1
h3k9me1k14ac0 = (3,len(common_celllines))
h3k9me1k14ac0 = np.zeros(h3k9me1k14ac0)

h3k9me1k14ac0[0] = ccle_common_cellline_expression[no_repeat_genes.index(
    'EHMT1')]
h3k9me1k14ac0[1] = ccle_common_cellline_expression[no_repeat_genes.index(
    'EHMT2')]
h3k9me1k14ac0[2] = ccle_common_cellline_expression[no_repeat_genes.index(
    'SETDB1')]

h3k9me1k14ac0_averaged = np.mean(h3k9me1k14ac0, axis = 0)

#H3K9me2K14ac0 --> EHMT1, EHMT2, SUV39H1 and SUV39H2
h3k9me2k14ac0 = (4,len(common_celllines))
h3k9me2k14ac0 = np.zeros(h3k9me2k14ac0)

h3k9me2k14ac0[0] = ccle_common_cellline_expression[no_repeat_genes.index(
    'EHMT1')]
h3k9me2k14ac0[1] = ccle_common_cellline_expression[no_repeat_genes.index(
    'EHMT2')]
h3k9me2k14ac0[2] = ccle_common_cellline_expression[no_repeat_genes.index(
    'SUV39H1')]
h3k9me2k14ac0[3] = ccle_common_cellline_expression[no_repeat_genes.index(
    'SUV39H2')]

h3k9me2k14ac0_averaged = np.mean(h3k9me2k14ac0, axis = 0)

#H3K9me3K14ac0 --> SETDB1, SETDB2, SUV39H1 and SUV39H2
h3k9me3k14ac0 = (4,len(common_celllines))
h3k9me3k14ac0 = np.zeros(h3k9me3k14ac0)

h3k9me3k14ac0[0] = ccle_common_cellline_expression[no_repeat_genes.index(
    'SETDB1')]
h3k9me3k14ac0[1] = ccle_common_cellline_expression[no_repeat_genes.index(
    'SETDB2')]
h3k9me3k14ac0[2] = ccle_common_cellline_expression[no_repeat_genes.index(
    'SUV39H1')]
h3k9me3k14ac0[3] = ccle_common_cellline_expression[no_repeat_genes.index(
    'SUV39H2')]

h3k9me3k14ac0_averaged = np.mean(h3k9me3k14ac0, axis = 0)

#H3K79me1 --> DOT1L

h3k79me1 = ccle_common_cellline_expression[no_repeat_genes.index('DOT1L')]

#H3K79me2 --> DOT1L

h3k79me2 = ccle_common_cellline_expression[no_repeat_genes.index('DOT1L')]


#May need to use this for more of the data
#Not used as of now
#Figure out way to correspond outlier to both data sets
def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]


h3_markers_list =list(dict.fromkeys(h3_markers))

h3k4me1_correlation = pearsonr(
    h3k4me1_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K4me1")])
h3k4me1_correlation = h3k4me1_correlation[0]

h3k9me1k14ac0_correlation = pearsonr(
    h3k9me1k14ac0_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K9me1K14ac0")])
h3k9me1k14ac0_correlation = h3k9me1k14ac0_correlation[0]

h3k9me2k14ac0_correlation = pearsonr(
    h3k9me2k14ac0_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K9me2K14ac0")])
h3k9me2k14ac0_correlation = h3k9me2k14ac0_correlation[0]

h3k9me3k14ac0_correlation = pearsonr(
    h3k9me3k14ac0_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K9me3K14ac0")])
h3k9me3k14ac0_correlation = h3k9me3k14ac0_correlation[0]

h3k79me1_correlation = pearsonr(
    h3k79me1, h3_common_cellline_expression[h3_markers_list.index("H3K79me1")])
h3k79me1_correlation = h3k79me1_correlation[0]

h3k79me2_correlation = pearsonr(
    h3k79me2, h3_common_cellline_expression[h3_markers_list.index("H3K79me2")])
h3k79me2_correlation = h3k79me2_correlation[0]

#New calculation with removed outlier at row 352
h3k9me2k14ac0_averaged = np.delete(h3k9me2k14ac0_averaged, 382)
deleted_outlier = np.delete(
    h3_common_cellline_expression[h3_markers_list.index("H3K9me2K14ac0")], 382)

h3k9me2k14ac0_correlation = pearsonr(h3k9me2k14ac0_averaged,deleted_outlier)
h3k9me2k14ac0_correlation = h3k9me2k14ac0_correlation[0]


#All correlations
markers = ('H3K4me1', 'H3K9me1K14ac0', 'H3K9me2K14ac0',
           'H3K9me3K14ac0', 'H3K79me1', 'H3K79me2')
y_pos = np.arange(len(markers))
correlation = [h3k4me1_correlation, h3k9me1k14ac0_correlation, h3k9me2k14ac0_correlation,
    h3k9me3k14ac0_correlation, h3k79me1_correlation, h3k79me2_correlation]

plt.bar(y_pos, correlation, align='center', alpha=0.5)
plt.xticks(y_pos, markers)
plt.xticks(rotation=45)
plt.xlabel('Histone Markers')
plt.ylabel('Correlation')
plt.title('Histone Marker Correlation Values')

plt.show()

#Scatter Plot
ccle_expression_axis = h3k9me2k14ac0_averaged
h3_expression_axis = deleted_outlier

plt.scatter(ccle_expression_axis, h3_expression_axis, label = 'Celllines')
plt.xlabel('CCLE Expression Data')
plt.ylabel('H3 Expression Data')
plt.legend(loc=4)
plt.title('H3K9me2K14ac0 Expression Comparison')

plt.show()

#ccle_expression_axis = h3k9me3k14ac0_averaged
#h3_expression_axis = h3_common_cellline_expression[h3_markers_list.index("H3K9me3K14ac0")]
#
#plt.scatter(ccle_expression_axis, h3_expression_axis, label = 'Celllines')
#plt.xlabel('CCLE Expression Data')
#plt.ylabel('H3 Expression Data')
#plt.legend(loc=4)
#plt.title('H3K9me3K14ac0 Expression Comparison')
#
#plt.show()
"""
