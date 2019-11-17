#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import re
import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from scipy.stats.stats import pearsonr


def standardizeCellLineLabels(CellLineList):
    """
    """
    tmp = []
    for CCL in CellLineList:
        CCL = str(CCL)
        cleanedCCL = CCL.split(',')[-1]
        cleanedCCL = cleanedCCL.replace('-', '')
        cleanedCCL = cleanedCCL.replace(' ', '')
        cleanedCCL = cleanedCCL.replace('.', '')
        cleanedCCL = cleanedCCL.upper()
        tmp.append(cleanedCCL)
    return tmp


def getCommonCellLineSet(geneCCLList, proteomicsCCLList):
    """
    """
    commonCCL = set(geneCCLList).intersection(set(proteomicsCCLList))
    commonCCL.add("GENES")
    return commonCCL


def medianImpute(df):
    """
    """
    df[df == np.inf] = np.nan
    MedianImputer = SimpleImputer(missing_values=np.nan, strategy='median')
    MedianImputer.fit(df)
    transformedDF = MedianImputer.transform(df)
    matrix = pd.DataFrame(
        transformedDF, index=df.index, columns=df.columns)
    return matrix


def makeProteomicsDF(df, commonCCL):
    """
    """
    df = df[df['Cell Line'].isin(list(commonCCL))]
    df = df.transpose()
    df.columns = df.iloc[0]
    df.drop(df.index[0], inplace=True)
    proteomicsMatrix = medianImpute(df)
    proteomicsMatrix = proteomicsMatrix.transpose()

    return proteomicsMatrix


def makeFinalMerge(geneDF, proteomicsDF):
    """
    """
    merged = pd.merge(geneDF, proteomicsDF,
                      how='inner', left_index=True, right_index=True)
    return merged


def constructGeneExpressionArray(geneExp, commonCCL, ptmModifierGenes):
    """
    """
    geneExp = geneExp.loc[:, list(commonCCL)]

    geneExpGenes = list(geneExp['GENES'])
    commonGenes = set(geneExpGenes).intersection(set(ptmModifierGenes))
    geneExp = geneExp.reset_index()
    geneExp = geneExp.drop('index', axis=1)
    geneExp = geneExp[geneExp['GENES'].isin(list(commonGenes))]

    geneExp = geneExp.rename(columns={'GENES': 'Cell line'})
    geneExp = geneExp.set_index('Cell line')
    geneExp = geneExp.transpose()

    geneMatrix = medianImpute(geneExp)
    geneMatrix = geneMatrix.reset_index()
    geneMatrix = geneMatrix.set_index('index')

    return geneMatrix


# Read in data
geneExp = pd.read_csv(
    '~/Data/RNASeq/CCLE/CCLE_934_TPM_EBI.tsv', delimiter='\t')
gcpExp = pd.read_csv(
    '~/Data/Proteomics/GCP/GCP_proteomics_remapped.csv')
leroy = pd.read_excel(
    '~/Data/Proteomics/LeRoy/LeRoy_et_al.xlsx', sheet_name='average')

# standardize cell line labels
geneExpCCL = list(geneExp.columns)
geneCellLineList = standardizeCellLineLabels(geneExpCCL)
LeRoyCCL = list(leroy['Cell Line'])
LeRoyCellLineList = standardizeCellLineLabels(LeRoyCCL)
gcpCCL = list(gcpExp['Cell Line'])
gcpCellLineList = standardizeCellLineLabels(gcpCCL)

geneExp.columns = geneCellLineList
leroy['Cell Line'] = LeRoyCellLineList
gcpExp['Cell Line'] = gcpCellLineList

# Get common cell lines between gene expression and proteomics datasets
commonCCL_LeRoy = getCommonCellLineSet(geneCellLineList, LeRoyCCL)
commonCCL_GCP = getCommonCellLineSet(geneCellLineList, gcpCCL)

# From this set, get histone modifying genes only.
#ptmModifiers = pd.read_csv('~/Data/Mappings/HistoneGeneMaps/EpiFactors/knownEpiFactors.csv')
ptmModifiers = pd.read_csv(
    '~/Data/Mappings/HistoneGeneMaps/EpiFactors/MeAcOnly.csv')
ptmModifierGenes = list(ptmModifiers['HGNC'])
ptmModifierGenes = [x.strip(' ') for x in ptmModifierGenes]

# Make gene expression array using the ptm modifier list
geneExpLeRoy = constructGeneExpressionArray(
    geneExp, commonCCL_LeRoy, ptmModifierGenes)

geneExpGCP = constructGeneExpressionArray(
    geneExp, commonCCL_GCP, ptmModifierGenes)

# Make the proteomics dataframe
gcpMatrix = makeProteomicsDF(gcpExp, commonCCL_GCP)
LeRoyMatrix = makeProteomicsDF(leroy, commonCCL_LeRoy)

CCLEMerged = makeFinalMerge(geneExpGCP, gcpMatrix)
LeRoyMerged = makeFinalMerge(geneExpLeRoy, LeRoyMatrix)

CCLEMerged.to_csv('~/Data/Regression/eGEMM/CCLERegression.csv', index=True)
LeRoyMerged.to_csv('~/Data/Regression/eGEMM/LeRoyRegression.csv', index=True)

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
