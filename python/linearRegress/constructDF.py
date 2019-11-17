"""
constructDF

Make the dataframes that will be used for your linear regression analyses
@author: Scott Campit
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


def makeDF(ExpressionFile):
    """
    """
    if ExpressionFile == "~/Data/Regression/eGEMM/MeAcErasers/CCLERegression.csv":
        df = pd.read_csv(ExpressionFile, index_col=0)

        Comp = pd.read_csv(
            '~/Data/FBA/eGEMM/CCLE_Comp_ATP.csv',
            index_col=0)

        CompModel = pd.merge(df, Comp, how='inner',
                             left_index=True, right_index=True)

    elif ExpressionFile == "~/Data/Regression/eGEMM/MeAcErasers/LeRoyRegression.csv":
        df = pd.read_csv(ExpressionFile, index_col=0)

        Comp = pd.read_csv(
            '~/Data/FBA/eGEMM/LeRoy_Comp_ATP.csv', index_col=0)

        CompModel = pd.merge(df, Comp, how='inner',
                             left_index=True, right_index=True)

    return CompModel


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
    '~/Data/Mappings/HistoneGeneMaps/EpiFactors/KDMT_HDAC_only.csv')
#ptmModifiers = pd.read_csv('~/Data/Mappings/Histone#GeneMaps/EpiFactors/MeAcOnly.csv')
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

CCLEMerged.to_csv(
    '~/Data/Regression/eGEMM/MeAcErasers/CCLERegression.csv', index=True)
LeRoyMerged.to_csv(
    '~/Data/Regression/eGEMM/MeAcErasers/LeRoyRegression.csv', index=True)

CCLECompModel = makeDF(
    "~/Data/Regression/eGEMM/MeAcErasers/CCLERegression.csv")
CCLECompModel.to_csv('~/Data/Regression/eGEMM/MeAcErasers/CCLECompModel.csv')

LeRoyCompModel = makeDF(
    "~/Data/Regression/eGEMM/MeAcErasers/LeRoyRegression.csv")
LeRoyCompModel.to_csv(
    '~/Data/Regression/eGEMM/MeAcErasers/LeRoyCompModel.csv')
