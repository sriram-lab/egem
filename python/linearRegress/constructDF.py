"""
constructDF

Make the dataframes that will be used for your linear regression analyses
@author: Scott Campit
"""
import pandas as pd

def makeDF(ExpressionFile):
    """
    """
    if ExpressionFile == "~/Data/Regression/eGEMM/CCLE.csv":
        df = pd.read_csv(ExpressionFile, index_col = 0)

        Comp = pd.read_csv(
            '~/Data/FBA/eGEMM/CCLE_Comp_ATP.csv',
            index_col=0)

        CompModel = pd.merge(df, Comp, how='inner',
                                 left_index=True, right_index=True)

    elif ExpressionFile == "~/Data/Regression/eGEMM/LeRoy.csv":
        df = pd.read_csv(ExpressionFile, index_col=0)

        Comp = pd.read_csv(
            '~/Data/FBA/eGEMM/LeRoy_Comp_ATP.csv', index_col=0)

        CompModel = pd.merge(df, Comp, how='inner',
                                  left_index=True, right_index=True)

    return CompModel

CCLECompModel = makeDF("~/Data/Regression/eGEMM/CCLE.csv")
CCLECompModel.to_csv('~/Data/Regression/eGEMM/CCLECompModel.csv')

LeRoyCompModel = makeDF("~/Data/Regression/eGEMM/LeRoy.csv")
LeRoyCompModel.to_csv('~/Data/Regression/eGEMM/LeRoyCompModel.csv')
