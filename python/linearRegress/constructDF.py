"""


@author: Scott Campit
"""
import pandas as pd

df = pd.read_csv('expressionDFNoMelt.csv', index_col=0)
CCLEComp = pd.read_excel(
    './../../data/eGEM/FluxData.xlsx', sheet_name='CCLE_Comp', index_col='Cell Lines')
CCLEFVA = pd.read_excel(
    './../../data/eGEM/FluxData.xlsx', sheet_name='CCLE_FVA', index_col='Cell Lines')
#LeRoyComp = pd.read_excel(
#    './../../data/eGEM/FluxData.xlsx', sheet_name='LeRoy_Comp', index_col='Cell Lines')
#LeRoyFVA = pd.read_excel(
#    './../../data/eGEM/FluxData.xlsx', sheet_name='LeRoy_FVA', index_col='Cell Lines')

CCLECompModel = pd.merge(df, CCLEComp, how='inner',
                         left_index=True, right_index=True)
CCLEFVAModel = pd.merge(df, CCLEFVA, how='inner',
                        left_index=True, right_index=True)
#LeRoyCompModel = pd.merge(df, LeRoyComp, how='inner',
#                          left_index=True, right_index=True)
#LeRoyFVAModel = pd.merge(df, LeRoyFVA, how='inner',
#                         left_index=True, right_index=True)

CCLECompModel.to_csv('CCLECompModel.csv')
CCLEFVAModel.to_csv('CCLEFVAModel.csv')
#LeRoyCompModel.to_csv('LeRoyCompModel.csv')
#LeRoyFVAModel.to_csv('LeRoyFVAModel.csv')
