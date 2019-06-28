#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 09:37:33 2019

@author: marcdimeo
"""

"""
This script is intended to serve as test/method to group cell lines according to their tissues
which will that be used for tissue specific analysis and histone correlation
"""


import pandas as pd
import numpy as np

"""
TISSUE SPECIFIC ANALYSIS
"""

#Take cell lines and match them with tissues and average out results
#Cell line --> Tissue correspondence is found in CCLE_GCP.csv\

#Change directory so it finds this file in 'data' and not in 'python' where I moved it.
df = pd.read_csv('CCLE_GCP.csv')


temp_df = pd.DataFrame(df['CellLineName'].str.split('_',1).tolist(),columns=['CellLine', 'Tissue'])

i = 0 
for tissue in temp_df.Tissue:
    tissue = tissue.split('_')
    tissue = " ".join(tissue)
    temp_df.Tissue[i] = tissue
    i = i+1
    
del df['CellLineName']

df.insert(0, 'CellLine', temp_df.CellLine , True)
df.insert(2, 'Tissue', temp_df.Tissue, True)

#Fill NaN values
df = df.fillna(method='ffill')

matrix = df.to_numpy()
        