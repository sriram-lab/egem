#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 09:37:33 2019

@author: marcdimeo
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

cellline_df = df['CellLineName'].copy()

cellline_dic = {}
for celllines in cellline_df:
    split = celllines.split('_')
    cellline_dic[split[0]] = " ".join(split[1:])

cellline_list =[]
tissue_list = []
for keys in cellline_dic:
    cellline_list.append(keys)
    value = cellline_dic[keys]
    tissue_list.append(value)

df.index.name = "Cell Line"

temp_df = {'Cell Line' : cellline_list, 'Tissue' : tissue_list }
temp_df = pd.DataFrame(temp_df)

del df['CellLineName']

