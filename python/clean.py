"""
Cleaning data from CCLE data obtained from Ghandi et al., 2019
@author: Scott Campit
"""

import pandas as pd

# Basic data clean up
df = pd.read_csv(r'/mnt/c/Users/scampit/Desktop/MeGEM/data/CCLE_GCP.csv')
df = df.drop('BroadID', axis=1)

media = pd.read_excel(r'/mnt/c/Users/scampit/Desktop/MeGEM/data/summary.xlsx', sheet_name='Cell Line Annotations', usecols = ['CCLE_ID', 'Growth.Medium'])

df = pd.merge(df, media, left_on='CellLineName', right_on='CCLE_ID')
tmp = df['Growth.Medium'].nunique()
print(tmp)

#df['CellLineName'] = df['CellLineName'].str.split('_').str[0]

# output everything for matlab

#df['CellLineName'].to_csv(r'/mnt/c/Users/scampit/Desktop/MeGEM/matlab/ccle_names.txt', header=False, index=False)

#_ = df.pop('CellLineName')

#h3marks = df.columns.tolist()

#with open(r'/mnt/c/Users/scampit/Desktop/MeGEM/matlab/ccle_h3marks2.txt', 'w') as f:
#    for i in h3marks:
#        print(i)
#        f.write("%s\n" % i)

#df.to_csv(r'/mnt/c/Users/scampit/Desktop/MeGEM/matlab/h3_relval.txt', header=False, index=False)
