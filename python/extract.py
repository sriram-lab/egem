"""
Extracting data from CCLE data obtained from Ghandi et al., 2019

The files that are needed as the input for the script:
  * CCLE_GCP.csv - obtained from the supplementary data under a different name, but it is the file that contains the global chromatin proteomics dataset

The files created from this script:
  * ccle_names.txt - the names of the CCLE cell lines corresponding to the H3 markers
  * h3marks.txt - the names of the H3 markers
  * h3_relval.txt - the values of the relative proteomics intensities obtained from Ghandi et al., 2019.
@author: Scott Campit
"""

import pandas as pd

def extract():
    """
    extract will perform the actual extraction
    """

    # Basic data clean up
    df = pd.read_csv(r'/mnt/c/Users/scampit/Desktop/MeGEM/data/CCLE_GCP.csv')
    df = df.drop('BroadID', axis=1)

    media = pd.read_excel(r'/mnt/c/Users/scampit/Desktop/MeGEM/data/summary.xlsx', sheet_name='Cell Line Annotations', usecols = ['CCLE_ID', 'Growth.Medium'])

    df = pd.merge(df, media, left_on='CellLineName', right_on='CCLE_ID')
    tmp = df['Growth.Medium'].nunique()
    print(tmp)

    df['CellLineName'] = df['CellLineName'].str.split('_').str[0]

    # output everything for matlab

    df['CellLineName'].to_csv(r'/mnt/c/Users/scampit/Desktop/MeGEM/matlab/ccle_names.txt', header=False, index=False)

    _ = df.pop('CellLineName')

    h3marks = df.columns.tolist()

    with open(r'/mnt/c/Users/scampit/Desktop/MeGEM/matlab/ccle_h3marks2.txt', 'w') as f:
        for i in h3marks:
            print(i)
            f.write("%s\n" % i)

    df.to_csv(r'/mnt/c/Users/scampit/Desktop/MeGEM/matlab/h3_relval.txt', header=False, index=False)

def iterred(sheetnam=''):
    """
    iterred will iterate over the files outputted from the eGEM metabolic model and construct a list of dataframes that can be used for visualizations.
    """
    df = pd.read_excel('/mnt/c/Users/scampit/Desktop/MeGEM/matlab/tables/eGEMn.xlsx', sheet_name=sheetnam)
    df = df[df.columns.drop(list(df.filter(regex='.1')))]
    df = df.drop('Unnamed: 22', axis=1)
    df = df.head(50)
    return df
