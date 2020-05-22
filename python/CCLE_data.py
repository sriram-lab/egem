#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 12:10:43 2019

@author: marcdimeo
"""

"""
This document serves as a way to find and manipulate all CCLE data into a dataset
which can be used for the histone correlation calculations
"""

<<<<<<< HEAD
import mygene
import pandas as pd
import numpy as np
from sklearn import preprocessing
import scipy.io as spio
import os
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns





"""
FUNCTIONS
"""
#Used to find common celllines
def common_celllines(a,b):
    a_set = set(a)
    b_set = set(b)
    common = a_set.intersection(b_set) 
    common = list(common)
    return common 

#Used to create correlation matrix
def Matrix(gene_list, histone_list, data1, data2):
    name = (len(gene_list), len(histone_list))
    name = np.zeros(name)
    i = 0 
    for gene in gene_list:
        j = 0
        for histone in histone_list:
            data1.loc[gene] = data1.loc[gene].fillna(method='ffill')
            data2.loc[histone] = data2.loc[histone].fillna(method='ffill')
            correlation = pearsonr(data1.loc[gene], data2.loc[histone])
            correlation = correlation[0]
            name[i][j] = correlation
            j = j+1
        i = i+1
    return name

#Used to make heatmaps
def Heatmap(data, size, colour, title):
    plt.figure(figsize=size)
    plt.xlabel('Histone Markers')
    plt.title(title)
    ax = sns.heatmap(data, cmap = colour)
    ax.set(xlabel = 'Histone Markers')
    plt.show

"""
GETTING EXPRESSION DATA
"""
data = pd.read_csv("./../../../CCLE_RNAseq_rsem_genes_tpm_20180929.txt", delim_whitespace=True)
data = data.replace(0, np.nan)
 
del data['transcript_ids']

gene_id_list = []

for gene_id in data['gene_id']:
    gene_id = gene_id.split('.')
    gene_id_list.append(gene_id[0])
    
data['gene_id'] = gene_id_list

mg = mygene.MyGeneInfo()
ginfo = mg.querymany(data['gene_id'], scopes='ensembl.gene', fields='symbol', species='human')

#removed repeated element
ginfo.pop(30544)

"""
CHANGE GENE IDs TO SYMBOLS
"""

gene = []
for entry in ginfo:
    if entry.get('symbol') != None:
        gene.append(entry['symbol'])
    else:
        gene.append("delete")

data['gene_id'] = gene

data.rename(columns={'gene_id' : 'Gene'}, inplace = True)

"""
FORMATTING DATA FRAME
"""

#get list of celllines and remove tissue type
celllines = list(data)
celllines.pop(0)

cellline_list =[]
tissue_list = []
for cell in celllines:
    temp = cell.split('_')
    cell = temp[0]
    cellline_list.append(cell)
    tissue = temp[1:]
    tissue = " ".join(tissue)
    tissue_list.append(tissue)

d = {'Celllines': cellline_list, 'Tissue' : tissue_list}
cell_df = pd.DataFrame(d)

df = data
df.columns = ['Gene'] + cellline_list

temp_df = df.copy()

del temp_df["Gene"]

x = temp_df.values.astype(float)
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
df_normalized = pd.DataFrame(x_scaled)

df_normalized.insert(0, 'Gene', df.Gene , True)
df_normalized.columns = ['Gene'] + cellline_list

#fill in NaN values
df_normalized = df_normalized.fillna(method='ffill')


"""
SEARCHING FOR MATCHING GENE ALOGORITHM
"""


search = "10013(HDAC6) 10014(HDAC5) 3065(HDAC1) 3066(HDAC2) 51564(HDAC7) 55869(HDAC8) 79885(HDAC11) 83933(HDAC10) 8841(HDAC3) 9734(HDAC9) 9759(HDAC4) 10524(KAT5) 10724(OGA) 11143(KAT7) 124359(CDYL2) 138474(TAF1L) 1387(CREBBP) 2033(EP300) 203611(CDY2B) 23522(KAT6B) 253175(CDY1B) 2648(KAT2A) 55140(ELP3) 6872(TAF1) 7994(KAT6A) 8202(NCOA3) 84148(KAT8) 8520(HAT1) 8648(NCOA1) 8850(KAT2B) 9085(CDY1) 9329(GTF3C4) 9425(CDYL) 9426(CDY2A) 9575(CLOCK) 10919(EHMT2) 11105(PRDM7) 2145(EZH1) 2146(EZH2) 23067(SETD1B) 29072(SETD2) 387893(KMT5A) 4297(KMT2A) 51111(KMT5B) 54904(NSD3) 55870(ASH1L) 55904(KMT2E) 56979(PRDM9) 58508(KMT2C) 6419(SETMAR) 64324(NSD1) 6839(SUV39H1) 7468(NSD2) 7799(PRDM2) 79723(SUV39H2) 79813(EHMT1) 8085(KMT2D) 80854(SETD7) 83852(SETDB2) 84193(SETD3) 84444(DOT1L) 84787(KMT5C) 93166(PRDM6) 9739(SETD1A) 9757(KMT2B) 9869(SETDB1) 22992(KDM2A) 23133(PHF8) 79697(RIOX1) 79831(KDM8) 84678(KDM2B)"

#Method 1

#search = search.split(" ")

#Method 2

search = search.split(" ")
search1 = []
s = False

for items in search:
    word  = ""
    for letters in items:
        if letters == "(":
            s = True
        if s == True:
            word = word + letters
        if letters == ")":
            s = False
    search1.append(word)

search = []
for items in search1: 
    items = items[1:len(items)-1]
    search.append(items)
        

for element in search:
    if element in gene:
        print(element + ":", "YES")
    else:
        print(element + ":","No")
        
"""
H3 RELVAL DATA
"""

os.chdir('/Users/marcdimeo/Desktop/University of Michigan Research/methylation-gem/matlab/new_var')

h3_celllines = np.array(spio.loadmat('correlation_value', squeeze_me=True)["h3_ccle_names_python"])
h3_celllines = h3_celllines.tolist()
h3_markers =   list(dict.fromkeys(spio.loadmat('correlation_value', squeeze_me=True)["h3_marks_python"]))
h3_expression = spio.loadmat('h3_relval', squeeze_me=True)["h3_relval"]

h3_relval = pd.DataFrame(h3_expression)
h3_relval = h3_relval.T
h3_relval.insert(0, 'Histone', h3_markers, True)
h3_relval.columns = ['Histone'] + h3_celllines

"""
NORMALIZE H3 DATA
"""

temp_df = h3_relval.copy()

del temp_df["Histone"]

x = temp_df.values.astype(float)
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
h3_normalized = pd.DataFrame(x_scaled)

h3_normalized.insert(0, 'Histone', h3_relval.Histone , True)
h3_normalized.columns = ['Histone'] + h3_celllines

h3_normalized = h3_normalized.fillna(method='ffill')



"""
REORGANIZING DATA
"""
#Finding common cell lines
h3_ccle_cellline = common_celllines(cellline_list, h3_celllines) 
h3_ccle_cellline.sort()
     
ccle_h3_df = df_normalized[['Gene']+ h3_ccle_cellline]
ccle_h3_df = ccle_h3_df.set_index('Gene')

h3_ccle_df = h3_normalized[['Histone']+h3_ccle_cellline]
h3_ccle_df= h3_ccle_df.set_index('Histone')  

h3_ccle_df= h3_ccle_df.fillna(method='ffill')
ccle_h3_df = ccle_h3_df.fillna(method='ffill')

#H3 and CCLE
h3_ccle_matrix = Matrix(search, h3_markers, ccle_h3_df, h3_ccle_df)
    
h3_ccle_matrix = pd.DataFrame(h3_ccle_matrix)
h3_ccle_matrix.columns = h3_markers
h3_ccle_matrix.insert(0, 'Genes', search , True)
h3_ccle_matrix = h3_ccle_matrix.set_index('Genes')

"""
LEROY DATA
"""
os.chdir('/Users/marcdimeo/Desktop/University of Michigan Research/methylation-gem/matlab/vars')

leroy_celllines =  spio.loadmat('supplementary_software_code', squeeze_me=True)["acetlevellist"]
leroy_celllines = leroy_celllines.tolist()
leroy_markers = spio.loadmat('methylation_proteomics_validation_data', squeeze_me=True)["acet_meth_list_rowlab"]
leroy_markers = leroy_markers.tolist()
leroy_expression=  spio.loadmat('hist_proteomics', squeeze_me=True)["acet_meth_listval"]

i = 0
for markers in leroy_markers:
    leroy_markers[i] = "H3" + markers
    i=i+1
    
i = 0
for markers in leroy_markers:
    if 'un' in markers:
        leroy_markers[i] = markers.replace('un', 'ac0')
    if 'ac' in markers:
        leroy_markers[i] = markers.replace('ac', 'ac1')
    i = i+1

leroy_expression = pd.DataFrame(leroy_expression)
leroy_expression.insert(0, 'Histones', leroy_markers , True)
leroy_expression = leroy_expression.set_index('Histones')
leroy_expression.columns = leroy_celllines

leroy_ccle_cellline = common_celllines(cellline_list, leroy_celllines) 
leroy_ccle_cellline.sort()

leroy_ccle_df = leroy_expression[leroy_ccle_cellline]
ccle_leroy_df = df_normalized[['Gene']+ leroy_ccle_cellline]
ccle_leroy_df = ccle_leroy_df.set_index('Gene')

        
leroy_ccle_matrix = Matrix(search, leroy_markers, ccle_leroy_df, leroy_ccle_df)   
    
leroy_ccle_matrix = pd.DataFrame(leroy_ccle_matrix)
leroy_ccle_matrix.columns = leroy_markers
leroy_ccle_matrix.insert(0, 'Genes', search , True)
leroy_ccle_matrix = leroy_ccle_matrix.set_index('Genes')


"""
HEATMAP
"""

Heatmap(h3_ccle_matrix, (10,5), 'Blues', 'H3 and CCLE Correlation Plot')
Heatmap(leroy_ccle_matrix, (10,5), 'Blues', 'LeRoy and CCLE Correlation Plot')

=======
import GEOparse
import mygene

"""
GETTING GENE NAMES
"""

file = open("GPL15308.txt", 'r')

gene_ids = []

for line in file:
    if line.startswith('!') or line.startswith('#') or line.startswith('I') or line.startswith('^') :
        pass
    else:
        ids = line.split("_")[0]
        gene_ids.append(ids)
file.close()

mg = mygene.MyGeneInfo()

gene_information = mg.querymany(gene_ids, scopes='entrezgene', fields='symbol', species='human')


genes = []

for line in gene_information:
    if 'symbol' not in line:
        pass
    else:
        gene = line['symbol']
        genes.append(gene)

"""
GETTING EXPRESSION DATA
"""

gse = GEOparse.get_GEO(filepath="./GSE36133_family.soft.gz")
>>>>>>> 48693a251744125a2ccc70193ab5dc8576deb6cb


