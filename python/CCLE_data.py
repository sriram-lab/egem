#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 12:10:43 2019

@author: Marc Di Meo
"""

"""
This document serves as a way to find and manipulate all CCLE data into a dataset
which can be used for the histone correlation calculations
"""

import mygene
import pandas as pd
import numpy as np
from sklearn import preprocessing
import scipy.io as spio
import os
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.io as pio
#import plotly.graph_objects as go
import plotly.figure_factory as ff
from scipy.spatial.distance import pdist, squareform



"""
FUNCTIONS
"""

def common(a,b):
    """ This functions purpose is to find common elements of lists by converting
    each list to a set"""
    a_set = set(a)
    b_set = set(b)
    common = a_set.intersection(b_set) 
    common = list(common)
    return common 

#Used to create correlation matrix
def Matrix(gene_list, histone_list, data1, data2):
    """This functions purpose is to create a numpy matrix which is filled with 
    correlation values based on two data sets and corresponding lists"""
    Rmat = (len(gene_list), len(histone_list))
    Rmat = np.zeros(Rmat)
    Pmat = Rmat
    i = 0 
    for gene in gene_list:
        j = 0
        for histone in histone_list:
            data1.loc[gene] = data1.loc[gene].fillna(method='ffill')
            data2.loc[histone] = data2.loc[histone].fillna(method='ffill')
            correlation = pearsonr(data1.loc[gene], data2.loc[histone])
            Rmat[i][j] = correlation[0]
            Pmat[i][j]
            j = j+1
        i = i+1
    return Rmat, Pmat

def Heatmap(data, ylabels, size, colour, title):
    """Using dataframes based on the Matrix fucntion, this function will create 
    a heatmap"""
    plt.figure(figsize=size)
    plt.xlabel('Histone Markers')
    plt.title(title)
    ax = sns.heatmap(data, cmap = colour, square = True, yticklabels = ylabels)
    ax.set(xlabel = 'Histone Markers')
    plt.show
    
def pearson_dfs(gene_list, histone_list, data1, data2):
    """This functions purpose is to create two dataframes. The first dataframe will contain all the r-values between the two sets of data. The second will contain all the p-values for each r-value.
    
    This function is very specific. The dataframe created will contain histone markers as the columns and gene names as the index. 
    
    This function assumes you already took the intersection between the two arrays. If you do not do that, then the matrix dimensions will not agree, resulting in an error.
    
    The format for inputting matrix should have columns correspond to the celllines.
    """
    Rmat = (len(gene_list), len(histone_list))
    Rmat = np.zeros(Rmat)
    Pmat = Rmat[:]
    i = 0 
    for gene in gene_list:
        j = 0
        for histone in histone_list:
            data1.loc[gene] = data1.loc[gene].fillna(method='ffill')
            data2.loc[histone] = data2.loc[histone].fillna(method='ffill')
            correlation = pearsonr(data1.loc[gene], data2.loc[histone])
            Rmat[i][j] = correlation[0]
            Pmat[i][j]
            j = j+1
        i = i+1
        
    Rdf = pd.DataFrame(Rmat)
    Rdf.columns = histone_list
    Rdf.insert(0, 'Genes', gene_list, True)
    Rdf = Rdf.set_index('Genes')
    
    Pdf = pd.DataFrame(Pmat)
    Pdf.columns = histone_list
    Pdf.insert(0, 'Genes', gene_list, True)
    Pdf = Pdf.set_index('Genes')
    return Rdf, Pdf



    
def Remove(duplicate): 
    """This function will remove duplicated elements from lists"""
    final_list = [] 
    for num in duplicate: 
        if num not in final_list: 
            final_list.append(num) 
    return final_list

def Clustermap(data, size, colour, method, metric):
    plt.figure(figsize=size)
    ax = sns.clustermap(data, cmap = colour, method = method, metric= metric)
    plt.show
    
def PlotlyHeat(df, size, title, xaxis, yaxis):
    pio.renderers.default = "chrome"
    fig = go.Figure(
    data=(go.Heatmap(z = df, x =xaxis, y = yaxis, colorscale = "rdbu")),
    layout_title_text=title)
    
    if size == None:
        fig.update_layout(
                autosize=True)
    else:
        fig.update_layout(
                autosize=False,
                width = size[0],
                height = size[1])

    fig.show()
    


"""
"""

"""
"""

def TissueAnalysis(dictionary, common_cellline, data1, data2, name):
    df1 = data1.copy()
    df2 = data2.copy()
    for celllines in common_cellline:
        if dictionary[celllines] != name:
            del df1[celllines]
            del df2[celllines]
    matrix = Matrix(search, h3_markers, df1, df2)
    matrix = pd.DataFrame(matrix)
    matrix.columns = h3_markers
    matrix.insert(0, 'Genes', search , True)
    matrix = matrix.set_index('Genes')
    return matrix

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
#query = []
#for entry in ginfo:
#   query.append(entry['query'])
#   
#query = Remove(query)
   


"""
CHANGE GENE IDs TO SYMBOLS
"""

gene = []
for entry in ginfo:
    if entry.get('symbol') != None:
        gene.append(entry['symbol'])
    else:
        gene.append("delete")

#gene = Remove(gene)
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
h3_ccle_cellline = common(cellline_list, h3_celllines) 
h3_ccle_cellline.sort()
     
ccle_h3_df = df_normalized[['Gene']+ h3_ccle_cellline]
ccle_h3_df = ccle_h3_df.set_index('Gene')

h3_ccle_df = h3_normalized[['Histone']+h3_ccle_cellline]
h3_ccle_df= h3_ccle_df.set_index('Histone')  

h3_ccle_df= h3_ccle_df.fillna(method='ffill')
ccle_h3_df = ccle_h3_df.fillna(method='ffill')

#H3 and CCLE
h3_ccle_matrix = Matrix(search, h3_markers, ccle_h3_df, h3_ccle_df)[0]
    
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

leroy_ccle_cellline = common(cellline_list, leroy_celllines) 
leroy_ccle_cellline.sort()

leroy_ccle_df = leroy_expression[leroy_ccle_cellline]
ccle_leroy_df = df_normalized[['Gene']+ leroy_ccle_cellline]
ccle_leroy_df = ccle_leroy_df.set_index('Gene')

        
leroy_ccle_matrix = Matrix(search, leroy_markers, ccle_leroy_df, leroy_ccle_df)[0]   
    
leroy_ccle_matrix = pd.DataFrame(leroy_ccle_matrix)
leroy_ccle_matrix.columns = leroy_markers
leroy_ccle_matrix.insert(0, 'Genes', search , True)
leroy_ccle_matrix = leroy_ccle_matrix.set_index('Genes')

"""
RECON1
"""
recon1_list = []
recon1_genes = pd.read_excel(r'/Users/marcdimeo/Desktop/University of Michigan Research/methylation-gem/data/RECON1_genes.xlsx')
for genes in recon1_genes['Genes']:
    genes = genes.split('\'')
    recon1_list.append(genes[1])
    
i = 0
for genes in recon1_list:
    genes = genes.split('_')
    recon1_list[i] = genes[0]
    i = i+1
    
recon1info = mg.querymany(recon1_list, scopes='entrezgene', fields='symbol', species='human')

recon1_list = []
for entry in recon1info:
    if entry.get('symbol') != None:
        recon1_list.append(entry['symbol'])
    else:
        pass
    
recon1_list = Remove(recon1_list)
recon1_list = common(gene, recon1_list)

"""
CCLE(RECON1) and LEROY
"""
leroy_ccle_r1 = Matrix(recon1_list, leroy_markers, ccle_leroy_df, leroy_ccle_df)[0]

leroy_ccle_r1 = pd.DataFrame(leroy_ccle_r1)
leroy_ccle_r1.columns = leroy_markers
leroy_ccle_r1.insert(0, 'Genes', recon1_list , True)
leroy_ccle_r1 = leroy_ccle_r1.set_index('Genes')

"""
TISSUE ANALYSIS: Only Run if you have to it takes a very long time approx 2 hours
"""

#cellline_tissue = cell_df.values
#tissue_dict ={}
#for data in cellline_tissue:
#    tissue_dict[data[0]] = data[1]
#
#
#h3_ccle_lung = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df,'LUNG')
#h3_ccle_ovary = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df,'OVARY')
#h3_ccle_li = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df,'LARGE INTESTINE')
#h3_ccle_cns = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df,'CENTRAL NERVOUS SYSTEM')
#h3_ccle_hlt = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df,'HAEMATOPOIETIC AND LYMPHOID TISSUE')
#h3_ccle_pancreas = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df,'PANCREAS')
#h3_ccle_uat =  TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df, 'UPPER AERODIGESTIVE TRACT')
#h3_ccle_breast =  TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df, 'BREAST')
#h3_ccle_prostate =  TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df,'PROSTATE')
#h3_ccle_stomach =  TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df, 'STOMACH')
#h3_ccle_endometrium =   TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df,'ENDOMETRIUM')
#h3_ccle_bone =  TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df, 'BONE')
#h3_ccle_skin = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df, 'SKIN')
#h3_ccle_liver =  TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df, 'LIVER')
#h3_ccle_fibroblast =  TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df, 'FIBROBLAST')
#h3_ccle_st =  TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df,'SOFT TISSUE')
#h3_ccle_bt=  TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df,  'BILIARY TRACT')
#h3_ccle_ag = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df, 'AUTONOMIC GANGLIA')
#h3_ccle_pleura = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df, 'PLEURA')
#h3_ccle_ut = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df, 'URINARY TRACT')
#h3_ccle_kidney = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df,'KIDNEY')
#h3_ccle_oesophagus =  TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df,'OESOPHAGUS')
#h3_ccle_thyroid = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df, 'THYROID')
#h3_ccle_sg = TissueAnalysis(tissue_dict, h3_ccle_cellline, ccle_h3_df, h3_ccle_df, 'SALIVARY GLAND')

"""
FULL HISTONE ANALYSIS
"""

#h3k4 = h3_ccle_df.loc['H3K4me0':'H3K4ac1']
#h3k9 = h3_ccle_df.loc['H3K9me0K14ac0':'H3K9ac1K14ac1']
#h3k18 = h3_ccle_df.loc['H3K18ac0K23ac0':'H3K18ac0K23ub1']
#h3k27 = h3_ccle_df.loc['H3K27me0K36me0':'H3.3K27me0K36me0']
#h3k56 = h3_ccle_df.loc['H3K56me0':'H3K56me1']
#h3k79 = h3_ccle_df.loc['H3K79me0':]
#
#h3k4 = h3k4.values
#h3k9 = h3k9.values
#h3k18 = h3k18.values
#h3k27 = h3k27.values
#h3k56 = h3k56.values
#h3k79 = h3k79.values
#
#h3k4 = np.mean(h3k4, axis = 0)
#h3k9 = np.mean(h3k9, axis = 0)
#h3k18= np.mean(h3k18, axis = 0)
#h3k27 = np.mean(h3k27, axis = 0)
#h3k56 = np.mean(h3k56, axis = 0)
#h3k79 = np.mean(h3k79, axis = 0)
#
#h3_histone_matrix = (6, len(h3_ccle_cellline)+1)
#h3_histone_matrix = np.zeros(h3_histone_matrix)
#
#h3_histone_list = ['H3K4', 'H3K9', 'H3K18', 'H3K27', 'H3K56', 'H3K79']
#
#h3_histone_matrix[0] = h3k4
#h3_histone_matrix[1] = h3k9
#h3_histone_matrix[2] = h3k18
#h3_histone_matrix[3] = h3k27
#h3_histone_matrix[4] = h3k56
#h3_histone_matrix[5] = h3k79
#
#h3_histone_matrix = pd.DataFrame(h3_histone_matrix)
##h3_histone_matrix.insert(0, 'Histones', h3_histone_list, True)
##h3_histone_matrix=h3_histone_matrix.set_index('Histones')
#h3_histone_matrix.columns = 
#
#A = list(h3_ccle_df.columns)

"""
GRAPHING
"""

Heatmap(h3_ccle_matrix, search, (12,12), 'RdBu', 'H3 and CCLE Correlation Plot')
Heatmap(leroy_ccle_matrix, search, (12,12), 'RdBu', 'LeRoy and CCLE Correlation Plot')
#Heatmap(leroy_ccle_r1, (10,5), 'Blues', 'LeRoy and CCLE Correlation Plot with Recon1 Genes') 
#Clustermap(leroy_ccle_r1, (10,5),'Blues', method = 'single' ,metric = 'correlation')
#Heatmap(h3_ccle_oesophagus, search, (12,12), 'Blues', 'Oesophagus Correlation Data')
#Heatmap(h3_ccle_sg, search, (12,12), 'Blues', 'Salivary Gland Correlation Data')
#PlotlyHeat(leroy_ccle_r1, (1000,10000), 'LeRoy and CCLE Data with all Recon1 Genes',leroy_markers, recon1_list)

""" 
Testing
"""
os.chdir('/Users/marcdimeo/Desktop/University of Michigan Research/methylation-gem/python')

import vis
leroy_ccle_r,leroy_ccle_p = pearson_dfs(recon1_list, leroy_markers, ccle_leroy_df, leroy_ccle_df)

vis.hierarchal_clustergram(leroy_ccle_r)

"""
PLOTLY
"""


"""
SAVE AS CSV FIL
"""
#export_csv = h3_ccle_matrix.to_csv(r'/Users/marcdimeo/Desktop/University of Michigan Research/methylation-gem/data/h3_ccle_correlation_matrix.csv', index = None, header=True)
#export_csv = leroy_ccle_matrix.to_csv(r'/Users/marcdimeo/Desktop/University of Michigan Research/methylation-gem/data/leroy_ccle_correlation_matrix.csv', index = None, header=True)


