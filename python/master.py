#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 11:28:01 2019

@author: marcdimeo
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
import plotly.figure_factory as ff
from scipy.spatial.distance import pdist, squareform
#import plotly.graph_objects as go

"""
FUNCTIONS
"""

def read_csv_df(file):
    """Read in a csv file and convert it to a dataframe. 
    
    This function will remove expression values of zero and replace it with NaN
    """
    data = pd.read_csv(file, delim_whitespace=True)
    data = data.replace(0, np.nan)
    return data

def change_filepath(path):
    """Change file path to look for other files
    """
    os.chdir(path)    

def column_remove(df, list_columns):    
    """Remove any unneeded columns of dataframes being created by the read_csv_df function
    """
    
    for columns in list_columns:
        del df[columns]
    return df

def remove(duplicate): 
    """This function will remove duplicated elements from a list
    """
    final_list = [] 
    for num in duplicate: 
        if num not in final_list: 
            final_list.append(num) 
    return final_list

def gene_converstion(gene_list, scopes, fields, species):
    """This serves as a way to convert gene ids into gene symbols which is much easier to work with. 
    
    It is important to note that the gene_list may require manipulation as mygene works differently for each scope you use. The gene_list may require some editting.
    
    The output is a list of dictionary which contain ids and gene symbols. Further manipulation is required regarding removing any duplpicate ids.
    """
    mg = mygene.MyGeneInfo()
    ginfo = mg.querymany(gene_list, scopes=scopes, fields=fields, species=species)

    return ginfo

def gene_symbols(ginfo_list):
    """ This function creates a list of strictly just gene symbols which should be used after gene_conversation. 

    It is important to note that you pop() any duplicate ids from ginfo so that the size can stay consistent and this can be placed back into the dataframe.
    
    """
    gene = []
    for entry in ginfo_list:
        if entry.get('symbol') != None:
            gene.append(entry['symbol'])
        else:
            gene.append("delete")
    return gene

def tissue_dict(df):
    """Using a df with number indexing and not gene symbol indexing (so not pearson_dfs) this will create a dictionary that will match celllines to their tissue type.
    
    This has tissue types split with an underscore so it will remove the underscore. This is only really required for tissue analysis.
    """
    celllines = list(df)
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
    cellline_tissue = cell_df.values
    tissue_dict ={}
    for data in cellline_tissue:
       tissue_dict[data[0]] = data[1]
    return tissue_dict, cellline_list

def normalize(df, delete_column, cell_list):
    """Takes in a dataframe with columns being celllines (with the exception of a column being genes or gene_ids) and normalizes the expression values.
    
    This will delete the desired column and then normalize the values and fill in the NaN values. It will then format the dataframe to back to how it originally looked. 
    """

    temp_df = df.copy()
    
    delete = temp_df[delete_column]
    del temp_df[delete_column]
    
    x = temp_df.values.astype(float)
    min_max_scaler = preprocessing.MinMaxScaler()
    x_scaled = min_max_scaler.fit_transform(x)
    df_normalized = pd.DataFrame(x_scaled)

    df_normalized.insert(0, 'Gene', delete , True)
    df_normalized.columns = ['Gene'] + cell_list
    

    #fill in NaN values
    df_normalized = df_normalized.fillna(method='ffill')
    return df_normalized

def common(a,b):
    """ This functions purpose is to find common elements of lists by converting each list to a set.
    
    The inputs of this function should both be lists. 
    """
    a_set = set(a)
    b_set = set(b)
    common = a_set.intersection(b_set) 
    common = list(common)
    return common 

#Used to create correlation matrix
def matrix(gene_list, histone_list, data1, data2):
    """This functions purpose is to create two numpy matrices. The first value is a matrix of the R-values between the two data sets, and the second matrix is the corresponding p-values. 
    
    The gene-list and histone-list allow you to input big expression value data sets and then specific what genes and histones you want to use. 

    This function is very specific. The dataframe created will contain histone markers as the columns and gene names as the index. 
    
    This function assumes you already took the intersection between the two arrays. If you do not do that, then the matrix dimensions will not agree, resulting in an error.
    
    The format for inputting matrix should have columns correspond to the celllines.
    
    """
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
    
def pearson_dfs(gene_list, histone_list, data1, data2):
    """This functions purpose is to create two dataframes. The first dataframe will contain all the r-values between the two sets of data. The second will contain all the p-values for each r-value.
    
    This function is very similar to the Matrix function however it will create a dataframe with columns being histones and rows being genes. 
   
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


def heatmap(data, ylabels, size, colour, title):
    """Using dataframes based on the pearson_dfs fucntion, this function will create a simple seaborn heatmap. By using a dataframe the xlabels will be automatically created from the dataframes columns
    """
    plt.figure(figsize=size)
    plt.xlabel('Histone Markers')
    plt.title(title)
    ax = sns.heatmap(data, cmap = colour, square = True, yticklabels = ylabels)
    ax.set(xlabel = 'Histone Markers')
    plt.show



    
def plotlyheat(df, size, title, xaxis, yaxis):
    """ This function will take in a dataframe and create a plotly heatmap based on that dataframe. This is a good tool for gene lists that are very big and are hard to visulize on seaborn. 
    
    The dataframe used should be one made with the pearson_dfs. The size arugment is only necessary for really big gene lists. If the list is small the graph will be created using autosizing. 
    
    The resulting image will open up in google chrome. 
    """
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
    
def tissue_analysis(dictionary, common_cellline, histone_list, gene_list, data1, data2, name ):
    """ This function is only necessary if you would like to perform some basic tissue analysis on data. The dictionary input should be a dictionary which matches cell lines to their specific tissue. 
    
    Common_cellline should be a list of the common_celllines between the two data sets. 
    
    Both data inputs should be dataframes created from the pearson_dfs function and the name will correspond to the type of tissue wishing to be analyzed. The output is a dataframe which can be graphed/
    
    This may not give you the best idea of tissue analysis as most datasets vary in the amount of celllines corresponding to each tissue. 
    """
    
    df1 = data1.copy()
    df2 = data2.copy()
    for celllines in common_cellline:
        if dictionary[celllines] != name:
            del df1[celllines]
            del df2[celllines]
    matrix = pearson_dfs(gene_list, histone_list, df1, df2)[0]
    return matrix

"""
CODE TO GET CCLE NORMALIZED
"""

df1 = read_csv_df("./../../../CCLE_RNAseq_rsem_genes_tpm_20180929.txt")
df1 = column_remove(df1, ['transcript_ids'])

gene_id_list = []

for gene_id in df1['gene_id']:
    gene_id = gene_id.split('.')
    gene_id_list.append(gene_id[0])
    
df1['gene_id'] = gene_id_list

ginfo = gene_converstion(df1['gene_id'], 'ensembl.gene', 'symbol', 'human')

#Remove any repeated gene_ids from ginfo. In this case there is a repeat at 30544

ginfo.pop(30544)

gene = gene_symbols(ginfo)

#Format the my df to my liking based on the gene names
df1['gene_id'] = gene
df1.rename(columns={'gene_id' : 'Gene'}, inplace = True)

#Retrieve the tissue dictionary in case you need it.
tissue_dictionary, cellline_list = tissue_dict(df1)

#Normalize the expression results
df1_normalized = normalize(df1, 'Gene', cellline_list)

"""
CODE TO GET H3 NORMALIZED
"""

change_filepath('/Users/marcdimeo/Desktop/University of Michigan Research/methylation-gem/matlab/new_var')

df2_celllines = np.array(spio.loadmat('correlation_value', squeeze_me=True)["h3_ccle_names_python"])
df2_celllines = df2_celllines.tolist()
df2_markers =   list(dict.fromkeys(spio.loadmat('correlation_value', squeeze_me=True)["h3_marks_python"]))
df2_expression = spio.loadmat('h3_relval', squeeze_me=True)["h3_relval"]

df2 = pd.DataFrame(df2_expression)
df2 = df2.T
df2.insert(0, 'Histone', df2_markers, True)
df2.columns = ['Histone'] + df2_celllines

df2_normalized = normalize(df2, 'Histone', df2_celllines)





