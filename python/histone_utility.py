#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 11:28:01 2019

@author: Marc Di Meo
"""

import mygene
import pandas as pd
import numpy as np
from sklearn import preprocessing
import os
from scipy.stats.stats import pearsonr


def read_csv_df(csv_file):
    """Read in a csv file and convert it to a dataframe. 
    
    This function will remove expression values of zero and replace it with NaN
    
    INPUT
        file: path to desired file
        
    OUTPUT
        data: data with empty values as NaN    
        
    """
    data = pd.read_csv(csv_file, delim_whitespace=True)
    data = data.replace(0, np.nan)
    return data


def change_filepath(path):
    """Change file path to look for other files
    
    INPUT 
        path: new directory you would like to look for file
        
    """
    
    os.chdir(path)  

def remove_columns(data, list_columns):    
    """Remove unnecssary dataframe columns
    
    INPUT
        data: dataframe with unnecessary columns
        list_columns: list of columns to remove from the dataframe
    
    OUTPUT
        data: dataframe without the columns you wished to remove
        
    """
    
    for columns in list_columns:
        del data[columns]
    return data


def remove_duplicate_entries(items_to_remove): 
    """This function will remove duplicated elements from a list
    
    INPUT
        list_of_items: list which you would like to remove any duplicated items
    
    OUTPUT
        final_list: list with no repeating elements
        
    """
    final_list = [] 
    for num in items_to_remove: 
        if num not in items_to_remove: 
            final_list.append(num) 
    return final_list


def gene_conversion(gene_list, scopes, fields, species):
    """This serves as a way to convert gene ids into gene  
    symbols which is much easier to work with. 
    
    It is important to note that the gene_list may require manipulation 
    as mygene works differently for each scope you use. 
    
    The gene_list may require some editing.
    
    The output is a list of dictionary which contain ids and gene symbols. 
    Further manipulation is required regarding removing any duplicate IDs.
    
    INPUT
    
        gene_list: list of genes which you would like to convert into desired 
                   id or symbol
        scopes: the current format type that the gene_list is written in
        fields: the resulting gene id or symbol desired
        species: type of species you are looking at
    
    OUTPUT
        ginfo: list of dictionaries containing gene ids, 
               symbols, scores, and query
               
    """
    mg = mygene.MyGeneInfo()
    ginfo = mg.querymany(gene_list, scopes=scopes, fields=fields, 
                         species=species)

    return ginfo


def gene_symbols(ginfo_list):
    """ This function creates a list of strictly just gene symbols 
    which should be used after gene_conversation. 

    It is important to note that you pop() any duplicate IDs from ginfo so that
    the size can stay consistent and this can be placed back into the dataframe.

    INPUT
        ginfo_list: list of dictionaries containing gene ids, 
                    symbols, scores, and query
    
    OUTPUT
        gene: list of just the gene symbols from the list of dictionaries 
              inputted (any gene that could not be converted 
              will be given the symbol 'delete')
              
    """
    gene = []
    for entry in ginfo_list:
        if entry.get('symbol') != None:
            gene.append(entry['symbol'])
        else:
            gene.append("delete")
    return gene


def tissue_dict(df):
    """Using a df with number indexing and not gene symbol 
    indexing (so not pearson_dfs) this will create a dictionary that will 
    match celllines to their tissue type.
    
    This has tissue types split with an underscore so it will remove 
    the underscore. This is only really required for tissue analysis.
   
    INPUT
        df: dataframe of standard indexing, with columns being celllines and 
            corresponding tissue (for the data used in the code, the tissue 
            and cellline is seperated by an underscore)
    
    OUTPUT
        tissue_dict: dictionary of cell line and corresponding tissue
        cellline_list: list of cell lines in the dataframe
        
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
    """Takes in a dataframe with columns being cell lines 
    (with the exception of a column being genes or gene_ids) and normalizes 
    the expression values.
    
    This will delete the desired column and then normalize the values and 
    fill in the NaN values. It will then format the dataframe to back 
    to how it originally looked. 
    
    INPUT
    
        df: dataframe with columns of cell lines (with the exception of a 
            column being genes or gene_ids)
        delete_column: column that does not contain expression values
        cell_list: list of cell lines within dataframe
    
    OUTPUT
        df_normalized: dataframe of normalized expression values and 
                       NaN values filled   
                       
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
    

    # Fill in NaN values.
    df_normalized = df_normalized.fillna(method='ffill')
    return df_normalized


def common(a,b):
    """ This functions purpose is to find common elements of 
    lists by converting each list to a set.
    
    The inputs of this function should both be lists. 
    
    INPUT
        a: list which you wish to find common elements
        b: list which you wish to find common elements
        
    OUTPUT
        common: list containing common elements from a and b
        
    """
    a_set = set(a)
    b_set = set(b)
    common = a_set.intersection(b_set) 
    common = list(common)
    return common 
    
def pearson_dfs(gene_list, histone_list, data1, data2, dtype):
    """This functions purpose is to create two dataframes. 
    
    The first dataframe will contain all the r-values between the two sets of data. 
    The second will contain all the p-values for each r-value.
    
    This function is very similar to the Matrix function however it will 
    create a dataframe with columns being histones and rows being genes. 
   
    This function is very specific. The dataframe created will contain histone
    markers as the columns and gene names as the index. 
    
    This function assumes you already took the intersection between the two arrays. 
    If you do not do that, then the matrix dimensions will not agree, 
    resulting in an error.
    
    The format for inputting matrix should have columns 
    correspond to the celllines.
    
    INPUT
        gene_list: list of genes which correlation study will be performed
        histone_list: list of histones which correlation study will be performed
        data1: dataframe containing genes which you would like to compare
        data2: dataframe containing histones which you would like to compare
        dtype: numpy or df, numpy returns a matrix, and df returns a dataframe
            
    OUTPUT
    
        Rmat: R-value matrix
        Pmat: p-value matrix
        Rdf: R-value dataframe
        Pdf: p-value dataframe
        
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
        
    if dtype=='numpy':
        return Rmat, Pmat
    
    if dtype=="df"  :  
        Rdf = pd.DataFrame(Rmat)
        Rdf.columns = histone_list
        Rdf.insert(0, 'Genes', gene_list, True)
        Rdf = Rdf.set_index('Genes')
    
        Pdf = pd.DataFrame(Pmat)
        Pdf.columns = histone_list
        Pdf.insert(0, 'Genes', gene_list, True)
        Pdf = Pdf.set_index('Genes')
        
        return Rdf, Pdf
    
    
def tissue_analysis(dictionary, common_cellline, histone_list, gene_list, 
                    data1, data2, tissue_type ):
    """ This function is only necessary if you would like to perform some basic 
    tissue analysis on data. The dictionary input should be a dictionary which 
    matches cell lines to their specific tissue. 
    
    Common_cellline should be a list of the common_celllines 
    between the two data sets. 
    
    Both data inputs should be dataframes created from the pearson_dfs
    function and the name will correspond to the type of tissue wishing 
    to be analyzed. The output is a dataframe which can be graphed.
    
    This may not give you the best idea of tissue analysis as most datasets 
    vary in the amount of celllines corresponding to each tissue. 
    
    
    INPUT
        dictionary: dictionary of cell lines and corresponding tissue
        common_cellline: list of cell lines both data sets have in common
        histone_list: list of histone markers wishing to be used
        gene_list: list of genes wishing to be used
        data1: the first data set being compared
        data2: the second data set being compared
        tissue_type: type of tissue wishing to isolate
        
    OUTPUT:
        matrix: new matrix with only cell lines corresponding to tissue_type
        
    """
    
    df1 = data1.copy()
    df2 = data2.copy()
    for celllines in common_cellline:
        if dictionary[celllines] != tissue_type:
            del df1[celllines]
            del df2[celllines]
    matrix = pearson_dfs(gene_list, histone_list, df1, df2, 'numpy')[0]
    return matrix
