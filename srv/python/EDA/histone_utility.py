#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 11:28:01 2019

@author: Marc Di Meo & Scott Campit
"""

import mygene
import pandas as pd
import numpy as np
from sklearn import preprocessing
import os
from scipy.stats.stats import pearsonr


def gene_symbols(ginfo_list):
    """ This function creates a list of strictly just gene symbols
    which should be used after gene_conversation.

    It is important to note that you pop() any duplicate IDs from ginfo so that
    the size can stay consistent and this can be placed back into the dataframe.

    params:
        ginfo_list: list of dictionaries containing gene ids, symbols, scores, and query

    return:
        gene: list of just the gene symbols from the list of dictionaries inputted (any gene that could not be
              converted will be given the symbol 'delete')

    """
    gene = []
    for entry in ginfo_list:
        if entry.get('symbol') != None:
            gene.append(entry['symbol'])
        else:
            gene.append("delete")
    return gene


def tissue_dict(df):
    """Using a dataframe with standard indexing, produce a dictionary that will
    match celllines to their tissue type.

    In the data being used for this function the cell line is written in the
    following way: 'cellline_tissue_type'. Therefore this function will split
    the cell line up, removing the underscores from the cell line and creating
    a dictionary the following way: 'cellline : tissue type'.
    INPUT
        df: dataframe of standard indexing, with columns being celllines and
            corresponding tissue (for the data used in the code, the tissue
            and cellline is seperated by an underscore)

    OUTPUT
        tissue_dict: dictionary of cell line and corresponding tissue
        cellline_list: list of cell lines in the dataframe
        cell_df: dataframe of cell line and corresponding tissue

    """
    celllines = list(df)
    celllines.pop(0)

    cellline_list = []
    tissue_list = []
    for cell in celllines:
        temp = cell.split('_')
        cell = temp[0]
        cellline_list.append(cell)
        tissue = temp[1:]
        tissue = " ".join(tissue)
        tissue_list.append(tissue)

    d = {'Celllines': cellline_list, 'Tissue': tissue_list}
    cell_df = pd.DataFrame(d)
    cellline_tissue = cell_df.values
    tissue_dict = {}
    for data in cellline_tissue:
        tissue_dict[data[0]] = data[1]
    return tissue_dict, cellline_list, cell_df


def normalize(df, delete_column, cell_list):
    """Takes in a dataframe with columns being cell lines
    (with the exception of a column being genes or gene_ids) and min-max
    normalizes the expression values.The expected output is an array of
    expression values between 0 and 1.

    This will remove the desired column which does not contain expression
    values (Histones, Genes, etc) and then normalize the values,
    filling in the NaN values. It will then output a dataframe of the same
    size of the dataframe which was used in the function with the deleted
    column inserted back into the dataframe.

    INPUT

        df: dataframe with the columns representing cell lines and one column which
            contains gene names/histone markers
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

    df_normalized.insert(0, delete_column, delete, True)
    df_normalized.columns = [delete_column] + cell_list

    # Fill in NaN values.
    df_normalized = df_normalized.fillna(method='ffill')
    return df_normalized


def common(list_a, list_b):
    """ This functions purpose is to find common elements of
    lists by converting each list to a set.

    The inputs of this function should both be lists.

    INPUT
        a: list which you wish to find common elements
        b: list which you wish to find common elements

    OUTPUT
        common: list containing common elements from a and b

    """
    a_set = set(list_a)
    b_set = set(list_b)
    common = a_set.intersection(b_set)
    common = list(common)
    return common


def pearson_dfs(gene_list, histone_list, gene_data, histone_data, dtype):
    """This functions purpose  is to compute the column-wise pearson
    correlation coefficient and p-value for each element in an array

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
        gene_data: dataframe containing genes which you would like to compare
        histone_data: dataframe containing histones which you would like to compare
        dtype: 'numpy' or 'df', numpy returns a matrix, and df returns a dataframe

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

    gene_data = gene_data.fillna(method='ffill')
    histone_data = histone_data.fillna(method='ffill')
    for gene in gene_list:
        j = 0
        for histone in histone_list:
            #            gene_data.loc[gene] = gene_data.loc[gene].fillna(method='ffill')
            #            histone_data.loc[histone] = histone_data.loc[histone].fillna(method='ffill')
            correlation = pearsonr(
                gene_data.loc[gene], histone_data.loc[histone])
            Rmat[i][j] = correlation[0]
            Pmat[i][j]
            j = j+1
        i = i+1

    if dtype == 'numpy':
        return Rmat, Pmat

    if dtype == "df":
        Rdf = pd.DataFrame(Rmat)
        Rdf.columns = histone_list
        Rdf.insert(0, 'Genes', gene_list, True)
        Rdf = Rdf.set_index('Genes')

        Pdf = pd.DataFrame(Pmat)
        Pdf.columns = histone_list
        Pdf.insert(0, 'Genes', gene_list, True)
        Pdf = Pdf.set_index('Genes')

        return Rdf, Pdf


def tissue_analysis(cell_df, common_cellline, histone_list, gene_list,
                    data1, data2, tissue_type):
    """ This function is only necessary if you would like to perform some basic
    tissue analysis on data. The dictionary input should be a dictionary which
    matches cell lines to their specific tissue.

    Common_cellline should be a list of the common_celllines
    between the two data sets.

    Both data inputs should be dataframes created from the pearson_dfs
    function and the tissue_type will correspond to the type of tissue wishing
    to be analyzed. The output is a dataframe which can be graphed.

    This may not give you the best idea of tissue analysis as most datasets
    vary in the amount of celllines corresponding to each tissue.


    INPUT
        cell_df: dataframe created from tissue_dict function
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

    cell_df = cell_df.T
    tissuelist = cell_df.loc['Tissue'].tolist()
    cell_df.columns = tissuelist
    cell_df = cell_df.drop(['Tissue'])
    tissuelist = list(dict.fromkeys(tissuelist))
    tissuelist.remove(tissue_type)
    tissuedf = cell_df.drop(tissuelist, axis=1)
    celllinelist = tissuedf.loc['Celllines'].tolist()
    celllinelist = common(celllinelist, common_cellline)

    tissue_matrix1 = pd.DataFrame(df1[celllinelist])
    tissue_matrix2 = pd.DataFrame(df2[celllinelist])

    matrix = pearson_dfs(gene_list, histone_list,
                         tissue_matrix1, tissue_matrix2, 'numpy')[0]
    return matrix
