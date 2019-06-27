#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 11:38:29 2019

@author: Marc Di Meo
"""
import numpy as np
import scipy.io as spio
import os
from scipy.stats.stats import pearsonr
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


"""
LOAD IN FILES
"""

#CCLE Data
os.chdir('/Users/marcdimeo/Desktop/University of Michigan Research/methylation-gem/matlab/vars')

ccle_genes = spio.loadmat('supplementary_software_code', squeeze_me=True)["ccleids_met"]
ccle_celllines = spio.loadmat('supplementary_software_code', squeeze_me=True)["celllinenames_ccle1"]
ccle_expression = spio.loadmat('supplementary_software_code', squeeze_me=True)["ccle_expression_metz"]

#H3 Data
os.chdir('/Users/marcdimeo/Desktop/University of Michigan Research/methylation-gem/matlab/new_var')

h3_celllines = spio.loadmat('correlation_value', squeeze_me=True)["h3_ccle_names_python"]
h3_markers =  spio.loadmat('correlation_value', squeeze_me=True)["h3_marks_python"]
h3_expression = spio.loadmat('h3_relval', squeeze_me=True)["h3_relval"]

#Common H3 and CCLE cellines
#n value
common_celllines = spio.loadmat('correlation_value', squeeze_me=True)["common_celllines"]
n = len(common_celllines)

"""
FIND REPEATED GENES
"""
#Convert to list so it is easier to work with
list_ccle_genes = np.ndarray.tolist(ccle_genes)

#Find repeats
no_repeat_genes= list(dict.fromkeys(ccle_genes))

#Create matrix with genes and eventually all their occurance positions
matrix_genes_occurence = np.transpose(np.array([no_repeat_genes, np.zeros(len(no_repeat_genes))], dtype=object))

#Find position of their occurance
def positions(list_data,gene):
    answer = [i for i,d in enumerate(list_data) if d==gene]
    return answer

i=0
for genes in no_repeat_genes:
    matrix_genes_occurence[i][1]=list(positions(list_ccle_genes, genes))
    i = i+1

        

"""
FIND AVERAGE OF REPEATED GENES
"""

expression_gene_no_repeats = (len(no_repeat_genes)),(len(ccle_celllines))
expression_gene_no_repeats = np.zeros(expression_gene_no_repeats)

j = 0 
k = 0 
for genes in no_repeat_genes:  
    
    for positions in matrix_genes_occurence[j][1]:


        if positions == (matrix_genes_occurence[j][1])[0]:
            temp_matrix = ccle_expression[positions]
        else:
            temp_matrix = np.vstack((temp_matrix,ccle_expression[positions]))
        
        if len(matrix_genes_occurence[j][1]) == 1:
            averaged = temp_matrix
        else:
            averaged = np.mean(temp_matrix,axis = 0)
        expression_gene_no_repeats[k] = averaged
    k = k+1
    j = j+1
    
"""
NEW EXPRESSION CHARTS WITH COMMON CELLLINES - DELETE UNCOMMON CELLLINES FROM EXPRESSION CHARTS
"""
#Finding the positon of cellline in CCLE data
ccle_common_celllines= np.transpose(np.array([common_celllines, np.zeros(len(common_celllines))], dtype=object))

list_ccle_common_celllines =  ccle_celllines.tolist()

i = 0
for cellline in common_celllines:
    if cellline in ccle_celllines:
        ccle_common_celllines[i][0] = cellline
        ccle_common_celllines[i][1] = list_ccle_common_celllines.index(cellline)
    i = i+1

#Finding the positon of cellline in H3 data
h3_common_celllines = np.transpose(np.array([common_celllines, np.zeros(len(common_celllines))], dtype=object))

list_h3_common_celllines = h3_celllines.tolist()

i = 0
for cellline in common_celllines:
    if cellline in h3_celllines:
        h3_common_celllines[i][0] = cellline
        h3_common_celllines[i][1] = list_h3_common_celllines.index(cellline)
    i = i+1

#New expression chart in the right order for CCLE    
ccle_common_cellline_expression = (len(common_celllines),len(no_repeat_genes))
ccle_common_cellline_expression = np.zeros(ccle_common_cellline_expression)

i=0
temp_transposed = np.transpose(expression_gene_no_repeats)
for cellline in common_celllines:
    ccle_common_cellline_expression[i] = temp_transposed[ccle_common_celllines[i][1]]
    i = i+1  

ccle_common_cellline_expression = np.transpose(ccle_common_cellline_expression)

#New expression chart in the right order for H3
h3_common_cellline_expression = (len(common_celllines),len(h3_markers),)
h3_common_cellline_expression = np.zeros(h3_common_cellline_expression)


i = 0
for cellline in common_celllines:
    h3_common_cellline_expression[i] = h3_expression[h3_common_celllines[i][1]]
    i = i+1  

h3_common_cellline_expression = np.transpose(h3_common_cellline_expression)

"""
HANDLING MISSING DATA
"""
#Replace nan with average --> this will be changed to KNN eventually
row_mean = np.nanmean(h3_common_cellline_expression, axis = 1) 
  
#Find indices where nan value is present 
inds = np.where(np.isnan(h3_common_cellline_expression)) 
  
#Replace inds with avg of column 
h3_common_cellline_expression[inds] = np.take(row_mean, inds[0])   

"""
SEARCHING FOR KEGG DATA
"""

search = "HDAC4, AHO3, BDMR, HA6116, HD4, HDAC-4, HDAC-A, HDACA"

search = search.split(", ")

for gene in search:
    if gene in no_repeat_genes:
        print(gene + ":", "YES")
    else:
        print(gene+ ":","No")
        
#found = [EHMT2, KMT5A --> SETD8, KMT2C --> MLL3, SUV39H1, SUV39H2, EHMT1, SETD7, SETDB2, DOT1L, SETDB1]
        
"""
HISTONE AVERAGING
"""

def histone_averaging(list_of_genes):
    matrix = (len(list_of_genes), len(common_celllines))
    matrix = np.zeros(matrix)
    
    i=0
    for genes in list_of_genes:
        matrix[i] = ccle_common_cellline_expression[no_repeat_genes.index(genes)]
        i = i+1
    
    result = np.mean(matrix, axis = 0)
    return result

h3k4me1_list = ['MLL3', 'SETD7']
h3k9me1k14ac0_list = ['SETDB1', 'EHMT1', 'EHMT2']
h3k9me2k14ac0_list = ['EHMT1', 'EHMT2', 'SUV39H1', 'SUV39H2']
h3k9me3k14ac0_list = ['SETDB1', 'SETDB2', 'SUV39H1', 'SUV39H2']
h3k79me1_list = ['DOT1L']
h3k79me2_list = ['DOT1L']

h3k4me1_averaged = histone_averaging(h3k4me1_list)
h3k9me1k14ac0_averaged = histone_averaging(h3k9me1k14ac0_list )
h3k9me2k14ac0_averaged = histone_averaging(h3k9me2k14ac0_list )
h3k9me3k14ac0_averaged = histone_averaging(h3k9me3k14ac0_list)
h3k79me1_averaged = histone_averaging(h3k79me1_list)
h3k79me2_averaged = histone_averaging(h3k79me2_list)


"""
REMOVING OUTLIERS
"""

#May need to use this for more of the data
#Not used as of now
#Figure out way to correspond outlier to both data sets
def reject_outliers(data, m=2):
    return data[abs(data - np.mean(data)) < m * np.std(data)]

"""
CORRELATION CALCULATION
"""

h3_markers_list =list(dict.fromkeys(h3_markers))

def correlation(x,y):
    result = pearsonr(x,y)
    result = result[0]
    return result


h3k4me1_correlation = correlation(h3k4me1_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K4me1")])
h3k9me1k14ac0_correlation = correlation(h3k9me1k14ac0_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K9me1K14ac0")])
h3k9me2k14ac0_correlation = correlation(h3k9me2k14ac0_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K9me2K14ac0")])
h3k9me3k14ac0_correlation = correlation(h3k9me3k14ac0_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K9me3K14ac0")])
h3k79me1_correlation = correlation(h3k79me1_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K79me1")])
h3k79me2_correlation = correlation(h3k79me2_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K79me2")])

#New calculation with removed outlier at row 352
h3k9me2k14ac0_averaged = np.delete(h3k9me2k14ac0_averaged, 382)
deleted_outlier = np.delete(h3_common_cellline_expression[h3_markers_list.index("H3K9me2K14ac0")], 382)

h3k9me2k14ac0_correlation = pearsonr(h3k9me2k14ac0_averaged,deleted_outlier)
h3k9me2k14ac0_correlation = h3k9me2k14ac0_correlation[0]

#Temporary individualized correlation
def individual_correlation(gene, histone):
    result = pearsonr(ccle_common_cellline_expression[no_repeat_genes.index(gene)],h3_common_cellline_expression[h3_markers_list.index(histone)])
    result = result[0]
    return result

MLL3 = individual_correlation('MLL3', 'H3K4me1')
SETD7 = individual_correlation('SETD7', 'H3K4me1')
SETDB1_w_h3k9me1 = individual_correlation('SETDB1', 'H3K9me1K14ac0')
EHMT1_w_h3k9me1 = individual_correlation('EHMT1', 'H3K9me1K14ac0')
EHMT2_w_h3k9me1 = individual_correlation('EHMT2', 'H3K9me1K14ac0')
EHMT1_w_h3k9me2 = individual_correlation('EHMT1', 'H3K9me2K14ac0')
EHMT2_w_h3k9me2 = individual_correlation('EHMT2', 'H3K9me2K14ac0')
SUV39H1_w_h3k9me2 = individual_correlation('SUV39H1', 'H3K9me2K14ac0')
SUV39H2_w_h3k9me2 = individual_correlation('SUV39H2', 'H3K9me2K14ac0')
SUV39H1_w_h3k9me2 = individual_correlation('SUV39H1', 'H3K9me2K14ac0')
SUV39H2_w_h3k9me2 = individual_correlation('SUV39H2', 'H3K9me2K14ac0')
SUV39H1_w_h3k9me3 = individual_correlation('SUV39H1', 'H3K9me3K14ac0')
SUV39H2_w_h3k9me3 = individual_correlation('SUV39H2', 'H3K9me3K14ac0')
SETDB1_w_h3k9me3 = individual_correlation('SETDB1', 'H3K9me3K14ac0')
SETDB2_w_h3k9me3 = individual_correlation('SETDB2', 'H3K9me3K14ac0')

"""
TESTING ALL POSSIBLE CORRELATIONS
"""
#for histone in h3_markers_list:
#    for gene in no_repeat_genes:
#        correlation = pearsonr(ccle_common_cellline_expression[no_repeat_genes.index(gene)],h3_common_cellline_expression[h3_markers_list.index(histone)])
#        correlation = correlation[0]
#        if correlation > 0.15:
#            print(gene, "-->", histone, ":", correlation)

"""
CONVERT ARRAYS INTO PANDAS 
"""

#For now just start with these and then convert as you go
#Eventually convert the whole thing to pandas but rightnow that seems a little more challenging than I thought
#ccle_common_cellline_expression = pd.DataFrame(ccle_common_cellline_expression)
#h3_common_cellline_expression = pd.DataFrame(h3_common_cellline_expression)

"""
CORRELATION OF LE ROY DATA AND H3
"""

os.chdir('/Users/marcdimeo/Desktop/University of Michigan Research/methylation-gem/matlab/vars')

leroy_celllines =  spio.loadmat('supplementary_software_code', squeeze_me=True)["acetlevellist"]
leroy_markers = spio.loadmat('methylation_proteomics_validation_data', squeeze_me=True)["acet_meth_list_rowlab"]
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
    
leroy_markers_list =list(dict.fromkeys(leroy_markers))
    
common_celllines_h3_leroy = list(set(leroy_celllines).intersection(h3_celllines))
common_celllines_h3_leroy.sort()

common_markers_h3_leroy = list(set(leroy_markers).intersection(h3_markers))
common_markers_h3_leroy.sort()

h3_leroy_common_marker_expression = (len(common_markers_h3_leroy), len(common_celllines_h3_leroy))
h3_leroy_common_marker_expression = np.zeros(h3_leroy_common_marker_expression)

i = 0
for markers in common_markers_h3_leroy:
    h3_leroy_common_marker_expression[i] = leroy_expression[leroy_markers_list.index(markers)]
    i = i+1   

h3_common_marker_expression = (len(common_markers_h3_leroy), 885)
h3_common_marker_expression = np.zeros(h3_common_marker_expression)

i=0
for markers in common_markers_h3_leroy:
    h3_common_marker_expression[i] = h3_common_cellline_expression[h3_markers_list.index(markers)]
    i = i+1   

h3_common_marker_expression = np.transpose(h3_common_marker_expression)

h3_common_marker = (len(common_celllines_h3_leroy), len(common_markers_h3_leroy))
h3_common_marker = np.zeros(h3_common_marker)


h3_cellline_list =list(dict.fromkeys(h3_celllines))


i = 0 
for celllines in common_celllines_h3_leroy:
    h3_common_marker[i] = h3_common_marker_expression[h3_cellline_list.index(celllines)]
    i = i+1
    
h3_common_marker = np.transpose(h3_common_marker)

#Correlation of Le Roy
leroy_correlation =[]

i = 0
for marker in h3_common_marker:
    leroy_correlation.append(correlation(h3_common_marker[i], h3_leroy_common_marker_expression[i]))
    i = i +1    
    
"""
CORRELATION BETWEEN CCLE AND LE ROY
"""

ccle_cellline_list = list(dict.fromkeys(ccle_celllines))

common_celllines_ccle_leroy = list(set(leroy_celllines).intersection(ccle_celllines))
common_celllines_ccle_leroy.sort()

ccle_leroy_expression = (len(common_celllines_ccle_leroy),len(no_repeat_genes))
ccle_leroy_expression = np.zeros(ccle_leroy_expression)

temp_flipped = np.transpose(expression_gene_no_repeats)

i = 0 
for celllines in common_celllines_ccle_leroy:
    ccle_leroy_expression[i] = temp_flipped[ccle_cellline_list.index(celllines)]
    i = i+1
    
ccle_leroy_expression = np.transpose(ccle_leroy_expression)

#def individual_correlation(gene, histone):
#    result = pearsonr(ccle_leroy_expression[no_repeat_genes.index(gene)], h3_leroy_common_marker_expression[h3_markers_list.index(histone)])
#    result = result[0]
#    return result
#
#MLL3 = individual_correlation('MLL3', 'H3K4me1')
#SETD7 = individual_correlation('SETD7', 'H3K4me1')
#SETDB1_w_h3k9me1 = individual_correlation('SETDB1', 'H3K9me1K14ac0')
#EHMT1_w_h3k9me1 = individual_correlation('EHMT1', 'H3K9me1K14ac0')
#EHMT2_w_h3k9me1 = individual_correlation('EHMT2', 'H3K9me1K14ac0')
#EHMT1_w_h3k9me2 = individual_correlation('EHMT1', 'H3K9me2K14ac0')
#EHMT2_w_h3k9me2 = individual_correlation('EHMT2', 'H3K9me2K14ac0')
#SUV39H1_w_h3k9me2 = individual_correlation('SUV39H1', 'H3K9me2K14ac0')
#SUV39H2_w_h3k9me2 = individual_correlation('SUV39H2', 'H3K9me2K14ac0')
#SUV39H1_w_h3k9me2 = individual_correlation('SUV39H1', 'H3K9me2K14ac0')
#SUV39H2_w_h3k9me2 = individual_correlation('SUV39H2', 'H3K9me2K14ac0')
#SUV39H1_w_h3k9me3 = individual_correlation('SUV39H1', 'H3K9me3K14ac0')
#SUV39H2_w_h3k9me3 = individual_correlation('SUV39H2', 'H3K9me3K14ac0')
#SETDB1_w_h3k9me3 = individual_correlation('SETDB1', 'H3K9me3K14ac0')
#SETDB2_w_h3k9me3 = individual_correlation('SETDB2', 'H3K9me3K14ac0')
    



"""
GRAPHING CORRELATION
"""

#All correlations --> Averaged
#markers = ('H3K4me1', 'H3K9me1K14ac0', 'H3K9me2K14ac0', 'H3K9me3K14ac0', 'H3K79me1', 'H3K79me2')
#y_pos = np.arange(len(markers))
#correlation = [h3k4me1_correlation, h3k9me1k14ac0_correlation, h3k9me2k14ac0_correlation, h3k9me3k14ac0_correlation, h3k79me1_correlation, h3k79me2_correlation]
#
#plt.bar(y_pos, correlation, align='center', alpha=0.5)
#plt.xticks(y_pos, markers)
#plt.xticks(rotation=45)
#plt.xlabel('Histone Markers')
#plt.ylabel('Correlation')
#plt.title('Histone Marker Correlation Values')
#
#plt.show()

#Scatter Plot
#ccle_expression_axis = h3k9me2k14ac0_averaged
#h3_expression_axis = deleted_outlier 
#
#plt.scatter(ccle_expression_axis, h3_expression_axis, label = 'Celllines')
#plt.xlabel('CCLE Expression Data')
#plt.ylabel('H3 Expression Data')
#plt.legend(loc=4)
#plt.title('H3K9me2K14ac0 Expression Comparison')
#
#plt.show()

#All correlations --> Individual
#genes = ('MLL3', 'SETD7', 'SETDB1 w h3k9me1', 'EHMT1 w h3k9me1', 'EHMT2 w h3k9me1', 'EHMT1 w h3k9me2', 'EHMT2 w h3k9me2', 'SUV39H1 w h3k9me2', 'SUV39H2 w h3k9me2', 'SETDB1 w h3k9me3', 'SETDB2 w h3k9me3', 'SUV39H1 w h3k9me3', 'SUV39H2 w h3k9me3', 'DOT1L', 'DOT1L')
#y_pos = np.arange(len(genes))
#correlation = [MLL3, SETD7,SETDB1_w_h3k9me1,EHMT1_w_h3k9me1, EHMT2_w_h3k9me1, EHMT1_w_h3k9me2, EHMT2_w_h3k9me2, SUV39H1_w_h3k9me2, SUV39H2_w_h3k9me2, SUV39H1_w_h3k9me3, SUV39H2_w_h3k9me3, SETDB1_w_h3k9me3, SETDB2_w_h3k9me3, h3k79me1_correlation, h3k79me2_correlation]
#
#plt.bar(y_pos, correlation, align='center', alpha=0.5)
#plt.xticks(y_pos, genes)
#plt.xticks(rotation=90)
#plt.xlabel('Genes w Histone Markers')
#plt.ylabel('Correlation')
#plt.title('H3 Individualized Genes and CCLE Data')
#
#plt.show()

#Le Roy et al Correlation Plot
#markers = common_markers_h3_leroy
#y_pos = np.arange(len(markers))
#correlation = leroy_correlation
#
#plt.bar(y_pos, correlation, align='center', alpha=0.5)
#plt.xticks(y_pos, markers)
#plt.xticks(rotation=90)
#plt.xlabel('Histone Markers')
#plt.ylabel('Correlation')
#plt.title('H3 Data and Le Roy et al Data')
#
#plt.show()

def plot(x,y,title, x_axis, y_axis):
    x_plot = x
    y_pos = np.arange(len(x_plot))
    y_plot = y

    plt.bar(y_pos, y_plot, align='center', alpha=0.5)
    plt.xticks(y_pos, x_plot)
    plt.xticks(rotation=90)
    plt.xlabel(x_axis)
    plt.ylabel(y_axis)
    plt.title(title)

    plt.show()
    
"""
HEAT MAP
"""

gene_list = ['MLL3', 'SETD7','SETDB1', 'EHMT1', 'EHMT2', 'SUV39H1', 'SUV39H2', 'DOT1L', 'SETDB2']






heatmap_matrix = (len(gene_list), len(h3_markers_list))
heatmap_matrix = np.zeros(heatmap_matrix)

i = 0
for gene in gene_list:
    j = 0
    for histone in h3_markers_list:
        correlation = pearsonr(ccle_common_cellline_expression[no_repeat_genes.index(gene)],h3_common_cellline_expression[h3_markers_list.index(histone)])
        correlation = correlation[0]
        heatmap_matrix[i][j] = correlation
        j = j+1
    i = i+1

def heatmap(data, row_labels, col_labels, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (N, M).
    row_labels
        A list or array of length N with the labels for the rows.
    col_labels
        A list or array of length M with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # We want to show all ticks...
    ax.set_xticks(np.arange(data.shape[1]))
    ax.set_yticks(np.arange(data.shape[0]))
    # ... and label them with the respective list entries.
    ax.set_xticklabels(col_labels, fontsize = 5)
    ax.set_yticklabels(row_labels, fontsize = 5)

    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=False, bottom=True,
                   labeltop=False, labelbottom=True)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
             rotation_mode="anchor")

    plt.figure(figsize=(100,200))
   
    return im, cbar


#fig, ax = plt.subplots()
#
#im, cbar = heatmap(heatmap_matrix, gene_list, h3_markers_list, ax=ax,
#                   cmap="YlGn", cbarlabel="Correlation")
#
#fig.tight_layout()
#plt.show()
<<<<<<< HEAD
    


=======
>>>>>>> aba65861ab151789e4308118fc0845516a72ace9


    


 