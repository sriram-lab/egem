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

#search = ""
#
#search = search.split(", ")
#
#for gene in search:
#    if gene in no_repeat_genes:
#        print(gene + ":", "YES")
#    else:
#        print(gene+ ":","No")
        
#found = [EHMT2, KMT5A --> SETD8, KMT2C --> MLL3, SUV39H1, SUV39H2, EHMT1, SETD7, SETDB2, DOT1L, SETDB1]
        
"""
HISTONE AVERAGING
"""
#H3K4me1 --> MLL3 and SETD7
h3k4me1 = (2, len(common_celllines))
h3k4me1  = np.zeros(h3k4me1)

h3k4me1[0] = ccle_common_cellline_expression[no_repeat_genes.index('MLL3')]
h3k4me1[1] = ccle_common_cellline_expression[no_repeat_genes.index('SETD7')]

h3k4me1_averaged = np.mean(h3k4me1 ,axis = 0)


#H3K9me1K14ac0 --> EHMT1, EHMT2, and SETDB1
h3k9me1k14ac0 = (3,len(common_celllines))
h3k9me1k14ac0 = np.zeros(h3k9me1k14ac0)

h3k9me1k14ac0[0] = ccle_common_cellline_expression[no_repeat_genes.index('EHMT1')]
h3k9me1k14ac0[1] = ccle_common_cellline_expression[no_repeat_genes.index('EHMT2')]
h3k9me1k14ac0[2] = ccle_common_cellline_expression[no_repeat_genes.index('SETDB1')]

h3k9me1k14ac0_averaged = np.mean(h3k9me1k14ac0, axis = 0)

#H3K9me2K14ac0 --> EHMT1, EHMT2, SUV39H1 and SUV39H2
h3k9me2k14ac0 = (4,len(common_celllines))
h3k9me2k14ac0 = np.zeros(h3k9me2k14ac0)

h3k9me2k14ac0[0] = ccle_common_cellline_expression[no_repeat_genes.index('EHMT1')]
h3k9me2k14ac0[1] = ccle_common_cellline_expression[no_repeat_genes.index('EHMT2')]
h3k9me2k14ac0[2] = ccle_common_cellline_expression[no_repeat_genes.index('SUV39H1')]
h3k9me2k14ac0[3] = ccle_common_cellline_expression[no_repeat_genes.index('SUV39H2')]

h3k9me2k14ac0_averaged = np.mean(h3k9me2k14ac0, axis = 0)

#H3K9me3K14ac0 --> SETDB1, SETDB2, SUV39H1 and SUV39H2
h3k9me3k14ac0 = (4,len(common_celllines))
h3k9me3k14ac0 = np.zeros(h3k9me3k14ac0)

h3k9me3k14ac0[0] = ccle_common_cellline_expression[no_repeat_genes.index('SETDB1')]
h3k9me3k14ac0[1] = ccle_common_cellline_expression[no_repeat_genes.index('SETDB2')]
h3k9me3k14ac0[2] = ccle_common_cellline_expression[no_repeat_genes.index('SUV39H1')]
h3k9me3k14ac0[3] = ccle_common_cellline_expression[no_repeat_genes.index('SUV39H2')]

h3k9me3k14ac0_averaged = np.mean(h3k9me3k14ac0, axis = 0)

#H3K79me1 --> DOT1L

h3k79me1 = ccle_common_cellline_expression[no_repeat_genes.index('DOT1L')]

#H3K79me2 --> DOT1L

h3k79me2 = ccle_common_cellline_expression[no_repeat_genes.index('DOT1L')]

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

h3k4me1_correlation = pearsonr(h3k4me1_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K4me1")])
h3k4me1_correlation = h3k4me1_correlation[0]

h3k9me1k14ac0_correlation = pearsonr(h3k9me1k14ac0_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K9me1K14ac0")])
h3k9me1k14ac0_correlation = h3k9me1k14ac0_correlation[0]

h3k9me2k14ac0_correlation = pearsonr(h3k9me2k14ac0_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K9me2K14ac0")])
h3k9me2k14ac0_correlation = h3k9me2k14ac0_correlation[0]

h3k9me3k14ac0_correlation = pearsonr(h3k9me3k14ac0_averaged, h3_common_cellline_expression[h3_markers_list.index("H3K9me3K14ac0")])
h3k9me3k14ac0_correlation = h3k9me3k14ac0_correlation[0]

h3k79me1_correlation = pearsonr(h3k79me1, h3_common_cellline_expression[h3_markers_list.index("H3K79me1")])
h3k79me1_correlation = h3k79me1_correlation[0]

h3k79me2_correlation = pearsonr(h3k79me2, h3_common_cellline_expression[h3_markers_list.index("H3K79me2")])
h3k79me2_correlation = h3k79me2_correlation[0]

#New calculation with removed outlier at row 352
h3k9me2k14ac0_averaged = np.delete(h3k9me2k14ac0_averaged, 382)
deleted_outlier = np.delete(h3_common_cellline_expression[h3_markers_list.index("H3K9me2K14ac0")], 382)

h3k9me2k14ac0_correlation = pearsonr(h3k9me2k14ac0_averaged,deleted_outlier)
h3k9me2k14ac0_correlation = h3k9me2k14ac0_correlation[0]

"""
GRAPHING CORRELATION
"""

#All correlations
markers = ('H3K4me1', 'H3K9me1K14ac0', 'H3K9me2K14ac0', 'H3K9me3K14ac0', 'H3K79me1', 'H3K79me2')
y_pos = np.arange(len(markers))
correlation = [h3k4me1_correlation, h3k9me1k14ac0_correlation, h3k9me2k14ac0_correlation, h3k9me3k14ac0_correlation, h3k79me1_correlation, h3k79me2_correlation]

plt.bar(y_pos, correlation, align='center', alpha=0.5)
plt.xticks(y_pos, markers)
plt.xticks(rotation=45)
plt.xlabel('Histone Markers')
plt.ylabel('Correlation')
plt.title('Histone Marker Correlation Values')

plt.show()

#Scatter Plot
ccle_expression_axis = h3k9me2k14ac0_averaged
h3_expression_axis = deleted_outlier 

plt.scatter(ccle_expression_axis, h3_expression_axis, label = 'Celllines')
plt.xlabel('CCLE Expression Data')
plt.ylabel('H3 Expression Data')
plt.legend(loc=4)
plt.title('H3K9me2K14ac0 Expression Comparison')

plt.show()

#ccle_expression_axis = h3k9me3k14ac0_averaged
#h3_expression_axis = h3_common_cellline_expression[h3_markers_list.index("H3K9me3K14ac0")] 
#
#plt.scatter(ccle_expression_axis, h3_expression_axis, label = 'Celllines')
#plt.xlabel('CCLE Expression Data')
#plt.ylabel('H3 Expression Data')
#plt.legend(loc=4)
#plt.title('H3K9me3K14ac0 Expression Comparison')
#
#plt.show()


 