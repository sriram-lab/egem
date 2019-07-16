#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 09:37:33 2019

@author: marcdimeo
"""

"""
This script is intended to serve as test/method to group cell lines according to their tissues and histones
which will that be used for tissue specific analysis and histone correlation
"""


import pandas as pd
import numpy as np

"""
FORMATTING DATA --> TISSUE ANALYSIS
"""

#Take cell lines and match them with tissues and average out results
#Cell line --> Tissue correspondence is found in CCLE_GCP.csv\
#CCLE_GCP is really the h3_relval data????

df = pd.read_csv(r'./../data/CCLE_GCP.csv')

temp_df = pd.DataFrame(df['CellLineName'].str.split('_',1).tolist(),columns=['CellLine', 'Tissue'])

i = 0 
for tissue in temp_df.Tissue:
    tissue = tissue.split('_')
    tissue = " ".join(tissue)
    temp_df.Tissue[i] = tissue
    i = i+1
    
del df['CellLineName']

df.insert(1, 'CellLine', temp_df.CellLine , True)
df.insert(2, 'Tissue', temp_df.Tissue, True)

#Fill NaN values
df = df.fillna(method='ffill')

#Convert values into arrays and lists
cellline_tissue  = temp_df.values

histones = list(df.columns.values)
histones = histones[3:]

expression_df = df.copy()
del expression_df['BroadID']
del expression_df['CellLine']
del expression_df['Tissue']

expression_df = expression_df.values

tissue_dict ={}
for data in cellline_tissue:
    tissue_dict[data[0]] = data[1]

"""
AVERAGING OUT EXPRESSION VALUES --> TISSUE ANALYSIS
"""

lung =[]
ovary = []
large_intestine = []
central_nervous_system = []
haematopoietic_and_lymphoid_tissue =[]
pancreas = []
upper_aerodigestive_tract = []
breast = []
prostate =[]
stomach = []
endometrium = []
bone = []
skin =[]
liver = []
fibroblast = []
soft_tissue = []
biliary_tract = []
autonomic_ganglia =[]
pleura = []
urinary_tract = []
kidney = []
oesophagus = []
thyroid = []
salivary_gland = []

i = 0
for celllines in tissue_dict:
    if tissue_dict[celllines] == "LUNG":
        lung.append(expression_df[i])
    if tissue_dict[celllines] == "OVARY":
        ovary.append(expression_df[i])
    if tissue_dict[celllines] == "LARGE INTESTINE":
        large_intestine.append(expression_df[i])
    if tissue_dict[celllines] == "CENTRAL NERVOUS SYSTEM":
        central_nervous_system.append(expression_df[i])
    if tissue_dict[celllines] == "HAEMATOPOIETIC AND LYMPHOID TISSUE":
        haematopoietic_and_lymphoid_tissue.append(expression_df[i])
    if tissue_dict[celllines] == "PANCREAS":
        pancreas.append(expression_df[i])
    if tissue_dict[celllines] == "UPPER AERODIGESTIVE TRACT":
        upper_aerodigestive_tract.append(expression_df[i])
    if tissue_dict[celllines] == "BREAST":
        breast.append(expression_df[i])
    if tissue_dict[celllines] == "PROSTATE":
        prostate.append(expression_df[i])
    if tissue_dict[celllines] == "STOMACH":
        stomach.append(expression_df[i])
    if tissue_dict[celllines] == "ENDOMETRIUM":
        endometrium.append(expression_df[i])
    if tissue_dict[celllines] == "BONE":
        bone.append(expression_df[i])
    if tissue_dict[celllines] == "SKIN":
        skin.append(expression_df[i])
    if tissue_dict[celllines] == "LIVER":
        liver.append(expression_df[i])
    if tissue_dict[celllines] == "FIBROBLAST":
        fibroblast.append(expression_df[i])
    if tissue_dict[celllines] == "SOFT TISSUE":
        soft_tissue.append(expression_df[i])
    if tissue_dict[celllines] == "BILIARY TRACT":
        biliary_tract.append(expression_df[i])
    if tissue_dict[celllines] == "AUTONOMIC GANGLIA":
        autonomic_ganglia.append(expression_df[i])
    if tissue_dict[celllines] == "PLEURA":
        pleura.append(expression_df[i])
    if tissue_dict[celllines] == "URINARY TRACT":
        urinary_tract.append(expression_df[i])
    if tissue_dict[celllines] == "KIDNEY":
        kidney.append(expression_df[i])
    if tissue_dict[celllines] == "OESOPHAGUS":
        oesophagus.append(expression_df[i])
    if tissue_dict[celllines] == "THYROID":
        thyroid.append(expression_df[i])
    if tissue_dict[celllines] == "SALIVARY GLAND":
        salivary_gland.append(expression_df[i])
    i = i + 1
   
#Convert to numpy arrays     
lung = np.array(lung)
ovary = np.array(ovary)
large_intestine = np.array(large_intestine)
central_nervous_system = np.array(central_nervous_system)
haematopoietic_and_lymphoid_tissue =np.array(haematopoietic_and_lymphoid_tissue)
pancreas = np.array(pancreas)
upper_aerodigestive_tract = np.array(upper_aerodigestive_tract)
breast = np.array(breast)
prostate =np.array(prostate)
stomach = np.array(stomach)
endometrium = np.array(endometrium)
bone = np.array(bone)
skin =np.array(skin)
liver = np.array(liver)
fibroblast = np.array(fibroblast)
soft_tissue = np.array(soft_tissue)
biliary_tract = np.array(biliary_tract)
autonomic_ganglia = np.array(autonomic_ganglia)
pleura = np.array(pleura)
urinary_tract = np.array(urinary_tract)
kidney = np.array(kidney)
oesophagus = np.array(oesophagus)
thyroid = np.array(thyroid)
salivary_gland = np.array(salivary_gland)
 
#Average
lung = np.mean(lung,axis = 0)
ovary = np.mean(ovary,axis = 0)
large_intestine = np.mean(large_intestine,axis = 0)
central_nervous_system = np.mean(central_nervous_system,axis = 0)
haematopoietic_and_lymphoid_tissue = np.mean(haematopoietic_and_lymphoid_tissue,axis = 0)
pancreas = np.mean(pancreas,axis = 0)
upper_aerodigestive_tract = np.mean(upper_aerodigestive_tract,axis = 0)
breast = np.mean(breast,axis = 0)
prostate = np.mean(prostate,axis = 0)
stomach = np.mean(stomach,axis = 0)
endometrium = np.mean(endometrium,axis = 0)
bone = np.mean(bone,axis = 0)
skin = np.mean(skin,axis = 0)
liver = np.mean(liver,axis = 0)
fibroblast = np.mean(fibroblast,axis = 0)
soft_tissue = np.mean(soft_tissue,axis = 0)
biliary_tract = np.mean(biliary_tract,axis = 0)
autonomic_ganglia = np.mean(autonomic_ganglia,axis = 0)
pleura = np.mean(pleura,axis = 0)
urinary_tract = np.mean(urinary_tract,axis = 0)
kidney = np.mean(kidney,axis = 0)
oesophagus = np.mean(oesophagus,axis = 0)
thyroid = np.mean(thyroid,axis = 0)
salivary_gland = np.mean(salivary_gland,axis = 0)  

"""
HISTONE ANALYSI
"""
h3k4 = df.iloc[:,3:7]
h3k9 = df.iloc[:,7:17]
h3k18 = df.iloc[:,17:22]
h3k27 = df.iloc[:,22:40]
h3k56 = df.iloc[:,40:42]
h3k79 = df.iloc[:,42:]

h3k4 = h3k4.values
h3k9 = h3k9.values
h3k18 = h3k18.values
h3k27 = h3k27.values
h3k56 = h3k56.values
h3k79 = h3k79.values

h3k4 = np.mean(h3k4, axis = 1)
h3k9 = np.mean(h3k9, axis = 1)
h3k18= np.mean(h3k18, axis = 1)
h3k27 = np.mean(h3k27, axis = 1)
h3k56 = np.mean(h3k56, axis = 1)
h3k79 = np.mean(h3k79, axis = 1)

