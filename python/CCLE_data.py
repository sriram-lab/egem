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

import GEOparse
import mygene
<<<<<<< HEAD
import pandas as pd

"""
GETTING GENE NAMES
""" 
=======

"""
GETTING GENE NAMES
"""
>>>>>>> aba65861ab151789e4308118fc0845516a72ace9

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

<<<<<<< HEAD
#gse = GEOparse.get_GEO(filepath="./GSE36133_family.soft.gz")
gse = GEOparse.get_GEO(geo="GSE36133", destdir="./")
=======
gse = GEOparse.get_GEO(filepath="./GSE36133_family.soft.gz")

>>>>>>> aba65861ab151789e4308118fc0845516a72ace9

