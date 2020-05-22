import pandas as pd
from biothings_client import get_client
import mygene
import csv
from itertools import chain


with open('./../data/CCLE/CCLEExp_geneSymbols.txt', 'r') as file:
    readLines = csv.reader(file)
    geneSymbols = list(readLines)

geneSymbols = list(chain.from_iterable(geneSymbols))

mg = get_client('gene')
output = mg.querymany(geneSymbols, scopes='symbol',
                      fields='entrezgene', species='human', as_dataframe=True)

output.to_csv('./../data/CCLE/GeneSymbol_to_Entrez.csv')
