# -*- coding: utf-8 -*-
"""
Created on Mon Dec  3 11:49:56 2018

SabioRK Specific Analysis For Single Organisms and Mammalia:
    1. single_org outputs heatmaps for the minimum, median, and maximum kinetic
    parameters associated with a specific metabolite as a function of the
    coding gene and the KEGG annotated pathway

    2. multi_org outputs box plots, where each point is either the minimum,
    median, or maximum kinetic parameter associated with a specific metabolite
    as a function of the kinetic parameter and the KEGG annotated pathway

@author: Scott Campit
"""

import pandas as pd
import numpy as np
import sys, os, mygene, plac
import seaborn as sns
import matplotlib.pyplot as plt

def dir_construct():

    # Makes new directories corresponding to each organism
    for fil in os.listdir():
        if fil.endswith('.txt'):
            name = fil.split('.')[0]
            dirs = ['kcat', 'km', 'sense']
            root_path = name+'/figures/'
            for folder in dirs:
                os.makedirs(os.path.join(root_path, folder))

#dir_construct()

def gene_symbol_map(df, organism):
    """
    """
    # mygene API will map UniProt IDs to gene symbols
    mg = mygene.MyGeneInfo()
    #df = pd.read_excel(organism)
    uniprot = df["UniprotID"]

    # mapping uniprot ID to gene symbol
    out = mg.querymany(uniprot, species=organism, scopes='uniprot', fields='symbol', as_dataframe=True)
    df = pd.merge(df, out, how='inner', left_on='UniprotID', right_index=True)

    # clean df to contain relevant info
    df = df[['symbol', 'Pathway', 'parameter.type','parameter.associatedSpecies', \
             'parameter.startValue', 'Substrate', \
             'Product']].set_index('symbol').drop_duplicates(keep='first')

    # Expand Substrate Scope
    s = df['Substrate'].str.split(';').apply(pd.Series, 1).stack()
    s.index = s.index.droplevel(-1)
    s.name = 'Substrate'

    del df['Substrate']
    df = df.join(s).drop_duplicates(keep='first')

    return df

def single_org():
    """
    single_org outputs heatmaps for the minimum, median, and maximum kinetic
    parameters associated with a specific metabolite as a function of the
    coding gene and the KEGG annotated pathway.

    """

    # mygene API will map UniProt IDs to gene symbols
    mg = mygene.MyGeneInfo()

    # Iterate with all organisms in directory
    for fil in os.listdir():
        if fil.endswith('.txt'):
            # name will be used as the output variable name
            name = fil.split('.')[0]

            df = pd.read_csv(fil, sep='\t')
            uniprot = df["UniprotID"]

            # mapping uniprot ID to gene symbol
            out = mg.querymany(uniprot, scopes='uniprot', fields='symbol', as_dataframe=True)
            df = pd.merge(df, out, how='inner', left_on='UniprotID', right_index=True)

            # clean df to contain relevant info
            df = df[['symbol', 'Pathway', 'parameter.type','parameter.associatedSpecies', \
                     'parameter.startValue', 'Substrate', \
                     'Product']].set_index('symbol').drop_duplicates(keep='first')

            # Expand Substrate Scope
            s = df['Substrate'].str.split(';').apply(pd.Series, 1).stack()
            s.index = s.index.droplevel(-1)
            s.name = 'Substrate'
            del df['Substrate']
            df = df.join(s).drop_duplicates(keep='first')

            # Extract specific kinetic parameters
            km = df[df['parameter.type'] == 'Km']
            km = km[np.isfinite(km['parameter.startValue'])].drop(columns = ['parameter.type', \
                    'Product'])
            kcat = df[df['parameter.type'] == 'kcat']
            kcat = kcat[np.isfinite(kcat['parameter.startValue'])].drop(columns = ['parameter.type', \
                    'Product'])
            sense = df[df['parameter.type'] == 'kcat/Km']
            sense = sense[np.isfinite(sense['parameter.startValue'])].drop(columns = ['parameter.type', \
                    'Product'])

            # Km heatmaps for a single organism
            query = km["Substrate"].drop_duplicates(keep='first').tolist()
            for metabolite in query:
                try:
                    # Km analysis
                    km_met = km[km['parameter.associatedSpecies'] == metabolite]
                    km_met = km_met[km_met['Substrate'] == metabolite]
                    km_met = km_met.groupby([km_met.index, 'Pathway'])['parameter.startValue'].median()
                    km_met = km_met.reset_index()
                    km_met["Km"] = -1*np.log2(km_met['parameter.startValue'])
                    km_met = km_met.pivot("symbol", "Pathway", "Km")
                    mask = km_met.isnull()

                    figs = {}
                    axs = {}

                    #for idx, plot in enumerate(km_met):
                    figs[metabolite] = plt.figure()
                    axs[metabolite] = figs[metabolite].add_subplot(111)
                    sns.heatmap(data=km_met, cmap="YlGnBu", linewidths=0.5, mask=mask)
                    plt.title(metabolite+" median Km values corresponding to each pathway in "+name)
                    plt.savefig("./"+name+"/figures/km/"+metabolite+" Median Km in "+name+".svg", dpi=600)

                except ValueError:
                    pass

            # kcat heatmaps for a single organism
            query = kcat["Substrate"].drop_duplicates(keep='first').tolist()
            for metabolite in query:
                try:
                    # kcat analysis
                    kcat_met = kcat[kcat['parameter.associatedSpecies'] == metabolite]
                    kcat_met = kcat_met[kcat_met['Substrate'] == metabolite]
                    kcat_met = kcat_met.groupby([kcat_met.index, 'Pathway'])['parameter.startValue'].median()
                    kcat_met = kcat_met.reset_index()
                    kcat_met["kcat"] = -1*np.log2(kcat_met['parameter.startValue'])
                    kcat_met = kcat_met.pivot("symbol", "Pathway", "kcat")
                    mask = kcat_met.isnull()

                    figs = {}
                    axs = {}

                    #for idx, plot in enumerate(kcat_met):
                    figs[metabolite] = plt.figure()
                    axs[metabolite] = figs[metabolite].add_subplot(111)
                    sns.heatmap(data=kcat_met, cmap="YlGnBu", linewidths=0.5, mask=mask)
                    plt.title(metabolite+" median kcat values corresponding to each pathway in "+name)
                    plt.savefig("./"+name+"/figures/kcat/"+metabolite+" Median kcat in "+name+".svg", dpi=600)

                except ValueError:
                    pass

            # kcat/Km heatmaps for a single organism
            query = sense["Substrate"].drop_duplicates(keep='first').tolist()

            for metabolite in query:
                try:
                    # kcat analysis
                    sense_met = sense[sense['parameter.associatedSpecies'] == metabolite]
                    sense_met = sense_met[sense_met['Substrate'] == metabolite]
                    sense_met = sense_met.groupby([sense_met.index, 'Pathway'])['parameter.startValue'].median()
                    sense_met = sense_met.reset_index()
                    sense_met["sense"] = -1*np.log2(sense_met['parameter.startValue'])
                    sense_met = sense_met.pivot("symbol", "Pathway", "sense")
                    mask = sense_met.isnull()

                    figs = {}
                    axs = {}

                    #for idx, plot in enumerate(sense_met):
                    figs[metabolite] = plt.figure()
                    axs[metabolite] = figs[metabolite].add_subplot(111)
                    sns.heatmap(data=sense_met, cmap="YlGnBu", linewidths=0.5, mask=mask)
                    plt.title(metabolite+" median kcat/Km values corresponding to each pathway in "+name)
                    plt.savefig("./"+name+"/figures/sense/"+metabolite+" Median kcat/Km in "+name+".svg", dpi=600)

                except ValueError:
                    pass

#single_org()

if __name__ == '__main__':
    plac.call(main)
