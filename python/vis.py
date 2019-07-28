"""
Creating visualizations.
@author: Scott Campit and Marc Di Meo
"""

import numpy as np
import pandas as pd

from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns

import chart_studio.plotly as py
import plotly.graph_objs as go
import plotly.figure_factory as ff
import plotly.io as pio

from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage

#import extract



def dynamic_surf(dfs):
    """
    """

    for idx, df in enumerate(dfs):
        df = df.set_index('Unnamed: 0')
        data = [
                go.surface(
                        z=np.array(dfs[0])
                )
                ]

        layout = go.Layout(
                autosize=True,
                width=500,
                height=500,
                margin=dict(
                        l=50,
                        r=50,
                        b=50,
                        t=50
                        )
                )
        fig = go.Figure(data=data, layout=layout)
        py.iplot(fig, filename='test')


def static_surf(dfs):
    """
    #surf creates the surface plot
    """

    fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(12, 10))

    for idx, df in enumerate(dfs):
        df = df.set_index('Unnamed: 0')
        x = range(len(df.index))
        y = range(len(df.columns))
        mx, my = np.meshgrid(x, y, indexing='ij')
        z = df.round(decimals=2)
        z = z.fillna(0)
        z = z.astype(int)
        surf = ax.plot_surface(mx, my, z, cmap=cm.coolwarm,
                               linewidth=0, alpha=(idx+0.01)/5)

    ax.get_proj = lambda: np.dot(
        Axes3D.get_proj(ax), np.diag([1.0, 1.0, 1.0, 1.0]))
    ax.set_xticklabels(df.index)
    ax.set_yticklabels(df.columns)
    ax.set_xlabel('Medium components')
    #ax.set_ylabel('Demand reactions')
    #ax.set_zlabel('Metabolic Flux')
    ax.view_init(60, 35)
    fig.colorbar(surf, shrink=0.5, aspect=5)
    fig.tight_layout()
    fig.savefig('test.png')




def heatmap(data, ylabels, size, colour, title):
    """Using dataframes based on the pearson_dfs function, this function will create a simple seaborn heatmap. By using a dataframe the xlabels will be automatically created from the dataframes columns
    
    INPUT
    
        data: dataframe based on the pearson_dfs function
            
        ylabels: y tick labels, for this histone correlation this would be the list of gene names used in pearson_r calcualtion
            
        size: (width, height)
            
        colour: colour scheme used to represent data (recommend 'rdbu')
            
        title: title of the graph
    
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
    
    INPUT
        
        df: dataframe based on the pearson_dfs function
            
        size: (width, hieght)
            
        title: title of the graph
            
        xaxis: x tick labels,for this histone correlation this would be the list of histone markers used in pearson_r calcualtion
            
        yaxis: y tick labels,for this histone correlation this would be the list of gene names used in pearson_r calcualtion
        
    
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




# Create list of dataframes
#dfs = []
#sheetnam = ['0.0001', '0.001', '0.01', '0.1', '1']
#for val in sheetnam:
#    df = extract.iterred(val)
#    dfs.append(df)

#static_surf(dfs)
#dynamic_surf(dfs)


def hierarchal_clustergram(data='path/to/data/file.txt'):
    """
    """
    # Set data parameters
    df = pd.DataFrame(data)
    histone_mark_names = df.columns.tolist()
    gene_names = df.index.tolist()

    # Create upper dendrogram that maps the distances for histone markers
    fig = ff.create_dendrogram(
        df.values.transpose(), orientation='bottom', labels=histone_mark_names)
    for i in range(len(fig['data'])):
        fig['data'][i]['yaxis'] = 'y2'

    # Create a side dendrogram that computes the distance between gene symbols
    dendro_side = ff.create_dendrogram(
        df, orientation='right', labels=gene_names)
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'

    # Add the side dendrogram to the figure
    for val in dendro_side['data']:
        fig.add_trace(val)

     #Create Heatmap
    dendro_leaves = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves = list(dendro_leaves)

    dist = pdist(df.values)
    heatmap_arr = squareform(dist)
    heatmap_arr = heatmap_arr[dendro_leaves, :]
    heatmap_arr = heatmap_arr[:, dendro_leaves]

    heatmap= [
            go.Heatmap(
            x=dendro_leaves,
            y=dendro_leaves,
            z=df.values,
            colorscale='RdBu'
                )
            ]

    heatmap[0]['x'] = fig['layout']['xaxis']['tickvals']
    heatmap[0]['y'] = dendro_side['layout']['yaxis']['tickvals']

    # Add Heatmap Data to Figure
    for val in heatmap:
        fig.add_trace(df.values)

    fig.update_layout({'width': 800, 'height': 800,
                                                 'showlegend': False, 'hovermode': 'closest',
                                                 })
                              # Edit xaxis
    fig.update_layout(xaxis={'domain': [.15, 1],
                                                       'mirror': False,
                                                       'showgrid': False,
                                                       'showline': False,
                                                       'zeroline': False,
                                                       'ticks':""})
                              # Edit xaxis2
    fig.update_layout(xaxis2={'domain': [0, .15],
                                                        'mirror': False,
                                                        'showgrid': False,
                                                        'showline': False,
                                                        'zeroline': False,
                                                        'showticklabels': False,
                                                        'ticks':""})

                              # Edit yaxis
    fig.update_layout(yaxis={'domain': [0, .85],
                                                       'mirror': False,
                                                       'showgrid': False,
                                                       'showline': False,
                                                       'zeroline': False,
                                                       'showticklabels': False,
                                                       'ticks': ""
                                                       })
                              # Edit yaxis2
    fig.update_layout(yaxis2={'domain': [.825, .975],
                                                        'mirror': False,
                                                        'showgrid': False,
                                                        'showline': False,
                                                        'zeroline': False,
                                                        'showticklabels': False,
                                                        'ticks':""})
    fig.write_image('test.svg')
    fig.show()
