"""
Creating visualizations.
@author: Scott Campit
"""

import numpy as np
import pandas as pd

from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
from matplotlib import cm
import seaborn as sns

import plotly.plotly as py
import plotly.graph_objs as go

import extract


def dynamic_surf(dfs):
    """
    """

    #for idx, df in enumerate(dfs):
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


# Create list of dataframes
dfs = []
sheetnam = ['0.0001', '0.001', '0.01', '0.1', '1']
for val in sheetnam:
    df = extract.iterred(val)
    dfs.append(df)

#static_surf(dfs)
dynamic_surf(dfs)
