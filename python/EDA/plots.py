"""
plots.py

TODO:
  * Create boxplots for each column

@author: Scott Campit
"""

import dash
import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output
import dash_bootstrap_components as dbc

import plotly
import pandas as pd
import plotly.graph_objects as go
from sklearn import preprocessing
import numpy as np

def gcp_heatmap(df, categoricals, labels):
    """
    gcp_heatmap create the heatmap for the global chromatin profiles.

    :param df:           Pandas dataframe containing CCLE GCP profiles
    :param categoricals: Pandas dataframe containing strings of categoricals
    :param labels:       Pandas dataframe containing label encodings for categorical dataframe
    :return gcp:         Plotly graph object containing data for GCP.
    """

    # Create each plot individually
    main = go.Heatmap(z=df.values, x=df.columns, y=df.index,
                      colorscale='RdBu', reversescale=True)
    side = go.Heatmap(z=labels.values, x=list(labels.columns), y=df.index.tolist(),
                      colorscale='Portland', showscale=False, text=categoricals.values)
    gcp = plotly.subplots.make_subplots(rows=1, cols=2,
                                        shared_xaxes=True,
                                        shared_yaxes=True,
                                        column_widths=[0.10, 0.90],
                                        horizontal_spacing=0.01)
    gcp.append_trace(side, row=1, col=1)
    gcp.append_trace(main, row=1, col=2)

    #fig = go.Figure(side)
    gcp.update_xaxes(tickangle=90)
    return gcp

def gcp_boxplot(df):
    """
    gcp_boxplot creates boxplots for the CCLE GCP data
    
    """
    df = pd.melt(df)
    gcp_bxplt = go.Figure(
        go.Box(x=df['variable'], y=df['value'], boxpoints='all', jitter=0.3)
    )
    return gcp_bxplt
    
def gcp_clustergram(df):
    """

    """


    return gcp_cluster