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
import plotly.express as px
import re
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
    gcp.update_layout(
        title_text="CCLE Global Chromatin Profiles",
        #xaxis_title="Histone markers",
        yaxis_title="Z-score distribution",
        font=dict(
            family="Courier New, monospace",
            size=18,
            color="#7f7f7f"
        )
    )
    gcp.layout.template = 'plotly_dark'

    return gcp

def gcp_boxplot(df):
    """
    gcp_boxplot creates boxplots for the CCLE GCP data
    
    """
    # Format to two columns: variable name and value associated with variable
    df = pd.melt(df)

    # Get groups of histone markers for colors
    pattern = '[a-z]+[a-z]+\d'
    df['Markers'] = df['variable'].replace(to_replace=pattern,
                                          value='',
                                          regex=True)
    # Make box plot
    gcp_bxplt = go.Figure(
        px.box(df,
               x='variable',
               y='value',
               points=False,
               color='Markers',
               color_discrete_sequence=px.colors.sequential.Agsunset,
               template='plotly_dark')
    )

    # Add categorical scatterplot / stripplot
    gcp_bxplt.add_trace(px.strip(df,
                                 x='variable',
                                 y='value',
                                 color='Markers').data[0])

    gcp_bxplt.update_traces(marker=dict(
                                        size=5,
                                        opacity=0.5
                            )
    )
    gcp_bxplt.update_xaxes(tickangle=-45)

    gcp_bxplt.update_layout(
        #title_text="CCLE Global Chromatin Profiles",
        xaxis_title="Histone markers",
        yaxis_title="Z-score distribution",
        font=dict(
            family="Courier New, monospace",
            size=18,
            color="#7f7f7f"
        )
    )

    return gcp_bxplt

def gcp_clustergram(df):
    """

    """


    return gcp_cluster