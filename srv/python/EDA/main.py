"""
Exploratory data analysis for the CCLE dataset
TODO:
  * Hierarchical clustering for GCP
  * Hierarchical clustering for Microarray
  * Hierarchical clustering for RNASeq
  * Barplots for GCP
  * Barplots for Microarray
  * Barplots for RNASeq
  * Add interactivity
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

import plots

# Load CCLE datasets to map
gcp = '~/Data/Expression/Proteomics/GCP/GCP_single.csv'
#rnaseq = '~/Data/Expression/RNASeq/CCLE/CCLE_934_TPM_EBI.tsv'
#microarray = '~/Data/Expression/Microarray/CCLE/'

# Preprocess GCP data
gcpdf = pd.read_csv(gcp)
gcpdf = gcpdf.dropna()
gcpdf = gcpdf.set_index(['CellLine'])

labels = {
    'tissue':  gcpdf['Lineage'].unique(),
    'medium':  gcpdf['Media'].unique(),
    'culture': gcpdf['Culture'].unique(),
    'transforms': ['minmax', 'zscore']
}


# Labels for heatmap
categoricals = gcpdf[['Lineage', 'Media', 'Culture']]
gcpdf = gcpdf.drop(['Lineage', 'Media', 'Culture'], axis=1)
df = categoricals.apply(preprocessing.LabelEncoder().fit_transform)
df.index = gcpdf.index

# Create Dash application
app = dash.Dash(__name__,
                external_stylesheets=[dbc.themes.SLATE],
                meta_tags=[
                    {"name": "viewport",
                     "content": "width=device-width, initial-scale=1"}
                ])

# Update dash layout
app.layout = html.Div([
    # Heatmap
    html.Div(dcc.Graph(
        figure=plots.gcp_heatmap(gcpdf, categoricals, df),
        config=dict(responsive=True),
        style={
            'height': '100%',
        },
    )),
    html.Div(dcc.Graph(
        figure=plots.gcp_boxplot(gcpdf)
    ))
])

# Label encode categoricals
#lineage = gcpdf.pop('Lineage')
#le = preprocessing.LabelEncoder()
#le_lineage = le.fit_transform(lineage)

#culture = gcpdf.pop('Culture')
#le = preprocessing.LabelEncoder()
#le_culture = le.fit_transform(culture)

#media = gcpdf.pop('Media')
#le = preprocessing.LabelEncoder()
#le_media = le.fit_transform(media)

#position = np.zeros(le_lineage.shape)
#idx = list(range(0, le_lineage.shape[0]))

if __name__ == '__main__':
    app.run_server(debug=True)