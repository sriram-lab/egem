"""
cluster.py
@author: Scott Campit
"""

import pandas as pd

import dash
import dash_bio as dashbio
import dash_html_components as html
import dash_core_components as dcc
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output

app = dash.Dash(__name__,
    external_stylesheets=[dbc.themes.BOOTSTRAP]
)

# Meta data
#metaFile = '~/Data/Mappings/CCLE/CCLE_Summary.xlsx'
#meta = pd.read_excel(metaFile, sheet_name='MatchedCells')
#geneID = list(meta['CellLine'])
#tissue = list(meta['Lineage'])
#medium = list(meta['Media'])
#culture = list(meta['Culture'])

# CCLE Expression filepaths
muArrayFile = '~/Data/Expression/Microarray/CCLE/CCLE_Microarray.xlsx'   # Z-score
gcpFile = '~/Data/Expression/Proteomics/GCP/GCP_proteomics_remapped.csv' # Z-score
proteomicsFile = None #'~/Data/Expression/Proteomics/CCLE/'
rppaFile = None
rnaseqFile = None
muRNAFile = None
methylationFile = None

# Expression dataframes
#microarray = pd.read_excel(muArrayFile)
#gcp = pd.read_csv(gcpFile)
proteomics = None
rppa = None
rnaseq = None
muRNA = None
methylation = None

opts = ["Microarray", "Global Chromatin Profiles"]
app.layout = html.Div([
    html.H1("CCLE Dataset"),
    dcc.Dropdown(
        id='clustergram-data',
        options=[
            {'label':o, 'value':o} for o in opts
        ],
        value="Microarray"
    ),
    html.Div(
        id='clustergram'
    )
])

@app.callback(
    Output('clustergram', 'children'),
    [Input('clustergram-data', 'value')]
)
def update_clustergram(type="Microarray"):
    """

    """
    if type is "Microarray":
        df = pd.read_excel(muArrayFile, index_col=0)
        df = df.fillna(0)
    else:
        df = pd.read_csv(gcpFile, index_col=0)
        df = df.fillna(0)
    ccl = df.index.to_list()
    feat = list(df.columns.values)
    return {
            dcc.Graph(
                figure=dashbio.Clustergram(
                    data=df.values
                    #row_labels=ccl,
                    #column_labels=feat,
                    #height=800,
                    #width=700
                    #col_group_marker=[
                    #    {'group': 1, 'annotation': 'largest column cluster', 'color': '#EF553B'}
                    #],
                    #row_group_marker=[
                    #    {'group': 2, 'annotation': 'cluster 2', 'color': '#AB63FA'},
                    #    {'group': 1, 'annotation': '', 'color': '#19D3F3'}
                    #]
                )
        )
    }

if __name__ == '__main__':
    app.run_server(debug=True)