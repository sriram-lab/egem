"""

"""

dbc.Row([
    dbc.Col(
        html.Div(
            dcc.Dropdown(
                id='tissue-dropdown',
                options=[{'label': t} for t in labels.get('tissue')],
                value=labels.get('tissue')[0]
            )
        )
    ),
    dbc.Col(
        html.Div(
            dcc.Dropdown(
                id='medium-dropdown',
                options=[{'label': m} for m in labels.get('medium')],
                value=labels.get('medium')[0]
            )
        )
    ),
    dbc.Col(
        html.Div(
            dcc.Dropdown(
                id='morphology-dropdown',
                options=[{'label': c} for c in labels.get('culture')],
                value=labels.get('culture')[0]
            )
        )
    ),
    dbc.Col(
        html.Div(
            dcc.Dropdown(
                id='transformation-dropdown',
                options=[{'label': r} for r in labels.get('culture')],
                value=labels.get('culture')[1]
            )
        )
    )
]),