import pandas as pd

media = pd.read_excel('./../data/CCLE/CCLE_Summary.xlsx', sheet_name='Media')
matchedCells = pd.read_excel(
    './../data/CCLE/CCLE_Summary.xlsx', sheet_name='MatchedCells')

print(matchedCells)

matchedMedia = pd.merge(media, matchedCells, how='inner',
                        left_on='Cell Line', right_on='CellLine')

matchedMedia = matchedMedia.Culture
matchedMedia.to_csv('matchedCulture.txt', index=False)
