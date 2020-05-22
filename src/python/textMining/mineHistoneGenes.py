"""
"""

from urllib.request import urlopen
from bs4 import BeautifulSoup
import re
import pandas as pd

pages = set()


def getGenesFromGO():
    global pages
    html = urlopen("http://epifactors.autosome.ru/genes")
    bso = BeautifulSoup(html, features="html.parser")

    table = bso.find('table')
    rows = table.find_all('tr')

    columnNames = ["HGNC_Symbol", "HGNC_ID", "Name", "UniProt",
                   "Function", "Complex", "Substrate", "Product"]
    for tr in rows:
        td = tr.find_all('td')
        row = [i.text for i in td]
        row = [i.replace('\n', '') for i in row]
        row = [i.replace('(details)', '') for i in row]
        del row[3, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15]
        print(row)
    #for link in bso.find_all('td'):
    #    print(link)


#getGenesFromGO()


def pdHTMLdf(html):
    array = pd.read_html(html, flavor='bs4')
    df = array[0]
    return df


#html = "http://epifactors.autosome.ru/genes"
#df = pdHTMLdf(html)
#df.to_csv('epifactors.csv')

#df = pd.read_csv('./epifactors.csv')
#df = df[~df['Product'].str.contains('#')]
#df.to_csv('cleanEpifactors.csv', index=False)
