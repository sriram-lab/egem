"""

"""
import re
import os
import numpy as np
import pandas as pd
from openpyxl import load_workbook

media = pd.read_excel(r'./../data/summary.xlsx', sheet_name='Cell Line Annotations', usecols = ['CCLE_ID', 'Growth.Medium', ''])
