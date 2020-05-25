"""
"""

import numpy as np
import pandas as pd
from openpyxl import Workbook
from openpyxl import load_workbook
from itertools import chain
import string
import logging

def parse_pdf(fileName):
    """
    """
    errors = []
    if fileName.split('.')[-1] is not '.pdf':
        raise Exception as e:
        None
    else:
        None

    return medium_df

def parse_html(fileName):
    """
    """
    return None

def parse_medium_formulation(fileName):
    """
    """
    errors = []
    try:
        if fileName.split('.')[-1] is '.pdf':
            None
        elif fileName.split('.')[-1] is '.html':
            None
    except Exception as e:
        errors.append([fileName, e])
        logging.error("Exception occurred", exc_info=True)
        continue

    return None
