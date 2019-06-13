"""
"""

import sys
import numpy as np
import scipy
import pandas as pd

import process
import random_forest
import validator
import visualizations
import save

# Get the frequency for each label in each dataset
df, df1, header, canc, targ, data, classes, orig_data, orig_classes, excl_targ = process.preprocess(datapath=sys.argv[1], fil=sys.argv[2], targ=sys.argv[3], exclude=sys.argv[4])
