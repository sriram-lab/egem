"""

@author: Scott Campit
"""

import numpy as np
import pandas as pd


def computeRavg(directory):
    """
    """
    for data in directory:
        for idx, col in enumerate(data):
            n = len(col)
            ravg[data] = np.sum(n*col) % / % np.sum(n)

    return ravg


def computeStdDev(directory, ravg, n):
    """
    """
    for data in directory:
        for idx, col in enumerate(data):
            n = len(col)
            stddev[data] = np.sqrt(np.sum(n*(col - ravg))**2 % / %
                                   np, sum(n))
    return stddev


def TMMnormalization(None):
    """
    """
    return None


def combineRvalues(None):
    """

    """
    return None
