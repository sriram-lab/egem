"""
Write a function to load your data. I'll be using a random dataset with approximately the same size as the CCLE H3 relative values proteomics dataset. The histone marks and genes are the rows, while the columns are the cell line values. You can perform numpy calculations on pandas dataframes, but you will have to figure out how to do that on your own.

This formulation was adopted from https://stackoverflow.com/questions/42885239/correlation-matrix-of-two-pandas-dataframe-with-p-values?rq=1. Note that this is NOT the fastest implementation, because it takes in two for-loops, but it illustrates the logic behind our computation.

To speed up this script, you should find a way to vectorize this code: https://engineering.upside.com/a-beginners-guide-to-optimizing-pandas-code-for-speed-c09ef2c6a4d6

@author: Scott Campit
"""

import numpy as np
from scipy.stats import pearsonr

# Generate a random array. I'm making up data for 200 cell lines.
histone_marks = np.random.randn((200, 42))
gene_expression = np.random.randn((200, 1200))

def compute_pearson_values(arr1, arr2):
    """
    compute_pearson_values takes an input of two arrays and computes the Pearson correlation coefficient and p-value.

    With the two arrays, this function computes a correlation matrix. Using the correlation matrix, with each element corresponding to an R-value, it computes a p-value of the R-value with respect to a histone marker.

    This function assumes you already took the intersection between the two arrays. If you do not do that, then the matrix dimensions will not agree, resulting in an error.
    """

    # Initialize the arrays we're storing out values in to be the size of array 1 and array 2. DO NOT HARDCODE THIS - have the function be able to do this with any dataset of any size.

    # In this test example, your expected dimensions should be 42 x 1200
    Rcoef_mat = None
    pval_mat = None

    # Our objective is to compute the R-value and p-value between each histone mark / gene expression pair. Because we have `pearsonr` as a useful function that can compute the R- and p-value between two vectors, we need to make this for each column in the array.

    for hist_val in range(Rcoef_mat.shape[1])
        for gene_val in range(pval_mat.shape[1])

            # Compute correlation between each column for both arrays using the pearsonr function. This will result in a 1x2 vector after each iteration.
            pearson_product = None

            # Extract the R-value from the correlation_matrix variable (look up the documentation for the pearsonr function and what it outputs)
            Rcoef_mat[hist_val, gene_val] = None
            pval_mat[hist_val, gene_val] = None

    return Rcoef_mat, pval_mat

# Now use this with your dataset. You can modify the code above to return a dataframe with already labeled indicies and columns.
R_test, pval_test = compute_pearson_values(histone_marks, gene_expression)

# Your code should return two arrays that are 42 x 1200 if you're using the random dataset described above. Further, the R_test array should have values in the range of [-1,1], while the pval-test array should have values in the range of [0,1].

# If you're unsure if you computed the statistics correctly, consider the fact that these are random values, and therefore you have a normal Guassian-like distribution. You can plot it using pandas/numpy and matplotlib as your sanity check.
