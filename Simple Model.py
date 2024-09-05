"""
Created on 29 Aug 2024
@author: Gila Ovadia

Apply classification model to data.
XGBoost classifier
"""

import Import_and_clean
import pandas as pd
import numpy as np

# region transpose

def transpose(df):
    df = df[['read_id', 'CpG_site', 'mod_qual']]
    df = df[lambda x: ~x['CpG_site'.isna()]]
    df = df.pivot(index='read_id', columns='CpG_site', values='mod_qual').reset_index().add_prefix('CpG_')

    return df

# transpose the data
print("Transpose")
train1 = transpose(train1)
train2 = transpose(train2)
test1 = transpose(test1)
test2 = transpose(test2)
print("Transposed")

