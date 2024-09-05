"""
Created on 27Aug2024:
@author: <Gila Ovadia>

This code imports genomic data (nanopore sequences after modkit procedure) and cleans it for use in modeling
"""

# region library
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pd.set_option('display.max_columns', 1000)

#endregion


def clean(df):
    # get the max value of methylation between m and h
    df = df.groupby(['read_id', 'forward_read_position', 'ref_position', 'chrom', 'mod_strand', 'ref_strand', 'ref_mod_strand',
                'fw_soft_clipped_start', 'fw_soft_clipped_end', 'read_length', 'base_qual']).agg({'mod_qual': 'max'}).reset_index()

    # merge with CpG sites
    df['ref_position_adj'] = np.where(df['ref_mod_strand'] == '+', df['ref_position'] + 1, df['ref_position'])
    df = df.merge(cpg, left_on='ref_position_adj', right_on='ref_pos', how='left')
    df.drop(columns=['chr', 'ref_pos'], inplace=True)

    return(df)

if __name__ == '__main__':

    # region import
    ## Import genomic data
    print("Importing Data")
    train1 = pd.read_csv(rf"C:\Users\gilao\Desktop\גנטי\Data\Genomic\train_test\set1\modkit\HG001_train", sep='\t')
    train2 = pd.read_csv(rf"C:\Users\gilao\Desktop\גנטי\Data\Genomic\train_test\set1\modkit\HG002_train", sep='\t')
    test1 = pd.read_csv(rf"C:\Users\gilao\Desktop\גנטי\Data\Genomic\train_test\set1\modkit\HG001_test", sep='\t')
    test2 = pd.read_csv(rf"C:\Users\gilao\Desktop\גנטי\Data\Genomic\train_test\set1\modkit\HG002_test", sep='\t')
    print("Data Imported")

    # import CpG data
    cpg = pd.read_csv(rf"C:\Users\gilao\Desktop\גנטי\Data\Genomic\Refrence Data\CpG_HG38\CpG_chr21.csv")
    cpg.rename(columns={'V1': 'chr', 'V2': 'ref_pos', 'V3': 'CpG_site'}, inplace=True)
    # endregion

    # region clean
    print("Cleaning Data")
    train1 = clean(train1)
    train2 = clean(train2)
    test1 = clean(test1)
    test2 = clean(test2)
    print("Data Clea")
    #endregion
