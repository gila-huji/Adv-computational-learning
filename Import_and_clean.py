import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

pd.set_option('display.max_columns', 1000)


# region import & clean
train1 = pd.read_csv(rf"C:\Users\gilao\Desktop\גנטי\Data\Genomic\train_test\set1\modkit\HG001_train", sep='\t')
cpg = pd.read_csv(rf"C:\Users\gilao\Desktop\גנטי\Data\Genomic\Refrence Data\CpG_HG38\CpG_chr21.csv")
cpg.rename(columns={'V1': 'chr', 'V2': 'ref_pos', 'V3': 'CpG_site'}, inplace=True)

def clean(df):
    # get the max value of methylation between m and h
    df = df.groupby(['read_id', 'forward_read_position', 'ref_position', 'chrom', 'mod_strand', 'ref_strand', 'ref_mod_strand',
                'fw_soft_clipped_start', 'fw_soft_clipped_end', 'read_length', 'base_qual']).agg({'mod_qual': 'max'}).reset_index()

    # merge with CpG sites
    df['ref_position_adj'] = np.where(df['ref_mod_strand'] == '+', df['ref_position'] + 1,
                                          df['ref_position'])
    df = df.merge(cpg, left_on='ref_position_adj', right_on='ref_pos', how='left')

    #   df.index = df.read_id
 #   df = df.drop(['read_id'], axis=1)
    return(df)

train1 = clean(train1)


train1['ref_position_adj'] = np.where(train1['ref_mod_strand'] == '+', train1['ref_position'] + 1, train1['ref_position'])
train1 = train1.merge(cpg, left_on='ref_position_adj', right_on='ref_pos', how='left')
train1 = train1.merge(cpg, left_on='ref_position', right_on='ref_pos', how='left', suffixes= ['_x', '_y'])

train1['CpG_site_new'] = np.where(train1['CpG_site_x'].isna(), train1['CpG_site_y'], train1['CpG_site_x'])
train1['CpG_site_new'].isna().sum()


train1[train1['CpG_site_new'].isna()]

x = train1[lambda x: x['read_id'] == '0000057b-40df-4801-8def-7737053fcac0']


train1.shape
train1['ref_position_adj']

train1[22]

train1.head()

train1['ref_position'].nunique()
train1['ref_position'].min()
train1['ref_position'].max()


train1[lambda x: x['ref_position'] == 25007657]