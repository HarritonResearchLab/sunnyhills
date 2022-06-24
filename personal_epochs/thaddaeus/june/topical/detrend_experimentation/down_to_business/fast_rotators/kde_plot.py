import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns

def make_kde_plot(): 
    
    fast_df = pd.read_csv('./personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/bad_detrend_catalog.csv') 
    fast_df = fast_df.drop(columns=['complex fast rotator', 'long gap', 'weird', 'comment on'])
    full_df = pd.read_csv('./routines/alpha_tls/data/download_log.csv')

    fast_df = fast_df.merge(full_df, left_on='clean fast rotator', right_on='TIC_ID')

    fast_df = fast_df['top_ls_power']
    full_df = full_df['top_ls_power']

    fig, ax = plt.subplots()
    sns.kdeplot(data=full_df, axes=ax, label='Full')
    sns.kdeplot(data=fast_df, axes=ax, label='Fast Subset')
    ax.legend()
    plt.savefig('./personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/fast_rotators/periods_dist.png', dpi=200)

make_kde_plot()