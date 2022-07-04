'''
TIC_ID,GDR2_ID,GDR2_RA,GDR2_DEC,TLC,mas_parrallax,distance,
gaia_g_mag,absolute_mag,bp-rp,probability_age_less_than_50,age_myr,age_err_lower,age_err_upper,
cluster_name,first_ls_period,second_ls_period,third_ls_period,
first_ls_power,second_ls_power,third_ls_power,95_fap_level,raw_counts,no_flare_counts,
BiWeight-w=0.5_counts,LombScargle-n=2_counts
'''

def update_local_catalog(tic_ids):
    
    from tqdm import tqdm 
    import pandas as pd
    import numpy as np
    
    df = ''
    
    first_ls_period = []
    second_ls_period = []
    third_ls_period = []
    first_ls_power = []
    second_ls_power = []
    third_ls_power = []
    fap_95_level = []

    raw_counts = []
    no_flare_counts = []
    wotan_counts = []
    lomb_scargle_counts = []

    for i in tqdm(df.index):  
        
