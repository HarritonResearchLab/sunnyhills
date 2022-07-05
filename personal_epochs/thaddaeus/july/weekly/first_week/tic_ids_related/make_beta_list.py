import pandas as pd
import numpy as np


def save_gaia_for_routine():
    df = pd.read_csv('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/data/archive/full_key.csv')

    probs = np.array(df['probability'])

    for i in [.68, .75, .90, .95, .99]: 
        print(len(np.where(probs>=i)[0]))

    r'''
    10215
    9385
    7454
    6564
    5419
    '''

    gaia_ids = np.array(df['GDR2_ID'])
    routine_ids = gaia_ids[np.where(probs>=.68)[0]]

    df = pd.DataFrame(routine_ids.reshape(-1,1), columns=['GDR2_ID'])
    print(df)
    df.to_csv('/ar1/PROJ/fjuhsd/shared/tessodyssey/routines/real/tls_beta_run/routine_ids.csv', index=False)

save_gaia_for_routine()

def make_base_catalog(): 
    '''
    TIC_ID,GDR2_ID,GDR2_RA,GDR2_DEC,TLC,mas_parrallax,distance,
    gaia_g_mag,absolute_mag,bp-rp,probability_age_less_than_50,age_myr,age_err_lower,age_err_upper,
    first_ls_period,second_ls_period,third_ls_period,
    first_ls_power,second_ls_power,third_ls_power,95_fap_level,raw_counts,no_flare_counts,
    BiWeight-w=0.5_counts,LombScargle-n=2_counts
    '''

    from sunnyhills.misc import return_kerr_cluster
    from tqdm import tqdm 

    df = pd.read_csv('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/data/archive/full_key.csv')

    import shutil 

    shutil.copyfile('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/data/archive/full_key.csv', 
                    '/ar1/PROJ/fjuhsd/shared/tessodyssey/data/full_base_catalog_key.csv')

#make_base_catalog()