import pandas as pd
import numpy as np

def make_list():
    all_gaia = np.array(pd.read_csv('/ar1/PROJ/fjuhsd/shared/tessodyssey/catalog-related/kerr_table_one_gaia_only.csv')['GDR2_ID']).astype(str)

    gaia_with_tic = np.array(pd.read_csv('/ar1/PROJ/fjuhsd/shared/tessodyssey/catalog-related/kerr_table_one_gaia_with_tic.csv')['GDR2_ID'])

    gaia_with_tic = np.array([str(i.split('_')[-1]) for i in gaia_with_tic])

    missing_gaia = np.setdiff1d(all_gaia, gaia_with_tic)

    df = pd.DataFrame(missing_gaia.reshape(-1,1), columns=['GDR2_ID'])

    df.to_csv('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/personal_epochs/thaddaeus/july/topical/tls_beta_run work/missing_tics/kerr_table_one_gaia_missing_tic.csv', index=False)

#make_list()

def query_some():
    from tqdm import tqdm 
    from sunnyhills.misc import gaia_to_tic
    gaia_ids = pd.read_csv('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/personal_epochs/thaddaeus/july/topical/tls_beta_run work/missing_tics/kerr_table_one_gaia_missing_tic.csv')['GDR2_ID']

    with open('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/personal_epochs/thaddaeus/july/topical/tls_beta_run work/missing_tics/kerr_table_one_gaia_missing_tic_but_found.csv', 'a') as f: 
        f.write('GDR2_ID,TIC_ID\n')    
        for index in tqdm(range(len(gaia_ids))): 
            gaia_id = gaia_ids[index]
            f.write('GDR2_'+str(gaia_id)+',TIC_'+str(gaia_to_tic([gaia_id])[1][0])+'\n')

query_some() 

# current process number: [2] 4014867
