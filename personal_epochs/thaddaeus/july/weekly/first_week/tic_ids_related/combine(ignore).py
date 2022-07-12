def initial(): 
    import pandas as pd 
    import numpy as np

    from tqdm import tqdm 

    small_subset = pd.read_csv('routines/real/tls_beta_run/routine_ids.csv')
    routine_gaia = np.array(small_subset['GDR2_ID']).astype(str)
    all_gaia_ids = np.array(pd.read_csv('routines/real/tls_beta_run/all_table_one_ids.csv')['GDR2_ID']).astype(str)

    missing_gaia_ids = []

    for GDR2_ID in tqdm(all_gaia_ids): 
        if len(np.where(routine_gaia==GDR2_ID)[0])<1:
            missing_gaia_ids.append(GDR2_ID)

    df = pd.DataFrame(np.array(missing_gaia_ids).reshape(-1,1), 
                    columns=['GDR2_ID'])

    df.to_csv('missing_tics.csv', index=False)

def fix_redownloaded():

    import numpy as np
    import pandas as pd

    second_log = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/personal_epochs/thaddaeus/july/weekly/first_week/tic_ids_related/second_log.txt'
    second_tics = []
    second_gaias = []

    with open(second_log, 'r') as f: 
        for line in f: 
            if 'TIC_ID' in line: 
                line_list = line.split('|')
                second_tics.append(line_list[0])
                second_gaias.append(line_list[1])

    print(len(second_tics))

    second_df = pd.DataFrame(list(zip(second_gaias, second_tics)), columns=['GDR2_ID','TIC_ID'])

    first_df = pd.read_csv('/ar1/PROJ/fjuhsd/shared/tessodyssey/routines/real/tls_beta_run/routine_ids.csv')

    first_gaias = np.array(first_df['GDR2_ID'])

    print(np.intersect1d(first_gaias, np.array(second_gaias)))

fix_redownloaded()