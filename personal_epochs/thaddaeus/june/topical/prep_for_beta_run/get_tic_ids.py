def get_em(routine_kerr_key:str='routines/real/beta_tls/routine_kerr_key.csv'): 
    import pandas as pd
    import numpy as np
    from astrobase.services.identifiers import gaiadr2_to_tic
    from tqdm import tqdm 

    key = pd.read_csv(routine_kerr_key)

    tic_ids = []
    gaia_ids = np.array(key['Gaia_DR2_ID'])

    for gaia_id in tqdm(gaia_ids):
        try: 
            tic_ids.append('TIC_'+gaiadr2_to_tic(str(gaia_id), verbose=False))

            r'''
            
            [I 220630 11:48:09 mast:329] getting cached MAST query result for request: {'format': 'json', 'params': {'ra': 82.60814975935784, 'dec': 28.23914413244664, 'radius': 0.008333333333333333}, 'service': 'Mast.Catalogs.Tic.Cone', 'timeout': 10.0}
            
            '''

        except: 
            tic_ids.append('')
    
    key['TIC_ID'] = tic_ids

    key.to_csv(routine_kerr_key, index=False)

get_em()
