import numpy as np
import pandas as pd
from multiprocessing import Pool

masterdf=pd.DataFrame()
TIC_index=[]
GAIA_index=[]

def gaia_to_tic(gaia_ids):
    from astrobase.services.identifiers import gaiadr2_to_tic
    tic_ids = []
    tic_ids.append(str(gaiadr2_to_tic(gaia_ids)))
    '''
    for gaia_id in gaia_ids:
        print(gaia_id)
        try: 
            if gaia_id !=None:
                tic_ids.append(str(gaiadr2_to_tic(gaia_id)))

            else: 
                tic_ids.append('')
        except: 
            tic_ids.append('')
            continue 
    '''

    return tic_ids
print('loaded function')

df=pd.read_csv('kerr1.csv')
GAIA_index=df['GAIA']

pool=Pool(8)
TIC_index=pool.map(gaia_to_tic, GAIA_index)


print('Finishing up...')
masterdf['TIC_ID']=TIC_index
masterdf.to_csv('output')