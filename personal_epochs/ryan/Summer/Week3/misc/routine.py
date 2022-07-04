import numpy as np
import pandas as pd
masterdf=pd.DataFrame()
TIC_index=[]
def gaia_to_tic(gaia_ids):
    from astrobase.services.identifiers import gaiadr2_to_tic
    tic_ids = []
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


    return dict(zip(gaia_ids, tic_ids)), tic_ids



df=pd.read_csv('kerr1.csv')
GAIA_index=df['GAIA']
print(str(GAIA_index[0]))
#x=gaia_to_tic(GAIA_index[0])