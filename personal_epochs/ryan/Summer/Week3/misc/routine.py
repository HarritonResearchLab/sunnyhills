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


    return tic_ids



df=pd.read_csv('dummy.csv')
GAIA_index=df['GAIA']
x=gaia_to_tic(GAIA_index)
print(x)