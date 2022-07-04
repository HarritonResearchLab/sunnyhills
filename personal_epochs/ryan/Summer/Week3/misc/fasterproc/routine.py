import numpy as np
import pandas as pd
import os
import time
dirindex=os.listdir()
dirindex=dirindex.remove('routine.py')
print(dirindex)
input('Correct dir listing?')
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
print('loaded function')


for i in dirindex:
    df=pd.read_csv(i)
    print('csv read')
    GAIA_index=df['GAIA']
    x=gaia_to_tic(GAIA_index)
    TIC_index.append(x)
    print('output appended and sleeping for 1 minute...')
    time.sleep(60)

print('Finishing up...')
masterdf['TIC_ID']=TIC_index
masterdf.to_csv('output')