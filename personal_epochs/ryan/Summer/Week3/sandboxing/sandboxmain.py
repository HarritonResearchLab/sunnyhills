import numpy as np
import pandas as pd
from multiprocessing import Pool

masterdf=pd.DataFrame()
TIC_index=[]
GAIA_index=[]

def gaia_to_tic(gaia_ids):
    from astrobase.services.identifiers import gaiadr2_to_tic
    tic_ids = []
    tic_ids.append(gaiadr2_to_tic(gaia_ids))
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

pool=Pool(16)
TIC_index=pool.map(gaia_to_tic, GAIA_index)


print('Finishing up...')
masterdf['TIC_ID']=TIC_index
masterdf.to_csv('output')


'''
multiprocessing.pool.RemoteTraceback:

Traceback (most recent call last):
  File "/home/fjuhsd/miniconda3/envs/py37/lib/python3.7/multiprocessing/pool.py", line 121, in worker
    result = (True, func(*args, **kwds))
  File "/home/fjuhsd/miniconda3/envs/py37/lib/python3.7/multiprocessing/pool.py", line 44, in mapstar
    return list(map(*args))
  File "sandboxmain.py", line 12, in gaia_to_tic
    tic_ids.append(str(gaiadr2_to_tic(gaia_ids)))
  File "/home/fjuhsd/miniconda3/envs/py37/lib/python3.7/site-packages/astrobase/services/identifiers.py", line 355, in gaiadr2_to_tic
    return str(matched_tic_id.item())
ValueError: can only convert an array of size 1 to a Python scalar


The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "sandboxmain.py", line 34, in <module>
    TIC_index=pool.map(gaia_to_tic, GAIA_index)
  File "/home/fjuhsd/miniconda3/envs/py37/lib/python3.7/multiprocessing/pool.py", line 268, in map
    return self._map_async(func, iterable, mapstar, chunksize).get()
  File "/home/fjuhsd/miniconda3/envs/py37/lib/python3.7/multiprocessing/pool.py", line 657, in get
    raise self._value
ValueError: can only convert an array of size 1 to a Python scalar
'''