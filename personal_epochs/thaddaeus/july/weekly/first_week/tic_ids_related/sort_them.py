import pandas as pd 
import numpy as np

fp = '/ar1/PROJ/fjuhsd/shared/tessodyssey/routines/real/tls_beta_run/routine_ids.csv'

df = pd.read_csv(fp)
gaia = np.array(df['GDR2_ID'])
tic = np.array(df['TIC_ID']).astype(int)

idx = np.argsort(gaia)

df['GDR2_ID'] = gaia[idx]
df['TIC_ID'] = tic[idx]

df.to_csv(fp, index=False)

print(df)