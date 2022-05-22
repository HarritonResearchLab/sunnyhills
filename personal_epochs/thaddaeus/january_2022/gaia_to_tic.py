import numpy as np
import pandas as pd 
from astrobase.services.identifiers import gaiadr2_to_tic

#
df = pd.read_csv('Dataset_nospaces.csv')
#thirty_gaia = np.array(df['Distance'])

all_distances = np.array(df['Distance'])

sorted_indices = np.argsort(all_distances)

gaia_ids = np.array(df['GDR2_ID'])[sorted_indices][0:50]

distances = all_distances[sorted_indices][0:50]

tic_ids = np.array([])

for gaia_id in gaia_ids: 
    tic_id = gaiadr2_to_tic(str(gaia_id))
    tic_ids = np.append(tic_ids, tic_id)

tic_df = pd.DataFrame(list(zip(gaia_ids, tic_ids, distances)), columns=['GDR2_ID', 'TIC_ID', 'Distance'])

tic_df.to_csv('closest_50.csv', index=False)

df.to_csv('full_key.csv', index=False)
