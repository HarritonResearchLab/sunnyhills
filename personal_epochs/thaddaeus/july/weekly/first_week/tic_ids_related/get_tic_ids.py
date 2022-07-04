from sunnyhills.misc import gaia_to_tic

import pandas as pd 
import numpy as np

all_gaia_ids = np.array(pd.read_csv('/ar1/PROJ/fjuhsd/shared/tessodyssey/routines/real/tls_beta_run/routine_ids.csv')['GDR2_ID'])

_, tic_ids = gaia_to_tic(gaia_ids=all_gaia_ids)

zipped = list(zip(tic_ids, all_gaia_ids))

df = pd.DataFrame(zipped, columns=['TIC_ID', 'GDR2_ID'])
df.to_csv('/ar1/PROJ/fjuhsd/shared/tessodyssey/routines/real/tls_beta_run/routine_ids.csv', index=False)

