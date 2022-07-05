from sunnyhills.misc import gaia_to_tic

import pandas as pd 
import numpy as np

all_gaia_ids = np.array(pd.read_csv('/ar1/PROJ/fjuhsd/shared/tessodyssey/routines/real/tls_beta_run/routine_ids.csv')['GDR2_ID'])

print(all)

with open('/ar1/PROJ/fjuhsd/shared/tessodyssey/routines/real/tls_beta_run/routine_ids_temp.csv', 'w') as f: 
    f.write('GDR2_ID,TIC_ID\n')
    for gaia_id in all_gaia_ids: 
        tic_id = gaia_to_tic(gaia_ids=[gaia_id])[1][0]
        f.write(str(gaia_id)+','+str(tic_id)+'\n')

