import pandas as pd
import numpy as np 

bad_ids = []
with open('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/bad_ids_log.txt', 'r') as f: 
    for line in f: 
        if 'pdcsap_flux' in line: 
            if 'TIC' and 'SPOC' in line: 
                line = line.split(' SECTOR=')[0]
                tic_id = line.split('LABEL=')[-1].replace('"','')
                if tic_id not in bad_ids: 
                    bad_ids.append(tic_id)

df = pd.DataFrame(bad_ids, columns=['TIC_ID'])
df.to_csv('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/bad_ids.txt', index=False)