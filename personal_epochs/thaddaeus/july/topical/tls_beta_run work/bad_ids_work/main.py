def download_a_couple(bad_ids):

    import lightkurve as lk 
    import pandas as pd
    from sunnyhills.pipeline_functions import remove_flares, remove_extreme_dips
    import matplotlib.pyplot as plt 
    import numpy as np 

    # get the light curve
    for tic_id in bad_ids: 
        try:
            lcc = lk.search_lightcurve(tic_id.replace('_', ' ')).download_all() 
            raw_list = [_l for _l in lcc if _l.meta['ORIGIN']=='NASA/Ames']
            raw_time = np.array([])
            raw_flux = np.array([])

            for lc in raw_list: 
                time = lc.time.value
                flux = lc.flux.value

                nan_mask = np.isnan(flux)

                time = time[~nan_mask]
                flux = np.array(flux[~nan_mask], dtype=float) 

                qual = lc.quality.value

                # remove non-zero quality flags
                sel = (qual[~nan_mask] == 0)

                time = time[sel]
                flux = flux[sel]

                # normalize around 1
                flux /= np.nanmedian(flux)

                raw_time = np.concatenate((raw_time, time))
                raw_flux = np.concatenate((raw_flux, flux)) 


            (no_flare_raw_time, no_flare_raw_flux), (_, _) = remove_flares(raw_time, raw_flux)
            no_flare_raw_time, no_flare_raw_flux = remove_extreme_dips(no_flare_raw_time, no_flare_raw_flux)


            fig, ax = plt.subplots(figsize=(20,2))
            ax.scatter(no_flare_raw_time, no_flare_raw_flux)
            ax.set(title='np.nanmedian(np.diff(time))='+str(np.nanmedian(np.diff(no_flare_raw_time))*24*60*60))
            plt.savefig('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/personal_epochs/thaddaeus/july/topical/tls_beta_run work/bad_ids_work/plots/'+tic_id+'.png', dpi=150)
            plt.clf()
            plt.close()

            print(np.median(np.diff(no_flare_raw_time))*24*60*60)
        
        except: 
            continue 

import pandas as pd 
download_a_couple(list(pd.read_csv('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/bad_ids.txt')['TIC_ID']))

    
    