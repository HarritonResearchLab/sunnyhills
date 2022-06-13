def download_pipeline(ids:str, download_dir:str, log_file:str): 
    '''
    Arguments: 
        ids : list of tic ids 
        download_dir: directory to download files to 
        log_file: file to append results to 
    '''

    import warnings
    from tqdm import tqdm 
    import pandas as pd
    import numpy as np 
    import os 
    import re
    import lightCurve as lk

    from sunnyhills.pipeline_functions import download_and_preprocess 
    from sunnyhills.false_alarm_checks import lombscargle

    warnings.filterwarnings("ignore")


    completed_ids = np.array([i.replace('.csv', '') for i in os.listdir(download_dir) if i.split('.')[-1]=='csv'])
    completed_ids = np.array([int(i.split('_')[-1]) for i in completed_ids])
    #print(completed_ids)
    ids = np.setdiff1d(ids, completed_ids)
 
    np.random.shuffle(ids)

    log_cols = ['tic_id','clean_counts','no_flare_counts','raw_counts','top_ls_period','top_ls_power','fap_95']
    log_cols += ['cdpp', 'with_flare_std', 'no_flare_std', 'no_flare_diff_std']
    log_cols += ['raw_amplitude_no_flare', 'raw_amplitude_with_flare']
    log_cols += ['first_harmonic_period', 'sub_harmonic_period']
    log_cols += [i+'_period' for i in ['second', 'third', 'fourth', 'fifth']]
    log_cols += [i+'_power' for i in ['second', 'third', 'fourth', 'fifth']]

    if not os.path.exists(log_file): 
        with open(log_file, 'w') as f:
            f.write('\n') 
    
    lines = []

    warnings.filterwarnings("ignore")

    for id in tqdm(ids): 
        try: 
            id = str(id)
            id = id.replace('TIC', '')
            tic_id = 'TIC '+id.replace('_', " ")
            tic_id = re.sub(' +', ' ', tic_id)
            lc_df, counts = download_and_preprocess(tic_id, download_dir)

            # counts = [clean_num_obs, no_flare_num_obs, raw_num_obs]

            time, flux = (lc_df['time'], lc_df['flux'])
            no_flare_mask = lc_df['no_flare_mask']

            ls_results = lombscargle(time=time[no_flare_mask], flux=[no_flare_mask])
            
            ls_period = ls_results[2]
            first_harmonic = ls_period / 2 
            sub_harmonic = ls_period*2

            other_four_periods = list(ls_results[1][1:5])
            other_four_powers = list(ls_results[0][1:5])

            ls_power = ls_results[3]
            fap_95 = ls_results[4][1]

            lc = lk.LightCurve(time=time[no_flare_mask], flux=flux[no_flare_mask])

            cdpp = lc.estimate_cdpp()

            with_flare_std = np.std(flux)
            std = np.std(flux[no_flare_mask])
            diff_std = np.std(np.diff(flux[no_flare_mask]))

            raw_flux_amplitude_no_flare = np.diff(np.percentile(flux[no_flare_mask], q=[90,10]))[0]
            raw_flux_flares_amplitude = np.diff(np.percentile(flux, q=[90,10]))[0]

            log_list = [tic_id, counts[0], counts[1], counts[2], round(ls_period, 3), ls_power, fap_95]

            log_list += [cdpp, with_flare_std, std, diff_std]
            log_list += [raw_flux_amplitude_no_flare, raw_flux_flares_amplitude]
            log_list += [first_harmonic, sub_harmonic]
            log_list += other_four_periods
            log_list += other_four_powers

            log_list = [str(i) for i in log_list]

            line = ','.join(log_list)+'\n'
            lines.append(line)
        
        except: 
            continue 
    
    with open(log_file, 'a') as f: 
        for line in lines: 
            f.write(line)
            
import os 

#ids = os.listdir(r'C:\Users\Research\Documents\GitHub\sunnyhills\data\current\processed\two_min_lightcurves')

import numpy as np
import pandas as pd 
ids = np.array(pd.read_csv('./data/current/current_key.csv')['tic_id'])

download_dir = './data/current/processed/two_min_lightcurves'
log_file = r'C:\Users\Research\Documents\GitHub\sunnyhills\data\current\download_log.txt'
download_pipeline(ids, download_dir, log_file)