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

    from sunnyhills.pipeline_functions import download_and_preprocess 
    from sunnyhills.false_alarm_checks import lombscargle

    warnings.filterwarnings("ignore")

    '''
    completed_ids = np.array([i.replace('.csv', '') for i in os.listdir(download_dir) if i.split('.')[-1]=='csv'])
    completed_ids = np.array([int(i.split('_')[-1]) for i in completed_ids])
    print(completed_ids)
    ids = np.setdiff1d(ids, completed_ids)

    print(ids)
    '''

    np.random.shuffle(ids)

    if not os.path.exists(log_file): 
        with open(log_file, 'w') as f:
            f.write('tic_id,clean_counts,no_flare_counts,raw_counts,ls_period,ls_power,fap_95'+'\n') 
    lines = []

    for id in tqdm(ids): 
        try: 
            id = str(id)
            id = id.replace('TIC', '')
            tic_id = 'TIC '+id.replace('_', " ")
            tic_id = re.sub(' +', ' ', tic_id)
            lc_df, counts = download_and_preprocess(tic_id, download_dir)

            # counts = [clean_num_obs, no_flare_num_obs, raw_num_obs]

            lc_df = lc_df[['no_flare_raw_time', 'no_flare_raw_flux']].dropna()

            no_flare_time, no_flare_flux = (np.array(i) for i in [lc_df['no_flare_raw_time'], lc_df['no_flare_raw_flux']])

            ls_results = lombscargle(time=no_flare_time, flux=no_flare_flux)
            
            ls_period = ls_results[2]
            ls_power = ls_results[3]
            fap_95 = ls_results[4][1]

            ls_sig = False 
            if ls_power>fap_95: 
                ls_sig = True

            log_list =[tic_id, counts[0], counts[1], counts[2], round(ls_period, 3), ls_power, fap_95]
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

ids = ['TIC_316295827','TIC_399020708','TIC_369743688','TIC_326538701']

download_dir = './data/current/processed/two_min_lightcurves'
log_file = r'C:\Users\Research\Documents\GitHub\sunnyhills\data\current\download_log.txt'
download_pipeline(ids, download_dir, log_file)