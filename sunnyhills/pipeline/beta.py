import numpy as np
import pandas as pd 
import os 
from tqdm import tqdm 


def tls_gamma_routine(tic_ids,cache_dir:str='/ar1/PROJ/fjuhsd/shared/tessodyssey/routines/real/tls_gamma_run/cache_dir/',
                      data_root='/ar1/PROJ/fjuhsd/shared/tessodyssey/data/processed/', 
                      detrend_folders = ['BiWeight-w=0.5', 'LombScargle-n=2'], 
                      piped_log:str='/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/beta_log.txt', 
                      results_file:str='/ar1/PROJ/fjuhsd/shared/tessodyssey/routines/real/tls_gamma_run/results/tls_only_key.txt'): 
    
    r'''
    
    Notes 
    -----
    This is a bare-bones routine maximally optimized for speed. 

    Also uses background process logging, so run it like you run alpha. --> I would like to write it all into one script we call *pure alpha* after Ray Dalio's fund (and because it wouldn't have a secondary beta component)
    '''
    
    import pickle 
    from sunnyhills.pipeline_functions import run_tls

    result_lines = []
    if os.path.exists(piped_log):
        with open(piped_log, 'r') as f:
            for line in piped_log: 
                if '|||START' in line: 
                    result_lines.append(line.split('|||START|||')[-1].split('|||END|||')[0])

        with open(piped_log, 'w') as f1: 
            f1.truncate()

    result_keys = ['TIC_ID', 'data_segment', 'detrend_method', 'SDE', 'period', 'T0', 'duration', 'depth', 'rp_rs', 'snr']
    if not os.path.exists(results_file):
        with open(results_file, 'w') as f: 
            f.write(','.join(result_keys)+'\n')

    if len(result_lines)>0: 
        with open(results_file, 'a') as f: 
            for line in result_lines: 
                f.write(line+'\n')

    extant_tic_ids = []
    with open(results_file, 'r') as f: 
        for line in f: 
            extant_tic_ids.append(line.split(',')[0])

    tic_ids = np.setdiff1d(tic_ids, extant_tic_ids)

    for tic_id in tqdm(tic_ids):   
        try:  
            for detrend_folder in detrend_folders: 

                data_path = data_root+detrend_folder+'/'+tic_id+'.csv'
                if os.path.exists(data_path):
                    data = pd.read_csv(data_path) 
                    data = data.dropna(subset=['clean_time','clean_flux'])
                    clean_time = np.array(data['clean_time'])
                    clean_flux = np.array(data['clean_flux'])
                    
                    # split by time 
                    # then run on Wotan and Lomb-Scargle Detrended

                    diff = np.diff(clean_time)
                    break_indices = np.where(diff>15)[0]+1
                    groups = len(break_indices)+1

                    if groups>1: 
                        time_segments = np.split(clean_time, break_indices)
                        flux_segments = np.split(clean_flux, break_indices)

                    else: 
                        time_segments = [clean_time]
                        flux_segments = [clean_flux]

                    for counter in range(groups):
                        time_segment = time_segments[counter]
                        flux_segment = flux_segments[counter]            

                        model_cache_path = cache_dir+tic_id+'_m:'+detrend_folder+'_s:'+str(counter)+'_tls-model.pickle' 
                        results_cache_path = cache_dir+tic_id+'_m:'+detrend_folder+'_s:'+str(counter)+'_tls-results.pickle'

                        if os.path.exists(results_cache_path) and os.path.exists(model_cache_path): 
                            with open(results_cache_path, 'rb') as file: 
                                tls_results = pickle.load(file)
                                
                            with open(model_cache_path, 'rb') as file: 
                                tls_model = pickle.load(file)
                                
                        else: 
                            tls_results, tls_model = run_tls(tic_id=tic_id, time=time_segment, flux=flux_segment, 
                                                                model_cache_path=model_cache_path, results_cache_path=results_cache_path, catalog_params=False)
                        
                        results_list = [tic_id, counter, detrend_folder, tls_results.SDE, tls_results.period, tls_results.T0,
                                    tls_results.duration, tls_results.depth, tls_results.rp_rs, tls_results.snr]
                    
                        results_list = [str(i) for i in results_list]
            
                        result_line = '|||START|||'+','.join(results_list)+'|||END|||'
                        print('\n'+result_line+'\n')
        
        except Exception as e: 
            print(e)
            continue 

tic_ids = []
with open('/ar1/PROJ/fjuhsd/shared/tessodyssey/data/temp_alpha_output.csv') as f: 
    for line in f: 
        tic_ids.append(line.split(',')[0])
tls_gamma_routine(tic_ids=tic_ids)
# Current PID: [18] 795511