def alpha_routine(key:str, data_dir:str, download_log_file:str, output_log_file:str, plots_dir:str=None): 
    r'''
    Parameters

    Returns

    Notes
        - Only runs BLS
    '''
    
    from sunnyhills.pipeline_functions import run_bls
    from sunnyhills.plotting import bls_validation_mosaic
    from misc import even_odd_transit, lombscargle
    import numpy as np
    import pandas as pd 
    from tqdm import tqdm 
    import os

    key_df = pd.read_csv(key)
    tic_ids = np.array(key_df['tic_id'])

    download_log = pd.read_csv(download_log_file)

    if plots_dir!=None: 
        if plots_dir[-1]!='/': 
            plots_dir+='/'

    if data_dir[-1]!='/':
        data_dir+='/'

    if type(tic_ids[0])==int: 
        tic_ids = np.array(['TIC_'+str(i) for i in tic_ids]) 
    elif type(tic_ids[0])==str: 
        if '_' not in tic_ids[0]: 
            tic_ids = np.array([i.replace('TIC', 'TIC_') for i in tic_ids])

    logged_ids = []
    results_dicts = []
    flag_dicts = []

    for tic_id in tqdm(tic_ids): 
        data = pd.read_csv(data_dir+tic_id+'.csv')
        clean_time = np.array(data['clean_time'])
        clean_flux = np.array(data['clean_flux'])

        clean_mask = np.isfinite(clean_time)
        clean_time=clean_time[clean_mask]
        clean_flux = clean_flux[clean_mask]

        # run bls
        bls_tuple = run_bls(clean_time, clean_flux)

        bls_best_params, bls_results = (bls_tuple[0], 
                                        bls_tuple[1])

        # best_params = [index, period, t0, duration, depth]
                                        
        bls_model, bls_in_transit, bls_stats = (bls_tuple[2], 
                                                bls_tuple[3],   
                                                bls_tuple[4])
        index = bls_best_params[0]
        bls_per = bls_best_params[1]
        results_key = ['period', 't0', 'duration', 'depth', 'depth_err', 'snr', 'log_likelihood']
        results_vals = bls_best_params[1:] + [bls_results.depth_err[index], bls_results.power[index], bls_results.log_likelihood[index]]
        results_dict = dict(zip(results_key, results_vals))

        false_alarm_dict = {}

        # Even-Odd Transit Check
        even_odd_check = even_odd_transit(routine='bls', stats=bls_stats)

        false_alarm_dict['even_odd_transit'] = even_odd_check

        # lomb scargle check 
        lombscargle_flag = False 
        ls_index = np.where(download_log['tic_id']==tic_id)[0]
        if len(ls_index)!=0: 
            ls_index = ls_index[0]
            if  0.95 < bls_per/download_log['ls_period'][ls_index] < 1.05:
                if download_log['ls_power'][ls_index] > download_log['fap_95'][ls_index]: 
                    lombscargle_flag = True 

        false_alarm_dict['lomb_scargle_check'] = lombscargle_flag

        logged_ids.append(tic_id)
        results_dicts.append(results_dict)
        flag_dicts.append(false_alarm_dict)

        if plots_dir!=None: 
            plot_path = plots_dir+tic_id+'.png'
            
            raw_data = data[['no_flare_raw_time', 'no_flare_raw_flux']].dropna()
            no_flare_raw_time = np.array(raw_data['no_flare_raw_time'])
            no_flare_raw_flux = np.array(raw_data['no_flare_raw_flux'])

            trend_data = data[['trend_time', 'trend_flux']].dropna()
            trend_time = np.array(trend_data['trend_time'])
            trend_flux = np.array(trend_data['trend_flux'])

            bls_validation_mosaic(tic_id=tic_id, clean_time=clean_time, clean_flux=clean_flux, trend_time=trend_time, trend_flux=trend_flux, raw_time=no_flare_raw_time, raw_flux=no_flare_raw_flux,
                                   best_params=bls_best_params, bls_results=bls_results, bls_model=bls_model, in_transit=bls_in_transit, bls_stats=bls_stats, path=plot_path)

            '''
            (tic_id: str, clean_time: Any, clean_flux: Any, raw_time: Any, raw_flux: Any, best_params: list, 
            bls_results: Any, bls_model: Any, in_transit: Any, bls_stats: Any, path: str = None, dpi: int = 150) -> None
            
            '''


    if not os.path.exists(output_log_file): 
        with open(output_log_file, 'w') as f: 
            header = ['tic_id']+list(results_dicts[0].keys())+list(flag_dicts[0].keys())
            f.write(','.join([str(i) for i in header])+'\n')

    with open(output_log_file, 'a') as f: 
        for index in range(len(logged_ids)): 
            line = [logged_ids[index]] + list(results_dicts[index].values())+list(flag_dicts[index].values())
            f.write(','.join([str(i) for i in line])+'\n')

def beta_routine(key:str, data_dir:str, download_log:str=None, output_log:str=None, plot_dir:str=None, 
                 tic_ids=None, cache_dir:str='/ar1/PROJ/fjuhsd/shared/github/sunnyhills/routines/alpha_tls/cache_dir', 
                 detrend_plot_dir:str='/ar1/PROJ/fjuhsd/shared/github/sunnyhills/routines/alpha_tls/plots/detrend_plots'):
    
    r'''
    Notes
    -----
        - This is the first TLS pipeline executor routine; only runs TLS
        - It starts by downloading light curves and logging data for any objects in the key lacking data in the data directory. 
        - In the process it creates detrend validation plots for downloaded data. 
        - It then runs TLS on them if there is no cached TLS model/results for them, and caches these. It then makes tls validation plot
        - It then saves tls run information with false alarm checks/info to output_log
    '''
    
    from sunnyhills.pipeline_functions import run_tls 
    import numpy as np
    import pandas as pd 
    from tqdm import tqdm 
    import os
    import pickle 
    from sunnyhills.plotting import tls_validation_mosaic, plot_detrend_validation
    from sunnyhills.false_alarm_checks import tls_even_odd, transit_outliers_fap_test, check_lombscargle
    from sunnyhills.pipeline_functions import download_pipeline

    if plot_dir!=None: 
        if plot_dir[-1]!='/': 
            plot_dir+='/'

    if data_dir[-1]!='/':
        data_dir+='/'

    if cache_dir[-1]!='/': 
        cache_dir+='/'

    if detrend_plot_dir!=None: 
        if detrend_plot_dir[-1]!='/': 
            detrend_plot_dir+='/'

    if tic_ids==None: 
        key_df = pd.read_csv(key)
        tic_ids = np.array(key_df['TIC_ID'])

        extant_tic_ids = [i.replace('.csv','') for i in os.listdir(data_dir) if i!='.gitkeep']

        tic_ids = np.setdiff1d(ar1=tic_ids, ar2=extant_tic_ids)

    # DOWNLOAD AND PLOTS FOR MISSING IDS # 
    download_pipeline(tic_ids=tic_ids, download_dir=data_dir, download_log=download_log)

    download_log = pd.read_csv(download_log)

    result_keys = ['SDE', 'period', 'T0', 'duration', 'depth', 'rp_rs', 'snr']
    flags_appended_to_key = False  
    result_lines = []

    for tic_id in tqdm(tic_ids): 
        #try: 
        data_path = data_dir+tic_id+'.csv'
        if os.path.exists(data_path):
            data = pd.read_csv(data_path) 
            clean_time = np.array(data['clean_time'])
            
            clean_flux = np.array(data['clean_flux'])

            if detrend_plot_dir!=None: 
                plot_detrend_validation(tic_id=tic_id, data_dir=data_path, plot_dir=detrend_plot_dir)
        
            if os.path.exists(cache_dir+tic_id+'_tls-model.pickle'): 
                pickle_results = cache_dir+tic_id+'_tls-results.pickle'
                with open(pickle_results, 'rb') as file: 
                    tls_results = pickle.load(file)

                pickle_model = cache_dir+tic_id+'_tls-model.pickle'
                with open(pickle_model, 'rb') as file: 
                    tls_model = pickle.load(file)
            
            else: 

                tls_results, tls_model = run_tls(tic_id=tic_id, time=clean_time, flux=clean_flux, cache_dir=cache_dir) 

            result_list = [tic_id]+[tls_results[key] for key in result_keys]
            
            tls_validation_mosaic(tic_id=tic_id, data=data_path, tls_results=tls_results, tls_model=tls_model, plot_dir=plot_dir)
            
            # FALSE ALARM CHECKS # 

            lombscargle_dict = check_lombscargle(tic_id=tic_id, tls_results=tls_results, download_log=download_log) 
            even_odd_dict = tls_even_odd(tls_results=tls_results)
            transit_outliers_dict = transit_outliers_fap_test(tls_results=tls_results)

            # ADDING RESULTS TO LIST # 
            if not flags_appended_to_key: 
                result_keys += list(lombscargle_dict.keys()) + list(even_odd_dict.keys()) + list(transit_outliers_dict.keys())
                flags_appended_to_key = True 

            result_list += list(lombscargle_dict.values()) + list(even_odd_dict.values()) + list(transit_outliers_dict.values())

            result_line = ','.join([str(i) for i in result_list])

            result_lines.append(result_line)
            
        #except: 
        #    continue 
        
    with open(output_log, 'a') as f: 
        f.write(','.join(result_keys)+'\n')
        
        for line in result_lines: 
            f.write(line+'\n')


key = 'data/current/current_key.csv'
data_dir = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/routines/alpha_tls/data/two_min_lightcurves'
download_log = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/routines/alpha_tls/data/download_log.txt'
output_log = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/routines/alpha_tls/routine-output.txt'
plot_dir = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/routines/alpha_tls/plots/tls_validation'


import os 
full_set = [i.replace('.csv', '') for i in os.listdir('data/current/processed/two_min_lightcurves') if i!='.gitkeep']

single_id = ['TIC_190885165']

beta_routine(tic_ids=single_id, key=key, data_dir=data_dir, download_log=download_log, output_log=output_log, plot_dir=plot_dir) 