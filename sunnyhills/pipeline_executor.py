def alpha_routine(key:str, data_dir:str, download_log_file:str, output_log_file:str, plots_dir:str=None): 
    '''
    Parameters

    Returns

    Notes
        - Only runs BLS
    '''
    
    from sunnyhills.pipeline_functions import run_bls
    from sunnyhills.plotting import bls_validation_mosaic
    from false_alarm_checks import even_odd_transit, lombscargle
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
        results_vals = bls_best_params[1:] + [bls_results.depth_err, bls_results.power[index], bls_results.log_likelihood]
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
            no_flare_raw_time = np.array(data['no_flare_raw_time'])
            no_flare_raw_flux = np.array(data['no_flare_raw_flux'])
            print(bls_model)
            bls_validation_mosaic(tic_id=tic_id, clean_time=clean_time, clean_flux=clean_flux, raw_time=no_flare_raw_time, raw_flux=no_flare_raw_flux,
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

key = r'C:\Users\Research\Documents\GitHub\sunnyhills\personal_epochs\thaddaeus\june\topical\misc\fake_key.csv'
download_log_file = r'C:\Users\Research\Documents\GitHub\sunnyhills\data\current\download_log.txt'
output_log_file = './personal_epochs/thaddaeus/june/topical/pipeline_alpha_dev/output.log'
plots_dir = './personal_epochs/thaddaeus/june/topical/pipeline_alpha_dev'
alpha_routine(key=key, data_dir='./data/current/processed/two_min_lightcurves', download_log_file=download_log_file, output_log_file=output_log_file, plots_dir=plots_dir)