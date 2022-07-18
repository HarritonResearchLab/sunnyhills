from sunnyhills.pipeline_functions import download
from sunnyhills.plotting import gen_cutout


def bls_alpha_routine(key:str, data_dir:str, download_log_file:str, output_log_file:str, plots_dir:str=None): 
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


def tls_alpha_routine(key:str, data_dir:str, download_log:str=None, output_log:str=None, plot_dir:str=None, 
                 tic_ids=None, cache_dir:str='/ar1/PROJ/fjuhsd/shared/github/sunnyhills/routines/alpha_tls/cache_dir'):
    
    r'''
    Notes
    -----
        - This is the first TLS pipeline executor routine; only runs TLS
        - It starts by downloading light curves and logging data for any objects in the key lacking data in the data directory. 
        - In the process it creates detrend validation plots for downloaded data. 
        - It then runs TLS on them if there is no cached TLS model/results for them, and caches these. It then makes tls validation plot
        - It then saves tls run information with false alarm checks/info to output_log
    '''

    import numpy as np
    import pandas as pd 
    from tqdm import tqdm 
    import os
    import pickle 
    from sunnyhills.plotting import tls_validation_mosaic, plot_detrend_validation		
    from sunnyhills.plotting import ls_subplots, transit_plots, phased_aliase_plots,gen_cutout	
    from sunnyhills.misc import merge_plots
    from sunnyhills.false_alarm_checks import tls_even_odd, transit_outliers_fap_test, check_lombscargle		
    from sunnyhills.pipeline_functions import download_pipeline, run_tls

    if plot_dir!=None: 
        if plot_dir[-1]!='/': 
            plot_dir+='/'

    if data_dir[-1]!='/':
        data_dir+='/'

    if cache_dir[-1]!='/': 
        cache_dir+='/'

    detrend_plot_dir = plot_dir+'detrend_plots/'
    ls_subplots_dir = plot_dir+'ls_subplots/'
    tls_validation_dir = plot_dir+'tls_validation/'
    transits_dir = plot_dir+'individual_transits/' 
    harmonics_dir = plot_dir + 'harmonics_dir/'
    cutout_dir = plot_dir+'cutouts/'
    combined_dir = plot_dir+'combined_plots/'

    plotting_dirs = [detrend_plot_dir, ls_subplots_dir, tls_validation_dir, transits_dir, harmonics_dir, cutout_dir,combined_dir]
    for dir in plotting_dirs: 
        if not os.path.exists(dir): 
            os.mkdir(dir)

    '''
    if type(tic_ids)==type(None): 
        key_df = pd.read_csv(key)
        tic_ids = np.array(key_df['TIC_ID'])

        extant_tic_ids = [i.replace('.csv','') for i in os.listdir(data_dir) if i!='.gitkeep']

        tic_ids_to_download = np.setdiff1d(ar1=tic_ids, ar2=extant_tic_ids)

    # DOWNLOAD AND PLOTS FOR MISSING IDS # 
    download_pipeline(tic_ids=tic_ids_to_download, download_dir=data_dir, download_log=download_log)
    '''

    tic_ids = [i.split('_tls')[0] for i in os.listdir(cache_dir) if i!='.gitkeep']
    tic_ids = list(set(tic_ids))

    download_log = pd.read_csv(download_log)

    result_keys = ['SDE', 'period', 'T0', 'duration', 'depth', 'rp_rs', 'snr']
    tls_result_keys = ['SDE', 'period', 'T0', 'duration', 'depth', 'rp_rs', 'snr']
    flags_appended_to_key = False  
    result_lines = []


    def check_in_dir(tic_id, dir): 
        import os 

        names = [i.split('.')[0] for i in os.listdir(dir) if i!='.gitkeep']

        return tic_id in names 

    for tic_id in tqdm(tic_ids):    
        try: 
            data_path = data_dir+tic_id+'.csv'
            if os.path.exists(data_path):
                data = pd.read_csv(data_path) 
                clean_time = np.array(data['clean_time'])
                
                clean_flux = np.array(data['clean_flux'])
                
                if not check_in_dir(tic_id, dir=detrend_plot_dir):
                    plot_detrend_validation(tic_id=tic_id, data_dir=data_dir, plot_dir=detrend_plot_dir)
  
                if os.path.exists(cache_dir+tic_id+'_tls-model.pickle'): 
                    pickle_results = cache_dir+tic_id+'_tls-results.pickle'
                    with open(pickle_results, 'rb') as file: 
                        tls_results = pickle.load(file)
                        
                    pickle_model = cache_dir+tic_id+'_tls-model.pickle'
                    with open(pickle_model, 'rb') as file: 
                        tls_model = pickle.load(file)
                        lombscargle_dict = check_lombscargle(tic_id=tic_id, tls_results=tls_results, download_log=download_log)  
                        
                        

                else: 
                    tls_results, tls_model = run_tls(tic_id=tic_id, time=clean_time, flux=clean_flux, cache_dir=cache_dir)

                all_false_alarms_dict = {}
  
                
          
                even_odd_dict = tls_even_odd(tls_results=tls_results)
                transit_outliers_dict = transit_outliers_fap_test(tls_results=tls_results)
                
                for alarm in [lombscargle_dict, even_odd_dict, transit_outliers_dict]: 
                    all_false_alarms_dict.update(alarm)

                if not check_in_dir(tic_id, dir=tls_validation_dir):
                  tls_validation_mosaic(tic_id=tic_id, data=data_path, tls_results=tls_results, tls_model=tls_model, plot_dir=tls_validation_dir, false_alarms_dictionary=all_false_alarms_dict)
               
                if not check_in_dir(tic_id,dir=transits_dir):
                    transit_plots(transits_dir,tic_id,clean_time,clean_flux,tls_results)

                if not check_in_dir(tic_id, dir=ls_subplots_dir):
                    if len(clean_time) > 10000:
                        ls_subplots(tic_id,ls_subplots_dir,clean_time[:10000],clean_flux[:10000])
                    else:
                        ls_subplots(tic_id,ls_subplots_dir,clean_time,clean_flux) 
                
                if not check_in_dir(tic_id,dir=cutout_dir):
                    gen_cutout(tic_id,cutout_dir)

                if not check_in_dir(tic_id, dir=harmonics_dir):
                    phased_aliase_plots(tic_id=tic_id, time=clean_time, flux=clean_flux, plot_path=harmonics_dir+tic_id+'.pdf',tls_results=tls_results)

                result_list = [tic_id]+[tls_results[key] for key in tls_result_keys]
                
                # FALSE ALARM CHECKS # 
                
                # ADDING RESULTS TO LIST # 
                if not flags_appended_to_key: 
                    result_keys += list(lombscargle_dict.keys()) + list(even_odd_dict.keys()) + list(transit_outliers_dict.keys())
                    flags_appended_to_key = True 

                result_list += list(lombscargle_dict.values()) + list(even_odd_dict.values()) + list(transit_outliers_dict.values())

                result_line = ','.join([str(i) for i in result_list])

                result_lines.append(result_line)
                
                if not check_in_dir(tic_id, dir=combined_dir):
                  merge_plots(tic_id=tic_id,plot_dir=plot_dir,export_dir=combined_dir)


            
        except: 
            continue 
    
    with open(output_log, 'a') as f: 
      f.write(','.join(['TIC_ID']+result_keys)+'\n')
        
      for line in result_lines: 
        f.write(line+'\n')

r'''
data_dir = './routines/real/alpha_tls/data/two_min_lightcurves/'
download_log = './routines/real/alpha_tls/data/download_log.csv'
key = './data/current/current_key.csv'
plot_dir = './routines/real/alpha_tls/plots/'
cache_dir = './routines/real/alpha_tls/cache_dir/' 
output_log = './routines/real/alpha_tls/output_log.txt'
# NEED TO IMPLEMENT HARMONICS DIR 

tls_alpha_routine(key=['TIC_316295827'], data_dir=data_dir, plot_dir=plot_dir, cache_dir=cache_dir, download_log=download_log, output_log=output_log)
'''

def tls_beta_routine_thaddaeus(key:str, working_dir:str, processed_to_search:list): 
    r'''
    arguments 
    ---------
    key : str
        Path to "key" file that has column called "TIC_ID" for doing everything 
    
    working_dir : str
        Path to working directory. All sub directories (such as data, plots, etc.) will be handled inside this function internally. 
    
    processed_to_search : list 
        List of the directory names of the different detrended data to use. E.g. processed_to_search=['LombScargle-n=2','BiWeight-d=0.5'] will run routine on those two detrending pipelines. 

    notes
    -----
    user should run /ar1/PROJ/fjuhsd/shared/tessodyssey/routines/initialize_directory.py before running this routine executor! 

    '''

    import os 
    import pandas as pd 
    import numpy as np
    from tqdm import tqdm 

    key = pd.read_csv(key)

    if working_dir[-1]!='/': 
        working_dir+='/'

    tic_ids = np.array(key['TIC_ID'])

    #try: 
    for tic_id in tqdm(tic_ids): 
        pass 

        # 1. RUN TLS # 

        # 2. PLOTS # 

        # 3. LOGGING # 

def tls_beta_routine(key:str, data_dir:str, download_log:str=None, output_log:str=None, plot_dir:str=None, 
                 tic_ids=None, cache_dir:str='/ar1/PROJ/fjuhsd/shared/github/sunnyhills/routines/alpha_tls/cache_dir'):
    
    r'''
    Notes
    -----
        - This is the first TLS pipeline executor routine; only runs TLS
        - It starts by downloading light curves and logging data for any objects in the key lacking data in the data directory. 
        - In the process it creates detrend validation plots for downloaded data. 
        - It then runs TLS on them if there is no cached TLS model/results for them, and caches these. It then makes tls validation plot
        - It then saves tls run information with false alarm checks/info to output_log
    '''

    import numpy as np
    import pandas as pd 
    from tqdm import tqdm 
    import os
    import pickle 
    from sunnyhills.plotting import tls_validation_mosaic, plot_detrend_validation		
    from sunnyhills.plotting import ls_subplots, transit_plots, phased_aliase_plots,gen_cutout	
    from sunnyhills.misc import merge_plots,find_data_gap
    from sunnyhills.false_alarm_checks import tls_even_odd, transit_outliers_fap_test, check_lombscargle		
    from sunnyhills.pipeline_functions import download_pipeline, run_tls


    if data_dir[-1]!='/':
        data_dir+='/'

    if cache_dir[-1]!='/': 
        cache_dir+='/'

    tic_ids = [i.split('_tls')[0] for i in os.listdir(cache_dir) if i!='.gitkeep']
    tic_ids = list(set(tic_ids))

    download_log = pd.read_csv(download_log)

    result_keys = ['SDE', 'period', 'T0', 'duration', 'depth', 'rp_rs', 'snr']
    tls_result_keys = ['SDE', 'period', 'T0', 'duration', 'depth', 'rp_rs', 'snr']
    flags_appended_to_key = False  
    result_lines = []


    def check_in_dir(tic_id, dir): 
        import os 

        names = [i.split('.')[0] for i in os.listdir(dir) if i!='.gitkeep']

        return tic_id in names 

    for tic_id in tqdm(tic_ids):    
        try: 
            data_path = data_dir+tic_id+'.csv'
            if os.path.exists(data_path):
                data = pd.read_csv(data_path) 
                data = data.dropna(subset=['clean_time','clean_flux'])
                clean_time = np.array(data['clean_time'])
                
                clean_flux = np.array(data['clean_flux'])

                gap_index = find_data_gap(np.array(clean_time))
                
                #if not check_in_dir(tic_id, dir=plotting_dirs[0]):
                    #plot_detrend_validation(tic_id=tic_id, data_dir=data_dir, plot_dir=plotting_dirs[0])
  
                if os.path.exists(cache_dir+tic_id+'_tls-model.pickle'): 
                    pickle_results = cache_dir+tic_id+'_tls-results.pickle'
                    with open(pickle_results, 'rb') as file: 
                        tls_results = pickle.load(file)
                        
                    pickle_model = cache_dir+tic_id+'_tls-model.pickle'
                    with open(pickle_model, 'rb') as file: 
                        tls_model = pickle.load(file)
                        lombscargle_dict = check_lombscargle(tic_id=tic_id, tls_results=tls_results, download_log=download_log)  
                        
                  

                else: 
                    tls_results, tls_model = run_tls(tic_id=tic_id, time=clean_time, flux=clean_flux, cache_dir=cache_dir)

                print(gap_index)
                if gap_index!= None:
                  run_routine(clean_time[:gap_index],clean_flux[:gap_index],tic_id,plot_dir,tls_results,tls_model,data_path)
                  run_routine(clean_time[gap_index:-1],clean_flux[:gap_index:-1],tic_id,plot_dir+'/second_plot/',tls_results,tls_model,data_path)

                else:
                  run_routine(clean_time,clean_flux,tic_id,plot_dir,tls_results,tls_model,data_path)
       
        except: 
            continue 
    
    #with open(output_log, 'a') as f: 
      #f.write(','.join(['TIC_ID']+result_keys)+'\n')
        
      #for line in result_lines: 
        #f.write(line+'\n')
        
def run_routine(clean_time,clean_flux,tic_id,plot_dir,tls_results,tls_model,data_path):
    import numpy as np
    import pandas as pd 
    from tqdm import tqdm 
    import os
    import pickle 
    from sunnyhills.plotting import tls_validation_mosaic, plot_detrend_validation		
    from sunnyhills.plotting import ls_subplots, transit_plots, phased_aliase_plots,gen_cutout	
    from sunnyhills.misc import merge_plots
    from sunnyhills.false_alarm_checks import tls_even_odd, transit_outliers_fap_test, check_lombscargle		
    from sunnyhills.pipeline_functions import download_pipeline, run_tls


    if plot_dir!=None: 
        if plot_dir[-1]!='/': 
            plot_dir+='/'

    for dir in plotting_dirs: 
      if not os.path.exists(dir): 
        os.mkdir(dir)

    detrend_plot_dir = plot_dir+'detrend_plots/'
    ls_subplots_dir = plot_dir+'ls_subplots/'
    tls_validation_dir = plot_dir+'tls_validation/'
    transits_dir = plot_dir+'individual_transits/' 
    harmonics_dir = plot_dir + 'harmonics_dir/'
    cutout_dir = plot_dir+'cutouts/'
    
    combined_dir = plot_dir+'combined_plots/'

    plotting_dirs = [detrend_plot_dir, ls_subplots_dir, tls_validation_dir, transits_dir, harmonics_dir, cutout_dir,combined_dir]
    for dir in plotting_dirs: 
        if not os.path.exists(dir): 
            os.mkdir(dir)
    def check_in_dir(tic_id, dir): 
        import os 

        names = [i.split('.')[0] for i in os.listdir(dir) if i!='.gitkeep']

        return tic_id in names 
    all_false_alarms_dict = {}
    
    lombscargle_dict = check_lombscargle(tic_id=tic_id, tls_results=tls_results, download_log=download_log)      
    even_odd_dict = tls_even_odd(tls_results=tls_results)
    transit_outliers_dict = transit_outliers_fap_test(tls_results=tls_results)
                
    for alarm in [lombscargle_dict, even_odd_dict, transit_outliers_dict]: 
      all_false_alarms_dict.update(alarm)

    if not check_in_dir(tic_id, dir=plotting_dirs[2]):
      tls_validation_mosaic(tic_id=tic_id, data=data_path, tls_results=tls_results, tls_model=tls_model, plot_dir=plotting_dirs[2], false_alarms_dictionary=all_false_alarms_dict)
                
      #if not check_in_dir(tic_id,dir=plotting_dirs[3]):
        #transit_plots(plotting_dirs[3],tic_id,clean_time,clean_flux,tls_results)
        
    if not check_in_dir(tic_id, dir=plotting_dirs[0]):
          plot_detrend_validation(tic_id=tic_id, data_dir=data_path, plot_dir=plotting_dirs[0])  

    if not check_in_dir(tic_id, dir=plotting_dirs[1]):
      ls_subplots(tic_id,plotting_dirs[1],clean_time,clean_flux)

    
    if not check_in_dir(tic_id,dir=plotting_dirs[5]):
      gen_cutout(tic_id,plotting_dirs[5])

    if not check_in_dir(tic_id, dir=plotting_dirs[4]):
          phased_aliase_plots(tic_id=tic_id, time=clean_time, flux=clean_flux, plot_path=plotting_dirs[4]+tic_id+'.pdf',tls_results=tls_results)

    #result_list = [tic_id]+[tls_results[key] for key in tls_result_keys]
    
      # FALSE ALARM CHECKS # 
    
      # ADDING RESULTS TO LIST # 
    '''
    if not flags_appended_to_key: 
        result_keys += list(lombscargle_dict.keys()) + list(even_odd_dict.keys()) + list(transit_outliers_dict.keys())
          flags_appended_to_key = True 

      result_list += list(lombscargle_dict.values()) + list(even_odd_dict.values()) + list(transit_outliers_dict.values())

      result_line = ','.join([str(i) for i in result_list])

      result_lines.append(result_line)
                  
                  if not check_in_dir(tic_id, dir=plotting_dirs[6]):
                    merge_plots(tic_id=tic_id,plot_dir=plot_dir,export_dir=plotting_dirs[6])

                  plotting_dirs = [dir.replace(plot_dir,'') for dir in plotting_dirs]
                  plotting_dirs = [plot_dir+'second_plot/'+dir for dir in plotting_dirs]
                  print(plotting_dirs[0])
                  for dir in plotting_dirs: 
                    if not os.path.exists(dir): 
                      os.mkdir(dir)
      '''

def run_tls(ls_path:str='/ar1/PROJ/fjuhsd/shared/data/processed/BiWeight-w=0.5/',biweight_path:str='/ar1/PROJ/fjuhsd/shared/data/processed/LombScargle-n=2'):
    from sunnyhills.pipeline_functions import modified_run_tls
    import os
    #TODO: Implement Ryan's multiprocessing so this process can be sped up
    
    # just using a for loop for now but will switch to multiprocessing later
    #NOTE TO SELF DONT RUN THE WHOLE LOOP AT ONCE! TEST ON ONE INDIVIDUAL CSV BEFORE RUNNING THIS PROCESS!
    
    for path in os.listdir(ls_path):
        modified_run_tls(ls_path+path)

    for path in os.listdir(biweight_path):
        modified_run_tls(path)