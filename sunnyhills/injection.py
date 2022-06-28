import numpy as np
import pandas as pd 

def inject(tic_id:str, time:np.array, flux:np.array, per:float=None, rp:float=None, t0:float=None,core_fraction:float=0.5):
    r'''
    
    Notes
    -----
        t0 : if t0 is not defined/None, then the algo will choose random time less than 0.05*max(time) to be t0
    '''
    
    import batman 
    import random 
    from sunnyhills.pipeline_functions import query_tls_vizier
    import multiprocessing

    if t0 is None:     
        percentiles = np.percentile(time, [0,5])
        t0 = random.uniform(percentiles[0], percentiles[1])

    R_J = 0.10045 # times R_S
    ab, mass, mass_min, mass_max, radius, radius_min, radius_max = query_tls_vizier(tic_id=tic_id)
    ab = list(ab)
    
    if rp is None: 
        rp = 0.5*R_J     
    
    if per is None: 
        per = np.random.uniform(1,15)

    core_count = int(core_fraction*multiprocessing.cpu_count())

    params = batman.TransitParams()
    params.n_threads = core_count
    params.t0 = t0                     #time of inferior conjunction
    params.per = per                      #orbital period
    params.rp = rp                   #planet radius (in units of stellar radii)
    params.a = 15.                       #semi-major axis (in units of stellar radii) --> fix this sometime!
    params.inc = 90.                     #orbital inclination (in degrees)
    params.ecc = 0.                      #eccentricity
    params.w = 90.                       #longitude of periastron (in degrees)
    if ab is not None: 
        params.u = ab  
        
    else: 
        ab = [] #limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"       #limb darkening model

    m = batman.TransitModel(params, time)    #initializes model
    transit_flux = m.light_curve(params)-1   

    flux = flux+transit_flux

    return time, flux, (per, rp, t0)

def inject_pipeline(ids:list, original_data_dir:str, new_data_dir:str, injection_key:str): 
    r'''
    injection_key : path for csv file to save injection info (i.e. log of all injected periods for all TIC_IDs)
    '''
    
    from sunnyhills.injection import inject 
    import os 
    import wotan 
    from tqdm import tqdm 
    from sunnyhills.pipeline_functions import remove_flares

    for dir in [original_data_dir, new_data_dir]: 
        if dir[-1]!='/': 
            dir+'/'

    injected_tic_ids = []
    injected_periods = []
    injected_t0s = []

    lower_sigma = 10
    dtrdict = {'method':'biweight',
                     'window_length':0.25,
                     'cval':5.0,
                     "break_tolerance":1.0}

    for tic_id in tqdm(ids): 
        path = original_data_dir+tic_id+'.csv'
        if os.path.exists(path): 
            data = pd.read_csv(path)

            mask = np.isfinite(data['raw_flux'])
            raw_time, raw_flux = (np.array(data[i])[mask] for i in ['raw_time', 'raw_flux'])

            indices = np.where(np.diff(raw_time)>75)[0]
            if len(indices)>0: 
                break_index = np.where(np.diff(raw_time)>75)[0][0]+1 # diff returns index of last item before break, so add one to this
                raw_time, raw_flux = (i[0:break_index] for i in [raw_time, raw_flux])

            time, flux, (per, rp, t0) = inject(tic_id=tic_id, time=raw_time, flux=raw_flux)
            injected_tic_ids.append(tic_id)
            injected_periods.append(per)
            injected_t0s.append(t0)

            clean_time = np.array([])
            clean_flux = np.array([])
            trend_time = np.array([])
            trend_flux = np.array([])
            no_flare_raw_time = np.array([])
            no_flare_raw_flux = np.array([])
            raw_time = np.array([])
            raw_flux = np.array([])

            raw_time = np.concatenate((raw_time, time))
            raw_flux = np.concatenate((raw_flux, flux))

            continue_lower_cut = True 
            while continue_lower_cut: 
                below_lower_cut = np.where(flux<(np.median(flux)-lower_sigma*np.std(flux)))[0]
                if len(below_lower_cut)>0: 
                    time = np.delete(time, below_lower_cut)
                    flux = np.delete(flux, below_lower_cut)

                else: 
                    continue_lower_cut=False

            (cleaned_time_temp, cleaned_flux_temp), (_, _) = remove_flares(time, flux)

            no_flare_raw_time = np.concatenate((no_flare_raw_time, cleaned_time_temp))
            no_flare_raw_flux = np.concatenate((no_flare_raw_flux, cleaned_flux_temp))
            
            detrended_flux_temp, trend_flux_temp = wotan.flatten(
                cleaned_time_temp, cleaned_flux_temp, return_trend=True,
                method=dtrdict['method'],
                break_tolerance=dtrdict['break_tolerance'],
                window_length=dtrdict['window_length'],
                cval=dtrdict['cval']
            )

            (cleaned_time_temp, detrended_flux_temp, trend_flux_temp), (_, _, _) = remove_flares(cleaned_time_temp, detrended_flux_temp, trend_flux_temp)

            clean_time = np.concatenate((clean_time, cleaned_time_temp))
            clean_flux = np.concatenate((clean_flux, detrended_flux_temp))
            trend_time = np.concatenate((trend_time, cleaned_time_temp))
            trend_flux = np.concatenate((trend_flux, trend_flux_temp))

            cols = [clean_time, clean_flux, trend_time, trend_flux, no_flare_raw_time, no_flare_raw_flux, raw_time, raw_flux]
            cols = [pd.Series(i) for i in cols]

            col_names = ['clean_time', 'clean_flux', 'trend_time', 'trend_flux', 'no_flare_raw_time', 'no_flare_raw_flux', 'raw_time', 'raw_flux']
        
            dictionary = {}
            for i in range(len(cols)):
                dictionary.update({col_names[i]:cols[i]})

            out_df = pd.DataFrame(dictionary)

            outpath = new_data_dir + tic_id + '.csv'

            out_df.to_csv(outpath, index=False)

        injection_df = pd.DataFrame(np.array([injected_tic_ids, injected_periods, injected_t0s]).T, columns=['TIC_ID', 'PER', 'T0'])
        injection_df.to_csv(injection_key)

'''
import os
orig_dir = './routines/alpha_tls/data/two_min_lightcurves/'
injection_key = './routines/first_bulk_injected/injection_key.csv'
ids = [i.replace('.csv','') for i in os.listdir(orig_dir) if i!='.gitkeep']
inject_pipeline(ids=ids, original_data_dir=orig_dir, new_data_dir='./routines/first_bulk_injected/data/', injection_key=injection_key)
'''

def recover_injected_routine(injection_key:str, data_dir:str, plot_dir:str, report_path:str, cache_dir:str=None):
    from sunnyhills.pipeline_functions import run_tls
    from sunnyhills.plotting import tls_validation_mosaic
    from tqdm import tqdm 
    import os 
    import pickle 

    key = pd.read_csv(injection_key)

    tic_ids = np.array(key['TIC_ID'])

    out_ids = []
    out_periods = []
    out_SDEs = []

    for tic_id in tqdm(tic_ids): 
        try: 

            data_path = data_dir+tic_id+'.csv' 
            data = pd.read_csv(data_path)  
            mask = np.isfinite(data['clean_flux'])
            clean_time, clean_flux = (np.array(data[i])[mask] for i in ['clean_time', 'clean_flux'])
            
            tls_model_cache_file = cache_dir+tic_id+'_tls-model.pickle'
            tls_results_cache_file = cache_dir+tic_id+'_tls-results.pickle'

            if os.path.exists(tls_model_cache_file) and os.path.exists(tls_results_cache_file): 
                pickle_results = cache_dir+tic_id+'_tls-results.pickle'
                with open(pickle_results, 'rb') as file: 
                    tls_results = pickle.load(file)

                pickle_model = cache_dir+tic_id+'_tls-model.pickle'
                with open(pickle_model, 'rb') as file: 
                    tls_model = pickle.load(file) 
            else: 
                tls_results, tls_model = run_tls(tic_id=tic_id, time=clean_time, flux=clean_flux, cache_dir=cache_dir, tls_params={'min_per':0.5, 'max_per':15, 'minimum_n_transit':3, 'freq_factor':5,'core_fraction':.95})
            
            tls_validation_mosaic(tic_id=tic_id, data=data_path, tls_model=tls_model, tls_results=tls_results, plot_dir=plot_dir, plot_type='png')        
            out_ids.append(tic_id)
            out_periods.append(tls_results.period)
            out_SDEs.append(np.max(tls_results.power))
        
        except: 
            continue

    report = pd.DataFrame(np.array([out_ids, out_periods, out_SDEs]).T, columns=['TIC_ID', 'PER', 'SDE'])
    report.to_csv(report_path, index=False)

key = './routines/first_bulk_injected/injection_key.csv'
data_dir = './routines/first_bulk_injected/data/'
plots_dir = './routines/first_bulk_injected/plots/'
report_path = './routines/first_bulk_injected/results.csv'
cache_dir = './routines/first_bulk_injected/cache_dir/'
recover_injected_routine(injection_key=key, data_dir=data_dir, plot_dir=plots_dir, report_path=report_path, cache_dir=cache_dir)


