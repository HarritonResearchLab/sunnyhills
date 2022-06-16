# cSpell:ignore lightkurve, biweight, cval, dtrdict, detrended
# cSpell: disable

from tkinter.ttk import Progressbar
import numpy as np

## DATA PROCESSING RELATED ##

## DATA PROCESSING RELATED ##

def download(
    ticstr: str, 
    outdir: str = 'none', 
    logdir: str = 'none'): 
    
    ''' 
    Args:
        outdir: directory where lightcurves will be saved. If not set, data will not be saved. 
        ticstr: e.g., 'TIC 441420236'
        logdir: directory for log file
    Returns: 
        raw_list: list of light curve ojects that mean criteria but have not been processed (i.e. not detrended, normalized, or sigma-clipped) 
        data_found: boolean variable telling you if data were found or not 
    '''

    import numpy as np 
    import lightkurve as lk 
    import os 
    import pickle 
    import warnings

    # get the light curve
    data_found = False
    
    warnings.warn('WARNING: THIS SHOULD NO LONGER BE USED BY ITSELF! NEEDS TO BE FIXED...USE download_preprocess INSTEAD')

    '''
    lcset = lk.search_lightcurve(ticstr) # otherwise it'll fail for TIC IDs without 120 second cadence data.
    if len(lcset) > 0:
        lcc = lcset[(lcset.author=='SPOC') & (lcset.exptime.value==120)].download_all()
        data_found = True
    '''

    lcc = lk.search_lightcurve(ticstr.replace('_', ' ')).download_all() # FIX THIS! 
    if lcc != None: 
        data_found = True 

    # select only the two-minute cadence SPOC-reduced data; convert to a list.
    # note that this conversion approach works for any LightCurveCollection
    # returned by lightkurve -- no need to hand-pick the right ones.  the exact
    # condition below says "if the interval is between 119 and 121 seconds,
    # take it".

    if data_found: 
        raw_list = [_l for _l in lcc
                if
                _l.meta['ORIGIN']=='NASA/Ames'
                and
                np.isclose(
                    120,
                    np.nanmedian(np.diff(_l.remove_outliers().time.value))*24*60*60,
                    atol=1
                )
        ]

        raw_list = [_l for _l in raw_list if _l.meta['FLUX_ORIGIN']=='pdcsap_flux']
        
        if len(raw_list) == 0: 
            data_found = False 

        if data_found: 
            new_raw_list = []

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

                new_raw_list.append({'time':time, 'flux':flux})

            raw_list = new_raw_list 
    
    if data_found: 
        if outdir != 'none': 
            joined = {'raw_list':raw_list}
            outfile = outdir + '/' + ticstr.replace(' ', '_') + '_raw_lc.pickle'
            
            with open(outfile, 'wb') as handle:
                pickle.dump(joined, handle, protocol=pickle.HIGHEST_PROTOCOL)

            warnings.warn("FIX THIS TO CSV!!!")
            warnings.warn('idea for the future: just only use download and preprocess together')

    if not data_found: 
        raw_list = None

    '''
    if logdir != 'none': 
            logfile = logdir + '/' + ticstr.replace(" ",'_')+'_log.pickle'

            logfile = convert_to_absolute_path(logfile)
            if os.path.exists(logfile): 
                with open(logfile, 'rb') as f: 
                    content = pickle.load(f)
                    content['sectors']=len(raw_list)

            else: 
                content = {'sectors':len(raw_list)}

            with open(logfile, 'wb') as f: 
                    pickle.dump(content, f, protocol=pickle.HIGHEST_PROTOCOL)
    '''
    return raw_list, data_found

def preprocess(
    raw_list: list,
    ticstr: str = '',
    outdir: str = "none", 
    dtrdict: dict = {'method':'biweight',
                     'window_length':0.5,
                     'cval':5.0,
                     "break_tolerance":1.0}, 
    lower_sigma: int = 10
    ):

    """
    Args:
        raw_list: raw list of light curves; see download() 
        outdir: directory where lightcurves will be saved. If not set, data will not be saved.  --> FIX!
        ticstr: e.g., 'TIC 441420236'.  
        dtrdict: dictionary with keys "window_length", "method", "cval",
                 "break_tolerance", or anything else needed by wotan's `flatten`
                 call.  These are documented at
                 https://wotan.readthedocs.io/en/latest/Usage.html 
        lower_sigma: sigma value for "lower" (non-Gunther) level sigma clip. Default is 10.  
    Returns: 
        lc_list: list of light curve ojects that have met all criteria, been removed of outliers, normalized, and flattened. 
        trend_list: list of light curve objects with x = time, y = trend
        raw_list: list of the raw light curve objects 
    """

    import numpy as np 
    import wotan 
    import pandas as pd
    import warnings

    clean_time = np.array([])
    clean_flux = np.array([])
    trend_time = np.array([])
    trend_flux = np.array([])
    no_flare_raw_time = np.array([])
    no_flare_raw_flux = np.array([])
    raw_time = np.array([])
    raw_flux = np.array([])

    raw_num_obs = 0
    no_flare_num_obs = 0 # was not treated after detrending, just before! 
    clean_num_obs = 0

    for lc in raw_list:

        time = lc['time']
        flux = lc['flux']

        raw_num_obs+=len(time)

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
        
        no_flare_num_obs += len(no_flare_raw_time)

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

        clean_num_obs += len(clean_time)

    if outdir != 'none':  
        if outdir[-1]!='/':  
            outdir+='/'  
        
        outfile = outdir+ticstr.replace(' ','_')+'.csv'
        
        cols = [clean_time, clean_flux, trend_time, trend_flux, no_flare_raw_time, no_flare_raw_flux, raw_time, raw_flux]
        cols = [pd.Series(i) for i in cols]

        col_names = ['clean_time', 'clean_flux', 'trend_time', 'trend_flux', 'no_flare_raw_time', 'no_flare_raw_flux', 'raw_time', 'raw_flux']
    
        dictionary = {}
        for i in range(len(cols)):
            dictionary.update({col_names[i]:cols[i]})

        out_df = pd.DataFrame(dictionary)

        out_df.to_csv(outfile, index=False)

    return out_df, [clean_num_obs, no_flare_num_obs, raw_num_obs]

def download_and_preprocess(
    ticstr: str = '',
    outdir: str = 'none', 
    logdir: str = 'none', 
    dtrdict: dict = {'method':'biweight',
                     'window_length':0.5,
                     'cval':5.0,
                     "break_tolerance":1.0}
    ): 
    
    '''
    Args: 
        ticstr: e.g. 'TIC 441420236'
        outdir: dir to save light curve to. default is none
        logdir: dir to save log file to. default is none
        dtrdict: detrending dictionary 
        
    Returns: 
        lc_list: list of light curve ojects that have met all criteria, been removed of outliers, normalized, and flattened. 
        trend_list: list of light curve objects with x = time, y = trend
        raw_list: list of the raw light curve objects 
        data_found: if data was not found during download, returns tuple of None objects
    '''

    import numpy as np
    import warnings 

    raw_list, data_found = download(ticstr=ticstr, logdir=logdir) 

    if data_found: 
        lc_df, counts = preprocess(raw_list=raw_list, ticstr=ticstr, outdir=outdir, dtrdict=dtrdict)

    else: 
        lc_df, counts = (None, None)
    
    return lc_df, counts 

def remove_flares(time, flux, flux_err=np.array([]), sigma:int=3): 
    ''' 
    Args:
        time: array-like (ideally Pandas.series) object of time values 
        flux: array-like (ideally Pandas.series) object of flux values
        flux_err: optional array-like (ideally Pandas.series) object of flux err values
    Returns: 
        two tuples (or three if an array was passed into flux_err); first tuple gives arrays that have been "cleaned" of flairs, and second gives the removed values. Array order in each tuple is time, flux, (flux_err)
    Additional Info: 
        Idea was taken from https://iopscience.iop.org/article/10.3847/1538-3881/ab5d3a
    '''
    import pandas as pd
    import numpy as np


    if isinstance(time, pd.Series) and isinstance(flux, pd.Series): 
        pass
    else: 
        time = pd.Series(time)
        flux = pd.Series(flux)
        if len(flux_err)>0: 
            flux_err = pd.Series(flux_err)
            removed_flux_err = np.array([])

    removed_time, removed_flux = (np.array([]), np.array([]))

    continue_global = True 
    while continue_global:
        length = len(flux)

        # # We use three copies and extract the flares from the middle to prevent edge effects --> BREDALLL! LOL
        three_flux = pd.concat((flux, flux, flux))
        global_rolling = three_flux.rolling(1024, center=False) 
        global_medians = global_rolling.median()[length:2*length]
        global_stds = global_rolling.std()[length:2*length]
        cutoffs = global_medians+sigma*global_stds

        remove_indices = []
        for i in flux.index: 
            if flux[i]>cutoffs[i] and global_medians[i]!=np.nan: 
                remove_indices.append(i)

        if len(remove_indices)==0: 
            continue_global = False 

        else:     
            removed_time = np.concatenate((removed_time, time[remove_indices])) 
            removed_flux = np.concatenate((removed_flux, flux[remove_indices])) 

            time = time.drop(remove_indices)
            flux = flux.drop(remove_indices)

            if len(flux_err)>0: 
                removed_flux_err = np.concatenate((removed_flux_err, flux_err[remove_indices]))
                flux_err = flux_err.drop(remove_indices)
            
    continue_local = True 
    while continue_local: 
        local_rolling = flux.rolling(128, center=True) 
        local_medians = local_rolling.median() 
        local_stds = local_rolling.std() 
        cutoffs = local_medians+sigma*local_stds 

        remove_indices = [] 
        for i in flux.index: 
            if flux[i]>cutoffs[i] and local_medians[i]!=np.nan: 
                remove_indices.append(i) 

        if len(remove_indices)==0: 
            continue_local = False 

        else:     
            removed_time = np.concatenate((removed_time, time[remove_indices]))
            removed_flux = np.concatenate((removed_flux, flux[remove_indices]))

            time = time.drop(remove_indices)
            flux = flux.drop(remove_indices)

            if len(flux_err)>0: 
                removed_flux_err = np.concatenate((removed_flux_err, flux_err[remove_indices]))
                flux_err = flux_err.drop(remove_indices)

    if len(flux_err)>0: 
        return (time.to_numpy(), flux.to_numpy(), flux_err.to_numpy()), (removed_time, removed_flux, removed_flux_err)
    else: 
        return (time.to_numpy(), flux.to_numpy()), (removed_time, removed_flux)

def download_pipeline(tic_ids:str, download_dir:str, download_log:str): 
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
    import lightkurve as lk 
    import re
    from transitleastsquares import catalog_info

    from sunnyhills.pipeline_functions import download_and_preprocess 
    from sunnyhills.misc import lombscargle

    np.random.shuffle(tic_ids)

    cols = ['TIC_ID', 'clean_num_obs', 'no_flare_num_obs', 'raw_num_obs']
    cols += ['no_flare_sigma', 'no_flare_diff_sigma', 'cdpp']
    cols += ['top_ls_period', 'top_ls_power', 'fap_95', 'top_ls_period_sub_harmonic', 'top_ls_period_harmonic']
    period_cols = ['second_ls_period', 'third_ls_period', 'fourth_ls_period', 'fifth_ls_period']
    cols += period_cols 
    power_cols = ['second_ls_power', 'third_ls_power', 'fourth_ls_power', 'fifth_ls_power']
    cols += power_cols 
    tls_catalog_cols = ['ab', 'mass', 'mass_min', 'mass_max', 'radius', 'radius_min', 'radius_max']
    cols += tls_catalog_cols
    cols += ['num_sectors']

    if not os.path.exists(download_log): 
        with open(download_log, 'w') as f:
            f.write(','.join(cols)+'\n') 

    lines = []

    for tic_id in tqdm(tic_ids): 
        #try: 
        lc_df, counts = (None, None)

        data_path= download_dir + tic_id + '.csv'
        if os.path.exists(data_path): 
            lc_df = pd.read_csv(data_path)
            counts = [len(lc_df['clean_time'].dropna()), 
                      len(lc_df['no_flare_raw_time'].dropna()), 
                      len(lc_df['raw_time'])]

        else: 
            lc_df, counts = download_and_preprocess(tic_id, download_dir)

        # counts = [clean_num_obs, no_flare_num_obs, raw_num_obs]


        no_flare_time, no_flare_flux = (np.array(i) for i in [lc_df['no_flare_raw_time'], lc_df['no_flare_raw_flux']])
        mask = np.isfinite(no_flare_time)
        no_flare_time, no_flare_flux = (i[mask] for i in [no_flare_time, no_flare_flux])

        # noise steps 
        no_flare_lc = lk.LightCurve(time=no_flare_time, flux=no_flare_flux)
        cdpp = no_flare_lc.estimate_cdpp()

        no_flare_sigma = np.std(no_flare_flux)
        no_flare_diff_sigma = np.std(np.diff(no_flare_flux))

        ls_results = lombscargle(time=no_flare_time, flux=no_flare_flux)
        
        top_ls_period = ls_results[2]

        top_ls_period_sub_harmonic = top_ls_period/2
        top_ls_period_harmonic = top_ls_period*2

        top_ls_power = ls_results[3]
        fap_95 = ls_results[4][1]

        four_other_powers, four_other_periods = (ls_results[0][1:5], ls_results[1][1:5])

        line_list = [tic_id] + counts + [no_flare_sigma, no_flare_diff_sigma, cdpp]
        line_list += [top_ls_period, top_ls_power, fap_95, top_ls_period_sub_harmonic, top_ls_period_harmonic]
        line_list += list(four_other_periods) + list(four_other_powers)

        # catalog info from TLS #
        tls_catalog_info = catalog_info(TIC_ID=int(tic_id.replace('TIC_', ''))) 

        ab = tls_catalog_info[0]
        ab = '"('+str(ab[0])+','+str(ab[1])+')"'

        line_list += [ab] + list(tls_catalog_info[1:]) 

        diffs = np.diff(lc_df['raw_time'])
        start_indices = np.where(diffs>25)[0]
        num_sectors = len(start_indices)+1
        
        line_list.append(num_sectors)

        line = ','.join([str(i) for i in line_list])+'\n'
        lines.append(line)
        
        #except: 
        #    continue 
    
    with open(download_log, 'a') as f: 
        for line in lines: 
            f.write(line)

## PERIOD SEARCH ROUTINES ##

## BLS ##
def run_bls(time, flux, 
            bls_params: dict = {'min_per':0.5, 'max_per':15, 
                                'minimum_n_transit':2, 
                                'freq_factor':1,
                                'durations':[0.05, 0.0667, 
                                             0.0834, 0.1,
                                             0.1167, 
                                             0.1334,
                                             0.15, 
                                             0.1667, 
                                             0.1834, 0.2], 
                                'objective':'snr'}): 

    '''
    args: 
        time: array of time values
        flux: array of flux values
        bls_params: params for bls execution. see documentation
    returns: 
        best_params: [index, period, t0, duration, sig_diff] for highest power period (sig_diff is sig_diff between left/right depths)
        results: the BLS results array 
        bls_model: the BLS model  
        in_transit_mask: mask for the in_transit points. to get not in transit, do ~in_transit_mask
        stats: the stats on the best period/duration/t0 
    '''

    from astropy.timeseries import BoxLeastSquares
    import numpy as np
    import warnings
    import lightkurve as lk
    from sunnyhills.borrowed import in_transit_mask

    durations = np.array(bls_params['durations'])

    lc = lk.LightCurve(time,flux)
    lc = lc.flatten()
    bls_model = lc.to_periodogram(method='bls',minimum_period=bls_params['min_per'],maximum_period=bls_params['max_per'],duration=bls_params['durations'])


    period = bls_model.period_at_max_power.value
    t0 = bls_model.transit_time_at_max_power.value
    duration = bls_model.duration_at_max_power.value
    in_transit = in_transit_mask(time,period,2*duration,t0)
    stats = bls_model.compute_stats(period, duration, t0)
    depth = stats['depth'][0]

    best_params = [period, t0, duration,depth]

    return best_params, bls_model, in_transit,stats 
        
def iterative_bls_runner(time:np.array, flux:np.array,  
                     iterations: int=1, 
                     bls_params: dict = {'min_per':0.1, 'max_per':15, 
                                'minimum_n_transit':2, 
                                'freq_factor':1,
                                'durations':[0.05, 0.06666666666666668, 
                                             0.08333333333333334, 0.1,
                                             0.11666666666666668, 
                                             0.13333333333333336,
                                             0.15000000000000002, 
                                             0.16666666666666669, 
                                             0.18333333333333335, 0.2], 
                                'objective':'snr'}, 
                    compute_stats: bool = True): 

    '''
    Args:
        stitched_lc: stitched_lc, per usual  
        iterations: number of times to run the search, can be between 1 and 10 
        bls_params: per usual, dictionary of BLS parameters
        compute_stats: will be set to true by default? 

    Returns: 
        results_dict: dictionary of results from each iteration 
        models_dict: dicitonary of bls models from each iteration
        in_transits_dict: dictionary of in_transit arrays from each iteration
        stats_dict: dictionary of stats from each iteration if compute_stats==True
    '''

    from astropy.timeseries import BoxLeastSquares
    import numpy as np
    import matplotlib.pyplot as plt

    durations = np.array(bls_params['durations'])

    if iterations < 1:
        iterations = 1
    elif iterations > 10: 
        iterations = 10 

    results_dict = {} 
    models_dict = {} 
    in_transits_dict = {} 
    stats_dict = {} 

    iteration_names = ['first', 'second', 'third', 'fourth',
                       'fifth', 'sixth', 'seventh', 'eigth',
                       'ninth', 'tenth']

    for index in range(iterations):
        iter_name = iteration_names[index] 
        bls_model = BoxLeastSquares(t=time, y=flux)

        results = bls_model.autopower(durations, frequency_factor=bls_params['freq_factor'], 
                                minimum_period=bls_params['min_per'], 
                                maximum_period=bls_params['max_per'],
                                objective=bls_params['objective'])
        
        results_dict[iter_name] = results

        index = np.argmax(results.power)
        period = results.period[index]
        t0 = results.transit_time[index]
        duration = results.duration[index]
        
        in_transit = bls_model.transit_mask(time, period, 2*duration, t0)

        in_transits_dict[iter_name] = in_transit

        models_dict[iter_name] = bls_model

        if compute_stats: 
            stats = bls_model.compute_stats(period, duration, t0)
            stats_dict[iter_name] = stats

        time = time[~in_transit]
        flux = flux[~in_transit]

    if compute_stats: 
        return results_dict, models_dict, in_transits_dict, stats_dict

    else: 
        return results_dict, models_dict, in_transits_dict

## TLS ##
def run_tls(tic_id:str, time, flux, 
            cache_dir:str=None, 
            tls_params: dict = {'min_per':0.5, 'max_per':15, 
                                'minimum_n_transit':3, 
                                'freq_factor':1,
                                'core_fraction':0.66}): 

    r'''
    args: 
        tic_id: obvious
        time: time array
        flux: flux array
        tls_params: don't worry about it lol 
    returns: 
        
    '''

    import numpy as np
    from transitleastsquares import transitleastsquares
    from transitleastsquares import transit_mask, catalog_info
    from transitleastsquares.stats import intransit_stats
    import multiprocessing 
    import pickle 

    num_cores = int(tls_params['core_fraction']*multiprocessing.cpu_count())
    tls_model = transitleastsquares(time, flux, verbose=False)
    
    ab, mass, mass_min, mass_max, radius, radius_min, radius_max = catalog_info(TIC_ID=int(tic_id.replace('TIC_',''))) 
    
    tls_results = tls_model.power(period_min=tls_params['min_per'],period_max=tls_params['max_per'],
                              verbose=False, show_progress_bar=False, use_threads=num_cores, u=ab)
    
    if cache_dir!=None: 
        if cache_dir[-1]!='/': 
            cache_dir+='/'

        tls_model_cache_file = cache_dir+tic_id+'_tls-model.pickle'
        tls_results_cache_file = cache_dir+tic_id+'_tls-results.pickle'

        with open(tls_model_cache_file, 'wb') as file: 
            pickle.dump(tls_model, file, protocol=pickle.HIGHEST_PROTOCOL)

        with open(tls_results_cache_file, 'wb') as file: 
            pickle.dump(tls_results, file, protocol=pickle.HIGHEST_PROTOCOL)

   

    return tls_results, tls_model 

def iterative_tls_runner(time, flux, 
                     iterations: int=1, 
                     tls_params: dict = {'min_per':0.5, 'max_per':15, 
                                'minimum_n_transit':2, 
                                'freq_factor':1,
                                'durations':[0.05, 0.06666666666666668, 
                                             0.08333333333333334, 0.1,
                                             0.11666666666666668, 
                                             0.13333333333333336,
                                             0.15000000000000002, 
                                             0.16666666666666669, 
                                             0.18333333333333335, 0.2], 
                                'objective':'snr'}): 

    '''
    Args:
        stitched_lc: stitched_lc, per usual  
        iterations: number of times to run the search, can be between 1 and 10 
        tls_params: per usual, dictionary of tls parameters
        compute_stats: will be set to true by default? 
    Returns: 
        results_dict: dictionary of results from each iteration 
        models_dict: dicitonary of tls models from each iteration
        in_transits_dict: dictionary of in_transit arrays from each iteration
    '''

    from transitleastsquares import transitleastsquares
    from transitleastsquares import transit_mask
    import numpy as np
    import matplotlib.pyplot as plt


    durations = np.array(tls_params['durations'])

    if iterations < 1:
        iterations = 1
    elif iterations > 10: 
        iterations = 10 

    results_dict = {} 
    models_dict = {} 
    in_transits_dict = {} 
    stats_dict = {} 

    iteration_names = ['first', 'second', 'third', 'fourth',
                       'fifth', 'sixth', 'seventh', 'eigth',
                       'ninth', 'tenth']

    for index in range(iterations):
        iter_name = iteration_names[index] 
        tls_model = transitleastsquares(t=time, y=flux)
#tls uses power instead of autopower?
#Removed duration because model.power only takes one argument
        results = tls_model.power(frequency_factor=tls_params['freq_factor'], 
                                minimum_period=tls_params['min_per'], 
                                maximum_period=tls_params['max_per'], 
                                objective=tls_params['objective'])
  
        
        results_dict[iter_name] = results

        index = np.argmax(results.power)
        period = results.periods[index]
        t0 = results.T0
        duration = results.duration
        
        in_transit = transit_mask(time, period, 2*duration, t0)

        in_transits_dict[iter_name] = in_transit

        models_dict[iter_name] = tls_model

        return results_dict, models_dict, in_transits_dict
