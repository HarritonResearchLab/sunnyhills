import numpy as np

def gaia_to_tic(gaia_ids):
    from astrobase.services.identifiers import gaiadr2_to_tic
    import numpy as np
    tic_ids = []
    for id in gaia_ids:
        if id !=None:
            tic_ids.append(gaiadr2_to_tic(id))

        else: 
            tic_ids.append(None)

    return dict(zip(gaia_ids, tic_ids))

def skycoord_to_tic(ra,dec):
    from astropy.coordinates import SkyCoord
    from eleanor import mast
    
    tic_id = mast.tic_from_coords((ra,dec))
    return tic_id

def get_best_period(periodogram):
  import numpy as np
  return periodogram.period[np.argmax(periodogram.power)]

def rebin(x, y, num_bins:int=20): 
    '''
    arguments: 
        x: time values 
        y: corresponding flux values
        num_bins: number of output points; default is 20
    returns: 
        x_rebinned: rebinned x 
        y_rebinned: rebinned y 
    '''
    success = True
    if len(x)>0: 
        step = int(len(x)/num_bins)

        ranges = []

        for i in range(0, len(x)-step, step): 
            ranges.append([i,i+step])

        x_rebinned = np.array([np.mean(x[range[0]:range[1]]) for range in ranges])
        y_rebinned = np.array([np.mean(y[range[0]:range[1]]) for range in ranges])

    else: 
        x_rebinned, y_rebinned = (x,y)
        success=False

    return x_rebinned, y_rebinned, success 

def phase(time, flux, period:float, t0:float=None, duration:float=None, bls_model=None, model_name:str=None, tls_results=None, fraction:int=0.2): 
        '''
        Arguments: 
            time: array of time values
            flux: corresponding array of flux values
            period: period to fold to 
            t0: t0 from BLS/TLS 
            bls_model: bls model (if model_name='BLS', otherwise None)
            fraction: how far from median time to go; default is 0.2
        Returns: 
            return_list: [folded_x, folded_flux] by default for model_name in (BLS, LS, TLS)
                         Note: if model_name is BLS, you need to give bls_model the bls_model you used to fit, so it can calculate the model transit shape
                         Note: if model_name is TLS, then you need to give tls_results the results from model.power() from tls bc TLS already calcs fold/not fold
                         Note: in all cases cases, as noted above, they will return folded model time and flux in addition to (after) folded data time and flux
        '''
        
        import numpy as np
        import warnings

        return_list = []

        if model_name!=None: 
            if model_name=='BLS' or model_name=='TLS': 

                x = (time - t0 + 0.5*period) % period - 0.5*period
                m = np.abs(x) < fraction
                
                return_list = [x[m], flux[m]]

                if model_name=='BLS' and bls_model!=None: 
                    f = bls_model.model(x + t0, period, duration, t0) # for BLS
                    return_list.append([x,f])

                elif model_name=='TLS' and tls_results!=None: 
                    return_list.append([tls_results.model_folded_phase,tls_results.model_folded_flux])

            elif model_name=='LS': 
                from astropy.timeseries import LombScargle

                phased_dates = (np.mod(time, period))/period
                return_list.append(phased_dates)
                return_list.append(flux)
                
                t_fit = np.linspace(0, 1)
                ls = LombScargle(time, flux)
                y_fit = ls.model(t_fit, 1/period)

                return_list.append([t_fit, y_fit])

        return np.array(return_list, dtype=object).flatten()


## MISC ##

def even_odd_transit(routine:str, stats):
    '''
    Arguments
        routine: BLS or TLS 
        stats: necessary object for computing left/right stats (for bls it is what is returned by compute_stats)

    '''

    import numpy as np

    even_odd_flag = False

    if routine=='BLS': 

        depth_odd = stats['depth_odd']
        err_odd = depth_odd[1]
        depth_odd = depth_odd[0]
        depth_even = stats['depth_even']
        err_even = depth_even[1]
        depth_even = depth_even[0]

        depths = np.sort([[depth_odd-err_odd, depth_odd+err_odd], 
                   [depth_even-err_even, depth_even+err_even]])

        if min(depths[1]>max(depths[0])): 
            even_odd_flag = True

    return even_odd_flag

def even_odd_phase_folded(time, flux, results):
    """Return even odd phase folded transits"""

    all_even_indices_in_transit = np.array([])
    all_odd_indices_in_transit = np.array([])

    all_indices = np.arange(0,len(time))

    transit_times = results.transit_times
    transit_duration_in_days = results.duration

    period = results.period
    T0 = results.T0

    for i in range(len(transit_times)):
        mid_transit = transit_times[i]
        tmin = mid_transit -  transit_duration_in_days
        tmax = mid_transit + transit_duration_in_days
        if np.isnan(tmin) or np.isnan(tmax):
            idx_intransit = []
            mean_flux = np.nan
        else:
            idx_intransit = np.where(np.logical_and(time > tmin, time < tmax))[0]

        # Check if transit odd/even to collect the flux for the mean calculations
        if i % 2 == 0:  # even
            all_even_indices_in_transit = np.concatenate((all_even_indices_in_transit, 
                                                          idx_intransit))
        else:  # odd
            all_odd_indices_in_transit = np.concatenate((all_odd_indices_in_transit,
                                                         idx_intransit))

    def fold(time, period, T0):
        """Normal phase folding"""
        #return (time - T0) / period - np.floor((time - T0) / period)
        return (time) / period - np.floor((time) / period)

    all_even_indices_in_transit = all_even_indices_in_transit.astype(int)
    all_odd_indices_in_transit = all_odd_indices_in_transit.astype(int)

    even_transit_flux = flux[all_even_indices_in_transit] 
    even_transit_time = time[all_even_indices_in_transit]
    odd_transit_flux = flux[all_odd_indices_in_transit]
    odd_transit_time = time[all_odd_indices_in_transit]

    even_transit_time_folded = fold(even_transit_time, period, T0)
    odd_transit_time_folded = fold(odd_transit_time, period, T0)

    even_transit_time_folded[even_transit_time_folded>0.5] = even_transit_time_folded[even_transit_time_folded>0.5] - 1
    odd_transit_time_folded[odd_transit_time_folded>0.5] = odd_transit_time_folded[odd_transit_time_folded>0.5] - 1

    return np.abs(even_transit_time_folded), even_transit_flux, np.abs(odd_transit_time_folded), odd_transit_flux, all_even_indices_in_transit, all_odd_indices_in_transit


def lombscargle(time,flux,flux_err:np.array=None,min_per:float=.1,max_per:int=15,calc_fap:bool=True,probabilities:list=[.1,.05,.01]):
    import numpy as np
    from astropy.timeseries import LombScargle
    
    best_period = 0
    best_period_power = 0

    fap_levels = None
    periodogram = LombScargle(time,flux,flux_err)
    frequencies,powers = periodogram.autopower(minimum_frequency=1/max_per, maximum_frequency=1/min_per)

    periods = 1/frequencies

    if calc_fap:
        fap_levels = periodogram.false_alarm_probability(probabilities)

    sorted = np.argsort(powers)[::-1] #descending order
    powers = powers[sorted]
    periods = periods[sorted]

    if len(sorted)>0: 
        best_period = periods[0] 
        best_period_power = powers[0]

    return powers, periods, best_period, best_period_power, fap_levels

## BELOW FUNCTIONS ARE VERONICA'S FOR STARS WITH CONFIRMED PLANETS ##

def download_data(path:str,dir:str=''):
    from sunnyhills.pipeline_functions import download_and_preprocess
    import pandas as pd
    import numpy as np
    from lightkurve.utils import LightkurveError
    from tqdm import tqdm

    df = pd.read_csv(path)
    # drop initials nans
    df = df.dropna(subset=['pl_orbper','sy_pnum','pl_controv_flag'])

    df = df.where(np.logical_and(df['pl_orbper']>.5,df['pl_orbper']<15))
    df = df.where(df['pl_controv_flag']==0)
    df = df.where(df['sy_pnum']==1)

    # filter out the data where the given fields did not meet the conditions
    filtered_df = df.dropna(subset=['pl_orbper','pl_controv_flag','sy_pnum','tic_id'])
    
    
    for item in tqdm(filtered_df['tic_id']):
        try:
          lc = download_and_preprocess(item,dir)
        except (KeyError,LightkurveError):
            pass
     
def append_download_status(path_to_csv:str,dir:str,save_dir):
    import os
    import pandas as pd
    import numpy as np
    
    df = pd.read_csv(path_to_csv)

    downloaded_lcs = np.array(os.listdir(dir))
    new_arr = []
    for lc in downloaded_lcs:
      new_arr.append(os.path.splitext(lc)[0].replace('_',' '))
    downloaded_2_min = []
    downloaded_lcs = new_arr
    
    for item in df['tic_id']:
      if item in downloaded_lcs[:]:
        downloaded_2_min.append(True)    
      else:
        downloaded_2_min.append(False)    
    downloaded_2_min = pd.Series(downloaded_2_min,name="has_two_min")

    df = df.merge(downloaded_2_min,left_index=True,right_index=True)
    df = df.drop_duplicates(subset=['tic_id'])
    df.to_csv(save_dir+'/current_key.csv')

    return df

def download_and_append_status(path_to_csv:str,lc_dir,save_dir:str):
    download_data(path_to_csv,lc_dir)
    append_download_status(path_to_csv,lc_dir,save_dir)

    return np.abs((time-transit_time+hp) % period - hp) < .5*duration

def merge_pdf(pdfs_dir:str,pages_per_pdf:int=5):
    import os
    import PyPDF2

    pdfs = os.listdir(pdfs)

    