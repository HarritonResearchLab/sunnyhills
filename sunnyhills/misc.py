import numpy as np

def gaia_to_tic(gaia_ids):
    from astrobase.services.identifiers import gaiadr2_to_tic
    import numpy as np
    tic_ids = []
    for id in np.array(gaia_ids):
        if id !=None:
            tic_ids.append(gaiadr2_to_tic(id))
    return np.array(tic_ids)

def skycoord_to_tic(ra,dec):
    from astropy.coordinates import SkyCoord
    from eleanor import mast
    
    tic_id = mast.tic_from_coords((ra,dec))
    return tic_id

def lomb_scargle(time,flux,flux_err:np.array=None,min_per:float=.1,max_per:int=15,calc_fap:bool=True,probilities:list=[.1,.05,.01]):
    import numpy as np
    from astropy.timeseries import LombScargle
    
    best_period = 0
    best_period_power = 0

    fap_levels = None
    periodogram = LombScargle(time,flux,flux_err)
    frequencies,powers = periodogram.autopower(min_per,max_per)

    periods = 1/frequencies

    if calc_fap:
        fap_levels = periodogram.false_alarm_probability(probilities)

    sorted = np.argsort(powers)[::-1] #descending order
    powers = powers[sorted]
    periods = periods[sorted]

    if len(sorted)>0: 
        best_period = periods[0] 
        best_period_power = powers[0]

    return powers,periods,best_period,best_period_power,fap_levels

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
    
    step = int(len(x)/num_bins)

    ranges = []

    for i in range(0, len(x)-step, step): 
        ranges.append([i,i+step])

    x_rebinned = np.array([np.mean(x[range[0]:range[1]]) for range in ranges])
    y_rebinned = np.array([np.mean(y[range[0]:range[1]]) for range in ranges])

    return x_rebinned, y_rebinned 

def phase(time, flux, best_params, model:None, model_name:None, fraction:int=0.2): 
        '''
        Arguments: 
            time: array of time values
            flux: corresponding array of flux values
            period: period to fold to 
            t0: t0 from BLS/TLS 
            bls_model: bls model (FIX FOR TLS!)
            fraction: how far from median time to go; default is 0.2
        Returns: 
            return_list: [folded x, folded flux, mask] additional things depending on model
        '''
        
        import numpy as np
        import warnings

        period = best_params[1]
        t0 = best_params[2]
        duration = best_params[3]
        x = (time - t0 + 0.5*period) % period - 0.5*period
        m = np.abs(x) < 0.2
        
        return_list = [x[m], flux[m], m]

        if model_name=='BLS' and model!=None: 
            f = model.model(x + t0, period, duration, t0) # for BLS
            return_list.append(x)
            return_list.append(f)

        return return_list
