import numpy as np

def gaia_to_tic(gaia_ids):
    from astrobase.services.identifiers import gaiadr2_to_tic
    import numpy as np
    tic_ids = []
    for id in np.array(gaia_ids):
        if id !=None:
            tic_ids.append(gaiadr2_to_tic(id))
    return np.array(tic_ids)

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

def rebin_lightcurve(folded_x, folded_flux, factor:int=30): 
        '''
        Arguments: 
            folded_x: folded time values
            folded_flux: corresponding flux values
        Returns: 
            binned_x: re-binned folded time values
            binned_flux: re-binned flux
        '''

        binned_x, binned_flux = ([], [])
        step = int(len(folded_x)/factor)
        indices = list(range(0, len(folded_x), step))
        for index, i in enumerate(indices): 
            if index==len(indices)-1: 
                binned_x.append(np.mean(folded_x[i:]))
                binned_flux.append(np.mean(folded_flux[i:]))
            else: 
                binned_x.append(np.mean(folded_x[i:indices[index+1]]))
                binned_flux.append(np.mean(folded_flux[i:indices[index+1]]))

        return binned_x, binned_flux 

def phase(time, flux, best_params, bls_model, fraction:int=0.2): 
        '''
        Arguments: 
            time: array of time values
            flux: corresponding array of flux values
            period: period to fold to 
            t0: t0 from BLS/TLS 
            bls_model: bls model (FIX FOR TLS!)
            fraction: how far from median time to go; default is 0.2
        Returns: 
            folded time values, folded flux values, x and f (latter two for red bls model)
        '''
        period = best_params[1]
        t0 = best_params[2]
        duration = best_params[3]
        x = (time - t0 + 0.5*period) % period - 0.5*period
        m = np.abs(x) < 0.2
        f = bls_model.model(x + t0, period, duration, t0)

        return x[m], flux[m], x, f