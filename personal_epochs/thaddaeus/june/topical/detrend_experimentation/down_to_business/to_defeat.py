

def the_huns(tic_id:str):  
    import matplotlib.pyplot as plt 
    import numpy as np 
    import pandas as pd
    import wotan  
    from astropy.timeseries import LombScargle
    from sunnyhills.pipeline_functions import remove_flares
    import random 
    from matplotlib.gridspec import GridSpec
    from transitleastsquares import transitleastsquares
    from transitleastsquares import transit_mask, catalog_info
    from sunnyhills.plotting import tls_validation_mosaic
    import batman 

    data = pd.read_csv('./routines/alpha_tls/data/two_min_lightcurves/'+tic_id+'.csv')

    no_flare_mask = np.isfinite(data['no_flare_raw_time'])
    time = np.array(data['no_flare_raw_time'])[no_flare_mask]
    flux = np.array(data['no_flare_raw_flux'])[no_flare_mask]

    per = 2.67
    rp = 0.07

    percentiles = np.percentile(time, [0,10])

    t0 = random.uniform(percentiles[0], percentiles[1])
    params = batman.TransitParams()
    params.t0 = t0                     #time of inferior conjunction
    params.per = per                      #orbital period
    params.rp = rp                   #planet radius (in units of stellar radii)
    params.a = 15.                       #semi-major axis (in units of stellar radii)
    params.inc = 87.                     #orbital inclination (in degrees)
    params.ecc = 0.                      #eccentricity
    params.w = 90.                       #longitude of periastron (in degrees)
    params.u = [0.1, 0.3]                #limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"       #limb darkening model

    m = batman.TransitModel(params, time)    #initializes model
    transit_flux = m.light_curve(params)-1   

    flux = flux+transit_flux

    fap_levels = None
    periodogram = LombScargle(time,flux, nterms=3)
    frequencies, powers = periodogram.autopower(minimum_frequency=1/15, maximum_frequency=1/0.1, method='fastchi2')

    periods = 1/frequencies

    best_frequency = frequencies[np.argmax(powers)]
    y_fit = periodogram.model(time, best_frequency)

    detrend_flux, trend_flux = wotan.flatten(
            time, flux, return_trend=True,
            method='biweight',
            break_tolerance=1.0,
            window_length=0.25,
            cval=5.0)

    (clean_time, detrend_flux, _), (_, _, _) = remove_flares(time, detrend_flux, trend_flux)
    (ls_clean_time, ls_clean_flux, _), (_, _, _) = remove_flares(time, flux/y_fit, y_fit)

    ab, mass, mass_min, mass_max, radius, radius_min, radius_max = catalog_info(TIC_ID=int(tic_id.replace('TIC_','')))
    ls_tls_model = transitleastsquares(ls_clean_time, ls_clean_flux, verbose=False)
    ls_tls_results = ls_tls_model.power(period_min=0.5,period_max=15,
                              verbose=False, show_progress_bar=False, use_threads=25, u=ab)
    wotan_tls_model = transitleastsquares(time, detrend_flux, verbose=False)
    wotan_tls_results = wotan_tls_model.power(period_min=0.5,period_max=15,
                              verbose=False, show_progress_bar=False, use_threads=25, u=ab)

    flat = np.concatenate((ls_tls_results.power, wotan_tls_results.power))
    ylim = (0, 1.15*np.max(flat))

    ls_p = './personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/fast_rotators/plots/ls_'
    ls_path = ls_p+str(per)+tic_id+'_p='+str(per)+'rp='+str(rp)+'_.png'

    tls_validation_mosaic(tic_id, data=None, tls_results=ls_tls_results, tls_model=ls_tls_model, clean_time=ls_clean_time, clean_flux=ls_clean_flux, 
                          trend_time=time, trend_flux=flux/y_fit, raw_time=time, raw_flux=flux, plot_path=ls_path) 

    wotan_p = './personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/fast_rotators/plots/wotan_'
    wotan_path = wotan_p+str(per)+tic_id+'_p='+str(per)+'rp='+str(rp)+'_.png'

    tls_validation_mosaic(tic_id, data=None, tls_model=wotan_tls_model, tls_results=wotan_tls_results, clean_time=clean_time, clean_flux=detrend_flux, 
                          trend_time=time, trend_flux=trend_flux, raw_time=time, raw_flux=flux, plot_path=wotan_path) 

for id in ['TIC_17417151','TIC_49954689','TIC_134631370','TIC_144282456','TIC_167664935','TIC_197005132','TIC_217182385','TIC_370705545','TIC_449941032']:
    the_huns(id)







        