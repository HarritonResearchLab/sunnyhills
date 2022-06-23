def show_with_plot(per, rp, tic_id='TIC_49954689'): 
    import matplotlib.pyplot as plt 
    import numpy as np 
    import pandas as pd
    import wotan  
    from astropy.timeseries import LombScargle
    from sunnyhills.pipeline_functions import remove_flares
    import batman 
    import random 
    from matplotlib.gridspec import GridSpec
    from transitleastsquares import transitleastsquares
    from transitleastsquares import transit_mask, catalog_info
    import multiprocessing

    #import batman

    data = pd.read_csv('./routines/alpha_tls/data/two_min_lightcurves/'+tic_id+'.csv')

    no_flare_mask = np.isfinite(data['no_flare_raw_time'])
    time = np.array(data['no_flare_raw_time'])[no_flare_mask]
    flux = np.array(data['no_flare_raw_flux'])[no_flare_mask]

    # add transits 

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

    # lomb scargle

    best_period = 0
    best_period_power = 0

    fap_levels = None
    periodogram = LombScargle(time,flux, nterms=3)
    frequencies, powers = periodogram.autopower(minimum_frequency=1/5, maximum_frequency=1/0.1, method='fastchi2')

    periods = 1/frequencies

    best_frequency = frequencies[np.argmax(powers)]
    y_fit = periodogram.model(time, best_frequency)

    detrend_flux, trend_flux = wotan.flatten(
            time, flux, return_trend=True,
            method='biweight',
            break_tolerance=1.0,
            window_length=0.1,
            cval=5.0)

    (clean_time, detrend_flux, _), (_, _, _) = remove_flares(time, detrend_flux, trend_flux)
    (ls_clean_time, ls_clean_flux, _), (_, _, _) = remove_flares(time, flux/y_fit, y_fit)

    ## PLOTTING ## 

    plt.style.use('seaborn-darkgrid')
    plt.rcParams['font.family']='serif'
    fig = plt.figure(constrained_layout=True, figsize=(22,12))

    gs = GridSpec(4, 4, figure=fig)
    ax1 = fig.add_subplot(gs[0, :-2]) # raw with Lomb-Scargle Trend
    ax2 = fig.add_subplot(gs[0, -2:]) # raw with wotan trend
    ax3 = fig.add_subplot(gs[1, :-2]) # lomb scargle detrended
    ax4 = fig.add_subplot(gs[1, -2:]) # wotan detrended 
    ax5 = fig.add_subplot(gs[2, 0]) # lomb scargle detrended tls periodogram
    ax6 = fig.add_subplot(gs[2,1]) # lomb scargle detrended tls phase folded  
    ax7 = fig.add_subplot(gs[2,2]) # wotan detrended tls phase folded
    ax8 = fig.add_subplot(gs[2, 3]) # wotan detrended tls periodogram 
    ax9 = fig.add_subplot(gs[3,0]) # lomb scargle periodogram 

    ax1.scatter(time, flux, s=1)
    ax1.plot(time, y_fit, lw=1, c='red')

    ax1.set(xlabel='Time (d)', ylabel='No Flare Flux', title='LOMB-SCARGLE (nterms=3)')

    ax3.scatter(ls_clean_time, ls_clean_flux, s=1)
    ax3.set(xlabel='Time (d)', ylabel='Detrend Flux')

    ax9.plot(periods, powers)
    ax9.set(xlabel='Period (d)', ylabel='LombScargle Power')

    ax2.scatter(time, flux, s=1)
    ax2.plot(time, trend_flux, c='red', lw=1)
    ax2.set(xlabel='Time (d)', ylabel='No Flare Flux', title='WOTAN (window-length=0.1d)')

    ax4.scatter(clean_time, detrend_flux, s=1)
    ax4.set(xlabel='Time (d)', ylabel='Detrend Flux')

    ab, mass, mass_min, mass_max, radius, radius_min, radius_max = catalog_info(TIC_ID=int(tic_id.replace('TIC_','')))
    ls_tls_model = transitleastsquares(ls_clean_time, ls_clean_flux, verbose=False)
    ls_tls_results = ls_tls_model.power(period_min=0.5,period_max=15,
                              verbose=False, show_progress_bar=False, use_threads=25, u=ab)
    wotan_tls_model = transitleastsquares(time, detrend_flux, verbose=False)
    wotan_tls_results = wotan_tls_model.power(period_min=0.5,period_max=15,
                              verbose=False, show_progress_bar=False, use_threads=25, u=ab)

    flat = np.concatenate((ls_tls_results.power, wotan_tls_results.power))
    ylim = (0, 1.15*np.max(flat))

    ax5.plot(ls_tls_results.periods, ls_tls_results.power, lw=1)
    ax5.axhline(y=9, ls='--', c='black')
    ax5.set(xlabel='Period (d)', ylabel='SDE', ylim=ylim)

    ax6.scatter(ls_tls_results.folded_phase, ls_tls_results.folded_y, s=3, c='grey')
    ax6.plot(ls_tls_results.model_folded_phase, ls_tls_results.model_folded_model, color='red')
    ax6.set(xlim=(0.4,0.6), title='PERIOD='+str(ls_tls_results.period))

    ax7.plot(wotan_tls_results.periods, wotan_tls_results.power, lw=1)
    ax7.set(xlabel='Period (d)', ylabel='SDE', ylim=ylim)
    ax7.axhline(y=9, ls='--', c='black')

    ax8.scatter(wotan_tls_results.folded_phase, wotan_tls_results.folded_y, s=3, c='grey')
    ax8.plot(wotan_tls_results.model_folded_phase, wotan_tls_results.model_folded_model, color='red')
    ax8.set(xlim=(0.4,0.6), title='PERIOD='+str(wotan_tls_results.period))

    plt.savefig('./personal_epochs/thaddaeus/june/topical/detrend_experimentation/'+tic_id+'_rp='+str(rp)+'_p='+str(per)+'.png', dpi=200)

#show_with_plot(per=2.4, rp=0.12)
#show_with_plot(per=2.4, rp=0.13)
#show_with_plot(per=0.8, rp=0.12)
#show_with_plot(per=0.8, rp=0.13)
show_with_plot(per=5, rp=0.14)
show_with_plot(per=5, rp=0.15)