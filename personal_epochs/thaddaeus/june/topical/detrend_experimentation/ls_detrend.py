def show_with_plot(tic_id='TIC_49954689'): 
    import matplotlib.pyplot as plt 
    import numpy as np 
    import pandas as pd
    import wotan  
    from astropy.timeseries import LombScargle
    from sunnyhills.pipeline_functions import remove_flares

    #import batman

    data = pd.read_csv('./routines/alpha_tls/data/two_min_lightcurves/'+tic_id+'.csv')

    no_flare_mask = np.isfinite(data['no_flare_raw_time'])
    time = np.array(data['no_flare_raw_time'])[no_flare_mask]
    flux = np.array(data['no_flare_raw_flux'])[no_flare_mask]

    '''
    # add transits 
    synthetic_transit = batman.TransitParams()
    # make the fake transits
    synthetic_transit.per = np.random.rand()*18+1
    synthetic_transit.rp = (np.random.rand()*18+2)
    synthetic_transit.t0 = time_start
    synthetic_transit.a = 19
    synthetic_transit.inc = 90
    synthetic_transit.ecc = 0
    synthetic_transit.w = 90
    synthetic_transit.u = [.4,.4]
    synthetic_transit.limb_dark = "quadratic"

    transit = batman.TransitModel(synthetic_transit,time)
    signal = transit.light_curve(synthetic_transit)
    noise = np.random.normal(0, 10**-6 * noise_level, int(samples))
    '''


    
    
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

    fig, axs = plt.subplots(5,1, figsize=(16,8))

    ax = axs[0]
    ax.scatter(time, flux, s=1)
    ax.plot(time, y_fit, lw=1, c='red')

    ax = axs[1]

    ax.scatter(time, flux-y_fit, s=1)

    ax = axs[2]

    ax.plot(periods, powers)

    ax = axs[3]
    ax.scatter(time, flux, s=1)
    ax.plot(time, trend_flux, c='red', lw=1)

    ax = axs[4]
    ax.scatter(clean_time, detrend_flux, s=1)

    plt.show()

show_with_plot()