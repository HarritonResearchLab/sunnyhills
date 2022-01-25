def make_plot(tic_id, data_dir): 
    import matplotlib.pyplot as plt
    from astropy.timeseries import LombScargle
    import numpy as np
    import pandas as pd
    import lightkurve as lk

    lc_df = pd.read_csv(data_dir+'/TIC '+tic_id+'.CSV')

    ls_df = lc_df[['time', 'flux', 'flux_err']]

    ls_df = ls_df.dropna()

    raw_times = np.array(ls_df['time'])
    raw_fluxes = np.array(ls_df['flux'])
    raw_flux_errs = np.array(ls_df['flux_err'])
    
    light_curve = lk.LightCurve(time=raw_times, flux=raw_fluxes, flux_err=raw_flux_errs)
    light_curve = light_curve.remove_outliers(sigma=3).normalize()

    times = np.array(light_curve.time.value)
    fluxes = np.array(light_curve.flux.value)
    flux_errs = np.array(light_curve.flux_err)

    # FIX to use detrended / clipped lc

    ### MAKE PLOT
    plt.style.use('./other/aesthetics/science.mplstyle') # set style...should work on all our computers because relative


    mosaic = """
        AAB
        CCD
        EEF
        """

    fig = plt.figure(constrained_layout=True, figsize=(12,6))
    ax_dict = fig.subplot_mosaic(mosaic)

    # A is raw light curve
    ax_dict['A'].scatter(raw_times, raw_fluxes, s=3)
    ax_dict['A'].set(xlabel='Time (d)', ylabel='Raw Normalized Flux')

    # to axis one of the subplots, do it like this: ax_dict['A'] (e.g. for subplot A)

    # lomb-periodogram 
    # max frequency = observation window / 3
    # min frequency = ???

    # ls won't work when nans in array 
    
    window_length = np.max(times)-np.min(times)
    periodogram = LombScargle(times, fluxes, flux_errs)

    false_alarm_levels = [0.01]
    faps = periodogram.false_alarm_level(false_alarm_levels)

    frequencies, powers = periodogram.autopower(minimum_frequency=3/window_length, maximum_frequency=10)
    # max freq and min freq 
    periods = 1/frequencies 

    # add FAP and harmonic notation

    ax_dict['D'].plot(periods, powers)
    
    ax_dict['D'].set(xlabel='Period (d)', ylabel='Power')#, yscale='log')

    grey_colors = ['lightgray','darkgray','dimgrey']
    for i in range(len(false_alarm_levels)-1,-1,-1): # Plot them in reverse order so the highest confidence label is 
        #confidence_label = str(100*(1-false_alarm_levels[i]))+'% FAP'
        ax_dict['D'].axhline(y=(faps[i]), xmin=1/10, xmax=window_length/3, color = grey_colors[i],lw=1)
    
    # sigma clipped light curve 
    ax_dict['E'].scatter(times, fluxes, s=3)
    ax_dict['E'].set(xlabel='Time (d)', ylabel=r'$\sigma$'+'-clipped Normalized Flux')

    period = periods[np.argmax(powers)] # get period corresponding to highest power
    phased_dates = np.mod(times, period)/period # phase the dates

    ax_dict['F'].scatter(phased_dates, fluxes, s=3)
    ax_dict['F'].set(xlabel='Phase (P='+str(round(period,2))+' d)', ylabel='Normalized Flux (FIX)')


    # metadata to save after each iteration: 
    '''

    Note: end results for the entire iteration should be one .csv with id column and columns for the parameters below, which are populated by row/ID basis  

    LombScargle Related: 
    > top_period
    > top_period_power
    > fap_99 
    > top period is harmonic (True/False) # fix because is vs possesses is different 
    
    Flare Related: 
    > 

    '''


    plt.show()

    ## LC MOSAIC KEY FOR NOW 
    """
    A: raw light curve 
    B: flare frequency plot
    C: raw light curve with annotated flares
    D: periodogram with annotated harmonics # detrending is doing weird things
    E: iterative 3 sigma clipped light curve
    F: Phase folded light curve 
    """

lc_dir = './data/LightCurve_keys'
tic_id = '118807628'
make_plot(tic_id, lc_dir)