def test(tic_id:str='TIC_14087575', which:str='left'): 
    from sunnyhills.pipeline_functions import run_tls
    from sunnyhills.plotting import tls_validation_mosaic 
    from sunnyhills.misc import inject 
    import numpy as np 
    import pandas as pd
    import matplotlib.pyplot as plt 
    from matplotlib.gridspec import GridSpec
    from scipy.stats import binned_statistic

    data = pd.read_csv('./routines/alpha_tls/data/two_min_lightcurves/'+tic_id+'.csv')
    mask = np.isfinite(data['clean_flux'])
    time = np.array(data['clean_time'])[mask]
    flux = np.array(data['clean_flux'])[mask]

    time, flux, tls_params = inject(time, flux, per=1.3, rp=0.13)

    break_indices = np.where(np.diff(time)>50)[0]
    first_break = break_indices[0]
    last_break = break_indices[-1]+1

    left_time = time[0:first_break]
    left_flux = flux[0:first_break]

    right_time = time[last_break:]
    right_flux = flux[last_break:]

    if which == 'left': 
        time = left_time
        flux = left_flux 

    elif which == 'right': 
        time = right_time 
        flux = right_flux 

    tls_results, tls_model = run_tls(tic_id=tic_id, time=time, flux=flux, show_progress_bar=True)

    fig = plt.figure(constrained_layout=True, figsize=(12,12))

    gs = GridSpec(2, 2, figure=fig) 

    ax1 = fig.add_subplot(gs[0, 0:]) # detrended light curve

    ax1.scatter(time, flux, s=1)
    ax1.plot(tls_results.model_lightcurve_time, 
             tls_results.model_lightcurve_model, 
             alpha=0.5, color='red', zorder=1)

    ax2 = fig.add_subplot(gs[1, 0:1]) # periodogram  

    for n in range(2, 10):
        ax2.axvline(n*tls_results.period, alpha=0.4, lw=1, linestyle="dashed")
        ax2.axvline(tls_results.period / n, alpha=0.4, lw=1, linestyle="dashed")

    ax2.plot(tls_results.periods, tls_results.power, color='black', lw=0.5)

    ax2.set(xlim=(np.min(tls_results.periods), np.max(tls_results.periods)), 
            xlabel='Period (days)', ylabel='SDE')

    ax3 = fig.add_subplot(gs[1, 1:]) # phase folded 

    folded_phase = tls_results.folded_phase 
    folded_y = tls_results.folded_y
    ax3.scatter(folded_phase, folded_y, s=3, c='grey')
        
    ax3.plot(tls_results.model_folded_phase, tls_results.model_folded_model, color='red')

    mask = np.logical_and(folded_phase<0.53, folded_phase>0.47)

    binned_time = binned_statistic(folded_phase[mask], folded_phase[mask], bins=20)[0]
    binned_flux = binned_statistic(folded_phase[mask], folded_y[mask], bins=20)[0]

    ax3.scatter(binned_time, binned_flux, s=35, c='orange', edgecolor='black')

    ax3.set(xlim=(0.47, 0.53))

    plot_path = './personal_epochs/thaddaeus/june/weekly/last_week/long_gap_experiment/'+tic_id+'_'+which+'.png'
    plt.savefig(plot_path, dpi=150)
 
test(which='both')


