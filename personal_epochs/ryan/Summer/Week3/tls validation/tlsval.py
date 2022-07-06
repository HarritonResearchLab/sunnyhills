
def tls_validation_mosaic(tic_id:str, data, tls_model, tls_results, false_alarms_dictionary:dict,
                          plot_dir:str=None, plot_path:str=None, dpi:int=150, clean_time=None, clean_flux=None, 
                          trend_time=None, trend_flux=None, raw_time=None, raw_flux=None, plot_type:str='pdf', 
                          even_odd_option:str='separate', true_transit_times:np.array=None): 

    '''
    arguments: 
        tic_id: tic id 
        data: path to data csv 
        clean_time: detrended and flare removed time array
        clean_flux: flux values corresponding to clean_time arg
        trend_time: time values of the trend
        trend_flux: flux values of the trend 
        raw_time: array of "raw" time values, i.e. not detrended and potentially with flares
        raw_flux: flux array corresponding to raw_time arg
        best_params, bls_results, bls_model, in_transit, bls_stats: the items returned by the run_bls function in pipeline_functions.py (NOTE: in run_bls, stats must be set to be calculated!)
        path: if defined, plot will be saved as the provided path. Otherwise, it will be displayed
        dpi: dpi of saved plot
 

    returns: 
    
    notes 
    -----

        if data is defined, you don't need to define all the individual arrays! I just added the option to do the individual arrays for some testing I was doing at one point 

        false_alarms_dictionary: dictionary of false_alarm_name : False alarm value items from different false alarm tests. 

        Also, if plot_type is pdf, it will save the file as pdf, otherwise if it's png the plot will get saved as png

        Finally, even_odd_option can be 'separate' which plots odd transit with time [0,1] and even transit with time [1,2], or 'together' which plots them over each other. 
    '''

    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from scipy.stats import binned_statistic
    from sunnyhills.misc import even_odd_phase_folded, return_kerr_cluster
    from transitleastsquares import (
        transitleastsquares,
        cleaned_array,
        catalog_info,
        transit_mask
    )

    from sunnyhills.pipeline_functions import query_simbad
    from sunnyhills.misc import normalize 
    
    if data is not None: 
        df = pd.read_csv(data)
        clean_time, clean_flux = (np.array(df[i]) for i in ['clean_time', 'clean_flux'])
        clean_mask = np.isfinite(clean_time)
        clean_time, clean_flux = (i[clean_mask] for i in [clean_time, clean_flux])
        raw_time, raw_flux = (np.array(df[i]) for i in ['no_flare_raw_time','no_flare_raw_flux'])
        trend_time, trend_flux = (np.array(df[i]) for i in ['trend_time','trend_flux'])
        trend_mask = np.isfinite(trend_time)
        trend_time, trend_flux = (i[trend_mask] for i in [trend_time, trend_flux])

    in_transit = transit_mask(clean_time, tls_results.period, tls_results.duration, tls_results.T0)

    transit_times = np.array(tls_results.transit_times) 
    transit_depths = tls_results.transit_depths
    transits_not_nan_mask = np.isfinite(transit_depths)

    transit_times = transit_times[transits_not_nan_mask]
    transit_depths = transit_depths[transits_not_nan_mask]

    split_axes = False
    break_index = None 
    last_break = None

    diff = np.diff(clean_time)
    trend_diff = np.diff(trend_time)
    raw_diff = np.diff(raw_time)
    if np.max(diff>50): 
        split_axes = True 
        break_indices = np.where(diff>50)[0]
        first_clean_break = break_indices[0]
        last_clean_break = break_indices[-1]+1

        raw_idx = np.where(raw_diff>50)[0]
        first_raw_break = raw_idx[0]
        last_raw_break = raw_idx[-1]+1

        trend_idx = np.where(trend_diff>50)[0]
        first_trend_break = trend_idx[0]
        last_trend_break = trend_idx[-1]+1

    ## PLOTTING ### 

    #plt.style.use('seaborn-darkgrid')
    plt.rcParams['font.family']='serif'
    fig = plt.figure(constrained_layout=True, figsize=(21,12))

    gs = GridSpec(3, 5, figure=fig)

    def orient_split_axes(ax_1, ax_2, flux): 
        ax_1.spines['right'].set_visible(False)
        ax_2.spines['left'].set_visible(False)
        ax_1.yaxis.tick_left()
        ax_2.yaxis.tick_right()

        temp = flux[np.isfinite(flux)]

        ylim = [0.95*np.min(temp), 1.05*np.max(temp)]

        ax_1.set(ylim=ylim)

        ax_2.set(ylim=ylim)

        d = .015 # how big to make the diagonal lines in axes coordinates
        kwargs = dict(transform=ax_1.transAxes, color='k', clip_on=False)
        ax_1.plot((1-d,1+d),(-d,+d), **kwargs) # top-left diagonal
        ax_1.plot((1-d,1+d),(1-d,1+d), **kwargs) # bottom-left diagonal

        kwargs.update(transform=ax_2.transAxes) 
        ax_2.plot((-d,d),(-d,+d), **kwargs) 
        ax_2.plot((-d,d),(1-d,1+d), **kwargs)

    if split_axes: 
        ax1a = fig.add_subplot(gs[0, 0:2])
        ax1b = fig.add_subplot(gs[0, 2:-1])
        ax2a = fig.add_subplot(gs[1, 0:2])
        ax2b = fig.add_subplot(gs[1, 2:-1])
        
        ax1a.scatter(clean_time[0:first_clean_break], clean_flux[0:first_clean_break], s=1)

        #ax1a.plot(tls_results.model_lightcurve_time[0:break_index], 
        #          tls_results.model_lightcurve_model[0:break_index], 
        #          alpha=0.5, color='red', zorder=1)

        max_left = np.max(clean_time[0:first_clean_break])
        min_right = np.min(clean_time[last_clean_break:])
        left_transits = transit_times[transit_times<=max_left]
        right_transits = transit_times[transit_times>=min_right]

        ax1a.set(ylabel='Detrended Flux')
        ax1a.set_title('TIC: '+str(tic_id).replace('_','')+' PERIOD: '+str(round(tls_results.period, 5)), size='xx-large')
        
        ax1b.scatter(clean_time[last_clean_break:], clean_flux[last_clean_break:], s=1)
        
        for left in left_transits: 
            ax1a.axvline(x=left, color='red', alpha=0.4, lw=2)

        for right in right_transits: 
            ax1b.axvline(x=right, color='red', alpha=0.4, lw=2)

        if true_transit_times is not None: 
            left_true_transits = true_transit_times[true_transit_times<=max_left]
            right_true_transits = true_transit_times[true_transit_times>=min_right]

            for left in left_true_transits: 
                ax1a.axvline(x=left, color='orange', alpha=0.4, lw=2)

            for right in right_true_transits: 
                ax1b.axvline(x=right, color='orange', alpha=0.4, lw=2)



        '''
        ax1b.plot(tls_results.model_lightcurve_time[last_break:], 
                  tls_results.model_lightcurve_model[last_break:], 
                  alpha=0.5, color='red', zorder=1)
        '''

        # FIX THIS! IDK WHY I NEED TO ADD 10 TO BREAK INDEX!
#######################################################S ALREADY EQUALS 1!!!!#######################################################
        orient_split_axes(ax1a, ax1b, clean_flux)

        ax2a.scatter(raw_time[0:first_raw_break], raw_flux[0:first_raw_break], s=1) 
        ax2a.plot(trend_time[0:first_trend_break], trend_flux[0:first_trend_break], lw=2, c='r') 
        ax2a.set(ylabel='Flux')
        ax2b.scatter(raw_time[last_raw_break:], raw_flux[last_raw_break:], s=1) 
        ax2b.plot(trend_time[last_trend_break:], trend_flux[last_trend_break:], lw=2, c='r') 
        
        orient_split_axes(ax2a, ax2b, raw_flux) 

        for ax in [ax1a, ax1b, ax2a, ax2b]: 
            ax.set(xlabel='Time (days)') 

    else: 
        ax1 = fig.add_subplot(gs[0, 0:-1]) # detrended light curve
        ax2 = fig.add_subplot(gs[1, 0:-1]) # no flare with trend light curve
####################################################################################################################################
    ######################################################TOP PLOT S ALREADY = 1!!!!!###############################################
        # detrend light curve 
        ax1.scatter(clean_time, clean_flux, s=1)
        ax1.plot(tls_results.model_lightcurve_time, 
                  tls_results.model_lightcurve_model, 
                  alpha=0.5, color='red', zorder=1)

        for i in transit_times: 
            ax1.axvline(x=i, color='red', alpha=0.5, lw=2)
        
        ax1.set(ylabel='Detrended Flux')
        ax1.set_title('TIC: '+str(tic_id).replace('_','')+' PERIOD: '+str(round(tls_results.period, 5)), size='xx-large')
        
        if true_transit_times is not None: 
            for true_transit_time in true_transit_times: 
                ax1.axvline(x=true_transit_time, color='orange', alpha=0.4, lw=2)

        # raw and trend light curve # 
        ax2.scatter(raw_time, raw_flux, s=1)
        ax2.plot(trend_time, trend_flux, lw=0.5, c='r')
        ax2.set(ylabel='Flux')

        for ax in [ax1, ax2]: 
                ax.set(xlabel='Time (days)')
#####################################################################################################################################
    ax3 = fig.add_subplot(gs[2, 0]) # phase folded transit 
    ax4 = fig.add_subplot(gs[2, 1]) # left right transits
    ax5 = fig.add_subplot(gs[2, 2]) # depth diffs 
    ax6 = fig.add_subplot(gs[2, 3]) # periodogram 
    ax7 = fig.add_subplot(gs[:, 4]) # notes 

    # phase folded
    folded_phase = tls_results.folded_phase 
    folded_y = tls_results.folded_y
#######################################################################bottom left corner
    ##s=1?
    ax3.scatter(folded_phase, folded_y, s=3, c='grey')
        
    ax3.plot(tls_results.model_folded_phase, tls_results.model_folded_model, color='red')

    mask = np.logical_and(folded_phase<0.53, folded_phase>0.47)

    binned_time = binned_statistic(folded_phase[mask], folded_phase[mask], bins=20)[0]
    binned_flux = binned_statistic(folded_phase[mask], folded_y[mask], bins=20)[0]

    ax3.scatter(binned_time, binned_flux, s=35, c='orange', edgecolor='black')

    ax3.set(xlim=(0.47, 0.53))

    # transit depths (odd, even)
    yerr = tls_results.transit_depths_uncertainties 

    yerr = yerr[transits_not_nan_mask]

    ax4.errorbar(x=transit_times, y=transit_depths, yerr=yerr, fmt='o', color='red')
    transit_x = [clean_time.min(), clean_time.max()]
    transit_base = 2*[np.mean(transit_depths)]
    ax4.plot(transit_x, transit_base, color='black', linestyle='dashed')
    ax4.plot(transit_x, 2*[1], color='black')
    ax4.set(xlabel='Time (days)', ylabel='Flux')

    ax4.xaxis.set_major_locator(plt.NullLocator())
#####################################################################################################
    # even odd transits # 

    even_transit_time_folded, even_transit_flux, odd_transit_time_folded, odd_transit_flux, even_indices, odd_indices = even_odd_phase_folded(time=clean_time, flux=clean_flux, results=tls_results)    
    
    if even_odd_option == 'separate': 
        even_transit_time_folded = normalize(even_transit_time_folded, output_range=[1,2])
        # I'm really proud of this because I was having so many problems fixing left right ... sometimes they were negative, sometimes positive, sometimes mixed...sometimes adding the 1.1*max(time) variable would make them super spead out, etc. 
        odd_transit_time_folded = normalize(odd_transit_time_folded)

    ax5.scatter(even_transit_time_folded, even_transit_flux, label='Even')
    ax5.scatter(odd_transit_time_folded, odd_transit_flux, label='Odd')
    ax5.get_xaxis().set_ticks([])
    ax5.set(xlabel='Time (d)', ylabel='Detrended Flux')
    ax5.legend()

    # periodogram #

    for n in range(2, 10):
        ax6.axvline(n*tls_results.period, alpha=0.4, lw=1, linestyle="dashed")
        ax6.axvline(tls_results.period / n, alpha=0.4, lw=1, linestyle="dashed")

    ax6.plot(tls_results.periods, tls_results.power, color='black', lw=0.5)

    ax6.set(xlim=(np.min(tls_results.periods), np.max(tls_results.periods)), 
            xlabel='Period (days)', ylabel='SDE')

    labels = ['period', 'depth', 'T0', 
                    'SDE', 'snr', 'rp/rs', 'transit_count', 
                    'distinct_transit_count']

    values = [tls_results.period, tls_results.depth, tls_results.T0, 
                    tls_results.SDE, tls_results.snr, tls_results.rp_rs, tls_results.transit_count, 
                    tls_results.distinct_transit_count]

    simbad_header, simbad_values = query_simbad(tic_id=tic_id)

    labels+=simbad_header 

    values+=simbad_values

    if false_alarms_dictionary is not None: 
        labels += list(false_alarms_dictionary.keys())
        values += list(false_alarms_dictionary.values())
    
    labels.append('cluster_name')
    cluster_name=return_kerr_cluster(tic_id)
    values.append(cluster_name)
    text_info = []
    for label, value in zip(labels, values):
        if type(value) is str and '|' in value: 
            value = '\n'+value.replace('|','\n')
        
        elif type(value) is bool: 
            value = str(value)
        
        else: 
            try: 
                value = str(round(value, 3))
            except: 
                value = str(value) 
            
        text_info.append(label+'='+value)

    ax7.text(x=0.1, y=0.5, s='\n\n'.join(str(i).replace('_',' ') for i in text_info), fontsize='xx-large', va='center', transform=ax7.transAxes)
    ax7.tick_params(labelbottom=False, labelleft=False, axis='both', which='both', length=0)
    ax7.spines['right'].set_visible(False)

    if plot_dir is None and plot_path is None:
        plt.show()
        plt.clf()
        plt.close()
    
    else: 
        if plot_path is not None: 
            plt.savefig(plot_path, dpi=dpi)
            plt.clf()
            plt.close()
        elif plot_dir is not None: 
            if plot_dir[-1]!='/': 
                plot_dir += '/'
            
            if plot_type=='png' or plot_type == 'pdf': 
                plot_path = plot_dir + tic_id + '.'+plot_type 
            else: 
                plot_path = plot_dir + tic_id + '.png'

            plt.savefig(plot_path, dpi=dpi)
            plt.clf()
            plt.close()
