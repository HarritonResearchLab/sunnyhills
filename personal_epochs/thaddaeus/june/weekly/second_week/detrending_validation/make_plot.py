from matplotlib.pyplot import xlabel


def plot_detrend_validation(tic_id, data_dir:str, plot_dir:str=None): 
    
    r'''
    Arguments
    ---------
    `tic_id`: tic id
    `data_dir`: directory to light curve data files
    `plot_dir`: directory you want to save light curve files to 
    '''
    
    import pandas as pd
    import numpy as np
    import matplotlib 
    import matplotlib.pyplot as plt

    if data_dir[-1]!='/': 
        data_dir += '/'

    data = pd.read_csv(data_dir+tic_id+'.csv')

    raw_time, raw_flux = (np.array(data[i]) for i in ['raw_time', 'raw_flux'])

    no_flare_raw_time, no_flare_raw_flux = (np.array(data[i]) for i in ['no_flare_raw_time', 'no_flare_raw_flux'])
    mask = np.isfinite(no_flare_raw_time) 
    no_flare_raw_time, no_flare_raw_flux = (i[mask] for i in [no_flare_raw_time, no_flare_raw_flux])

    trend_time, trend_flux =  (np.array(data[i]) for i in ['trend_time', 'trend_flux'])
    mask = np.isfinite(trend_flux)
    trend_time, trend_flux = (i[mask] for i in [trend_time, trend_flux])

    clean_time, clean_flux = (np.array(data[i]) for i in ['clean_time', 'clean_flux'])
    mask = np.isfinite(clean_time)
    clean_time, clean_flux = (i[mask] for i in [clean_time, clean_flux])

    diffs = np.diff(raw_time)
    start_indices = [0]+list(np.where(diffs>25)[0]+1)
    inclusive_time_ranges = []

    for pos, index in enumerate(start_indices):
        if pos<len(start_indices)-1:
            inclusive_time_range = (raw_time[index], raw_time[start_indices[pos+1]-1])
        else:
            inclusive_time_range = (raw_time[index], raw_time[-1])

        inclusive_time_ranges.append(inclusive_time_range)

    num_sectors = len(inclusive_time_ranges)

    plt.style.use('seaborn-darkgrid')
    font = {'family' : 'serif', 'size' : 5}

    matplotlib.rc('font', **font)

    fig, axs = plt.subplots(3,num_sectors, figsize=(num_sectors*3.5, 4))

    if num_sectors>1: 

        # raw with flares #
        for i in range(num_sectors):

            inclusive_range = inclusive_time_ranges[i]

            # raw # 
            ax = axs[0, i]
            raw_mask = np.logical_and(raw_time>=inclusive_range[0], 
                                    raw_time<=inclusive_range[1])

            ax.scatter(raw_time[raw_mask], raw_flux[raw_mask], s=0.25, color='#408ee0')
            
            ax.set(xlabel='Time (d)', xlim=inclusive_range)
            if i==0: 
                ax.set(ylabel='Raw Flux', xlabel='Time (d)', xlim=inclusive_range)
            else: 
                ax.set(xlim=inclusive_range)
                plt.setp(ax.get_yticklabels(), visible=True)
                plt.setp(ax.get_xticklabels(), visible=True)

            # no flare and trend # 

            ax = axs[1, i]

            no_flare_mask = np.logical_and(no_flare_raw_time>=inclusive_range[0], 
                                        no_flare_raw_time<=inclusive_range[1])

            ax.scatter(no_flare_raw_time[no_flare_mask], no_flare_raw_flux[no_flare_mask], s=0.25, color='#408ee0')
            
            trend_mask = np.logical_and(trend_time>=inclusive_range[0], 
                                        trend_time<=inclusive_range[1])

            ax.plot(trend_time[trend_mask], trend_flux[trend_mask], c='red', lw=0.5)
            
            ax.set(xlabel='Time (d)', xlim=inclusive_range)
            if i==0: 
                ax.set(ylabel='No Flare Flux', xlabel='Time (d)', xlim=inclusive_range)
            else: 
                ax.set(xlim=inclusive_range)
                plt.setp(ax.get_yticklabels(), visible=True)
                plt.setp(ax.get_xticklabels(), visible=True)

            # clean # 

            ax = axs[2,i]
            clean_mask = np.logical_and(clean_time>=inclusive_range[0], 
                                        clean_time<=inclusive_range[1])

            ax.scatter(clean_time[clean_mask], clean_flux[clean_mask], s=0.25, color='#408ee0') 
            
            ax.set(xlabel='Time (d)', xlim=inclusive_range)
            if i==0: 
                ax.set(ylabel='Detrended Clean Flux')
            else: 
                ax.set(xlim=inclusive_range)
                plt.setp(ax.get_yticklabels(), visible=True)
                plt.setp(ax.get_xticklabels(), visible=True)

    else: 
            # raw with flares #
        for i in range(num_sectors):

            inclusive_range = inclusive_time_ranges[i]

            # raw # 
            ax = axs[0]
            raw_mask = np.logical_and(raw_time>=inclusive_range[0], 
                                    raw_time<=inclusive_range[1])

            ax.scatter(raw_time[raw_mask], raw_flux[raw_mask], s=0.25, color='#408ee0')
            
            ax.set(xlabel='Time (d)', xlim=inclusive_range)
            if i==0: 
                ax.set(ylabel='Raw Flux', xlabel='Time (d)', xlim=inclusive_range)
            else: 
                ax.set(xlim=inclusive_range)
                plt.setp(ax.get_yticklabels(), visible=True)
                plt.setp(ax.get_xticklabels(), visible=True)

            # no flare and trend # 

            ax = axs[1]

            no_flare_mask = np.logical_and(no_flare_raw_time>=inclusive_range[0], 
                                        no_flare_raw_time<=inclusive_range[1])

            ax.scatter(no_flare_raw_time[no_flare_mask], no_flare_raw_flux[no_flare_mask], s=0.25, color='#408ee0')
            
            trend_mask = np.logical_and(trend_time>=inclusive_range[0], 
                                        trend_time<=inclusive_range[1])

            ax.plot(trend_time[trend_mask], trend_flux[trend_mask], c='red', lw=0.5)
            
            ax.set(xlabel='Time (d)', xlim=inclusive_range)
            if i==0: 
                ax.set(ylabel='No Flare Flux', xlabel='Time (d)', xlim=inclusive_range)
            else: 
                ax.set(xlim=inclusive_range)
                plt.setp(ax.get_yticklabels(), visible=True)
                plt.setp(ax.get_xticklabels(), visible=True)

            # clean # 

            ax = axs[2]
            clean_mask = np.logical_and(clean_time>=inclusive_range[0], 
                                        clean_time<=inclusive_range[1])

            ax.scatter(clean_time[clean_mask], clean_flux[clean_mask], s=0.25, color='#408ee0') 
            
            ax.set(xlabel='Time (d)', xlim=inclusive_range)
            if i==0: 
                ax.set(ylabel='Detrended Clean Flux')
            else: 
                ax.set(xlim=inclusive_range)
                plt.setp(ax.get_yticklabels(), visible=True)
                plt.setp(ax.get_xticklabels(), visible=True) 
            
    #plt.tight_layout()
    plt.subplots_adjust(wspace=0.1, hspace=0.45)

    if plot_dir!=None: 
        if plot_dir[-1]!='/': 
            plot_dir+='/'

        plot_path = plot_dir+tic_id+'.png'
        plt.savefig(plot_path, dpi=250)

data_dir = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/data/current/processed/two_min_lightcurves'
plot_dir = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/personal_epochs/thaddaeus/june/weekly/second_week/detrending_validation/plots'
plot_detrend_validation(tic_id='TIC_405040121', data_dir=data_dir, plot_dir=plot_dir)