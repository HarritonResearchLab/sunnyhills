"""
Re-useable plotting functions can go in this file.  They can be called using
driver scripts.

Contents:
    plot_kerr21_XY
    plot_au_mic_detrending
"""
import os
from glob import glob
import numpy as np, matplotlib.pyplot as plt, pandas as pd
import matplotlib as mpl

from numpy import array as nparr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from astropy import units as u, constants as const
from astropy.coordinates import SkyCoord
from astropy.table import Table
import numpy as np
import warnings

#from sunnyhills.paths import DATADIR, EPOCHSDIR

def plot_kerr21_XY(outdir, colorkey=None):
    """
    Make a top-down plot like Figure 7 of Kerr+2021
    (https://ui.adsabs.harvard.edu/abs/2021ApJ...917...23K/abstract).

    Optionally color by stellar age, similar to Fig 24 of McBride+2021.
    (https://ui.adsabs.harvard.edu/abs/2021AJ....162..282M/abstract).

    Args:

        outdir (str): path to output directory

        colorkey (str or None): column string from Kerr+2021 table to color
        points by: for instance "Age", "bp-rp", "plx", or "P" (latter means
        "Probability of star being <50 Myr old).
    """

    assert isinstance(outdir, str)
    assert isinstance(colorkey, str) or (colorkey is None)

    #
    # get data, calculate galactic X,Y,Z positions (assuming well-measured
    # parallaxes).
    #
    kerrpath = os.path.join(DATADIR, "Kerr_2021_Table1.txt")
    df = Table.read(kerrpath, format='cds').to_pandas()

    from sunnyhills.physicalpositions import calculate_XYZ_given_RADECPLX
    x,y,z = calculate_XYZ_given_RADECPLX(df.RAdeg, df.DEdeg, df.plx)

    #
    # make the plot
    #
    fig, ax = plt.subplots(figsize=(4,4))

    x_sun, y_sun = -8122, 0
    ax.scatter(
        x_sun, y_sun, c='black', alpha=1, zorder=1, s=20, rasterized=True,
        linewidths=1, marker='x'
    )

    if colorkey is None:
        # By default, just show all the stars as the same color.  The
        # "rasterized=True" kwarg here is good if you save the plots as pdfs,
        # to not need to save the positions of too many points.
        ax.scatter(
            x, y, c='black', alpha=1, zorder=2, s=2, rasterized=True,
            linewidths=0, marker='.'
        )

    else:
        # Add a colorbar.
        color = df[colorkey]

        # Only show points for which color is defined (e.g., not all stars have
        # ages reported in the table).
        sel = ~pd.isnull(color)

        # Define the colormap.  See
        # https://matplotlib.org/stable/tutorials/colors/colormaps.html, there
        # are better choices, but this one is OK for comparing against
        # McBride+21.
        cmap = mpl.cm.get_cmap('rainbow')

        _p = ax.scatter(
            x[sel], y[sel], c=color[sel], alpha=1, zorder=3, s=2, rasterized=True,
            linewidths=0, marker='.', cmap=cmap
        )

        # For the colorbar, inset it into the main plot to keep the square
        # aspect ratio.
        axins1 = inset_axes(ax, width="3%", height="20%", loc='lower right',
                            borderpad=0.7)

        cb = fig.colorbar(_p, cax=axins1, orientation="vertical",
                          extend="neither")
        cb.ax.tick_params(labelsize='x-small')
        cb.ax.yaxis.set_ticks_position('left')
        cb.ax.yaxis.set_label_position('left')

        KEYLABELDICT = {
            'Age': 'Age [years]',
            'P': 'P$_{\mathrm{<50\ Myr}}$',
            'plx': 'Parallax [mas]',
        }
        if colorkey in KEYLABELDICT:
            cb.set_label(KEYLABELDICT[colorkey], fontsize='x-small')
        else:
            cb.set_label(colorkey, fontsize='x-small')

    ax.set_xlabel("X [pc]")
    ax.set_ylabel("Y [pc]")

    s = ''
    if colorkey:
        s += f"_{colorkey}"

    outpath = os.path.join(outdir, f'kerr21_XY{s}.png')
    fig.savefig(outpath, bbox_inches='tight', dpi=400)
    print(f"Made {outpath}")

def plot_star_detrending(
    outdir: str,
    ticstr: str,
    dtrdict: dict = {'method':'biweight',
                     'window_length':0.5,
                     'cval':5.0,
                     "break_tolerance":1.0},
    known_ephemeris: dict = None
    ) -> None:
    """
    Args:
        outdir: directory where plots will be written.

        ticstr: e.g., 'TIC 441420236'.  will be passed to lightkurve and used
        in naming plots.

        dtrdict: dictionary with keys "window_length", "method", "cval",
            "break_tolerance", or anything else needed by wotan's `flatten`
            call.  These are documented at
            https://wotan.readthedocs.io/en/latest/Usage.html

    Kwargs:
        known_ephemeris: optional dictionary with keys 't0', 'period', to
        overplot the known transit times.
    """

    import lightkurve as lk
    from wotan import flatten, slide_clip

    plotpath = os.path.join(
        outdir,
        f'{ticstr.replace(" ","_")}_{dtrdict["method"]}_'+
        f'window{dtrdict["window_length"]:.2f}_cval{dtrdict["cval"]:.1f}_'+
        f'breaktol{dtrdict["break_tolerance"]:.1f}.png'
    )
    if os.path.exists(plotpath):
        print(f"Found {plotpath}, skipping.")
        return 1

    # get the light curve
    lcc = lk.search_lightcurve(ticstr).download_all()

    # select only the two-minute cadence SPOC-reduced data; convert to a list.
    # note that this conversion approach works for any LightCurveCollection
    # returned by lightkurve -- no need to hand-pick the right ones.  the exact
    # condition below says "if the interval is between 119 and 121 seconds,
    # take it".
    lc_list = [_l for _l in lcc
         if
         _l.meta['ORIGIN']=='NASA/Ames'
         and
         np.isclose(
             120,
             np.nanmedian(np.diff(_l.remove_outliers().time.value))*24*60*60,
             atol=1
         )
    ]

    # how many sectors of 2-minute cadence TESS data (processed by the Ames
    # group) are available for this object?
    N_sectors = len(lc_list)

    fig, axs = plt.subplots(nrows=2*N_sectors, ncols=1, figsize=(12,4*N_sectors))

    titlestr = os.path.basename(plotpath).replace('.png','').replace('_',' ')
    axs[0].set_title(titlestr, fontsize='small')

    for ix, lc in enumerate(lc_list):

        time = lc.time.value
        flux = lc.pdcsap_flux.value
        qual = lc.quality.value

        # remove non-zero quality flags
        sel = (qual == 0)

        time = time[sel]
        flux = flux[sel]

        # normalize around 1
        flux /= np.nanmedian(flux)

        # remove outliers before local window detrending
        clipped_flux = slide_clip(time, flux, window_length=0.5, low=3,
                                  high=2, method='mad', center='median')

        # see https://wotan.readthedocs.io/en/latest/Usage.html for other
        # possible options.
        flat_flux, trend_flux = flatten(
            time, clipped_flux, return_trend=True,
            method=dtrdict['method'],
            break_tolerance=dtrdict['break_tolerance'],
            window_length=dtrdict['window_length'],
            cval=dtrdict['cval']
        )

        # top axis: the data and model
        axs[2*ix].scatter(time, flux, s=1, c='black', zorder=2)
        axs[2*ix].plot(time, trend_flux, lw=1, c='red',alpha=0.6, zorder=3)

        # bottom axis: the residual
        axs[2*ix+1].scatter(time, flat_flux, s=1, c='black', zorder=2)

    # assume 't0' in known_ephemeris dictionary is given in same time-system as
    # "time".  then make a bunch of vertical lines, but set the x-axis and
    # y-axis limits to hide them.
    for ax in axs:
        if not isinstance(known_ephemeris, dict):
            break
        t0 = known_ephemeris['t0']
        period = known_ephemeris['period']
        epochs = np.arange(-1000,1000,1)
        tra_times = t0 + period*epochs
        xlim, ylim = ax.get_xlim(), ax.get_ylim()
        ax.vlines(
            tra_times, ylim[0], ylim[1], colors='darkgray', alpha=0.8,
            linestyles='--', zorder=-2, linewidths=0.5
        )
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    fig.text(-0.01,0.5, 'Relative flux', va='center', rotation=90,
             fontsize='large')
    fig.text(0.5,-0.01, 'Time [BTJD]', va='center', ha='center',
             fontsize='large')

    fig.tight_layout()

    fig.savefig(plotpath, bbox_inches='tight', dpi=400)
    print(f"Made {plotpath}")

def bls_validation_mosaic(tic_id:str, clean_time:np.array, clean_flux:np.array, 
                          trend_time:np.array, trend_flux:np.array,
                          raw_time:np.array, raw_flux:np.array, 
                          best_params:list, bls_model, in_transit, bls_stats, 
                          path:str=None, dpi:int=150): 

    '''
    arguments: 
        tic_id: tic id 
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
    '''
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from sunnyhills.borrowed import tls_intransit_stats
    from sunnyhills.misc import phase, rebin
    import numpy as np

    plt.rcParams['font.family']='serif'
    plt.style.use('seaborn-darkgrid')
    plt.rcParams['font.family']='serif'
    fig = plt.figure(constrained_layout=True, figsize=(12,12))

    gs = GridSpec(4, 3, figure=fig)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, :])
    ax3 = fig.add_subplot(gs[2, :-1])
    ax4 = fig.add_subplot(gs[-1, :-1])
    ax5 = fig.add_subplot(gs[2:, -1])

    # raw and trend light curve
        
    ax1.scatter(raw_time, raw_flux, s=1)
    ax1.scatter(trend_time, trend_flux, s=1, c='r')
    ax1.set(ylabel='Flux')
    
    # detrend light curve
    ax2.scatter(clean_time, clean_flux, s=1)
    period = best_params[0]
    t0 = best_params[1]
    duration = best_params[2]  

    params = {
      'power': bls_model.max_power.value,
      'period': period,
      't0': t0,
      'duration': duration,
      'depth': best_params[3],
      'transit_time': np.argmax(bls_stats['transit_times'])
    }  


    phased_time, phased_flux, (x, f) = phase(time=clean_time, flux=clean_flux, period=period, 
                                           t0=t0, duration=duration, bls_model=bls_model, model_name='BLS')
    ax2.vlines(clean_time[in_transit], min(clean_flux), max(clean_flux), color='red', lw=0.05, alpha=0.4, zorder=0)
    ax2.set(ylabel='Detrended Flux')

    for ax in [ax1, ax2]: 
        ax.set(xlabel='Time (days)')

    # phase folded
    ax3.scatter(phased_time, phased_flux, s=3, c='grey')

    binned_x, binned_flux, success = rebin(phased_time, phased_flux)
    
 
    #x3.plot(x.value, f.value, color='red', alpha=0.5)
    if success: 
        ax3.scatter(binned_x, binned_flux, c='orange', s=40, edgecolor='black')

    ax3.set(xlim=(-0.2, 0.2))

    for ax in [ax3,ax4]: 
        ax.set(xlabel='Time from mid-transit (days)', ylabel='Detrended Flux')



    # periodogram 
    bls_model.plot(ax=ax4,xlabel='Period (d)',ylabel='Power',style='seaborn-darkgrid')

    ax5.tick_params(labelbottom=False, labelleft=False, axis='both', which='both', length=0)

    text_info = []
    for key_name,val in params.items(): 
      val = round(float(val),5)
      result = key_name
      if type(result)!=str: 
          result = str(round(result, 5))

      text_info.append(key_name+': '+ str(val)+'\n')

    ax5.text(x=0.1, y=0.5, s='\n'.join(str(i).replace('_','') for i in text_info), fontsize='large', va='center', transform=ax5.transAxes)

    ax1.set_title('TIC: '+str(tic_id).replace('_','')+' PERIOD: '+str(round(period, 5)), size='xx-large')

    if path==None:
        plt.show()
    else: 
        plt.savefig(path, dpi=dpi)

def tls_validation_mosaic(tic_id:str, data, tls_model, tls_results,
                          plot_dir:str=None, dpi:int=150): 

    '''
    arguments: 
        tic_id: tic id 
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
    '''

    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec
    from sunnyhills.misc import rebin 
    from sunnyhills.false_alarm_checks import even_odd_phase_folded
    from transitleastsquares import (
        transitleastsquares,
        cleaned_array,
        catalog_info,
        transit_mask
    )

    #plt.style.use('https://raw.githubusercontent.com/thissop/MAXI-J1535/main/code/misc/stolen_science.mplstyle?token=GHSAT0AAAAAABP54PQO2X2VXMNS256IWOBOYRNCFBA')
    #plt.style.use('https://raw.githubusercontent.com/thissop/MAXI-J1535/main/code/misc/stolen_science.mplstyle?token=GHSAT0AAAAAABVPXDLXKPLDOD6CCSMC3WMSYVDQ6XA')
    plt.style.use('seaborn-darkgrid')
    plt.rcParams['font.family']='serif'
    fig = plt.figure(constrained_layout=True, figsize=(12,12))

    gs = GridSpec(4, 3, figure=fig)
    ax1 = fig.add_subplot(gs[0, :]) # detrended light curve
    ax2 = fig.add_subplot(gs[1, :-1]) # no flare with trend light curve
    ax3 = fig.add_subplot(gs[2, :-1]) # phase folded transit 
    ax4 = fig.add_subplot(gs[-1, 0]) # left right transits
    ax5 = fig.add_subplot(gs[-1, -2]) # depth diffs 
    ax6 = fig.add_subplot(gs[-1,-1]) # periodogram 
    ax7 = fig.add_subplot(gs[1:-1, -1]) # notes 

    
    df = pd.read_csv(data)
    clean_time, clean_flux = (np.array(df[i]) for i in ['clean_time', 'clean_flux'])
    clean_mask = np.isfinite(clean_time)
    clean_time, clean_flux = (i[clean_mask] for i in [clean_time, clean_flux])
    in_transit = transit_mask(clean_time, tls_results.period, tls_results.duration, tls_results.T0)

    raw_time, raw_flux = (np.array(df[i]) for i in ['no_flare_raw_time','no_flare_raw_flux'])

    trend_time, trend_flux = (np.array(df[i]) for i in ['trend_time','trend_flux'])
    trend_mask = np.isfinite(trend_time)
    trend_time, trend_flux = (i[trend_mask] for i in [trend_time, trend_flux])

    # detrend light curve
    ax1.scatter(clean_time, clean_flux, s=1)
    
    in_transit = transit_mask(clean_time, tls_results.period, tls_results.duration, tls_results.T0)
    ax1.scatter(clean_time[in_transit], clean_flux[in_transit], color='red', s=2, zorder=0)
    ax1.scatter(clean_time[~in_transit], clean_flux[~in_transit], color='grey', alpha=0.5, s=2, zorder=0)
    ax1.plot(tls_results.model_lightcurve_time, tls_results.model_lightcurve_model, alpha=0.5, color='red', zorder=1)
    ax1.set(ylabel='Detrended Flux')

    # raw and trend light curve
        
    ax2.scatter(raw_time, raw_flux, s=1)
    ax2.plot(trend_time, trend_flux, lw=0.5, c='r')
    ax2.set(ylabel='Flux')
    
    lc_xlim = (min((min(clean_time), min(raw_time))), max((max(clean_time), max(raw_time))))
    for ax in [ax1, ax2]: 
        ax.set(xlabel='Time (days)', xlim=lc_xlim)

    # phase folded
    ax3.scatter(tls_results.folded_phase, tls_results.folded_y, s=3, c='grey')
    
    binned_x, binned_flux, success = rebin(tls_results.folded_phase, tls_results.folded_y)
    
    ax3.plot(tls_results.model_folded_phase, tls_results.model_folded_model, color='red')

    ax3.scatter(binned_x, binned_flux, s=9, c='orange', edgecolor='black')

    ax3.set(xlim=(0.47, 0.53))

    # transit depths (odd, even)

    # https://github.com/hippke/tls/blob/71da590d3e199264822db425ab9f9f633253986e/transitleastsquares/stats.py#L338

    transit_times = np.array(tls_results.transit_times) 
    transit_depths = tls_results.transit_depths
    yerr = tls_results.transit_depths_uncertainties 

    mask = np.isfinite(transit_depths)

    transit_times, transit_depths, yerr = (i[mask] for i in [transit_times, transit_depths, yerr])
    
    ax4.errorbar(x=transit_times, y=transit_depths, yerr=yerr, fmt='o', color='red')
    transit_x = [clean_time.min(), clean_time.max()]
    transit_base = 2*[np.min(transit_depths)]
    ax4.plot(transit_x, transit_base, color='black', linestyle='dashed')
    ax4.plot(transit_x, 2*[1], color='black')
    ax4.set(xlabel='Time (days)', ylabel='Flux')
    
    ax4.xaxis.set_major_locator(plt.NullLocator())

    # even odd transits
    #ax6.

    even_transit_time_folded, even_transit_flux, odd_transit_time_folded, odd_transit_flux, even_indices, odd_indices = even_odd_phase_folded(time=clean_time, flux=clean_flux, results=tls_results)    
    ax5.scatter(even_transit_time_folded, even_transit_flux, label='Even')
    #ax5.scatter(tls_results.folded_phase[even_indices], tls_results.folded_y[even_indices])
    max_even = np.max(even_transit_time_folded)
    shifted_odd_time = odd_transit_time_folded+max_even
    ax5.scatter(shifted_odd_time, odd_transit_flux, label='Odd')
    ax5.get_xaxis().set_ticks([])
    ax5.set(xlabel='Time (d)', ylabel='Detrended Flux')
    ax5.legend()

    # periodogram 
    #ax6.plot(bls_results.period, bls_results.power)
    #ax6.axvline(tls_results.period, alpha=0.4, lw=3)
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

    text_info = []
    for label, value in zip(labels, values):
        text_info.append(label+'='+str(round(value, 5)))

    ax7.text(x=0.1, y=0.5, s='\n\n'.join(str(i).replace('_',' ') for i in text_info), fontsize='large', va='center', transform=ax7.transAxes)
    ax7.tick_params(labelbottom=False, labelleft=False, axis='both', which='both', length=0)

    ax1.set_title('TIC: '+str(tic_id).replace('_','')+' PERIOD: '+str(round(tls_results.period, 5)), size='xx-large')

    if plot_dir==None:
        plt.show()
    else: 
        if plot_dir[-1]!='/': 
            plot_dir += '/'
        plot_path = plot_dir + tic_id + '.png'
        plt.savefig(plot_path, dpi=dpi)

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