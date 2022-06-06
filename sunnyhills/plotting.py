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

import numpy as np
def bls_validation_mosaic(tic_id:str, clean_time:np.array, detrend_flux:np.array, 
                          raw_time:np.array, raw_flux:np.array, 
                          best_params:list, bls_results, bls_model, in_transit, bls_stats, 
                          path:str=None, dpi:int=150): 

    '''
    arguments: 
        tic_id: tic id 
        clean_time: detrended and flare removed time array
        detrend_flux: flux values corresponding to clean_time arg
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

    import astropy.units as units
    from lightkurve.periodogram import Periodogram

    plt.style.use('https://raw.githubusercontent.com/thissop/MAXI-J1535/main/code/misc/stolen_science.mplstyle?token=GHSAT0AAAAAABP54PQO2X2VXMNS256IWOBOYRNCFBA')
    fig = plt.figure(constrained_layout=True, figsize=(12,12))

    gs = GridSpec(4, 3, figure=fig)
    ax1 = fig.add_subplot(gs[0, :])
    ax2 = fig.add_subplot(gs[1, :])
    ax3 = fig.add_subplot(gs[2, :-1])
    ax4 = fig.add_subplot(gs[-1, 0])
    ax5 = fig.add_subplot(gs[-1, -2])
    ax6 = fig.add_subplot(gs[2:, -1])

    
    # raw and trend light curve
    p = Periodogram(bls_results.period*units.microhertz,units.Quantity(bls_results.power))
    p.flatten()
    p.plot(ax=ax1,xlabel='period',ylabel='power',style='https://raw.githubusercontent.com/thissop/MAXI-J1535/main/code/misc/stolen_science.mplstyle?token=GHSAT0AAAAAABP54PQO2X2VXMNS256IWOBOYRNCFBA')
    '''    
    ax1.scatter(raw_time, raw_flux, s=1)
    ax1.set(ylabel='Flux')
    '''
    # detrend light curve
    ax2.scatter(clean_time, detrend_flux, s=1)
    index = np.argmax(bls_results.power)
    period = bls_results.period[index]
    t0 = bls_results.transit_time[index]
    duration = bls_results.duration[index]


    phased_time, phased_flux, x, f = phase(clean_time, detrend_flux, best_params, bls_model)

    ax2.vlines(clean_time[in_transit], min(detrend_flux), max(detrend_flux), color='red', lw=0.05, alpha=0.4, zorder=0)
    ax2.set(ylabel='Detrended Flux')

    for ax in [ax1, ax2]: 
        ax.set(xlabel='Time (days)')

    # phase folded
    ax3.scatter(phased_time, phased_flux, s=3, c='grey')

    binned_x, binned_flux = rebin(phased_time, phased_flux)
    
    ax3.plot(x, f, color='red', alpha=0.5)
    ax3.scatter(binned_x, binned_flux, c='orange', s=40, edgecolor='black')

    ax3.set(xlim=(-0.2, 0.2))

    for ax in [ax3, ax4, ax5]: 
        ax.set(xlabel='Time from mid-transit (days)', ylabel='Detrended Flux')

    # transit depths (odd, even)

    # https://github.com/hippke/tls/blob/71da590d3e199264822db425ab9f9f633253986e/transitleastsquares/stats.py#L338

    sig_diff = best_params[4]

    intransit_stats = tls_intransit_stats(clean_time, detrend_flux, 
                                          bls_stats['transit_times'], 
                                          best_params[3])

    odd_ = intransit_stats[4]
    odd_time, odd_flux = (odd_[0], odd_[1])
    odd_phased_time, odd_phased_flux, odd_x, odd_f = phase(odd_time, odd_flux, best_params, bls_model)
    odd_binned_time, odd_binned_flux = rebin(odd_phased_time, odd_phased_flux)

    even_ = intransit_stats[5]
    even_time, even_flux = (even_[0], even_[1])
    even_phased_time, even_phased_flux, even_x, even_f = phase(even_time, even_flux, best_params, bls_model)
    even_phased_time = even_phased_time + 2.5*np.max(odd_phased_time) # shift over to show together 
    
    odd_even_median = (np.max(odd_phased_time)+np.min(even_phased_time))/2

    even_binned_time, even_binned_flux = rebin(even_phased_time, even_phased_flux)

    ax4.scatter(odd_phased_time, odd_phased_flux, s=3, c='grey')
    ax4.scatter(even_phased_time, even_phased_flux, s=3, c='grey')
    ax4.scatter(odd_binned_time, odd_binned_flux, c='orange', s=40, edgecolor='black')
    ax4.scatter(even_binned_time, even_binned_flux, c='orange', s=40, edgecolor='black')

    ax4.axvline(x=odd_even_median, color='black', lw=0.5, label='Diff: '+str(round(sig_diff, 5))+r'$\sigma$')

    ax4.legend(loc='upper right', handlelength=0)

    ax4.set(xlabel='Odd (left) and Even (right) Folded Transits')
    ax4.xaxis.set_major_locator(plt.NullLocator())

    # periodogram 
    ax5.plot(bls_results.period, bls_results.power)
    ax5.set(xlabel='Period (d)', ylabel='Power')

    #ax6.axis('off')
    ax6.tick_params(labelbottom=False, labelleft=False, axis='both', which='both', length=0)

    index = np.argmax(bls_results.power)
    text_info = []
    for key_name in bls_results.keys(): 
        result = bls_results[key_name]
        if type(result)!=str: 
            result = str(round(result[index], 5))

        text_info.append(key_name+': '+result+'\n')

    ax6.text(x=0.1, y=0.5, s='\n'.join(str(i).replace('_','') for i in text_info), fontsize='large', va='center', transform=ax6.transAxes)

    ax1.set_title('TIC: '+str(tic_id).replace('_','')+' PERIOD: '+str(round(period, 5)), size='xx-large')

    if path==None:
        plt.show()
    else: 
        plt.savefig(path, dpi=dpi)