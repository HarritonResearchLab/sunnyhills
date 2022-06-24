import matplotlib.pylab as plt
import numpy as np
from matplotlib.gridspec import GridSpec
import pandas as pd 
from sunnyhills.pipeline_functions import run_tls 
from scipy.stats import binned_statistic
from sunnyhills.misc import even_odd_phase_folded

tic_id = 'TIC_17417151'

df = pd.read_csv('./routines/alpha_tls/data/two_min_lightcurves/'+tic_id+'.csv')

clean_time, clean_flux = (np.array(df[i]) for i in ['clean_time', 'clean_flux'])
clean_mask = np.isfinite(clean_time)
clean_time, clean_flux = (i[clean_mask] for i in [clean_time, clean_flux])

# delete this!! 

tls_results, tls_model = run_tls(tic_id, clean_time, clean_flux)

raw_mask = np.isfinite(df['no_flare_raw_time'])
raw_time, raw_flux = (np.array(df[i])[raw_mask] for i in ['no_flare_raw_time','no_flare_raw_flux'])

trend_time, trend_flux = (np.array(df[i]) for i in ['trend_time','trend_flux'])
trend_mask = np.isfinite(trend_time)
trend_time, trend_flux = (i[trend_mask] for i in [trend_time, trend_flux])

split_axes = False
break_index = None 

diff = np.diff(clean_time)
if np.max(diff>50): 
    split_axes = True 
    break_index = np.where(np.diff(clean_time)>50)[0][0]+1 # so it's [0:break_index), [break_index:]
    break_index = int(break_index)  

## PLOTTING ### 

#plt.style.use('seaborn-darkgrid')
plt.rcParams['font.family']='serif'
fig = plt.figure(constrained_layout=True, figsize=(18,12))

gs = GridSpec(3, 5, figure=fig)

def orient_split_axes(ax_1, ax_2, flux): 
    ax_1.spines['right'].set_visible(False)
    ax_2.spines['left'].set_visible(False)
    ax_1.yaxis.tick_left()
    ax_2.yaxis.tick_right()

    ylim = [0.95*np.min(flux), 1.05*np.max(flux)]

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
    
    ax1a.scatter(clean_time[0:break_index], clean_flux[0:break_index], s=1)
    ax1a.set(ylabel='Detrended Flux')
    ax1b.scatter(clean_time[break_index:], clean_flux[break_index:], s=1)
    
    orient_split_axes(ax1a, ax1b, clean_flux)

    ax2a.scatter(raw_time[0:break_index], raw_flux[0:break_index], s=1)
    ax2a.plot(trend_time[0:break_index], trend_flux[0:break_index], lw=0.5, c='r')
    ax2a.set(ylabel='Flux')
    ax2b.scatter(raw_time[break_index:], raw_flux[break_index:], s=1)
    ax2b.plot(trend_time[break_index:], trend_flux[break_index:], lw=0.5, c='r')
    
    orient_split_axes(ax2a, ax2b, raw_flux)

    for ax in [ax1a, ax1b, ax2a, ax2b]: 
        ax.set(xlabel='Time (days)')

else: 
    ax1 = fig.add_subplot(gs[0, 0:-1]) # detrended light curve
    ax2 = fig.add_subplot(gs[1, 0:-1]) # no flare with trend light curve
   
    # detrend light curve 
    ax1.scatter(clean_time, clean_flux, s=1)
    ax1.set(ylabel='Detrended Flux')
    
    # raw and trend light curve # 
    ax2.scatter(raw_time, raw_flux, s=1)
    ax2.plot(trend_time, trend_flux, lw=0.5, c='r')
    ax2.set(ylabel='Flux')

    for ax in [ax1, ax2]: 
            ax.set(xlabel='Time (days)')

ax3 = fig.add_subplot(gs[2, 0]) # phase folded transit 
ax4 = fig.add_subplot(gs[2, 1]) # left right transits
ax5 = fig.add_subplot(gs[2, 2]) # depth diffs 
ax6 = fig.add_subplot(gs[2, 3]) # periodogram 
ax7 = fig.add_subplot(gs[:, 4]) # notes 

# phase folded
folded_phase = tls_results.folded_phase 
folded_y = tls_results.folded_y
ax3.scatter(folded_phase, folded_y, s=3, c='grey')
    
ax3.plot(tls_results.model_folded_phase, tls_results.model_folded_model, color='red')

mask = np.logical_and(folded_phase<0.53, folded_phase>0.47)

binned_time = binned_statistic(folded_phase[mask], folded_phase[mask], bins=20)[0]
binned_flux = binned_statistic(folded_phase[mask], folded_y[mask], bins=20)[0]

ax3.scatter(binned_time, binned_flux, s=35, c='orange', edgecolor='black')

ax3.set(xlim=(0.47, 0.53))

# transit depths (odd, even)

transit_times = np.array(tls_results.transit_times) 
transit_depths = tls_results.transit_depths
yerr = tls_results.transit_depths_uncertainties 

mask = np.isfinite(transit_depths)

transit_times, transit_depths, yerr = (i[mask] for i in [transit_times, transit_depths, yerr])

ax4.errorbar(x=transit_times, y=transit_depths, yerr=yerr, fmt='o', color='red')
transit_x = [clean_time.min(), clean_time.max()]
transit_base = 2*[np.mean(transit_depths)]
ax4.plot(transit_x, transit_base, color='black', linestyle='dashed')
ax4.plot(transit_x, 2*[1], color='black')
ax4.set(xlabel='Time (days)', ylabel='Flux')

ax4.xaxis.set_major_locator(plt.NullLocator())

# even odd transits # 

even_transit_time_folded, even_transit_flux, odd_transit_time_folded, odd_transit_flux, even_indices, odd_indices = even_odd_phase_folded(time=clean_time, flux=clean_flux, results=tls_results)    
ax5.scatter(even_transit_time_folded, even_transit_flux, label='Even')
max_even = np.max(even_transit_time_folded)
shifted_odd_time = odd_transit_time_folded+max_even
ax5.scatter(shifted_odd_time, odd_transit_flux, label='Odd')
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

text_info = []
for label, value in zip(labels, values):
    text_info.append(label+'='+str(round(value, 5)))

ax7.text(x=0.1, y=0.5, s='\n\n'.join(str(i).replace('_',' ') for i in text_info), fontsize='xx-large', va='center', transform=ax7.transAxes)
ax7.tick_params(labelbottom=False, labelleft=False, axis='both', which='both', length=0)
ax7.spines['right'].set_visible(False)
ax1.set_title('TIC: '+str(tic_id).replace('_','')+' PERIOD: '+str(round(tls_results.period, 5)), size='xx-large')

plt.savefig('./personal_epochs/thaddaeus/june/topical/misc/img.png')