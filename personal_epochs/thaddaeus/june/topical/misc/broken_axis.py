import matplotlib.pylab as plt
import numpy as np
from matplotlib.gridspec import GridSpec
import pandas as pd 
from sunnyhills.pipeline_functions import run_tls 

tic_id = 'TIC_466265409'

df = pd.read_csv('./routines/alpha_tls/data/two_min_lightcurves/'+tic_id+'.csv')

clean_time, clean_flux = (np.array(df[i]) for i in ['clean_time', 'clean_flux'])
clean_mask = np.isfinite(clean_time)
clean_time, clean_flux = (i[clean_mask] for i in [clean_time, clean_flux])

raw_time, raw_flux = (np.array(df[i]) for i in ['no_flare_raw_time','no_flare_raw_flux'])

trend_time, trend_flux = (np.array(df[i]) for i in ['trend_time','trend_flux'])
trend_mask = np.isfinite(trend_time)
trend_time, trend_flux = (i[trend_mask] for i in [trend_time, trend_flux])

split_axes = False
break_index = None 

run_tls(tic_id, clean_time, clean_flux, verbose=True)

diff = np.diff(clean_time)
if np.max(diff>50): 
    split_axes = True 
    break_index = np.where(np.diff(clean_time)>50)[0][0]+1 # so it's [0:break_index), [break_index:]
    break_index = int(break_index)  

## PLOTTING ### 

#plt.style.use('seaborn-darkgrid')
plt.rcParams['font.family']='serif'
fig = plt.figure(constrained_layout=True, figsize=(12,12))

gs = GridSpec(4, 4, figure=fig)

if split_axes: 
    ax1a = fig.add_subplot(gs[0, 0:2])
    ax1b = fig.add_subplot(gs[0, 2:])
    # hide the spines between ax and ax2
    ax1a.spines['right'].set_visible(False)
    ax1b.spines['left'].set_visible(False)
    ax1a.yaxis.tick_left()
    #ax.tick_params(labeltop='off') # don't put tick labels at the top
    ax1b.yaxis.tick_right()

    # detrend light curve 

    ax1_ylim = [0.95*np.min(clean_flux), 1.05*np.max(clean_flux)]

    ax1a.scatter(clean_time[0:break_index], clean_flux[0:break_index], s=1)
    ax1a.set(ylabel='Detrended Flux', ylim=ax1_ylim)

    ax1b.scatter(clean_time[break_index:], clean_flux[break_index:], s=1)
    ax1b.set(ylim=ax1_ylim)

    d = .015 # how big to make the diagonal lines in axes coordinates
    kwargs = dict(transform=ax1a.transAxes, color='k', clip_on=False)
    ax1a.plot((1-d,1+d),(-d,+d), **kwargs) # top-left diagonal
    ax1a.plot((1-d,1+d),(1-d,1+d), **kwargs) # bottom-left diagonal

    kwargs.update(transform=ax1b.transAxes) 
    ax1b.plot((-d,d),(-d,+d), **kwargs) 
    ax1b.plot((-d,d),(1-d,1+d), **kwargs)

    for ax in [ax1a, ax1b]: 
        ax.set(xlabel='Time (days)')

else: 
    ax1 = fig.add_subplot(gs[0, :]) # detrended light curve
    
    # detrend light curve 

    ax1.scatter(clean_time, clean_flux, s=1)
    ax1.set(ylabel='Detrended Flux')

    for ax in [ax1, ax2]: 
        ax.set(xlabel='Time (days)', xlim=lc_xlim)

ax2 = fig.add_subplot(gs[1, :-1]) # no flare with trend light curve
ax3 = fig.add_subplot(gs[2, :-1]) # phase folded transit 
ax4 = fig.add_subplot(gs[-1, 0]) # left right transits
ax5 = fig.add_subplot(gs[-1, -2]) # depth diffs 
ax6 = fig.add_subplot(gs[-1,-1]) # periodogram 
ax7 = fig.add_subplot(gs[1:-1, -1]) # notes 

# raw and trend light curve # 
    
ax2.scatter(raw_time, raw_flux, s=1)
ax2.plot(trend_time, trend_flux, lw=0.5, c='r')
ax2.set(ylabel='Flux')

lc_xlim = (min((min(clean_time), min(raw_time))), max((max(clean_time), max(raw_time))))


plt.savefig('./personal_epochs/thaddaeus/june/topical/misc/img.png')