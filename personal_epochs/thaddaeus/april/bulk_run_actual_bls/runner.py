import sunnyhills
from sunnyhills.pipeline_functions import run_bls, download, preprocess
import pandas as pd
import numpy as np 
from tqdm import tqdm 
import matplotlib.pyplot as plt 
import pickle 
import os 

key_df = pd.read_csv('./data/current/current_key.csv')

mask = np.array(key_df['has_two_min'])
tic_ids = np.array(key_df['TIC_ID'])[mask]

lc_out_dir = './data/current/processed/two_min_lightcurves'
plot_dir = './personal_epochs/thaddaeus/april/bulk_run_actual_bls/first_results/plots'
for tic_id in tqdm(tic_ids): 
    tic_id = 'TIC_'+str(tic_id)

    run_analysis = False 
    already_downloaded = False
    for file in os.listdir(lc_out_dir): 
        if tic_id in file: 
            already_downloaded = True 
            run_analysis = True 
            break 
    
    lc_dict = {}
    if already_downloaded: 
        file_path = lc_out_dir+'/'+tic_id+'_lc.pickle' 
        with open(file_path, 'rb') as handle:
            lc_dict = pickle.load(handle)

        stitched_lc = lc_dict['stitched_lc']
        stitched_trend = lc_dict['stitched_trend']
        stitched_raw = lc_dict['stitched_raw']
        

    else: 

        raw_lc, data_found = download(ticstr=tic_id)
                    
        if data_found: 
            run_analysis = True 
            stitched_lc, stitched_trend, stitched_raw = preprocess(raw_list=raw_lc, ticstr=tic_id, outdir=lc_out_dir)

    if run_analysis: 
        results, bls_model, in_transit, stats = run_bls(stitched_lc)

        time = stitched_lc['time']
        flux = stitched_lc['flux']
        fig = plt.figure(constrained_layout=True, figsize=(7, 7))
        mosaic = """
        ab
        cd"""

        ax_dict = fig.subplot_mosaic(mosaic)
        
        period = results.period[np.argmax(results.power)]

        ax = ax_dict['a']
        ax.axvline(period, alpha=0.4, lw=3)
        for n in range(2, 10):
            ax.axvline(n*period, alpha=0.4, lw=1, linestyle="dashed")
            ax.axvline(period / n, alpha=0.4, lw=1, linestyle="dashed")
        ax.plot(results.period, results.power, color='black', lw=0.5)
        ax.set(xlim=(0, max(results.period)), ylabel='SDE BLS', xlabel='Period (days)')

        index = np.argmax(results.power)
        period = results.period[index]
        t0 = results.transit_time[index]
        duration = results.duration[index]

        ax = ax_dict['b']
        x = (time - t0 + 0.5*period) % period - 0.5*period
        m = np.abs(x) < 0.5
        ax.scatter(
            x[m],
            flux[m],
            color='blue',
            s=10,
            alpha=0.5,
            zorder=2)

        x = np.linspace(-0.13, 0.13, 1000)
        f = bls_model.model(x + t0, period, duration, t0)
        ax.plot(x, f, color='red')
        #ax.set_xlim(-0.13, 0.13)
        #plt.ylim(0.9985, 1.00025)
        ax.set_xlabel("Time from mid-transit (days)")
        ax.set_ylabel("Flux")

        ax = ax_dict['c']

        ax.scatter(stitched_raw['time'], stitched_raw['flux'], c='black', s=3, label='Raw')

        ax.scatter(stitched_trend['time'], stitched_trend['flux'], c='red', s=3, label='Trend')

        ax.legend()

        ax.set(xlabel='Time (mjd)', ylabel='Flux')

        ax = ax_dict['d']

        ax.scatter(stitched_lc['time'], stitched_lc['flux'], c='black', s=3, label='Detrended LC')
        ax.set(xlabel='Time (MJD)', ylabel='Flux')

        ax.legend()

        plt.show()
        #plt.savefig(plot_dir+'/'+tic_id+'.png', dpi=150)
        plt.clf()
        plt.close()

    
        
    if not already_found: 
        out_dict = {'stitched_lc':stitched_lc, 'stitched_trend':stitched_trend, 'stitched_raw':stitched_raw}
        file_path = lc_out_dir+'/'+tic_id+'_lc.pickle' 
        with open(file_path, 'wb') as handle:
            pickle.dump(out_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)