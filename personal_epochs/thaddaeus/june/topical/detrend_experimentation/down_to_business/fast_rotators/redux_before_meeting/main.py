def test(fast_df, data_dir, plot_dir): 
    import pandas as pd 
    import numpy as np
    from sunnyhills.injection import better_inject
    import matplotlib.pyplot as plt 
    from astropy.timeseries import LombScargle
    from sunnyhills.pipeline_functions import remove_flares, run_tls
    import wotan 
    from tqdm import tqdm 

    fast_df = pd.read_csv(fast_df)
    tic_ids = np.array(fast_df['TIC_ID'])
    periods = np.array(fast_df['top_ls_period'])  

    if data_dir[-1]!='/':
        data_dir += '/'

    if plot_dir[-1]=='/':
        data_dir += '/'       

    for index, tic_id in tqdm(enumerate(tic_ids)): 
        try: 
            data = pd.read_csv(data_dir+tic_id+'.csv')

            mask = np.isfinite(data['no_flare_raw_flux'])
            no_flare_time, no_flare_flux = [np.array(data[i])[mask] for i in ['no_flare_raw_time', 'no_flare_raw_flux']]  

            time, flux, (per, _, t0) = better_inject(tic_id=tic_id, time=no_flare_time, flux=no_flare_flux, rp=0.10045*0.5)

            list_of_trend_fluxes = []
            list_of_trend_times = []

            list_of_detrended_times = []
            list_of_detrended_fluxes = []

            # DETRENDING # 

            for i in range(1,4): 
                periodogram = LombScargle(time, flux, nterms=i)
                frequencies, powers = periodogram.autopower(minimum_frequency=1/15, maximum_frequency=1/0.1, method='fastchi2')
                best_frequency = frequencies[np.argmax(powers)]
                y_fit = periodogram.model(no_flare_time, best_frequency)

                (clean_time, clean_flux, _), (_, _, _) = remove_flares(time, flux/y_fit, y_fit)

                list_of_trend_times.append(time)
                list_of_trend_fluxes.append(y_fit)
            
                list_of_detrended_times.append(clean_time)
                list_of_detrended_fluxes.append(clean_flux)

            for i in (0.25, 0.5): 
                clean_flux, trend_flux = wotan.flatten(time, flux, return_trend=True,method='biweight',break_tolerance=1,window_length=i,cval=5)

                (clean_time, clean_flux, _), (_, _, _) = remove_flares(time, clean_flux, trend_flux)
                
                list_of_trend_times.append(time)
                list_of_trend_fluxes.append(trend_flux)
            
                list_of_detrended_times.append(clean_time)
                list_of_detrended_fluxes.append(clean_flux)

            # PLOT # 
            transits = np.arange(t0, np.max(time), per)

            fig, axs = plt.subplots(5,2, figsize=(16,15))

            titles = ['1-term LS', '2-term LS', '3-term LS', '0.25 d BW', '0.5 d BW']

            for i in range(0, 5): 
                tls_results, tls_model = run_tls(tic_id, list_of_detrended_times[i], list_of_detrended_fluxes[i])

                local_title = titles[i]+'\nTLS PERIOD = '+str(round(tls_results.period, 3))
                
                ax = axs[i, 0]

                ax.scatter(time, flux, s=1, c='grey')
                ax.plot(list_of_trend_times[i], list_of_trend_fluxes[i], color='cornflowerblue', lw=2)

                ax.set(xlabel='Time (d)', ylabel='Flux', title=local_title)
                for transit in transits: 
                    ax.axvline(x=transit, lw=2, color='orange', alpha=0.6)

                ax = axs[i, 1]
                ax.scatter(list_of_detrended_times[i], list_of_detrended_fluxes[i], color='grey', s=1)

                ax.plot(tls_results.model_lightcurve_time, 
                    tls_results.model_lightcurve_model, color='cornflowerblue', lw=1, zorder=1)

                ax.set(xlabel='Time (d)', ylabel='Flux', title=local_title)
                for transit in transits: 
                    ax.axvline(x=transit, lw=2, color='orange', alpha=0.6)
            
            fig.suptitle('STELLAR PERIOD ='+str(round(periods[index], 2))+'\nINJECTED PERIOD = '+str(round(per, 2)), size='xx-large')
            plt.tight_layout()
            plt.savefig(plot_dir+tic_id+'.png', dpi=150)
            plt.clf()
            plt.close()

        except: 
            continue

fast_df = './personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/fast_rotators/fast_df.csv'
data_dir = './data/current/processed/two_min_lightcurves/'
plot_dir = './personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/fast_rotators/redux_before_meeting/plots/'       
test(fast_df, data_dir, plot_dir)