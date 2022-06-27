import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
import seaborn as sns

def make_kde_plot(): 
    
    fast_df = pd.read_csv('./personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/bad_detrend_catalog.csv') 
    fast_df = fast_df.drop(columns=['complex fast rotator', 'long gap', 'weird', 'comment on'])
    full_df = pd.read_csv('./routines/alpha_tls/data/download_log.csv')

    fast_df = fast_df.merge(full_df, left_on='clean fast rotator', right_on='TIC_ID')

    fast_df = fast_df['top_ls_power']
    full_df = full_df['top_ls_power']

    fig, ax = plt.subplots()
    sns.histplot(data=full_df, axes=ax, label='Full', kde=True, stat='probability')
    sns.histplot(data=fast_df, axes=ax, label='Fast Subset', kde=True, stat='probability', color='indianred')
    ax.legend()
    plt.savefig('./personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/fast_rotators/periods_dist.png', dpi=200)

#make_kde_plot()

def save_periods(): 
    fast_df = pd.read_csv('./personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/bad_detrend_catalog.csv') 
    fast_df = fast_df.drop(columns=['complex fast rotator', 'long gap', 'weird', 'comment on'])
    full_df = pd.read_csv('./routines/alpha_tls/data/download_log.csv')

    fast_df = fast_df.merge(full_df, left_on='clean fast rotator', right_on='TIC_ID')

    df = pd.DataFrame()
    for i in ['TIC_ID', 'top_ls_period', 'top_ls_power', 'fap_95']: 
        df[i] = fast_df[i]

    print(df)

    df.to_csv('./personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/fast_rotators/fast_df.csv', index=False)

#save_periods()

def ls_detrend(): 
    from sunnyhills.pipeline_functions import remove_flares
    from sunnyhills.plotting import plot_detrend_validation
    from astropy.timeseries import LombScargle
    from sunnyhills.misc import inject 
    from tqdm import tqdm 

    df = pd.read_csv('./personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/fast_rotators/fast_df.csv')

    tic_ids = df['TIC_ID']

    for tic_id in tqdm(tic_ids): 
        initial_data = pd.read_csv('./routines/alpha_tls/data/two_min_lightcurves/'+tic_id+'.csv')
        
        raw_time, raw_flux = (np.array(initial_data[i]) for i in ['raw_time', 'raw_flux'])

        mask = np.isfinite(initial_data['no_flare_raw_flux'])

        no_flare_time, no_flare_flux = (np.array(initial_data[i]) for i in ['no_flare_raw_time', 'no_flare_raw_flux'])

        no_flare_time, no_flare_flux = (i[mask] for i in [no_flare_time, no_flare_flux])

        periodogram = LombScargle(no_flare_time, no_flare_flux, nterms=2)

        frequencies, powers = periodogram.autopower(minimum_frequency=1/15, maximum_frequency=1/0.1, method='fastchi2')

        periods = 1/frequencies

        best_frequency = frequencies[np.argmax(powers)]
        y_fit = periodogram.model(no_flare_time, best_frequency)

        (clean_time, clean_flux, _), (_, _, _) = remove_flares(no_flare_time, no_flare_flux/y_fit, y_fit)

        cols = [clean_time, clean_flux, no_flare_time, y_fit, no_flare_time, no_flare_flux, raw_time, raw_flux]
        cols = [pd.Series(i) for i in cols]

        col_names = ['clean_time', 'clean_flux', 'trend_time', 'trend_flux', 'no_flare_raw_time', 'no_flare_raw_flux', 'raw_time', 'raw_flux']
    
        dictionary = {}
        for i in range(len(cols)):
            dictionary.update({col_names[i]:cols[i]})

        out_df = pd.DataFrame(dictionary)

        data_dir = './personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/fast_rotators/data/'
        outfile = data_dir+tic_id+'.csv'
        out_df.to_csv(outfile, index=False)

        plot_dir = './personal_epochs/thaddaeus/june/topical/detrend_experimentation/down_to_business/fast_rotators/plots/detrend_testing/'
        plot_detrend_validation(tic_id=tic_id, data_dir=data_dir, plot_dir=plot_dir, plot_type='png')

ls_detrend()