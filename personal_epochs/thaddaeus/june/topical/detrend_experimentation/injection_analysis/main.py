import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
plt.style.use('https://gist.githubusercontent.com/thissop/44b6f15f8f65533e3908c2d2cdf1c362/raw/fab353d758a3f7b8ed11891e27ae4492a3c1b559/science.mplstyle')

def initial_plot(true:str='./routines/simulations/second_bulk_injected/injection_key.csv', 
                 predicted:str='./routines/simulations/second_bulk_injected/results.csv'): 
    pass 
    true_df = pd.read_csv(true)       
    predicted_df = pd.read_csv(predicted)

    true_df.to_csv(true, index=False)

    merged = true_df.merge(predicted_df, on='TIC_ID')

    mask = np.where(predicted_df['SDE']>9)[0]

    merged_sig = merged.iloc[mask]

    fig, axs = plt.subplots(1, 2, figsize=(12,5))

    sig_true, sig_tls = (np.array(i) for i in (merged_sig['INJECTED_PER'], merged_sig['TLS_PER']))

    recovered_sig = str(round(100*np.sum([1 for i in range(len(sig_true)) if 0.99<=sig_true[i]/sig_tls[i]<=1.1])/len(sig_true), 1))

    scat = axs[0].scatter(sig_true, sig_tls, c=merged_sig['SDE'], cmap='inferno', label=recovered_sig+' % Recovery Rate')
    plt.colorbar(scat, ax=axs[0])
    axs[0].axline((0,0), slope=1)
    axs[0].set(xlabel='True Period (d)', ylabel='TLS Period (d)', title='SDE > 9', xlim=(0.1, 17), ylim=(0.1, 17))
    axs[0].legend()

    injected_true, injected_per = (np.array(i) for i in (merged['INJECTED_PER'], merged['TLS_PER']))

    recovered_total = str(round(100*np.sum([1 for i in range(len(injected_true)) if 0.99<=injected_true[i]/injected_per[i]<=1.1])/len(injected_per), 1))

    scat = axs[1].scatter(merged['INJECTED_PER'], merged['TLS_PER'], c=merged['SDE'], cmap='inferno', label=recovered_total+' % Recovery Rate')
    plt.colorbar(scat, ax=axs[1])
    axs[1].axline((0,0), slope=1)
    axs[1].set(xlabel='True Period (d)', ylabel='TLS Period (d)', title='SDE > 0', xlim=(0.1, 17), ylim=(0.1, 17))
    axs[1].legend()

    plt.savefig('./personal_epochs/thaddaeus/june/topical/detrend_experimentation/injection_analysis/RJ=1_results.png')

initial_plot()

def high_sde_poor_predictions(true:str='./routines/simulations/second_bulk_injected/injection_key.csv', 
                              predicted:str='./routines/simulations/second_bulk_injected/results.csv'): 
    pass 
    true_df = pd.read_csv(true)       
    predicted_df = pd.read_csv(predicted)

    true_df.to_csv(true, index=False)

    merged = true_df.merge(predicted_df, on='TIC_ID')

    mask = np.where(predicted_df['SDE']>9)[0]

    mask = np.where(predicted_df['TLS_PER']<2)
    merged = merged.iloc[mask]

    merged.to_csv('./personal_epochs/thaddaeus/june/topical/detrend_experimentation/injection_analysis/RJ=1_to_look_into.csv', index=False)

high_sde_poor_predictions()

