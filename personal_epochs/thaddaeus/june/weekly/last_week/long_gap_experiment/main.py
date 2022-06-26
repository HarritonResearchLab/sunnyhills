def test(tic_id:str='TIC_14087575'): 
    from sunnyhills.pipeline_functions import run_tls
    from sunnyhills.misc import inject 
    import numpy as np 
    import pandas as pd
    import matplotlib.pyplot as plt 

    data = pd.read_csv('./routines/alpha_tls/data/two_min_lightcurves/'+tic_id+'.csv')
    mask = np.isfinite(data['clean_flux'])
    time = np.array(data['clean_time'])[mask]
    flux = np.array(data['clean_flux'])[mask]

    time, flux, tls_params = inject(time, flux, 4.3)

    break_indices = np.where(np.diff(time)>50)[0]
    first_break = break_indices[0]
    last_break = break_indices[-1]+1

    left_time = time[0:first_break]
    left_flux = flux[0:first_break]

    right_time = time[last_break:]
    right_flux = flux[last_break:]

    fig, axs = plt.subplots(2,1)

    axs[0].scatter(left_time, left_flux, s=1)

    axs[1].scatter(right_time, right_flux, s=1)

    plt.savefig('./personal_epochs/thaddaeus/june/weekly/last_week/long_gap_experiment/two_gaps.png', dpi=150)
    plt.clf()
    plt.close()

test()


