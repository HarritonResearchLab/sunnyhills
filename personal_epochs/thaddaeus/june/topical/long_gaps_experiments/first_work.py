def first_work(tic_id:str='TIC_14087575'): 
    import numpy as np
    import pandas as pd
    import pandas as pd 

    from sunnyhills.pipeline_functions import run_tls

    data = pd.read_csv('./data/current/processed/two_min_lightcurves/'+tic_id+'.csv')
    mask = np.isfinite(data['clean_time'])

    clean_time, clean_flux = (np.array(data[col])[mask] for col in ['clean_time', 'clean_flux'])

    max_per = 15
    num_transits = int((np.max(clean_time)-np.min(clean_time))/max_per)

    print(num_transits)

    tls_results, tls_model = run_tls(tic_id=tic_id, time=clean_time, flux=clean_flux, show_progress_bar=True, min_transits=num_transits)

first_work()