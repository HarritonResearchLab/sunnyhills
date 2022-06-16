def first_look(results:str='personal_epochs/thaddaeus/june/weekly/second_week/tls_routine/output_log.txt'): 
    import pandas as pd
    import numpy as np
    

    df = pd.read_csv(results).sort_values(by='SDE', ascending=False)
    print(df)

#first_look()

def make_best_plots(data_dir:str, plot_dir:str='/ar1/PROJ/fjuhsd/shared/github/sunnyhills/personal_epochs/thaddaeus/june/weekly/second_week/tls_routine/plots'):
    r'''
    Parameters

    Returns

    Notes
        - This is the first TLS pipeline executor routine; only runs TLS
    '''
    
    from sunnyhills.pipeline_functions import run_tls 
    from sunnyhills.false_alarm_checks import tls_even_odd
    from sunnyhills.plotting import tls_validation_mosaic
    from transitleastsquares import transit_mask
    import numpy as np
    import pandas as pd 
    from tqdm import tqdm 
    import pickle 
    import os

    results = 'personal_epochs/thaddaeus/june/weekly/second_week/tls_routine/output_log.txt'
    df = pd.read_csv(results).sort_values(by='SDE', ascending=False)

    ids = np.array(df['TIC_ID'][0:10])
    np.random.shuffle(ids)


    #for tic_id in tqdm(ids): 
    tic_id = 'TIC_190885165'
    print(tic_id)
    data_path = data_dir+tic_id+'.csv'
    if os.path.exists(data_path):
        data = pd.read_csv(data_path) 
        clean_time = np.array(data['clean_time'])
        clean_flux = np.array(data['clean_flux'])

        cache_dir = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/personal_epochs/thaddaeus/june/weekly/second_week/tls_routine/caching/'
        
        if os.path.exists(cache_dir+tic_id+'_tls-model.pickle'): 
            pickle_results = cache_dir+tic_id+'_tls-results.pickle'
            print(pickle_results)
            with open(pickle_results, 'rb') as file: 
                tls_results = pickle.load(file)

            pickle_model = cache_dir+tic_id+'_tls-model.pickle'
            with open(pickle_model, 'rb') as file: 
                tls_model = pickle.load(file)
        
        else: 

            tls_results, tls_model, transit_mask = run_tls(tic_id=tic_id, time=clean_time, flux=clean_flux, cache_dir=cache_dir) 

        eb_flag = tls_even_odd(tls_results=tls_results)

        flags = {}

        #tls_validation_mosaic(tic_id=tic_id, data=data_path, tls_results=tls_results, tls_model=tls_model, plot_dir=plot_dir)
        


data_dir = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/data/current/processed/two_min_lightcurves/'
make_best_plots(data_dir)