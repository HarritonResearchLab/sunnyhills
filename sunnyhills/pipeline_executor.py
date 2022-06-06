def alpha_routine(key:str, make_plots:bool=True): 
    '''
    Parameters

    Returns

    Notes
        - Runs BLS, if no signal found for BLS, runs TLS. Runs some FAP tests, saves results.
    '''
    
    from pipeline_functions import run_bls
    import numpy as np
    import pandas as pd 
    from tqdm import tqdm 

    key_df = pd.read_csv(key)
    tic_ids = np.array(key_df['tic_id'])
    if type(tic_ids[0])==int: 
        tic_ids = np.array(['TIC_'+str(i) for i in tic_ids]) 
    elif type(tic_ids[0])==str: 
        if '_' not in tic_id: 
            tic_ids = np.array([i.replace('TIC', 'TIC_') for i in tic_ids])

    for tic_id in tic_ids: 
        data = pd.read_csv('')
        clean_time = np.array(data['cleaned_time'])
        detrend_flux = np.array(data['detrended_flux'])

        # run bls
        bls_tuple = run_bls(clean_time, detrend_flux)

        bls_best_params, bls_results = (bls_tuple[0], 
                                        bls_tuple[1])
                                        
        bls_model, bls_in_transit, bls_stats = (bls_tuple[2], 
                                                bls_tuple[3],   
                                                bls_tuple[4])

        


