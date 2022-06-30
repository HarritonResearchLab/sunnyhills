def first_play(tic_id:str='TIC_312185848'): 
    import pandas as pd
    import numpy as np
    import pickle 
    from sunnyhills.false_alarm_checks import tls_even_odd 

    pickle_results = './routines/simulations/first_bulk_injected/cache_dir/'+tic_id+'_tls-results.pickle'
    with open(pickle_results, 'rb') as file: 
        tls_results = pickle.load(file)


    print(tls_results.depth_mean_even)

    print(tls_even_odd(tls_results))

first_play()