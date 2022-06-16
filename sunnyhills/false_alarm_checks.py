from sunnyhills.pipeline_functions import download


def tls_even_odd(tls_results): 
    
    '''
    Description 
    -----------
    returns true if mean odd and mean even transit depths are inconsistent at one sigma level 

    '''

    import numpy as np

    eb_flag = True 

    odd = tls_results.depth_mean_even 

    odd = [odd[0]-odd[1], odd[0]+odd[1]]

    even = tls_results.depth_mean_odd
    even = [even[0]-even[1], even[0]+even[1]]

    temp = np.sort([even, odd])

    if temp[0][1]>temp[1][0]: 
        eb_flag = False

    return {'Even-Odd Flag':eb_flag}

def transit_outliers_fap_test(tls_results): 
    '''
    Description
    -----------
    See section 3.3 from edi-vetter paper (jon zink et al. 2020)
    '''
    
    import numpy as np

    folded_phase = tls_results.folded_phase 
    folded_y = tls_results.folded_y
    
    folded_model_flux = tls_results.model_folded_model
    
    half_dur = tls_results.duration 
    transit_mask = np.logical_and(folded_phase<=0.5+half_dur, folded_phase>=0.5-half_dur)
    
    subtracted = folded_y-folded_model_flux 

    lc_sigma = np.std(subtracted[~transit_mask])
    transit_sigma = np.std(subtracted[transit_mask]) 

    return {'lc_sigma':lc_sigma, 'transit_sigma':transit_sigma} 

def check_lombscargle(tic_id, tls_results, download_log): 

    import pandas as pd 
    import numpy as np

    if isinstance(download_log, pd.DataFrame): 
        df = download_log  
    else: 
        df = pd.read_csv(download_log)

    index = np.where(df['TIC_ID']==tic_id)[0]
    top_ls_period = df['top_ls_period'][index]  
    print(top_ls_period)

    ls_flag = False
    if 0.99<tls_results.period/top_ls_period<1.01: 
        ls_flag = True 

    return {'ls_top_period_test':ls_flag}
