# CURRENT USE #

def tls_even_odd(tls_results): 
    
    '''
    Description 
    -----------
    returns true if mean odd and mean even transit depths are inconsistent at three sigma level 

    '''

    import numpy as np

    eb_flag = True 

    odd = tls_results.depth_mean_even 

    odd = [odd[0]-3*odd[1], odd[0]+3*odd[1]] # note 99% confidence because mean +/- 3 sigma

    even = tls_results.depth_mean_odd
    even = [even[0]-3*even[1], even[0]+3*even[1]] # note 99% confidence because mean +/- 3 sigma

    idx = np.argsort([even[0], odd[0]])

    temp = np.array([even, odd])[idx]

    if temp[0][1]>=temp[1][0]: 
        eb_flag = False

    even = [round(i, 3) for i in even]
    odd = [round(i, 3) for i in odd]

    return {'Even-Odd Flag':eb_flag, 'even':even, 'odd':odd} # save both parts of even and odd as columns in TLS results, e.g. even_low and even_high, odd_low and odd_high 

def check_lombscargle(tic_id, tls_results, alpha_catalog): 

    import pandas as pd 
    import numpy as np

    if isinstance(alpha_catalog, pd.DataFrame): 
        df = alpha_catalog  
    else: 
        df = pd.read_csv(alpha_catalog)

    index = np.where(df['TIC_ID']==tic_id)[0]
    top_ls_period = np.array(df['top_ls_period'])[index]
    if len(top_ls_period)>0: 
        top_ls_period = top_ls_period[0]   

    first_alias_ls_flag = False

    ls_flag = False
    if 0.99<tls_results.period/top_ls_period<1.01: 
        ls_flag = True 
    
    elif 0.99<(2*tls_results.period)/top_ls_period<1.01: 
        first_alias_ls_flag = True 

    elif 0.99<(0.5*tls_results.period)/top_ls_period<1.01: 
        sub_alias_ls_flag = True 

    return {'ls_top_period_test':ls_flag, 'ls_period':top_ls_period, 'n=2_harmonic_test':first_alias_ls_flag, 'n=0.5_harmonic_test':sub_alias_ls_flag}

# NOT IN CURRENT USE # 

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

    return {'lc_sigma':round(lc_sigma, 3), 'transit_sigma':round(transit_sigma, 3)} 
