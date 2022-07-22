# CURRENT USE #

def tls_even_odd(tls_results): 
    
    '''
    Description 
    -----------
    returns true if mean odd and mean even transit depths are inconsistent at three sigma level 

    '''

    import numpy as np

    sig_multiple = 1

    eb_flag = True 

    odd = tls_results.depth_mean_even 

    odd_one = [odd[0]-odd[1], odd[0]+odd[1]] 
    odd_three = [odd[0]-3*odd[1], odd[0]+3*odd[1]] # note 99% confidence because mean +/- 3 sigma

    even = tls_results.depth_mean_odd
    even_one = [even[0]-even[1], even[0]+even[1]] 
    even_three = [even[0]-3*even[1], even[0]+3*even[1]] # note 99% confidence because mean +/- 3 sigma

    idx = np.argsort([even_one[0], odd_one[0]])
    temp = np.array([even_one, odd_one])[idx]
    if temp[0][1]>=temp[1][0]:
        pass
    else: 
        pass 
        sig_multiple = (odd[0]-even[0])/(even[1]-odd[1])

    idx = np.argsort([even_three[0], odd_three[0]])
    temp = np.array([even_three, odd_three])[idx]

    if temp[0][1]>=temp[1][0]: 
        eb_flag = False

    even_three = [round(i, 3) for i in even_three]
    odd_three = [round(i, 3) for i in odd_three]

    return {'Even-Odd Flag':eb_flag, 'even':even_three, 'odd':odd_three} # save both parts of even and odd as columns in TLS results, e.g. even_low and even_high, odd_low and odd_high 

def check_lombscargle(tic_id, tls_results, catalog=None, on_the_fly:bool=False, raw_data_dir:str='/ar1/PROJ/fjuhsd/shared/tessodyssey/data/raw/'): 
    r'''
    note that this checks whether or not harmonics of tls period are within 1% of a single highest ls period (doesn't incorporate LS period harmonics or secondary/tertiary LS periods)
    '''

    import pandas as pd 
    import numpy as np
    from sunnyhills.misc import lombscargle

    if catalog is not None: 
        if isinstance(catalog, pd.DataFrame): 
            df = catalog  
        else: 
            df = pd.read_csv(catalog)

        index = np.where(df['TIC_ID']==tic_id)[0]
        top_ls_period = np.array(df['top_ls_period'])[index]
        if len(top_ls_period)>0: 
            top_ls_period = top_ls_period[0]   

    elif on_the_fly is True and raw_data_dir is not None: 

        data = pd.read_csv(raw_data_dir+tic_id+'.csv').dropna()
        time = np.array(data['no_flare_raw_time'])
        flux = np.array(data['no_flare_raw_flux'])

        top_ls_period = lombscargle(time, flux, calc_fap=False)[1][0]

    else: 
        raise Exception('Either catalog needs to be defined, or on_the_fly needs to be true')

    first_alias_ls_flag = False
    sub_alias_ls_flag = False

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
