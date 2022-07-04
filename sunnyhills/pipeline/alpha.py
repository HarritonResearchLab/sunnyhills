'''
TIC_ID,GDR2_ID,GDR2_RA,GDR2_DEC,TLC,mas_parrallax,distance,
gaia_g_mag,absolute_mag,bp-rp,probability_age_less_than_50,age_myr,age_err_lower,age_err_upper,
cluster_name,first_ls_period,second_ls_period,third_ls_period,
first_ls_power,second_ls_power,third_ls_power,95_fap_level,raw_counts,no_flare_counts,
BiWeight-w=0.5_counts,LombScargle-n=2_counts
'''

def create_local_catalog(tic_ids, base_catalog, working_dir:str='/ar1/PROJ/fjuhsd/shared/tessodyssey/'):
    
    from tqdm import tqdm 
    import pandas as pd
    import numpy as np
    import os 
    from sunnyhills.pipeline_functions import query_tls_vizier, query_simbad, remove_flares
    from astropy.timeseries import LombScargle
    
    if working_dir[-1]!='/':
        working_dir+='/'

    base_df = pd.read_csv(base_catalog)
    base_ids = list(base_df['TIC_ID'])

    downloaded_sectors = [] #
    sector_starts = [] #
    sector_ends = [] #

    first_ls_periods = [] 
    second_ls_periods = []  
    third_ls_periods = []  
    first_ls_powers = [] 
    second_ls_powers = [] 
    third_ls_powers = [] 
    fap_95_levels = [] 

    raw_counts = []
    no_flare_counts = []
    wotan_counts = [] # 
    lomb_scargle_counts = [] 

    ab_values = [] 
    mass_values = [] 
    mass_mins = [] 
    mass_maxs = []
    radii = []
    radius_mins = []
    radius_maxs = []

    OTYPES = [] 
    SP_TYPES = []
    MAIN_OTYPES = []

    data_dir = working_dir+'data/'
    raw_dir = data_dir + 'raw/'
    wotan_dir = data_dir+'processed/BiWeight-w=0.5/'
    LombScargle_dir = data_dir+'processed/LombScargle-n=2/'

    processed_data_col_names = ['clean_time', 'clean_flux', 'trend_time', 'trend_flux', 'no_flare_raw_time', 'no_flare_raw_flux', 'raw_time', 'raw_flux']

    for tic_id in tqdm(base_ids):  
        #try: 
        if tic_id in tic_ids:
            
            # 1. DOWNLOAD DATA 
            raw_path = raw_dir + tic_id + '.csv'
        
            raw_time = []
            raw_flux = []
            no_flare_raw_time = []
            no_flare_raw_flux = []

            raw_counts.append(len(raw_time))
            no_flare_counts.append(len(no_flare_raw_time))

            # 2. PREPROCESS DATA

            # 2.a Wotan Detrend
            wotan_path = wotan_dir+tic_id+'.csv'

            # 2.b Lomb-Scargle on Data
            periodogram = LombScargle(no_flare_raw_time, no_flare_raw_flux, nterms=2)

            frequencies, powers = periodogram.autopower(minimum_frequency=1/15, maximum_frequency=1/0.1, method='fastchi2')

            sort_idx = np.argsort(powers)[::-1]
            powers, frequencies = (i[sort_idx] for i in (powers, frequencies))
            periods = 1/frequencies 
            best_frequency = frequencies[0]
            fap_95_level = periodogram.false_alarm_level([0.05])

            first_ls_periods.append(periods[0])
            second_ls_periods.append(periods[1])
            third_ls_periods.append(periods[2])
            first_ls_powers.append(powers[0])
            second_ls_powers.append(powers[1])
            third_ls_powers.append(powers[2])
            fap_95_levels.append(fap_95_level)

            # 2.c Lomb-Scargle Detrend

            y_fit = periodogram.model(no_flare_raw_time, best_frequency)

            (ls_clean_time, ls_clean_flux), (_, _) = remove_flares(no_flare_raw_time, no_flare_raw_flux/y_fit, y_fit)

            lomb_scargle_counts.append(len(ls_clean_time))

            ls_data_cols = [ls_clean_time, ls_clean_flux, no_flare_raw_time, y_fit, no_flare_raw_time, no_flare_raw_flux, raw_time, raw_flux]
            ls_data_cols = [pd.Series(i) for i in ls_data_cols]
        
            ls_dictionary = {}
            for i in range(len(processed_data_col_names)):
                ls_dictionary.update({processed_data_col_names[i]:ls_data_cols[i]})

            ls_out_df = pd.DataFrame(ls_dictionary)

            ls_data_path = LombScargle_dir+tic_id+'.csv'
            ls_out_df.to_csv(ls_data_path, index=False)

            # 3. ADD CATALOG INFO 

            # Vizier Information # 
            ab, mass, mass_min, mass_max, radius, radius_min, radius_max = query_tls_vizier(tic_id)
            ab_values.append('"("'+str(ab[0])+','+str(ab[1])+')"')
            mass_values.append(mass)
            mass_mins.append(mass_min)
            mass_maxs.append(mass_max)
            radii.append(radius)
            radius_mins.append(radius_min)
            radius_maxs.append(radius_max) 

            # Simbad Information #
            _, simbad_results = query_simbad(tic_id)

            OTYPES.append(simbad_results[0])
            SP_TYPES.append(simbad_results[1])
            MAIN_OTYPES.append(simbad_results[2])

        else:  

            # Lomb-Scargle Information # 
            first_ls_periods.append(np.nan)
            second_ls_periods.append(np.nan)
            third_ls_periods.append(np.nan)
            first_ls_powers.append(np.nan)
            second_ls_powers.append(np.nan)
            third_ls_powers.append(np.nan)
            fap_95_levels.append(np.nan)

            # Detrend Information # 
            raw_counts.append(np.nan)
            no_flare_counts.append(np.nan)
            wotan_counts.append(np.nan)
            lomb_scargle_counts.append(np.nan)

            # Vizier Information #

            ab_values.append('')
            mass_values.append(np.nan)
            mass_mins.append(np.nan)
            mass_maxs.append(np.nan)
            radii.append(np.nan)
            radius_mins.append(np.nan)
            radius_maxs.append(np.nan) 

            # Simbad Information # 
            OTYPES.append('')
            SP_TYPES.append('')
            MAIN_OTYPES.append('')


