'''
TIC_ID,GDR2_ID,GDR2_RA,GDR2_DEC,TLC,mas_parrallax,distance,
gaia_g_mag,absolute_mag,bp-rp,probability_age_less_than_50,age_myr,age_err_lower,age_err_upper,
cluster_name,first_ls_period,second_ls_period,third_ls_period,
first_ls_power,second_ls_power,third_ls_power,95_fap_level,raw_counts,no_flare_counts,
BiWeight-w=0.5_counts,LombScargle-n=2_counts
'''

# FUTURE FIX: MAKE IT SO THAT IT DOESN'T OVERWRITE PREVIOUS FILES/ENTRIES #
# FUTURE FIX: MAKE IT APPEND TO A CSV LINE-BY-LINE MID-ITERATION # 

def create_base_catalog(tic_ids, out_catalog_path:str='/ar1/PROJ/fjuhsd/shared/tessodyssey/data/temp_alpha_output.csv',
                        working_dir:str='/ar1/PROJ/fjuhsd/shared/tessodyssey/'):
    
    from tqdm import tqdm 
    import pandas as pd
    import numpy as np
    import os 
    from sunnyhills.pipeline_functions import query_tls_vizier, query_simbad, remove_flares, better_download, better_preprocess
    from astropy.timeseries import LombScargle
    
    if working_dir[-1]!='/':
        working_dir+='/'

    all_downloaded_sectors = [] 
    all_sector_starts = [] 
    all_sector_ends = []

    first_ls_periods = [] 
    second_ls_periods = []  
    third_ls_periods = []  
    first_ls_powers = [] 
    second_ls_powers = [] 
    third_ls_powers = [] 
    fap_95_levels = [] 

    all_raw_counts = []
    all_no_flare_counts = []
    all_wotan_counts = [] # 
    all_lomb_scargle_counts = [] 

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

    completed_tic_ids = []
    for tic_id in tqdm(tic_ids):  
        #try: 
        if tic_id in tic_ids:
            
            # 1. DOWNLOAD DATA 
            raw_path = raw_dir + tic_id + '.csv'

            raw_df, downloaded_sectors, sector_starts, sector_ends, last_dates = better_download(tic_id=tic_id, save_directory=raw_dir)
            if raw_df is not None and downloaded_sectors is not None: 
                completed_tic_ids.append(tic_id)
                all_downloaded_sectors.append(downloaded_sectors) 
                all_sector_starts.append(sector_starts)
                all_sector_ends.append(sector_ends)

                raw_time, raw_flux = (np.array(raw_df[i]) for i in ['raw_time', 'raw_flux'])
                no_flare_mask = np.isfinite(raw_df['no_flare_raw_time'])
                no_flare_raw_time, no_flare_raw_flux = (np.array(raw_df[i])[no_flare_mask] for i in ['no_flare_raw_time', 'no_flare_raw_flux'])

                all_raw_counts.append(len(raw_time))
                all_no_flare_counts.append(len(no_flare_raw_time))

                # 2. PREPROCESS DATA

                # 2.a Wotan Detrend

                _, wotan_counts = better_preprocess(tic_id=tic_id, raw_data=raw_path, last_dates=last_dates, 
                                                    save_directory=wotan_dir)

                all_wotan_counts.append(wotan_counts)

                # 2.b Lomb-Scargle on Data
                periodogram = LombScargle(no_flare_raw_time, no_flare_raw_flux, nterms=2)

                frequencies, powers = periodogram.autopower(minimum_frequency=1/15, maximum_frequency=1/0.1, method='fastchi2')

                sort_idx = np.argsort(powers)[::-1]
                powers, frequencies = (i[sort_idx] for i in (powers, frequencies))
                periods = 1/frequencies 
                best_frequency = frequencies[0]

                first_ls_periods.append(periods[0])
                second_ls_periods.append(periods[1])
                third_ls_periods.append(periods[2])
                first_ls_powers.append(powers[0])
                second_ls_powers.append(powers[1])
                third_ls_powers.append(powers[2])

                # 2.c Lomb-Scargle Detrend

                y_fit = periodogram.model(no_flare_raw_time, best_frequency)

                (ls_clean_time, ls_clean_flux), (_, _) = remove_flares(no_flare_raw_time, no_flare_raw_flux/y_fit)

                all_lomb_scargle_counts.append(len(ls_clean_time))

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

            # Bookkeeping and Detrend Information #
            all_downloaded_sectors.append(np.nan) 
            all_sector_starts.append(np.nan)
            all_sector_ends.append(np.nan)
            all_raw_counts.append(np.nan)
            all_no_flare_counts.append(np.nan)
            all_wotan_counts.append(np.nan)

            # Lomb Scargle Information # 
            first_ls_periods.append(np.nan)
            second_ls_periods.append(np.nan)
            third_ls_periods.append(np.nan)
            first_ls_powers.append(np.nan)
            second_ls_powers.append(np.nan)
            third_ls_powers.append(np.nan)
            fap_95_levels.append(np.nan)
            all_lomb_scargle_counts.append(np.nan)

            # Vizier Information # 
            ab_values.append(np.nan)
            mass_values.append(np.nan)
            mass_mins.append(np.nan)
            mass_maxs.append(np.nan)
            radii.append(np.nan)
            radius_mins.append(np.nan)
            radius_maxs.append(np.nan)

            # Simbad Information #
            OTYPES.append(np.nan)
            SP_TYPES.append(np.nan)
            MAIN_OTYPES.append(np.nan)

    # Bookkeeping and Detrend Information #
    out_df = pd.DataFrame()
    out_df['TIC_ID'] = tic_ids
    out_df['downloaded_sector_nums'] = all_downloaded_sectors 
    out_df['sector_starts'] = all_sector_starts
    out_df['sector_ends'] = all_sector_ends
    out_df['raw_counts'] = all_raw_counts
    out_df['no_flare_counts'] = all_no_flare_counts
    out_df['wotan_counts'] = all_wotan_counts

    # Lomb Scargle Information # 
    out_df['first_ls_periods'] = first_ls_periods
    out_df['second_ls_periods'] = second_ls_periods
    out_df['third_ls_periods'] = third_ls_periods
    out_df['first_ls_powers'] = first_ls_powers
    out_df['second_ls_powers'] = second_ls_powers
    out_df['third_ls_powers'] = third_ls_powers
    out_df['lomb_scargle_counts'] = all_lomb_scargle_counts

    # Vizier Information # 
    out_df['ab_values'] = ab_values 
    out_df['mass_value'] = mass_values
    out_df['mass_min'] = mass_mins
    out_df['mass_max'] = mass_maxs
    out_df['radius'] = radii
    out_df['radius_min'] = radius_mins
    out_df['radius_max'] = radius_maxs

    # Simbad Information #
    out_df['OTYPES'] = OTYPES
    out_df['SP_TYPE'] = SP_TYPES
    out_df['MAIN_OTYPE'] = MAIN_OTYPES

    # SAVE UPDATED CATALOG # 
    out_df.to_csv(out_catalog_path, index=False)

create_base_catalog(['TIC_1232360', 'TIC_43149283'])