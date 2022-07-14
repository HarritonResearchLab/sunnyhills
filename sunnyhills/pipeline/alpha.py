'''
TIC_ID,GDR2_ID,GDR2_RA,GDR2_DEC,TLC,mas_parrallax,distance,
gaia_g_mag,absolute_mag,bp-rp,probability_age_less_than_50,age_myr,age_err_lower,age_err_upper,
cluster_name,first_ls_period,second_ls_period,third_ls_period,
first_ls_power,second_ls_power,third_ls_power,95_fap_level,raw_counts,no_flare_counts,
BiWeight-w=0.5_counts,LombScargle-n=2_counts
'''

# FUTURE FIX: MAKE IT SO THAT IT DOESN'T OVERWRITE PREVIOUS FILES/ENTRIES #
# FUTURE FIX: MAKE IT APPEND TO A CSV LINE-BY-LINE MID-ITERATION # 

from sunnyhills.pipeline_functions import download


def create_base_catalog(tic_ids,working_dir:str='/ar1/PROJ/fjuhsd/shared/tessodyssey/', piped_log:str='/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/log.txt'):
    
    r'''
    
    this only works if you run it the Thaddaeus/Bouma way as a subprocess

    i.e. python -u alpha.py &>> log.txt &


    the log_file needs to be the file that the routine is logged into 
    
    '''

    from tqdm import tqdm 
    import pandas as pd
    import numpy as np
    from sunnyhills.pipeline_functions import query_tls_vizier, query_simbad, remove_flares, better_download, better_preprocess
    from astropy.timeseries import LombScargle
    import warnings 
    import os 

    #warnings.filterwarnings("ignore")

    if working_dir[-1]!='/':
        working_dir+='/'

    data_dir = working_dir+'data/'
    raw_dir = data_dir + 'raw/'
    wotan_dir = data_dir+'processed/BiWeight-w=0.5/'
    LombScargle_dir = data_dir+'processed/LombScargle-n=2/'

    downloaded_log = data_dir + 'temp_alpha_output.csv'

    # reset tic ids, log results from last run # 
    if not os.path.exists(downloaded_log):
        log_col_names = ['TIC_ID','all_downloaded_sector','all_sector_start','all_sector_end','all_raw_count',
                     'all_no_flare_count','all_wotan_count','first_ls_period','second_ls_period','third_ls_period',
                     'first_ls_power','second_ls_power','third_ls_power','all_lomb_scargle_count','ab_value',
                     'mass_value','mass_min','mass_max','radius','radius_min','radius_max','OTYPE','SP_TYPE','MAIN_OTYPE']

        with open(downloaded_log, 'w') as f0: 
            f0.write(','.join(log_col_names)+'\n')

    if os.path.exists(piped_log):

        print('piped_log exists')

        result_lines = []

        with open(piped_log, 'r+') as f1: 
            for line in f1: 
                if '|||START_RESULT_LINE|||' in line: 
                    result_lines.append(line.split('|||START_RESULT_LINE|||')[-1].split("|||END_RESULT_LINE|||")[0].strip())
        
        with open(piped_log, 'w') as f1: 
            f1.truncate()

        print(len(result_lines))

        with open(downloaded_log, 'a') as f2: 
            for line in result_lines: 
                f2.write(line+'\n')

    extant_tic_ids = [] 
    if os.path.exists(downloaded_log): 
        with open(downloaded_log, 'r') as f: 
            for line in f: 
                tic_id = line.split(',')[0]
                if 'TIC_' in tic_id: 
                    extant_tic_ids.append(tic_id)

    tic_ids = list(np.setdiff1d(tic_ids, extant_tic_ids))

    for downloaded_tic_id in os.listdir(LombScargle_dir): 
        if downloaded_tic_id!='.gitkeep': 
            downloaded_tic_id = downloaded_tic_id.split('.')[0]
            if downloaded_tic_id not in tic_ids: 
                tic_ids.append(downloaded_tic_id)

    processed_data_col_names = ['clean_time', 'clean_flux', 'trend_time', 'trend_flux', 'no_flare_raw_time', 'no_flare_raw_flux', 'raw_time', 'raw_flux']

    for tic_id in tqdm(tic_ids):  

        #try: 
            
        # 1. DOWNLOAD DATA 
        raw_path = raw_dir + tic_id + '.csv'

        #sys.stderr = open(os.devnull, "w")
        raw_df, downloaded_sectors, sector_starts, sector_ends, last_dates = better_download(tic_id=tic_id, save_directory=raw_dir, verbose=True)
        #sys.stderr = sys.__stderr__
        
        if raw_df is not None: 
            out_list = []
            
            out_list.append(tic_id)
            out_list.append(downloaded_sectors) 
            out_list.append(sector_starts)
            out_list.append(sector_ends)

            raw_time, raw_flux = (np.array(raw_df[i]) for i in ['raw_time', 'raw_flux'])
            no_flare_mask = np.isfinite(raw_df['no_flare_raw_time'])
            no_flare_raw_time, no_flare_raw_flux = (np.array(raw_df[i])[no_flare_mask] for i in ['no_flare_raw_time', 'no_flare_raw_flux'])

            out_list.append(len(raw_time))
            out_list.append(len(no_flare_raw_time))

            # 2. PREPROCESS DATA

            # 2.a Wotan Detrend

            _, wotan_counts = better_preprocess(tic_id=tic_id, raw_data=raw_path, last_dates=last_dates, 
                                                save_directory=wotan_dir, window_length=0.3)

            out_list.append(wotan_counts)

            # 2.b Lomb-Scargle on Data
            periodogram = LombScargle(no_flare_raw_time, no_flare_raw_flux, nterms=2)

            frequencies, powers = periodogram.autopower(minimum_frequency=1/15, maximum_frequency=1/0.1, method='fastchi2')

            sort_idx = np.argsort(powers)[::-1]
            powers, frequencies = (i[sort_idx] for i in (powers, frequencies))
            periods = 1/frequencies 
            best_frequency = frequencies[0]

            out_list.append(periods[0])
            out_list.append(periods[1])
            out_list.append(periods[2])
            out_list.append(powers[0])
            out_list.append(powers[1])
            out_list.append(powers[2])

            # 2.c Lomb-Scargle Detrend

            y_fit = periodogram.model(no_flare_raw_time, best_frequency)

            (ls_clean_time, ls_clean_flux), (_, _) = remove_flares(no_flare_raw_time, no_flare_raw_flux/y_fit)

            out_list.append(len(ls_clean_time))

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
            out_list.append('"("'+str(ab[0])+','+str(ab[1])+')"')
            out_list.append(mass)
            out_list.append(mass_min)
            out_list.append(mass_max)
            out_list.append(radius)
            out_list.append(radius_min)
            out_list.append(radius_max) 

            # Simbad Information #
            _, simbad_results = query_simbad(tic_id)

            out_list.append(simbad_results[0])
            out_list.append(simbad_results[1])
            out_list.append(simbad_results[2])

            print('\n|||START_RESULT_LINE|||'+','.join([str(i) for i in out_list])+'|||END_RESULT_LINE|||\n')

            #sys.stderr = sys.__stderr__
        r'''
        except: 
            continue 

        finally: 
            pass 
        '''

import pandas as pd 
import numpy as np

tic_ids = np.array(pd.read_csv('/ar1/PROJ/fjuhsd/shared/tessodyssey/catalog-related/kerr_table_one_gaia_tic.csv')['TIC_ID'])
#tic_ids = np.setdiff1d(tic_ids, np.array(pd.read_csv('/ar1/PROJ/fjuhsd/shared/tessodyssey/data/temp_alpha_output.csv')['TIC_ID']))
tic_ids = list(tic_ids)

tic_ids = list(pd.read_csv('bad_ids.txt')['TIC_ID'])

#tic_ids = ['TIC_1232360','TIC_9966678']

create_base_catalog(tic_ids=tic_ids, piped_log='bad_ids_log.txt') 
# current PID: [1] 231907
#TIC_100398936
# randomly stopped ~1000 iterations in? does it get killed for using too much memory or something? 

# what will the sub functions do if data files already exist? 