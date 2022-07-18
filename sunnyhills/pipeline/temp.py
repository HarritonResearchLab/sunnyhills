import numpy as np
import pandas as pd 
import os 

def grab_lines(): 

    isfinite_IDs = []
    result_lines = []
    successful_IDs = []
    no_data_ids = []

    with open('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/log.txt', 'r') as f: 
        for line in f: 
            start_not_found = True 
            if  '|||HAS_DATA' in line: 
                line_list = line.split(',') 
                tic_id = line_list[2].split('|||')[0] 
                if line_list[1] == 'FALSE': 
                    no_data_ids.append(tic_id) 
                else: 
                    isfinite_IDs.append(tic_id)
               
            elif '|||START_RESULT_LINE|||' in line: 
                result_line = line.split('|||START_RESULT_LINE|||')[-1].split("|||END_RESULT_LINE|||")[0].strip()
                result_lines.append(result_line)

                successful_IDs.append(result_line.split(',')[0])

    isfinite_IDs = np.setdiff1d(isfinite_IDs, successful_IDs)

    with open('/ar1/PROJ/fjuhsd/shared/tessodyssey/data/temp_alpha_output.csv', 'a') as f2: 
        for line in result_lines: 
            f2.write(line+'\n')

    missing_ids_log = '/ar1/PROJ/fjuhsd/shared/tessodyssey/data/no_data_TICs.csv'

    no_data_ids = np.concatenate((no_data_ids, list(pd.read_csv(missing_ids_log)['TIC_ID'])))
    
    if len(no_data_ids)>0: 
        no_data_ids_df = pd.DataFrame(no_data_ids, columns=['TIC_ID'])
        no_data_ids_df.to_csv(missing_ids_log, index=False)

    with open('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/current_weird_isfinite.txt', 'w') as f: 
        f.write('TIC_ID'+'\n') 
        for i in isfinite_IDs: 
            f.write(i+'\n') 

def count_no_isfinite(): 
    
    is_finite_IDs = np.array(pd.read_csv('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/current_weird_isfinite.txt')['TIC_ID'])
    for tic_id in is_finite_IDs: 
        df = pd.read_csv('/ar1/PROJ/fjuhsd/shared/tessodyssey/data/raw/'+tic_id+'.csv')
        if len(df.index)>0: 
            print('non-zero')