import numpy as np
import pandas as pd
import os 
from tqdm import tqdm 

def in_term_check():


    TIC_IDs = []
    SDEs = []

    with open('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/beta_log.txt', 'r') as f: 
        for line in f: 
            if '|||START' in line: 
                line_list = line.split('|||START|||')[-1].split('|||END|||')[0].split(',')
                TIC_IDs.append(line_list[0])
                SDEs.append(float(line_list[3]))

    idx = np.argsort(SDEs)
    for i in idx: 
        print(TIC_IDs[i], SDEs[i])

def make_some_plots():
    from sunnyhills.plotting import tls_validation_mosaic
    import pickle 

    df = pd.read_csv('/ar1/PROJ/fjuhsd/shared/tessodyssey/routines/real/tls_gamma_run/results/tls_only_key.txt')
    df = df.dropna().sort_values(by='SDE', ascending=False)

    tic_ids = np.array(df['TIC_ID'])
    data_segments = np.array(df['data_segment'])
    detrend_methods = np.array(df['detrend_method'])
    SDEs = np.array(df['SDE'])
    #TIC_ID,data_segment,detrend_method,SDE

    cache_dir = '/ar1/PROJ/fjuhsd/shared/tessodyssey/routines/real/tls_gamma_run/cache_dir/'

    for index, tic_id in tqdm(enumerate(tic_ids)):
        try: 
            data_segment = data_segments[index]
            detrend_method = detrend_methods[index]
            SDE = SDEs[index]

            plot_path = '/ar1/PROJ/fjuhsd/shared/tessodyssey/routines/real/tls_gamma_run/plots/initial/'
            plot_path += str(SDE)+'_'+tic_id+'_m:'+detrend_method+'_s:'+str(data_segment)+'.png'

            if not os.path.exists(plot_path):

                data_path ='/ar1/PROJ/fjuhsd/shared/tessodyssey/data/processed/'    

                model_cache_path = cache_dir+tic_id+'_m:'+detrend_method+'_s:'+str(data_segment)+'_tls-model.pickle' 
                results_cache_path = cache_dir+tic_id+'_m:'+detrend_method+'_s:'+str(data_segment)+'_tls-results.pickle'

                if os.path.exists(results_cache_path) and os.path.exists(model_cache_path): 
                    with open(results_cache_path, 'rb') as file: 
                        tls_results = pickle.load(file)
                        
                    with open(model_cache_path, 'rb') as file: 
                        tls_model = pickle.load(file)


                if detrend_method == 'BiWeight-w=0.5':
                    data_path += 'BiWeight-w=0.5/'
                else: 
                    data_path += 'LombScargle-n=2/'

                data_path+=tic_id+'.csv'

                tls_validation_mosaic(tic_id=tic_id, data=data_path, tls_model=tls_model, tls_results=tls_results, false_alarms_dictionary={}, plot_path=plot_path)
        
        except Exception as e: 
            print(e)
            continue 

make_some_plots()