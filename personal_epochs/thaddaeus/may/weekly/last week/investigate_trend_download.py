import pandas as pd
import numpy as np 
from sunnyhills.pipeline_functions import download_and_preprocess 
import warnings 
warnings.filterwarnings("ignore")

from tqdm import tqdm 

for id in tqdm(['TIC 36352297', 'TIC 405010127', 'TIC 75266189', 'TIC 21744120']): 
    download_and_preprocess(id, './personal_epochs/thaddaeus/may/weekly/last week/quick_download_check')
    from sunnyhills.plotting import bls_validation_mosaic 

    from sunnyhills.pipeline_functions import run_bls

    data_csv = pd.read_csv('./personal_epochs/thaddaeus/may/weekly/last week/quick_download_check/'+id.replace(' ', '_')+'.csv')
    clean_time, clean_flux = (np.array(data_csv[i]) for i in ['cleaned_time','detrended_flux'])

    mask = np.isfinite(clean_flux)
    clean_time, clean_flux = (arr[mask] for arr in (clean_time, clean_flux))

    raw_time,raw_flux = (np.array(data_csv[i]) for i in ['raw_time','raw_flux'])
    mask = np.isfinite(raw_flux)
    raw_time, raw_flux = (arr[mask] for arr in (raw_time, raw_flux))

    bls_results, bls_model, in_transit, bls_stats = run_bls(clean_time, clean_flux)

    bls_validation_mosaic(id, clean_time, clean_flux, 
                          raw_time, raw_flux, 
                          bls_results, bls_model, in_transit, bls_stats, 
                          path='./personal_epochs/thaddaeus/may/weekly/last week/quick_download_check/'+id+'.png')