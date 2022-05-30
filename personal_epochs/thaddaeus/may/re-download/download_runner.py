# this took 7 hours and 34 minutes to download and process 430 light curves lol

import warnings
from tqdm import tqdm 
import pandas as pd
import numpy as np 

from sunnyhills.pipeline_functions import download_and_preprocess 

warnings.filterwarnings("ignore")

ids = np.array(pd.read_csv('./data/current/current_key.csv')['tid'])

import os 

completed_ids = np.array([i.replace('.csv', '') for i in os.listdir('./data/current/processed/two_min_lightcurves') if i!='.gitkeep'])
completed_ids = np.array([int(i.split('_')[-1]) for i in completed_ids])

ids = np.setdiff1d(ids, completed_ids)

for id in tqdm(ids): 
    try: 
        tic_id = 'TIC '+str(id).replace('_', " ")
        download_and_preprocess(tic_id, './data/current/processed/two_min_lightcurves')
        
    except Exception as e: 
        print(e)
        continue

    break 