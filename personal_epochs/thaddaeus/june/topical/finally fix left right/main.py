from sunnyhills.pipeline_functions import run_tls 
from sunnyhills.plotting import tls_validation_mosaic
import numpy as np
import pandas as pd

tic_id = 'TIC_239134011'
data_path = './routines/simulations/first_bulk_injected/data/'+tic_id+'.csv'
data = pd.read_csv(data_path)
mask = np.isfinite(data['clean_flux'])
time, flux = (np.array(data[i])[mask] for i in ['clean_time', 'clean_flux'])

tls_results, tls_model = run_tls(tic_id=tic_id, time=time, flux=flux, tls_params={'min_per':0.5, 'max_per':15, 'minimum_n_transit':3, 'freq_factor':5, 'core_fraction':0.9})

plot_dir = './personal_epochs/thaddaeus/june/topical/finally fix left right/'
tls_validation_mosaic(tic_id, data=data_path, tls_model=tls_model, tls_results=tls_results, plot_dir=plot_dir, plot_type='png')