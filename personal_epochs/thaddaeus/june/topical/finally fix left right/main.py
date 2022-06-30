from sunnyhills.false_alarm_checks import check_lombscargle, tls_even_odd, transit_outliers_fap_test
from sunnyhills.pipeline_functions import run_tls
from sunnyhills.plotting import tls_validation_mosaic
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt 


tic_id = 'TIC_239134011'
data_path = './routines/simulations/first_bulk_injected/data/'+tic_id+'.csv'
data = pd.read_csv(data_path)
mask = np.isfinite(data['clean_flux'])
time, flux = (np.array(data[i])[mask] for i in ['clean_time', 'clean_flux'])

cache_dir = './personal_epochs/thaddaeus/june/topical/finally fix left right/'

tls_results, tls_model = run_tls(tic_id=tic_id, time=time, flux=flux, cache_dir=cache_dir, 
                                 tls_params={'min_per':0.5, 
                                             'max_per':15,
                                             'minimum_n_transit':3,
                                             'freq_factor':5,
                                             'core_fraction':0.9})

fig, ax = plt.subplots()

per = tls_results.period 
t0 = tls_results.T0

from sunnyhills.plotting import phased_aliase_plots

phased_aliase_plots(tic_id, time, flux, tls_results, plot_path='./personal_epochs/thaddaeus/june/topical/finally fix left right/testing_phasing.png')

print(tls_results.keys())
print(tls_results.odd_even_mismatch)
quit()

plot_dir = './personal_epochs/thaddaeus/june/topical/finally fix left right/'
download_log = './routines/real/alpha_tls/data/download_log.csv'

all_false_alarms_dict = {}

lombscargle_dict = check_lombscargle(tic_id=tic_id, tls_results=tls_results, download_log=download_log) 
even_odd_dict = tls_even_odd(tls_results=tls_results)
transit_outliers_dict = transit_outliers_fap_test(tls_results=tls_results)

for alarm in [lombscargle_dict, even_odd_dict, transit_outliers_dict]: 
    all_false_alarms_dict.update(alarm)

tls_validation_mosaic(tic_id=tic_id, data=data_path, tls_results=tls_results, tls_model=tls_model, plot_dir=plot_dir, false_alarms_dictionary=all_false_alarms_dict, plot_type='png')


