import pandas as pd

import numpy as np

df = pd.read_csv('./data/current/processed/two_min_lightcurves/TIC_9966678.csv')

time = np.array(df['cleaned_time'])

not_nan = np.isfinite(time)

time = time[not_nan]




flux = np.array(df['detrended_flux'])

flux = flux[not_nan]    

print(time)