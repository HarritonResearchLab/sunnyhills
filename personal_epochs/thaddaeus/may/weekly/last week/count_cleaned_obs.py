import os 
import numpy as np
import pandas as pd

df = pd.read_csv('./data/current/lightcurve_key.csv')

df['0'] = [int(i.replace('TIC_','')) for i in df['0']]

export_df = pd.read_csv('./personal_epochs/veronica/may/toi_matching/toi_export.csv')
export_df['tid'] = [int(i) for i in export_df['tid']]

merged = export_df.merge(df, left_on='tid', right_on='0') 
print(merged)
print(merged['tfopwg_disp'])
