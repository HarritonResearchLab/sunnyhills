from matplotlib.pyplot import table
from sunnyhills.misc import download_and_append_status


def filter_data(path:str):
    '''
    filters out the data from the known exoplanets table
    '''
    import pandas as pd
    import numpy as np

    df = pd.read_csv(path)
    # drop initials nans
    df = df.dropna(subset=['pl_orbper','sy_pnum','pl_controv_flag'])

    df = df.where(np.logical_and(df['pl_orbper']>.5,df['pl_orbper']<15))
    df = df.where(df['pl_controv_flag']==0)
    df = df.where(df['sy_pnum']==1)
    df = df.where(df['has_two_min']==True)

    # filter out the data where the given fields did not meet the conditions
    filtered_df = df.dropna(subset=['pl_orbper','pl_controv_flag','sy_pnum','tic_id','has_two_min'])

    return filtered_df


def bls_check(table_path:str,lc_dir:str,display_mosaics:bool):
  import pandas as pd    
  import numpy as np
  import matplotlib.pyplot as plt
  from tqdm import tqdm
  from sunnyhills.plotting import bls_validation_mosaic
  from sunnyhills.misc import get_best_period
  from sunnyhills.pipeline_functions import run_bls


  error = np.array([])
  df = filter_data(table_path)
  modified_dir = lc_dir + '/'
  stellar_rads = np.array([])
  planet_rads = np.array([])
  graphed_error = np.array([])
  for item,real_per,stellar_rad,planet_rad in tqdm(zip(df['tic_id'][:5],df['pl_orbper'][:5],df['st_rad'][:5],df['pl_rade'][:5])):
    print(item)
    lc = pd.read_csv(modified_dir+item.replace(' ','_')+'.csv')
    lc = lc.dropna(subset=['clean_time','clean_flux'])
    time = lc['clean_time']
    flux = lc['clean_flux']
    best_params, results, bls_model, in_transit, stats = run_bls(time, flux)
    error = np.append(error,(real_per-get_best_period(results))/real_per)
    if display_mosaics:
      try:
        bls_validation_mosaic(item.replace('TIC',''), time, flux, time, flux, best_params, results, bls_model, in_transit, stats)
      except ValueError:
        print('trouble displaying validation mosaic')
    if stellar_rad != np.NaN and planet_rad != np.NaN:
      graphed_error = np.append(graphed_error,(real_per-get_best_period(results))/real_per)
      stellar_rads = np.append(stellar_rads,stellar_rad)    
      planet_rads = np.append(planet_rads,planet_rad*109)
  plt.figure(112)
  plt.scatter(np.absolute(graphed_error),np.divide(planet_rads,stellar_rads))
  plt.xlabel('ratio of planetary radius to stellar radius')
  plt.ylabel('error percentage')
  plt.figure(122)
  plt.hist(np.absolute(error),50)
  plt.xlabel('ratio of planetary radius to stellar radius')
  plt.ylabel('error percentage')
          
  error = 100*np.mean(np.absolute(error))
  print("average error for all the samples is: ", error)


  return error