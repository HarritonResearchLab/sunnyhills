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


def bls_check(table_path:str,lc_dir:str,display_mosaics:bool=False):  
  '''
  runs the bls pipeline on known light curves and then computes the average error of the dataset
  
  args:
  table_path: path to confirmed candidates file
  lc_dir: path to light curves to compare to
  display_mosaics: option to display mosaic plots 
  '''    
  import pandas as pd    
  import numpy as np
  import matplotlib.pyplot as plt
  from tqdm import tqdm
  from sunnyhills.pipeline_functions import run_bls

  error = np.array([])
  df = filter_data(table_path)
  modified_dir = lc_dir + '/'
  stellar_rads = np.array([])
  planet_rads = np.array([])
  graphed_error = np.array([])
  for item,real_per,stellar_rad,planet_rad in tqdm(zip(df['tic_id'],df['pl_orbper'],df['st_rad'],df['pl_rade'])):
    print(item)
    lc = pd.read_csv(modified_dir+item.replace(' ','_')+'.csv')
    lc = lc.dropna(subset=['clean_time','clean_flux','trend_time','trend_flux'])
    time = lc['clean_time']
    flux = lc['clean_flux']
    trend_time = lc['trend_time']
    trend_flux = lc['trend_flux']
    best_params, bls_model, in_transit,stats = run_bls(time, flux)
    error = np.append(error,(real_per-best_params[0])/real_per)
    if display_mosaics:
      bls_validation_mosaic(item.replace('TIC',''), time, flux, trend_time, trend_flux, time, flux, best_params, bls_model, in_transit,stats)
    if stellar_rad != np.NaN and planet_rad != np.NaN:
      graphed_error = np.append(graphed_error,(real_per-best_params[0])/real_per)
      stellar_rads = np.append(stellar_rads,stellar_rad)    
      planet_rads = np.append(planet_rads,planet_rad*109)
  plt.figure(112)
  plt.title('radius vs error')
  plt.scatter(np.divide(planet_rads,stellar_rads),np.log(100*np.absolute(graphed_error)))
  plt.xlabel('ratio of planetary radius to stellar radius')
  plt.ylabel('log error percentage')
  plt.savefig('data/bls_check/error_vs_radius_scatter.png')
  plt.figure(122)
  plt.title('error histogram')
  plt.hist(100*np.absolute(error),50)
  plt.xticks(np.arange(0, 350, 25))
  plt.xlabel('error percentage')
  plt.ylabel('frequency')

  error = 100*np.mean(np.absolute(error))
  print("average error for all the samples is: ", error)
  plt.legend([str(error)])
  plt.savefig('data/bls_check/percentage_histogram.png')

  return error

def bls_check_runner(display_validation_plots:bool=False):
  bls_check('personal_epochs/veronica/june/current_key.csv','data/current/processed/known_two_min',display_validation_plots)


def tls_check(table_path:str,lc_dir:str):
  from sunnyhills.pipeline_functions import run_tls    
  import pandas as pd    
  import numpy as np
  import matplotlib.pyplot as plt
  from tqdm import tqdm



  error = np.array([])
  df = filter_data(table_path)
  modified_dir = lc_dir + '/'
  stellar_rads = np.array([])
  planet_rads = np.array([])
  graphed_error = np.array([])
  for item,real_per,stellar_rad,planet_rad in tqdm(zip(df['tic_id'],df['pl_orbper'],df['st_rad'],df['pl_rade'])):
    item = item.replace(' ','_')
    print(item)
    lc = pd.read_csv(modified_dir+item+'.csv')
    lc = lc.dropna(subset=['clean_time','clean_flux','trend_time','trend_flux'])
    time = lc['clean_time']
    flux = lc['clean_flux']
    results,model = run_tls(item,time, flux)
    per = results.period
    error = np.append(error,(real_per-per)/real_per)
    if stellar_rad != np.NaN and planet_rad != np.NaN:
      graphed_error = np.append(graphed_error,(real_per-per)/real_per)
      stellar_rads = np.append(stellar_rads,stellar_rad)    
      planet_rads = np.append(planet_rads,planet_rad*109)
  plt.figure(112)
  plt.title('radius vs error')
  plt.scatter(np.divide(planet_rads,stellar_rads),np.log(100*np.absolute(graphed_error)))
  plt.xlabel('ratio of planetary radius to stellar radius')
  plt.ylabel('log error percentage')
  plt.savefig('data/tls_check/error_vs_radius_scatter.png')
  plt.figure(122)
  plt.title('error histogram')
  plt.hist(100*np.absolute(error),50)
  plt.xticks(np.arange(0, 350, 25))
  plt.xlabel('error percentage')
  plt.ylabel('frequency')

  error = 100*np.mean(np.absolute(error))
  print("average error for all the samples is: ", error)
  plt.legend([str(error)])
  plt.savefig('data/tls_check/percentage_histogram.png')

  return error

def tls_check_runner():
  tls_check('personal_epochs/veronica/june/current_key.csv','data/current/processed/known_two_min')