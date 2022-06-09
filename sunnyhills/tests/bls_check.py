def filter_data(path:str):
    '''
    filters out the data from the known exoplanets table
    '''
    from sunnyhills.pipeline_functions import download
    import pandas as pd
    import numpy as np
    import math
    from lightkurve.utils import LightkurveDeprecationWarning, LightkurveError

    df = pd.read_csv(path)
    # drop initials nans
    df = df.dropna(subset=['pl_orbper','sy_pnum','pl_controv_flag'])

    df = df.where(np.logical_and(df['pl_orbper']>.5,df['pl_orbper']<15))
    df = df.where(df['pl_controv_flag']==0)
    df = df.where(df['sy_pnum']==1)

    # filter out the data where the given fields did not meet the conditions
    filtered_df = df.dropna(subset=['pl_orbper','pl_controv_flag','sy_pnum','tic_id'])

    testable_ids = []
    periods = []
    for item,per in zip(filtered_df['tic_id'],filtered_df['pl_orbper']):
      try:
        result,found = download(item)
        if found==True:
          testable_ids.append(item)
          periods.append(per)
        else:
          pass
      except (KeyError,LightkurveError):
          pass
    return testable_ids,periods

def bls_check(download_dir:str,display_mosaics:bool=False):
    '''
    calculates the error between the periods calculated by the bls pipeline and the actual discovered orbital periods
    
    arguments:
    download_dir: the download directory for any tic ids that were not downloaded and preprocessed
    display_mosaics: whether to display the bls mosaic for each id parsed
    '''
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    from sunnyhills.pipeline_functions import download_and_preprocess,run_bls
    from sunnyhills.plotting import bls_validation_mosaic
    from sunnyhills.misc import get_best_period
    import warnings

    ids,periods = filter_data('/content/known_TESS_planets.csv')
    error = np.array([])
    not_found = np.array([])
    for tid,real_per in zip(ids,periods):
        lc = download_and_preprocess(tid)
        try:
          best_params, results, bls_model, in_transit, stats = run_bls(np.array(lc[0][0]), np.array(lc[0][1]))
          error = np.append(error,(real_per-get_best_period(results))/real_per)
          if display_mosaics:
            bls_validation_mosaic(tid, lc[0][0], lc[0][1], lc[0][0], lc[0][1], best_params, results, bls_model, in_transit, stats)
        except TypeError:
          not_found = np.append(not_found,tid)
    error = 100*np.mean(np.absolute(error))
    print("average error for all the samples is: ", error)
    warnings.warn('please contact Veronica with the results of the not_found file, if it is empty then that is good!')
    not_found.tofile(download_dir+'/search_errors.csv',sep=',')
    return error,not_found