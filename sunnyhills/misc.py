import numpy as np

def gaia_to_tic(gaia_ids):
    from astrobase.services.identifiers import gaiadr2_to_tic
    import numpy as np
    tic_ids = []
    for id in gaia_ids:
        if id !=None:
            tic_ids.append(gaiadr2_to_tic(id))

        else: 
            tic_ids.append(None)

    return dict(zip(gaia_ids, tic_ids))

def skycoord_to_tic(ra,dec):
    from astropy.coordinates import SkyCoord
    from eleanor import mast
    
    tic_id = mast.tic_from_coords((ra,dec))
    return tic_id

def get_best_period(periodogram):
  import numpy as np
  return periodogram.period[np.argmax(periodogram.power)]

def rebin(x, y, num_bins:int=20): 
    '''
    arguments: 
        x: time values 
        y: corresponding flux values
        num_bins: number of output points; default is 20
    returns: 
        x_rebinned: rebinned x 
        y_rebinned: rebinned y 
    '''
    success = True
    if len(x)>0: 
        step = int(len(x)/num_bins)

        ranges = []

        for i in range(0, len(x)-step, step): 
            ranges.append([i,i+step])

        x_rebinned = np.array([np.mean(x[range[0]:range[1]]) for range in ranges])
        y_rebinned = np.array([np.mean(y[range[0]:range[1]]) for range in ranges])

    else: 
        x_rebinned, y_rebinned = (x,y)
        success=False

    return x_rebinned, y_rebinned, success 

def phase(time, flux, period:float, t0:float=None, duration:float=None, bls_model=None, model_name:str=None, tls_results=None, fraction:int=0.2): 
        '''
        Arguments: 
            time: array of time values
            flux: corresponding array of flux values
            period: period to fold to 
            t0: t0 from BLS/TLS 
            bls_model: bls model (if model_name='BLS', otherwise None)
            fraction: how far from median time to go; default is 0.2
        Returns: 
            return_list: [folded_x, folded_flux] by default for model_name in (BLS, LS, TLS)
                         Note: if model_name is BLS, you need to give bls_model the bls_model you used to fit, so it can calculate the model transit shape
                         Note: if model_name is TLS, then you need to give tls_results the results from model.power() from tls bc TLS already calcs fold/not fold
                         Note: in all cases cases, as noted above, they will return folded model time and flux in addition to (after) folded data time and flux
        '''
        
        import numpy as np
        import warnings

        return_list = []

        if model_name!=None: 
            if model_name=='BLS' or model_name=='TLS': 

                x = (time - t0 + 0.5*period) % period - 0.5*period
                m = np.abs(x) < fraction
                
                return_list = [x[m], flux[m]]

                if model_name=='BLS' and bls_model!=None: 
                    f = bls_model.model(x + t0, period, duration, t0) # for BLS
                    return_list.append([x,f])

                elif model_name=='TLS' and tls_results!=None: 
                    return_list.append([tls_results.model_folded_phase,tls_results.model_folded_flux])

            elif model_name=='LS': 
                from astropy.timeseries import LombScargle

                phased_dates = (np.mod(time, period))/period
                return_list.append(phased_dates)
                return_list.append(flux)
                
                t_fit = np.linspace(0, 1)
                ls = LombScargle(time, flux)
                y_fit = ls.model(t_fit, 1/period)

                return_list.append([t_fit, y_fit])

        return np.array(return_list, dtype=object).flatten()

## BELOW FUNCTIONS ARE VERONICA'S FOR STARS WITH CONFIRMED PLANETS ##

def download_data(path:str,dir:str=''):
    from sunnyhills.pipeline_functions import download_and_preprocess
    import pandas as pd
    import numpy as np
    from lightkurve.utils import LightkurveError
    from tqdm import tqdm

    df = pd.read_csv(path)
    # drop initials nans
    df = df.dropna(subset=['pl_orbper','sy_pnum','pl_controv_flag'])

    df = df.where(np.logical_and(df['pl_orbper']>.5,df['pl_orbper']<15))
    df = df.where(df['pl_controv_flag']==0)
    df = df.where(df['sy_pnum']==1)

    # filter out the data where the given fields did not meet the conditions
    filtered_df = df.dropna(subset=['pl_orbper','pl_controv_flag','sy_pnum','tic_id'])
    
    
    for item in tqdm(filtered_df['tic_id']):
        try:
          lc = download_and_preprocess(item,dir)
        except (KeyError,LightkurveError):
            pass
     
def append_download_status(path_to_csv:str,dir:str,save_dir):
    import os
    import pandas as pd
    import numpy as np
    
    df = pd.read_csv(path_to_csv)

    downloaded_lcs = np.array(os.listdir(dir))
    new_arr = []
    for lc in downloaded_lcs:
      new_arr.append(os.path.splitext(lc)[0].replace('_',' '))
    downloaded_2_min = []
    downloaded_lcs = new_arr
    
    for item in df['tic_id']:
      if item in downloaded_lcs[:]:
        downloaded_2_min.append(True)    
      else:
        downloaded_2_min.append(False)    
    downloaded_2_min = pd.Series(downloaded_2_min,name="has_two_min")

    df = df.merge(downloaded_2_min,left_index=True,right_index=True)
    df = df.drop_duplicates(subset=['tic_id'])
    df.to_csv(save_dir+'/current_key.csv')

    return df

def download_and_append_status(path_to_csv:str,lc_dir,save_dir:str):
    download_data(path_to_csv,lc_dir)
    append_download_status(path_to_csv,lc_dir,save_dir)

    return np.abs((time-transit_time+hp) % period - hp) < .5*duration

def merge_pdf(pdfs_dir:str,pages_per_pdf:int=5):
    import os
    import PyPDF2

    pdfs = os.listdir(pdfs)

    