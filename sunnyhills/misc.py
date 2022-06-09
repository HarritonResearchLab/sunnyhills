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
    
    step = int(len(x)/num_bins)

    ranges = []

    for i in range(0, len(x)-step, step): 
        ranges.append([i,i+step])

    x_rebinned = np.array([np.mean(x[range[0]:range[1]]) for range in ranges])
    y_rebinned = np.array([np.mean(y[range[0]:range[1]]) for range in ranges])

    return x_rebinned, y_rebinned 

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

        return np.flatten(return_list)
