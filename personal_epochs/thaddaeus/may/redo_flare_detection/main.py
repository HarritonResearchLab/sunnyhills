import pandas as pd

df = pd.read_csv(r'C:\Users\Research\Documents\GitHub\sunnyhills\data\archive\LightCurve_keys\TIC 140045538.csv')

time = pd.Series(df['time'])
flux = pd.Series(df['pdcsap_flux'])
fluxerr = pd.Series(df['pdcsap_flux_err'])

def remove_flares(time, flux, flux_err=None, sigma:int=3): 
    ''' 
    Args:
        time: array-like (ideally Pandas.series) object of time values 
        flux: array-like (ideally Pandas.series) object of flux values
        flux_err: optional array-like (ideally Pandas.series) object of flux err values
    Returns: 
        two tuples (or three if an array was passed into flux_err); first tuple gives arrays that have been "cleaned" of flairs, and second gives the removed values. Array order in each tuple is time, flux, (flux_err)
    Additional Info: 
        Idea was taken from https://iopscience.iop.org/article/10.3847/1538-3881/ab5d3a
    '''
    import pandas as pd
    import numpy as np


    if isinstance(time, pd.Series) and isinstance(flux, pd.Series): 
        pass
    else: 
        time = pd.Series(time)
        flux = pd.Series(flux)
        if flux_err!=None: 
            flux_err = pd.Series(flux_err)
            removed_flux_err = np.array([])

    removed_time, removed_flux = (np.array([]), np.array([]))

    continue_global = True 
    while continue_global:
        length = len(flux)

        # # We use three copies and extract the flares from the middle to prevent edge effects --> BREDALLL! LOL
        three_flux = pd.concat((flux, flux, flux))
        global_rolling = three_flux.rolling(1024, center=False) 
        global_medians = global_rolling.median()[length:2*length]
        global_stds = global_rolling.std()[length:2*length]
        cutoffs = global_medians+sigma*global_stds

        remove_indices = []
        for i in flux.index: 
            if flux[i]>cutoffs[i] and global_medians[i]!=np.nan: 
                remove_indices.append(i)

        if len(remove_indices)==0: 
            continue_global = False 

        else:     
            removed_time = np.concatenate((removed_time, time[remove_indices])) 
            removed_flux = np.concatenate((removed_flux, flux[remove_indices])) 

            time = time.drop(remove_indices)
            flux = flux.drop(remove_indices)

            if flux_err!=None: 
                removed_flux_err = np.concatenate((removed_flux_err, flux_err[remove_indices]))
                flux_err = flux_err.drop(remove_indices)
            
    continue_local = True 
    while continue_local: 
        local_rolling = flux.rolling(128, center=True) 
        local_medians = local_rolling.median() 
        local_stds = local_rolling.std() 
        cutoffs = local_medians+sigma*local_stds 

        remove_indices = [] 
        for i in flux.index: 
            if flux[i]>cutoffs[i] and local_medians[i]!=np.nan: 
                remove_indices.append(i) 

        if len(remove_indices)==0: 
            continue_local = False 

        else:     
            removed_time = np.concatenate((removed_time, time[remove_indices]))
            removed_flux = np.concatenate((removed_flux, flux[remove_indices]))

            time = time.drop(remove_indices)
            flux = flux.drop(remove_indices)

            if flux_err!=None: 
                removed_flux_err = np.concatenate((removed_flux_err, flux_err[remove_indices]))
                flux_err = flux_err.drop(remove_indices)

    if flux_err!=None: 
        return (time.to_numpy(), flux.to_numpy(), flux_err.to_numpy()), (removed_time, removed_flux, removed_flux_err)
    else: 
        return (time.to_numpy(), flux.to_numpy()), (removed_time, removed_flux)


(clean_time, clean_flux), (flair_time, flair_flux) = remove_flares(time, flux)