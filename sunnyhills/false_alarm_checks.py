import numpy as np

def even_odd_transit(routine:str, stats):
    '''
    Arguments
        routine: BLS or TLS 
        stats: necessary object for computing left/right stats (for bls it is what is returned by compute_stats)

    '''

    import numpy as np

    even_odd_flag = False

    if routine=='BLS': 

        depth_odd = stats['depth_odd']
        err_odd = depth_odd[1]
        depth_odd = depth_odd[0]
        depth_even = stats['depth_even']
        err_even = depth_even[1]
        depth_even = depth_even[0]

        depths = np.sort([[depth_odd-err_odd, depth_odd+err_odd], 
                   [depth_even-err_even, depth_even+err_even]])

        if min(depths[1]>max(depths[0])): 
            even_odd_flag = True

    return even_odd_flag

def lombscargle(time,flux,flux_err:np.array=None,min_per:float=.1,max_per:int=15,calc_fap:bool=True,probabilities:list=[.1,.05,.01]):
    import numpy as np
    from astropy.timeseries import LombScargle
    
    best_period = 0
    best_period_power = 0

    fap_levels = None
    periodogram = LombScargle(time,flux,flux_err)
    frequencies,powers = periodogram.autopower(minimum_frequency=1/max_per, maximum_frequency=1/min_per)

    periods = 1/frequencies

    if calc_fap:
        fap_levels = periodogram.false_alarm_probability(probabilities)

    sorted = np.argsort(powers)[::-1] #descending order
    powers = powers[sorted]
    periods = periods[sorted]

    if len(sorted)>0: 
        best_period = periods[0] 
        best_period_power = powers[0]

    return powers,periods,best_period,best_period_power,fap_levels