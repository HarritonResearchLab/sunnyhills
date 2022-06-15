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

def even_odd_phase_folded(time, flux, results):
    """Return even odd phase folded transits"""

    all_even_indices_in_transit = np.array([])
    all_odd_indices_in_transit = np.array([])

    all_indices = np.arange(0,len(time))

    transit_times = results.transit_times
    transit_duration_in_days = results.duration

    period = results.period
    T0 = results.T0

    for i in range(len(transit_times)):
        mid_transit = transit_times[i]
        tmin = mid_transit -  transit_duration_in_days
        tmax = mid_transit + transit_duration_in_days
        if np.isnan(tmin) or np.isnan(tmax):
            idx_intransit = []
            mean_flux = np.nan
        else:
            idx_intransit = np.where(np.logical_and(time > tmin, time < tmax))[0]

        # Check if transit odd/even to collect the flux for the mean calculations
        if i % 2 == 0:  # even
            all_even_indices_in_transit = np.concatenate((all_even_indices_in_transit, 
                                                          idx_intransit))
        else:  # odd
            all_odd_indices_in_transit = np.concatenate((all_odd_indices_in_transit,
                                                         idx_intransit))

    def fold(time, period, T0):
        """Normal phase folding"""
        #return (time - T0) / period - np.floor((time - T0) / period)
        return (time) / period - np.floor((time) / period)

    all_even_indices_in_transit = all_even_indices_in_transit.astype(int)
    all_odd_indices_in_transit = all_odd_indices_in_transit.astype(int)

    even_transit_flux = flux[all_even_indices_in_transit] 
    even_transit_time = time[all_even_indices_in_transit]
    odd_transit_flux = flux[all_odd_indices_in_transit]
    odd_transit_time = time[all_odd_indices_in_transit]

    even_transit_time_folded = fold(even_transit_time, period, T0)
    odd_transit_time_folded = fold(odd_transit_time, period, T0)

    return even_transit_time_folded, even_transit_flux, odd_transit_time_folded, odd_transit_flux, all_even_indices_in_transit, all_odd_indices_in_transit

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