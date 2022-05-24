
def lomb_scargle(time,flux,flux_err:np.Array()=None,min_per:float=.1,max_per:int=15,calc_fap:boolean=True):
    import numpy as np
    from astropy.timeseries import LombScargle
    
    best_period = 0
    best_period_power = 0

    fap_levels = None
    periodogram = LombScargle(time,flux,flux_err)
    frequencies,powers = periodogram.autopower(min_per,max_per)

    periods = 1/frequencies

    if calc_fap:
        fap_levels = periodogram.false_alarm_probability(min_per,max_per)
        fap_99 = fap_levels[0]
        fap_mask = np.where(powers>fap_99)
        periods = periods[fap_mask]
        powers = powers[fap_mask]

    sorted = np.argsort(powers)[::-1] #descending order
    powers = powers[sorted]
    periods = periods[sorted]

    if len(sorted)>0: 
        best_period = periods[0] 
        best_period_power = powers[0]

    return powers,periods,best_period,best_period_power,fap_levels


    
