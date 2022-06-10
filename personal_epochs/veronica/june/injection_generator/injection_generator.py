def injection_generator():
    import batman
    import numpy as np
    from astropy.constants import R_earth
    
    time_start = 3.14
    data_duration = 5 #keep it short for now
    samples_per_day = 720
    
    samples = int(data_duration*samples_per_day)
    time =  np.linspace(time_start,time_start+data_duration,samples)
    
    synthetic_transit = batman.TransitParams()
    # make the fake transits
    synthetic_transit.per = np.random.rand()*18+1
    synthetic_transit.rp = (np.random.rand()*18+2)*R_earth.value
    