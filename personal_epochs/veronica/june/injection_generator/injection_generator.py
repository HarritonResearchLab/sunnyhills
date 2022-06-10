def injection_generator():
    import batman
    import numpy as np
    from astropy.constants import R_earth
    time_start = 3.14
    data_duration = 100 #keep it short for now
    samples_per_day = 720

    samples = int(data_duration*samples_per_day)
    time =  np.linspace(time_start,time_start+data_duration,samples)

    synthetic_transit = batman.TransitParams()
    # make the fake transits
    synthetic_transit.per = np.random.rand()*18+1
    synthetic_transit.rp = (np.random.rand()*18+2)
    synthetic_transit.t0 = time_start
    synthetic_transit.a = 19
    synthetic_transit.inc = 90
    synthetic_transit.ecc = 0
    synthetic_transit.w = 90
    synthetic_transit.u = [.4,.4]
    synthetic_transit.limb_dark = "quadratic"

    transit = batman.TransitModel(synthetic_transit,time)
    signal = transit.light_curve(synthetic_transit)
