def injection_generator(sectors:int =1, noise_level:int =50):
    '''
    Arguments:
    sectors: the amount of sectors to sample, light curve will be longer if more sectors are specified
    noise_level: flux noise level in ppm
    '''
    import batman
    import numpy as np
    from sunnyhills.pipeline_functions import run_bls

    time_start = 3.14
    data_duration = 27.4*sectors #keep it short for now
    samples_per_day = 720

    samples = int(data_duration*samples_per_day)
    time =  np.linspace(time_start,time_start+data_duration,samples)
    period = 0
    
    while not(period<.5 and period>15):
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
        noise = numpy.random.normal(0, 10**-6 * noise_level, int(samples))

        best_params, results, bls_model, in_transit, stats = run_bls(time,signal+noise)

        period = best_params[1]
    
    return results,best_params


