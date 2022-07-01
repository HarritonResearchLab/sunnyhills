def compare(tic_id:str='TIC_197005132', data_dir:str='./data/current/processed/two_min_lightcurves/'): 
    from sunnyhills.injection import inject
    from sunnyhills.pipeline_functions import query_tls_vizier
    import pandas as pd
    import numpy as np
    from pytransit import QuadraticModel

    import matplotlib.pyplot as plt 

    data_path = data_dir+tic_id+'.csv'

    ab, mass, mass_min, mass_max, radius, radius_min, radius_max = query_tls_vizier(tic_id) 

    ab = [ab[0], ab[1]]

    R_J = 0.10045 
    R_P = 0.5*R_J 

    rp_over_rs = R_P/radius 

    per = 2.93

    raw_data = pd.read_csv(data_path)
    mask = np.isfinite(raw_data['no_flare_raw_time'])
    raw_time, raw_flux = (np.array(raw_data[i])[mask] for i in ['no_flare_raw_time', 'no_flare_raw_flux'])

    diff = np.diff(raw_time)
    idx = np.where(diff>75)[0]
    if len(idx)>0: 
        idx = idx[0]

        raw_time = raw_time[0:idx]
        raw_flux = raw_flux[0:idx]

    t0 = np.percentile(raw_time, q=2.2)

    print('###TIC_ID###')
    print(tic_id)

    print('stellar radius: ', radius)

    scale_factor = 1/radius
    print('scale factor:', scale_factor)

    scaled = R_P*scale_factor

    print('original R_P:', str(R_P)+' R_SUN')
    print('scaled R_P:', scaled)

    raw_time_batman, raw_flux_batman, (_, _, _) = inject(tic_id, time=raw_time, flux=raw_flux, 
                                                         per=per, t0=t0, rp=scaled)
    r'''
    tm = QuadraticModel()
    tm.set_data(raw_time)

    raw_time_pytransit = raw_time 

    raw_flux_pytransit = tm.evaluate(k=0.4, ldc=ab, t0=t0, p=per, a=15, i=87)
    '''
    transits = np.arange(t0, np.max(raw_time), step=per)

    fig, ax = plt.subplots(figsize=(7,3))

    ax.scatter(raw_time_batman, raw_flux_batman, color='grey', s=2)
    ax.set(title='BATMAN', xlabel='Time (d)', ylabel='Flux')

    for transit in transits: 
        ax.axvline(x=transit, color='red', alpha=0.5, lw=2)
    
    plt.subplots_adjust(hspace=0.5)
    plt.tight_layout()
    plt.savefig('./personal_epochs/thaddaeus/june/topical/pytransit_vs_batman/'+tic_id+'.png', dpi=150)
    
for tic_id in ['TIC_43451370', 'TIC_197005132', 'TIC_436783866', 'TIC_266002317', 'TIC_232965288',
               'TIC_372918535', 'TIC_409372963']:
    compare(tic_id=tic_id)