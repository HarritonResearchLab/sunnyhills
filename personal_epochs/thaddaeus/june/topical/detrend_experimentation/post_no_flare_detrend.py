def detrend(time, flux, period): 

    if 0 < period <= 2: 
        # do two term ls detrend
        pass 
    elif period <=5: 
        # do wotan, experiment with per/10 or per/20 window length 
        pass 
    else: 
        # do wotan, detrend window = 0.5 
        pass  