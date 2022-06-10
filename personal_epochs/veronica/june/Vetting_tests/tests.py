# adapted from: https://arxiv.org/pdf/2001.11515.pdf

def tran_dur(bls_params):
    '''
    Checks if the Transit Duration is too long when compared to the period
    https://arxiv.org/pdf/2001.11515.pdf (see section 3.11)
    '''
    trans_dur_false_per = False
    
    fit_P=bls_params[1]
    fit_t0=bls_params[2]
    fit_tdur=bls_params[3]
    
    if fit_P<.5:
        trans_dur_false_per=True
    elif fit_tdur/fit_P>.1:
        trans_dur_false_per=True
    else:
        trans_dur_false_per=False
    bls_params.append(trans_dur_false_per)    

    return bls_params    
                   