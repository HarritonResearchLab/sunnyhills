

def gaia_to_tic(gaia_ids):
    from astrobase.services.identifiers import gaiadr2_to_tic
    import numpy as np
    tic_ids = []
    for id in np.array(gaia_ids):
        if id !=None:
            tic_ids.append(gaiadr2_to_tic(id))
    return np.array(tic_ids)

