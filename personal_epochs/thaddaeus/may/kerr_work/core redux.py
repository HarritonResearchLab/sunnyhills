def figure_coords_out(table:str=''):
    import pandas as pd 
    import numpy as np
    import matplotlib.pyplot as plt 
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    
    df = pd.read_csv(table)
    ras = np.array(df['RAdeg'])
    decs = np.array(df['DEdeg'])

    distances = 1000/np.array(df['plx'])

    ra = ras[0]
    dec = decs[0]
    c = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    print(c, c.galactic)

    c = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='icrs')
    gc = c.transform_to('galactic') #fk5 didn't work...
    
    gc_x = np.array([i.value for i in gc.l])

    gc_y = np.array([i.value for i in gc.b])

    #x = distances*np.cos(ras)*np.cos(decs)
    #y = distances*np.sin(ras)*np.cos(decs)

    #print(ras[0:10], x[0:10])

    # X=Dcoslcosb, Y=Dsinlcosb, Z=Dsinb --> this works! so cool!
    # https://galaxiesbook.org/chapters/A.-Coordinate-systems.html
    
    # https://iopscience.iop.org/article/10.3847/1538-4357/ac0251#apjac0251f7
    # https://docs.astropy.org/en/stable/coordinates/transforming.html
    # https://docs.astropy.org/en/stable/coordinates/index.html 

    plt.rcParams['font.family']='serif'

    fig, axs = plt.subplots(1,2)

    ax = axs[0]

    ax.scatter(ras, decs, c='black', s=1, label='Original ICRS')
    ax.legend()
    ax = axs[1]

    ax.scatter(gc_x, gc_y, c='black', s=1, label='Galactic (not galacteocentric)')
    
    ax.legend()

    plt.show()

figure_coords_out(r'C:\Users\Research\Documents\GitHub\sunnyhills\personal_epochs\thaddaeus\may\kerr_work\table_one.csv')