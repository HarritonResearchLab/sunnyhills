def make_plot(table_two:str=r'C:\Users\Research\Documents\GitHub\sunnyhills\personal_epochs\thaddaeus\may\kerr_work\table_two.csv'): 
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt

    table_two = pd.read_csv(table_two)

    ras = np.array(table_two['RAdeg'])
    decs = np.array(table_two['DEdeg'])

    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import numpy as np 

    gc = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='fk5').transform_to('galactic')
    
    ras = np.array([i.value for i in gc.l])
    decs = np.array([i.value for i in gc.b])
    
    # https://iopscience.iop.org/article/10.3847/1538-4357/ac0251#apjac0251f7
    # https://docs.astropy.org/en/stable/coordinates/transforming.html
    # https://docs.astropy.org/en/stable/coordinates/index.html 

    fig, ax = plt.subplots(figsize=(5, 5))

    ax.scatter(ras, decs, c='black', s=1)

    lims = [-370,370]

    #ax.set(xlim=lims, ylim=lims, xlabel='X', ylabel='Y')

    plt.show()

make_plot()

def test_coords(): 
    import astropy.units as u
    from astropy.coordinates import SkyCoord
    import numpy as np 

    gc = SkyCoord(ra=np.array([309.8679, 309.8679])*u.degree, dec=np.array([79.892, 79.892])*u.degree, frame='fk5').transform_to('galactic')
    print(gc.l.value, gc.b.value)
    #print(gc.)

#test_coords()