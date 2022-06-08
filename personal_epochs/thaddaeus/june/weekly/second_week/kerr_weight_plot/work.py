def make_plot(table_two:str=r'C:\Users\Research\Documents\GitHub\sunnyhills\personal_epochs\thaddaeus\may\kerr_work\table_two.csv'): 
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.units as u
    from astropy.coordinates import SkyCoord

    table_two = pd.read_csv(table_two)

    weights = np.array(table_two['Weight'])

    percents = [0,90,99,99.9]
    percentiles = np.percentile(weights, percents)

    ras = np.array(table_two['RAdeg'])
    decs = np.array(table_two['DEdeg'])

    distances = 1000/np.array(table_two['plx'])

    gc = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='fk5').transform_to('galactic')
    
    ras = np.array([i.value for i in gc.l])
    decs = np.array([i.value for i in gc.b])

    x = distances*np.cos(ras)*np.cos(decs)
    y = distances*np.sin(ras)*np.cos(decs)

    # X=Dcoslcosb, Y=Dsinlcosb, Z=Dsinb --> this works! so cool!
    # https://galaxiesbook.org/chapters/A.-Coordinate-systems.html
    
    # https://iopscience.iop.org/article/10.3847/1538-4357/ac0251#apjac0251f7
    # https://docs.astropy.org/en/stable/coordinates/transforming.html
    # https://docs.astropy.org/en/stable/coordinates/index.html 

    plt.rcParams['font.family']='serif'

    fig, axs = plt.subplots(2,2,figsize=(7, 7))

    pos = [[0,0], [0,1], [1,0], [1,1]]

    lims = [-425,425]
    for index, pair in enumerate(pos):
        mask = weights>percentiles[index]
        label = 'Percentile: '+str(percents[index])+'\nNum. Obs. '+str(np.sum(mask))
        axs[pair[0],pair[1]].scatter(x[mask],y[mask], c='black', s=1, label=label)
        axs[pair[0],pair[1]].set(xlim=lims, ylim=lims, xlabel='X', ylabel='Y')
        axs[pair[0],pair[1]].legend(fancybox=False, loc='lower right', fontsize='x-small', edgecolor='black')
    
    plt.show()

    table_one = pd.read_csv(r'C:\Users\Research\Documents\GitHub\sunnyhills\personal_epochs\thaddaeus\may\kerr_work\table_one.csv')

    ras = np.array(table_one['RAdeg'])
    decs = np.array(table_one['DEdeg'])

    distances = 1000/np.array(table_one['plx'])

    #gc = SkyCoord(ra=ras*u.degree, dec=decs*u.degree, frame='fk5').transform_to('galactic')
    
    #ras = np.array([i.value for i in gc.l])
    #decs = np.array([i.value for i in gc.b])

    x = distances*np.cos(ras)*np.cos(decs)
    y = distances*np.sin(ras)*np.cos(decs)

    mask = np.array(table_one['P'])>0.5

    plt.rcParams['font.family']='serif'
    fig, ax = plt.subplots(figsize=(4,4))

    ax.scatter(x[mask],y[mask], c='black', s=0.25, label='P>0.5 Num. Obs. = '+str(np.sum(mask)))
    ax.set(xlim=lims, ylim=lims, xlabel='X', ylabel='Y')
    ax.legend(fancybox=False, loc='lower right', fontsize='x-small', edgecolor='black')

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