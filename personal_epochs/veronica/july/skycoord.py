import pandas as pd
import numpy as np

test_data = pd.read_csv('/mnt/c/users/60002/Documents/Github/sunnyhills/data/current/current_key.csv')
ra = test_data['GDR2_RA']
dec = test_data['GDR2_DEC']

#plt.figure()
def skyplot(ra:np.array,dec:np.array,style:str='points'):
    import matplotlib.pyplot as plt
    from matplotlib import projections
    import numpy as np
    from scipy.stats import gaussian_kde

    plt.figure(figsize=(8,4.2))
    plt.subplot(111,projection='aitoff')
    if style == 'points':
        plt.grid(True)
        plt.plot(ra,dec,'o',markersize=2,alpha=.95,color='#6495ED')
    elif style == 'density':
        plt.grid(True)
        points = np.vstack([ra,dec])
        point_density = gaussian_kde(points)(points)
        plt.scatter(ra,dec,c=point_density,s=10,cmap='Blues')



    plt.title('Skyplot Map')
    plt.show()
    plt.savefig('test.png')
skyplot(np.array(ra),np.array(dec),'density')
