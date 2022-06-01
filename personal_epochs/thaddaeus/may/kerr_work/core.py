def merge_them(table_one:str='./personal_epochs/thaddaeus/may/kerr_work/table_one.csv',
               table_two:str='./personal_epochs/thaddaeus/may/kerr_work/table_two.csv'): 

    import pandas as pd
    import numpy as np

    df_one = pd.read_csv(table_one)
    #print(df_one.iloc[-1])
    df_two = pd.read_csv(table_two)
    
    ids_one = np.array(df_one['Gaia'])
    ids_two = np.array(df_two['Gaia'])

    _, _ , drop_indices =  np.intersect1d(ids_one, ids_two, return_indices=True) # 46 overlapping
    keep_indices = np.array(df_two.index) 
    mask = ~np.isin(keep_indices, drop_indices)
    
    df_two = df_two.iloc[mask]

    # query spectral class from gaia?

    combined = pd.concat((df_one, df_two))

    combined.to_csv(table_one.replace('table_one', 'combined'), index=False)

#merge_them()

def make_some_plots(combined=r'C:\Users\Research\Documents\GitHub\sunnyhills\personal_epochs\thaddaeus\may\kerr_work\combined.csv'):
    import pandas as pd
    import numpy as np 
    import matplotlib.pyplot as plt

    df = pd.read_csv(combined) 
    distances = np.array(df['distance_pc'])
    origins = np.array(df['table_origin'])

    fig, axs = plt.subplots(2,2, figsize=(6,6))
    plt.tight_layout()

    ax = axs[0,0]

    ax.hist(distances)
    y = 'Count'
    ax.set(xlabel='Distance (both tables)', ylabel=y)

    ax = axs[1,0]

    ax.hist(distances[np.where(origins==1)[0]])
    ax.set(xlabel='Distance (table one)', ylabel=y)

    ax = axs[1,1]
    ax.hist(distances[np.where(origins==2)[0]])
    ax.set(xlabel='Distance (table two)', ylabel=y)

    axs[0,1].axis('off')

    plt.savefig(r'C:\Users\Research\Documents\GitHub\sunnyhills\personal_epochs\thaddaeus\may\kerr_work\dist_hists.png', dpi=150)

make_some_plots()
