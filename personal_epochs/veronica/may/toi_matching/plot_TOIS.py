def merge_keys():
    import pandas as pd
    #TODO: change paths to relative paths
    curr = pd.read_csv('C:/Users/60002/Documents/GitHub/sunnyhills/data/current/lightcurve_key.csv')
    toi = pd.read_csv('C:/Users/60002/Documents/GitHub/sunnyhills/personal_epochs/veronica/may/toi_matching/toi_export.csv')

    merged = curr.merge(toi)

    merged.to_csv('merged_key.csv')

def plot_toi_percentages():
    import matplotlib.pyplot as plt
    import pandas as pd
    toi_keys = {
        'APC':0,
        'CP':0,
        'FA':0,
        'KP':0,
        'PC':0
    }

    df = pd.read_csv('merged_key.csv')
    df = df.dropna()
    for toi in df['tfopwg_disp']:
        toi_keys[str(toi)]+=1

    keys = list(toi_keys.keys())
    vals = list(toi_keys.values())

    plt.bar(keys,vals)
    plt.show()
    


    