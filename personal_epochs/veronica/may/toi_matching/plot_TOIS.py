def merge_keys():
    import pandas as pd
    #TODO: change paths to relative paths
    curr = pd.read_csv('C:/Users/60002/Documents/GitHub/sunnyhills/data/current/lightcurve_key.csv')
    toi = pd.read_csv('C:/Users/60002/Documents/GitHub/sunnyhills/personal_epochs/veronica/may/toi_matching/toi_export.csv')

    merged = curr.merge(toi)

    merged.to_csv('merged_key.csv')

def addAnnotations(x,y):
    import matplotlib.pyplot as plt
    for pos in range(len(x)):
        plt.text(pos,y[pos],y[pos],ha='center')

def plot_toi_percentages():
    import matplotlib.pyplot as plt
    import pandas as pd
    toi_keys = {
        'APC':0,
        'CP':0,
        'FA':0,
        'FP':0,
        'KP':0,
        'PC':0
    }

    old_df = pd.read_csv('C:/Users/60002/Documents/GitHub/sunnyhills/personal_epochs/veronica/may/toi_matching/toi_export.csv')
    df = old_df.dropna()
    planetary_keys =['ambiguous planetary candidate', 'Confirmed Planet','False Alarm','False Positive','Known Planet','Planetary Candidate']
    for toi in df['tfopwg_disp']:
        toi_keys[str(toi)]+=1

    keys = list(toi_keys.keys())
    vals = list(toi_keys.values())

    for key,val in zip(keys,vals):
        plt.bar(key,val)
    plt.legend(planetary_keys,title='keys')
    plt.xlabel('TOI Key')
    plt.ylabel('Count')
    plt.title('TOI Dispositions')
    addAnnotations(keys,vals)
    plt.show()
    
    return (len(old_df['tfopwg_disp'])-len(df['tfopwg_disp']))/len(old_df['tfopwg_disp'])*100
    

count = plot_toi_percentages()
print(count)
    