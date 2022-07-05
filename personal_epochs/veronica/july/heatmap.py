def heatmap(path_to_csv:str,plot_dir:str = None):
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

    injection_results = pd.read_csv(path_to_csv)
    recovery_proportion = np.divide(np.array(injection_results['injected_period']),np.array(injection_results['recovered_period']))
    injection_results = injection_results.join(pd.Series(recovery_proportion,name='recovery_proportion'))
    fig,ax = plt.subplots()

    ax =sns.heatmap(injection_results.pivot('R_P','rotation_period','recovery_proportion'),annot=True,cmap="YlGnBu")
    ax.invert_yaxis()
    ax.set_title('TLS Recovery Percentages')
    ax.set_xlabel('Stellar Rotation Period (s)')
    ax.set_ylabel('Planetary Radius ($R_{{{Jupiter}}}$)')

    if plot_dir!=None:
        if plot_dir[-1]!="/":
            plot_dir+="/"
    fig.savefig(plot_dir +"heatmap.png")

def heatmap_fn(path_to_csv:str,plot_dir:str=None):
    import pandas as pd
    import matplotlib.pyplot as plt
    import numpy as np
    import seaborn as sns

    injection_results = pd.read_csv(path_to_csv)
    recovery_pctg = []

    for injected_periods, recovered_periods in zip(injection_results['injected_period'],injection_results['recovered_period']):
        injected = injected_periods.replace("\"","").split(',') 
        recovered = recovered_periods.replace("\"","").split(',')
        injected = list(map(float,injected))
        recovered = list(map(float,recovered))
        error = np.abs(np.divide(np.array(injected),np.array(recovered)) -1)
        recovery_pctg.append(len(np.where(error<=.01)[0])/len(error))
        
    injection_results = injection_results.join(pd.Series(recovery_pctg,name='recovery_pctg'))

    fig,ax = plt.subplots()

    ax =sns.heatmap(injection_results.pivot('R_P','rotation_period','recovery_pctg'),annot=True,cmap="YlGnBu")
    ax.invert_yaxis()
    ax.set_title('TLS Recovery Percentages')
    ax.set_xlabel('Stellar Rotation Period (s)')
    ax.set_ylabel('Planetary Radius ($R_{{{Jupiter}}}$)')

    if plot_dir!=None:
        if plot_dir[-1]!="/":
            plot_dir+="/"
    fig.savefig(plot_dir +"heatmap.png")
heatmap('C:/Users/60002/Documents/GitHub/sunnyhills/personal_epochs/veronica/july/test_data.csv','C:/Users/60002/Documents/GitHub/sunnyhills/personal_epochs/veronica/july')