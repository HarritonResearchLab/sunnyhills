def collect_plots(): 

    import os
    import shutil 
    import pandas as pd 
    from tqdm import tqdm 

    plots_to_copy = list(pd.read_csv('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/for_email/plot_names.txt')['filename'])
    base = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/initial/'
    new_dir = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/for_email/'
    for plot in plots_to_copy: 
        shutil.copyfile(base+plot, new_dir+plot.replace('SDE__SDE__SDE__','SDE__'))

collect_plots()