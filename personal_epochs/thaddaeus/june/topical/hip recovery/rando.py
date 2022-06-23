tic_id = 'TIC_166527623'
base_dir = './personal_epochs/thaddaeus/june/topical/hip recovery/'
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt

def download(): 
    from sunnyhills.pipeline_functions import download_and_preprocess

    download_and_preprocess(ticstr='TIC 166527623', outdir='./personal_epochs/thaddaeus/june/topical/hip recovery/')

#download()

def first_play(): 
    from sunnyhills.pipeline_functions import run_tls
    from sunnyhills.plotting import tls_validation_mosaic
    import pandas as pd
    import numpy as np

    data_path = './personal_epochs/thaddaeus/june/topical/hip recovery/TIC_166527623.csv'
    data = pd.read_csv(data_path)
    mask = np.isfinite(data['clean_flux'])

    time, flux = (np.array(data[i])[mask] for i in ['clean_time', 'clean_flux'])

    second_mask = time>1700
    time, flux = (i[second_mask] for i in [time, flux])

    tls_results, tls_model = run_tls('TIC_166527623', time, flux, show_progress_bar=True)
    tls_validation_mosaic('second_half_TIC_166527623', data_path, tls_model, tls_results, plot_dir='./personal_epochs/thaddaeus/june/topical/hip recovery/')

first_play()

def first_plot(): 
    data_path = './personal_epochs/thaddaeus/june/topical/hip recovery/TIC_166527623.csv'
    data = pd.read_csv(data_path)
    mask = np.isfinite(data['clean_flux'])

    time, flux = (np.array(data[i])[mask] for i in ['clean_time', 'clean_flux'])    
    second_mask = time>1700
    time, flux = (i[second_mask] for i in [time, flux])

    plt.scatter(time, flux, s=1)

    plt.savefig(base_dir+'second_half.png', dpi=150)

#first_plot()