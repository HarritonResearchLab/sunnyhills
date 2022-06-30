from pdb import _rstr
from re import T


def tls_even_odd(tls_results): 
    
    '''
    Description 
    -----------
    returns true if mean odd and mean even transit depths are inconsistent at three sigma level 

    '''

    import numpy as np

    eb_flag = True 

    odd = tls_results.depth_mean_even 

    odd = [odd[0]-3*odd[1], odd[0]+3*odd[1]]

    even = tls_results.depth_mean_odd
    even = [even[0]-3*even[1], even[0]+3*even[1]]

    idx = np.argsort([even[0], odd[0]])

    temp = np.array([even, odd])[idx]

    if temp[0][1]>=temp[1][0]: 
        eb_flag = False

    even = [round(i, 3) for i in even]
    odd = [round(i, 3) for i in odd]

    return {'Even-Odd Flag':eb_flag, 'even':even, 'odd':odd}

def transit_outliers_fap_test(tls_results): 
    '''
    Description
    -----------
    See section 3.3 from edi-vetter paper (jon zink et al. 2020)
    '''
    
    import numpy as np

    folded_phase = tls_results.folded_phase 
    folded_y = tls_results.folded_y
    
    folded_model_flux = tls_results.model_folded_model
    
    half_dur = tls_results.duration 
    transit_mask = np.logical_and(folded_phase<=0.5+half_dur, folded_phase>=0.5-half_dur)
    
    subtracted = folded_y-folded_model_flux 

    lc_sigma = np.std(subtracted[~transit_mask])
    transit_sigma = np.std(subtracted[transit_mask]) 

    return {'lc_sigma':round(lc_sigma, 3), 'transit_sigma':round(transit_sigma, 3)} 

def check_lombscargle(tic_id, tls_results, download_log): 

    import pandas as pd 
    import numpy as np

    if isinstance(download_log, pd.DataFrame): 
        df = download_log  
    else: 
        df = pd.read_csv(download_log)

    index = np.where(df['TIC_ID']==tic_id)[0]
    top_ls_period = np.array(df['top_ls_period'])[index]
    if len(top_ls_period)>0: 
        top_ls_period = top_ls_period[0]   

    ls_flag = False
    if 0.99<tls_results.period/top_ls_period<1.01: 
        ls_flag = True 

    return {'ls_top_period_test':ls_flag, 'ls_top_period':top_ls_period}

def generate_cutout(tic_id:str='',large_size:int=30,small_size:int=5):
    import matplotlib.pyplot as plt
    import lightkurve as lk
    from matplotlib.patches import Rectangle
    import os
    import eleanor
    from PyPDF2 import PdfMerger

    os.mkdir('temp')

    data = eleanor.Source(tic=int(tic_id.replace('TIC', ''))) 
    data = eleanor.TargetData(data)
    vis = eleanor.Visualize(data)

    result = lk.search_tesscut(tic_id)
    small_tpf = result.download(cutout_size=small_size)
    large_tpf = result.download(cutout_size=large_size)
    fig, ax = plt.subplots()
    small_tpf.plot(ax=ax)
    rect_x = ax.get_xlim()[1]-ax.get_xlim()[0]
    rect_y = ax.get_ylim()[1]-ax.get_ylim()[0]
    plt.savefig('temp/small.pdf')
    fig1, ax1 = plt.subplots()
    large_tpf.plot(ax=ax1)
    ax1 = vis.plot_gaia_overlay(int(tic_id.replace('TIC ','')),large_tpf)
    ax1.add_patch(Rectangle((ax.get_xlim()[0],ax.get_ylim()[0]),rect_x,rect_y,fc='none',linewidth=1,ec='yellow',ls='--'))
    plt.savefig('temp/large.pdf')

    merger = PdfMerger()

    for pdf in ['temp/small.pdf','temp/large.pdf']:
      merger.append(pdf)
    merger.write(tic_id+'.pdf')
    merger.close()
    os.remove('temp/small.pdf')
    os.remove('temp/large.pdf')
    os.rmdir('temp/')

def centroid_offset_test(tic_id:str, per:float, t0:float, dur:float, plot_path:str=None): 
    import vetting 
    from vetting import centroid_test
    import lightkurve as lk

    #offset_detected_list, pvalues, fig = centroid_offset_test(tic_id='TIC 13023738', per=8.0635, t0=2459104.87232, dur=3.396/24)

    search_result = lk.search_targetpixelfile(tic_id.replace('_', ' '))

    tpfs = []

    if len(search_result)>0: 
        for i in range(len(search_result)): 
            tpf = search_result[i].download(quality_bitmask='hard')
            tpfs.append(tpf)
        length = len(tpfs)
        pers = length*[per]
        t0s = length*[t0]
        durs = length*[dur]
        
        r = centroid_test(tpf, pers, t0s, durs, aperture_mask='pipeline', plot=return_plots)

        offset_detected_list = r['centroid_offset_detected'][0]
        fig = r['figs'][0]
        pvalues = r['pvalues'][0]

        if plot_path is not None: 
            fig.savefig(plot_path, dpi=150)

        return offset_detected_list, pvalues, fig 

    else: 
        return None, None, None  