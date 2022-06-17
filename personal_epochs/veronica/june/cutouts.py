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
