def make_it():
    import os 
    import pandas as pd 

    plot_dir = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/personal_epochs/thaddaeus/july/topical/tls_gamma_run/vetting/initial/'
    tic_ids = list(set(['TIC_'+i.split('_')[2] for i in os.listdir(plot_dir) if i!= '.gitkeep']))

    filenames = [i for i in os.listdir(plot_dir) if i!='.gitkeep']

    files_to_keep = []
    SDEs = []

    for tic_id in tic_ids: 
        max_SDE = 0
        file_to_save = ''
        for file in filenames: 
            if tic_id in file: 
                SDE = float(file.split('_')[0])
                if SDE > max_SDE: 
                    max_SDE = SDE
                    file_to_save = file
        SDEs.append(max_SDE)
        files_to_keep.append(file_to_save)

    for file in os.listdir(plot_dir): 
        if file not in files_to_keep: 
            os.remove(plot_dir+file)        

    for file in files_to_keep: 
        os.rename(plot_dir+file, plot_dir+file.replace(':','_').replace('.','_'))

    df = pd.DataFrame(list(zip(SDEs, [i.replace(':','_').replace('.',',') for i in files_to_keep])), columns=['SDE','filename']) 
    df.to_csv('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/personal_epochs/thaddaeus/july/topical/tls_gamma_run/vetting/quick_key.csv', index=False)

make_it()