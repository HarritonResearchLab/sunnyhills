def as_subprocess(file_directory:str, file_name:str, verbose=False): 
    import subprocess
    import os 

    os.chdir(file_directory)

    if '.py' not in file_name: 
        file_name+='.py'

    if verbose: 
        subprocess.Popen(["python",file_name])

    else: 
        os.system('python '+file_name+' >/dev/null 2>&1 &')

