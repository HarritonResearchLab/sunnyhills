import os
import shutil 
from tqdm import tqdm 

base = '/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/initial/'
for f in tqdm(os.listdir(base)): 

    shutil.copyfile(base+f, base+f.replace(':','__').replace('_png','.png'))
    os.remove(base+f)