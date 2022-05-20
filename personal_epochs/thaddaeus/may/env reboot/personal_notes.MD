# personal notes

## setting up wsl env

* delete dir with files in it: ```rm -r dir_name```

* **making two seperate envs***
    * wotan land: ```pip install numpy==1.21.6 wotan``` then ```pip install pandas```
        * `ERROR: numba 0.55.1 has requirement numpy<1.22,>=1.18, but you'll have numpy 1.22.3 which is incompatible.`

    * tls land: first ```pip install numpy==1.21.6``` then ```pip install numpy pandas matplotlib seaborn scipy astropy lightkurve transitleastsquares```

    * pip install --force-reinstall numpy==1.19.3
        * tls import was messed up earlier
   
