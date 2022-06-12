![CI](https://github.com/HarritonResearchLab/sunnyhills/actions/workflows/main.yml/badge.svg?event=push)
![python versions](https://img.shields.io/badge/python-3.7-blue&logo=Python&style=Plastic)

# sunnyhills


## Environment Setup
1. Download [miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/windows.html) (and despite the warning the installed gives, make sure to select the "Add Anaconda to my PATH enviroment variable" option) 
2. If you don't already have windows build tools installed, install them from here: ```https://visualstudio.microsoft.com/visual-cpp-build-tools/``` 
3. From Anaconda Prompt (miniconda3)<sup>1</sup>, execute ```conda create --name py37 python=3.7```
4. Execute ```conda activate py37```
5. Execute ```conda install astropy h5py beautifulsoup4 html5lib bleach pyyaml pandas sortedcontainers pytz setuptools mpmath bottleneck ipython pip```
6. Execute ```pip install -U transitleastsquares```
7. Execute ```pip install -U wotan```
8. ```cd``` into the top level ```sunnyhills``` directory on your machine within the conda prompt (for me this is ```cd Documents``` then ```cd GitHub``` then ```cd sunnyhills```)
9. Execute the command ```conda install conda-build```
10. Execute ```conda develop .```
11. Execute ```pip install -e .```
    * Note that if you ever want to *uninstall* the sunnyhills repo, you can use ```pip uninstall .``` or ```conda develop -u .```
12. To start using this environment in VSCode, go to ```view``` then ```command palette``` then ```Python: Select Interpreter``` then click ```Enter interpreter path...``` then enter the path to the ```python.exe``` file within the environment folder oft this enviroment you just made. For me, this is something like ```C:\Users\Research\miniconda3\envs\py37\python.exe``` (note that you can find the environment under ```envs``` under ```miniconda3```)

Finally, you should now be able to use everything in this perfectly working environment! 

[1] I recommend pinning the Anaconda Prompt to your taskbar once you find it (for me it's located around ```C:\Users\Research\AppData\Roaming\Microsoft\Windows\Start Menu\Programs\Anaconda3 (64-bit)\Anaconda Prompt (miniconda3).lnk```). I wouldn't use the Anaconda Powershell Prompt because I didn't test this approach on it. 

* Also note that you might have to launch VSCode from Anaconda Prompt (execute ```code```) after activating the virtual environment. 

## Test
- hello from wh1!

## Miscellaneous Notes
* **Colab Setup:** 1). ```pip3``` install ```wotan```, ```transitleastsquares```, and ```lightkurve``` 2). install fortan  via  ```sudo apt-get install gfortran``` 3). Update ```matplotlib``` version via ```pip3 install 'matplotlib==3.3.1'``` 
* **Installing the `sunnyhills` python package**
    * Not really recommended way to do things anymore, just keeping it here for storage. 
    * `$ python setup.py develop`, when executed under the directory of this repository, will install the `sunnyhills` python package to
    your computer's working path.  It'll also let you edit code in `/sunnyhills/`,
    and then import it from anywhere else on your computer.
* **Shields**: pretty nice [website](https://shields.io/category/platform-support) for shields

## About this repository
![Visualization of this repo](./diagram.svg)
