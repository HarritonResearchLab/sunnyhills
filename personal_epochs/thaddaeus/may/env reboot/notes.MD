# notes

* install wsl. Do all code execution from wsl environment (e.g. cd into directories with the given files you are editing in vs code locally on your computer, and execute the command ```python filename.py```, where ```filename.py``` is the file you want to run)

## summary (read this first)

1. make sure ```sunnyhills``` repo is cloned locally on your computer
2. install windows terminal 
3. install wsl on your computer
4. link windows terminal with wsl 
5. open a new wsl session in windows terminal
6. execute the following command: ```pip install numpy pandas matplotlib seaborn sklearn scipy astropy transitleastsquares wotan, altaipony```
7. within wsl terminal session, navigate to 

## usage notes

* Whenever you want to edit your files, edit them in vscode *outside* wsl, and then run them within wsl by navigating to their directory within wsl and then running them like ```python file_name.py``` (this would run the file ```file_name.py``` in the current directory)
    * Note that vscode will probably go crazy with you not having installed the modules you are working with, this is because it does not initially look for the wsl interpreter. this is okay. It will run fine when you run in wsl, as shown as above. 
    * Note that with wsl you won't be able to use ```plt.show()```, because wsl doesn't support GUI's...instead ```plt.savefig()``` everything and open them up within your normal windows environment. 
* whenever you want to change directories in wsl terminal, use ```cd new_dir``` command, e.g. if there is a folder in my current directory called ```Downloads```, I could cd into it via ```cd Downloads```
    * handy ```ls``` command in terminal lists all files/folders in the current directory. 

## wsl notes

* install WSL2, and makes sure you have hyper-v on (I'm pretty sure you can select this by completely restarting your computer...microsoft will prompt you to choose...otherwise, it's probably a simple thing you can look up how to do on the internet)

## trouble shooting

* this is a little fuzzy to me, but you might need to install the visual studio (note *not* vscode) c/c++ build tools etc. outside of wsl...you can do this at this [link](https://visualstudio.microsoft.com/vs/features/cplusplus/)

* this is also a little fuzzy to me, but if any of the modules are giving you issues, install fortran at the command prompt in wsl with the command ```sudo apt-get install gfortran``` (sudo so you'll have to enter your wsl password)