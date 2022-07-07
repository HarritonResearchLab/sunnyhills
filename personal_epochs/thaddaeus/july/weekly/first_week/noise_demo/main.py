import numpy as np
import matplotlib.pyplot as plt 

x = np.linspace(0,1,500)
rms = [0.01,0.02,0.03,0.04,0.05]

fig, axs = plt.subplots(5,1, figsize=(4,15))

for i in range(5): 
    ax = axs[i]
    ax.scatter(x, 1+np.random.normal(0,rms[i], size=500), s=1, color='cornflowerblue')
    ax.set(xlabel='Time',ylabel='Flux (rms='+str(rms[i])+')')

plt.tight_layout()

plt.savefig('personal_epochs/thaddaeus/july/weekly/first_week/noise_demo/noise_demo.png', dpi=150)