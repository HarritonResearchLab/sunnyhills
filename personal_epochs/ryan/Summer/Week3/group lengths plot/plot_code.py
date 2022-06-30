import pandas as pd
import numpy as np
import seaborn as sea

df = pd.read_csv('https://raw.githubusercontent.com/HarritonResearchLab/sunnyhills/main/personal_epochs/ryan/Summer/Week2/final.csv?token=GHSAT0AAAAAABWBVZ7BGXXFIXLO6PQVL65MYV2C6WQ')
array = np.array(df['NUM_OBS'])
x=[90,95,99,99.9,99.99]
percentiles=np.array(x)
print(np.percentile(array, percentiles))
print("***********************")
print(np.std(array))
sea.kdeplot(data=array)
