import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
df = pd.read_csv('FINAL.csv')
foo = np.array(df['NUM_OBS'])
bar=np.log10(foo)
baz=np.array(bar)
#print(bar)
x = np.array(df['GDR2_RA'])
y = np.array(df['GDR2_DEC'])

plt.scatter(x, y, c=baz)
plt.show()
