import pandas as pd
import os
tick=1
match=[]
df = pd.read_csv('List.csv')
foo=df['TIC']

os.chdir('/Users/Ryan Axe/Documents/sunnyhills/data/current')
df2 = pd.read_csv('current_key.csv')
bar = df2['tid']
#print(bar)

#print(bar.pop(553))

for i in foo:
    #print('#',i)
    
    tick=0
    #tick was set to 2 to match the csv file but it needs to be set to 0 to match the ra and dec lists.
    #print(tick)
    for x in bar:

        if(i==x):
            #print('##',tick)
            match.append(tick)
        else:
            tick+=1
#print(match)

ra = df2['GDR2_RA']
dec = df2['GDR2_DEC']

for x in match:
    ra.pop(x)

for z in match:
    print(dec.pop(z))



