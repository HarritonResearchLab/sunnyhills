import os
import pandas as pd
import numpy as np
import time
tick=0

os.chdir("/Users/Ryan Axe/Documents/sunnyhills/data/current/processed/two_min_lightcurves--temp")
bar=os.listdir()
bar.remove(".gitkeep")

os.chdir('/Users/Ryan Axe/Documents/HRL/Summer')
with open('Compile.csv', 'w') as fb:
        for i in bar:
            os.chdir("/Users/Ryan Axe/Documents/sunnyhills/data/current/processed/two_min_lightcurves--temp")  
            #print(i)
            df = pd.read_csv(i)
            #print("read csv")
            df.dropna()
            #print("dropped na's")
            cleaned_time=df['cleaned_time']
            #print(len(cleaned_time))
            baz=len(cleaned_time)
            foo=str(baz)
            os.chdir('/Users/Ryan Axe/Documents/HRL/Summer')
            fb.write('\n')
            fb.write(foo+',')
            tick+=1
            print(tick)
        fb.close()
        
'''
        print("Stored column")
        os.chdir("/Users/Ryan Axe/Documents/HRL/Summer")
        print("Changed directory...")
        df_compile = pd.read_csv("Compile.csv")
        df_compile.dropna()
        print("Read Compile.csv")


        df_compile["num_obs"] = cleaned_time
        df_compile.dropna()

        print("Appended Column")
        df_compile.to_csv("Compile.csv")
'''

  
  

