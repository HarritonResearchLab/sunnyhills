import pandas as pd
import numpy as np
import os
os.remove("Export.csv")
input("please confirm file was deleted...")
#np.unique() to remove duplicate rows?
df = pd.read_csv('https://raw.githubusercontent.com/HarritonResearchLab/sunnyhills/main/personal_epochs/veronica/may/validation/tess%20nasa%20exo%20archive%20confirmed.csv')
df.dropna()
df.pop('loc_rowid')
df.pop('sy_pnum')
df.pop('disc_year')
df.pop('pl_controv_flag')
df.pop('pl_orbper')
df.pop('pl_radj')
heading = ['pl_name','hostname','sy_snum','st_spectype','rastr','ra','decstr','dec','sy_dist','sy_vmag','sy_kmag','sy_gaiamag','sy_tmag','releasedate']


for i in heading:
    df.drop_duplicates(subset=[i])


df.to_csv('Export.csv')

#All of this code down here was a waste of time because of df.drop_doplicates. It also didnt work...
'''
for i in heading:
  #Reads the specific column associated with the given header using df[]. Then an array of that column is stored to foo for further manipulation
  foo=np.array(df[i])
  #set is the same as np.unique() but it works with strings. Set() removes duplicates in a list
  bar=set(foo)
  baz=np.array(bar)
  #print(bar)
  df[i]=pd.Series(baz)
#Write each new column to a new csv file
#lol could have used df.drop_duplicates([subset])
df.dropna()
df.to_csv('Export.csv')
'''
