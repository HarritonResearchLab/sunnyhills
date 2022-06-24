def find_lc_timegroups(lctimes, mingap=4.0):
    '''Finds gaps in the provided time-series and indexes them into groups.
    This finds the gaps in the provided `lctimes` array, so we can figure out
    which times are for consecutive observations and which represent gaps
    between seasons or observing eras.
    Parameters
    ----------
    lctimes : array-like
        This contains the times to analyze for gaps; assumed to be some form of
        Julian date.
    mingap : float
        This defines how much the difference between consecutive measurements is
        allowed to be to consider them as parts of different timegroups. By
        default it is set to 4.0 days.
    Returns
    -------
    tuple
        A tuple of the form: `(ngroups, [slice(start_ind_1, end_ind_1), ...])`
        is returned.  This contains the number of groups as the first element,
        and a list of Python `slice` objects for each time-group found. These
        can be used directly to index into the array of times to quickly get
        measurements associated with each group.
    '''

    lc_time_diffs = np.diff(lctimes)
    group_start_indices = np.where(lc_time_diffs > mingap)[0]

    if len(group_start_indices) > 0:

        group_indices = []

        for i, gindex in enumerate(group_start_indices):

            if i == 0:
                group_indices.append(slice(0,gindex+1))
            else:
                group_indices.append(slice(group_start_indices[i-1]+1,gindex+1))

        # at the end, add the slice for the last group to the end of the times
        # array
        group_indices.append(slice(group_start_indices[-1]+1,len(lctimes)))

    # if there's no large gap in the LC, then there's only one group to worry
    # about
    else:
        group_indices = [slice(0,len(lctimes))]

    return len(group_indices), group_indices



import os
import pandas as pd
import numpy as np
outputdf=pd.DataFrame()
#Append these lists as new columns to dataframe
TIC_ID=[]
START_DATE=[]
END_DATE=[]
NUM_OBS=[]
counter=0
os.chdir("/Users/Ryan Axe/Documents/sunnyhills/data/current/processed/two_min_lightcurves")
bar=os.listdir()
bar.remove(".gitkeep")
os.chdir("/Users/Ryan Axe/Documents/HRL/Summer/Week2")

with open('Temp.csv', 'w') as fb:
    for i in bar:
        os.chdir("/Users/Ryan Axe/Documents/sunnyhills/data/current/processed/two_min_lightcurves")
        df = pd.read_csv(i)
        df=df.dropna()
        ct=df["clean_time"]
        array1=np.array(ct)
        foo=find_lc_timegroups(ct)
        groupnum=foo[0]
        for j in range(groupnum):
            #refering to slice as section because its referenced in function which causes error
            section=foo[1][j]
            #slice:
            x=[]
            x=array1[section]
            y=len(x)
            #append y to dataframe
            #append TICid
            #append min and max using x[0] and x[y-1]


            TIC_ID.append(i)
            END_DATE.append(x[y-1])
            NUM_OBS.append(y)
            #x is indexed based on the splice range. This means for a lot of x cases there is no x[0]
            START_DATE.append(x[0])


        os.chdir("/Users/Ryan Axe/Documents/HRL/Summer/Week2")
#        fb.write('\n')
#        fb.write(i+str(foo))
    outputdf['TIC_ID']=TIC_ID
    outputdf['START_DATE']=START_DATE
    outputdf['END_DATE']=END_DATE
    outputdf['NUM_OBS']=NUM_OBS



    os.chdir("/Users/Ryan Axe/Documents/HRL/Summer/Week2")
    fb.close()
    outputdf.to_csv("final.csv", encoding='utf-8', index=False)
