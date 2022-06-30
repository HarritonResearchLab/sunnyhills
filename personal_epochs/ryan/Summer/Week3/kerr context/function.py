def function(gaiaid: str):
    import pandas as pd
    import numpy as np
    names=['Cepheus Flare','Pleiades','Taurus-Orion','Ophiuchus Southeast','Fornax-Horologium','CMa North','Aquila East','Cepheus Far North','Vela-CG7','ASCC 123','Cepheus-Cygnus','Lyra','Cerberus','Carina-Musca','Perseus','Perseus','Taurus-Orion II','Greater Taurus','IC 2391	101','NGC 2451A','Chamaeleon','Sco-Cen	7394','Taurus-Orion III','Vela-CG4','Taurus-Orion IV','Monoceros Southwest','Greater Orion']
    #note for some reason when making this csv with excel pandas returned a dtype warning. To fix this just specify the dtype as unicode
    df1=pd.read_csv('table1.csv',dtype='unicode')
    df2=pd.read_csv('table2.csv',dtype='unicode')

    id_table1=np.array(df1['GAIA ID'])
    EOM_table1=np.array(df1['EOM'])
    TLC_table1=np.array(df1['TLC'])
    id_table2=np.array(df2['GAIA ID'])
    EOM_table2=np.array(df2['EOM'])
    TLC_table2=np.array(df2['TLC'])
    breakcheck='start'
    stellar_name=0
    x=0
#use str() input again to make sure gaiaid is infact a string. Might be redundant...
    match_id1=np.where(str(gaiaid)==id_table1)[0]
    match_id2=np.where(str(gaiaid)==id_table2)[0]
#Use .any() to remove python warnings... Function still produces accurate results.
        
        
    if(EOM_table2[match_id2].any()=='-1' and TLC_table2[match_id2].any()!='-1'):
        x=TLC_table2[match_id2]
        stellar_name=names[int(x)]
        breakcheck='if_loop1'
    elif(EOM_table1[match_id1].any()=='-1' and TLC_table1[match_id1].any()!='-1'):
        x=TLC_table1[match_id1]
        stellar_name=names[int(x)]
        breakcheck='elif_loop1'
###########################################################################################################
    elif(EOM_table2[match_id2].any()=='-1' and TLC_table1[match_id1].any()=='-1'):
        stellar_name='Field Star'
        breakcheck='elif_loop2'
    elif(EOM_table1[match_id1].any()=='-1' and TLC_table1[match_id1].any()=='-1'):
        stellar_name='Field Star'
        breakcheck='elif_loop3'
############################################################################################################
    else:
        stellar_name='NOT FOUND'
        breakcheck='else_loop1'
    print(stellar_name)
#can return and print breakcheck for debugging reasons...
    return stellar_name, breakcheck
