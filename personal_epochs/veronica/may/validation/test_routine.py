def bls_check():
    import pandas as pd

    df = pd.read_csv('C:/Users/60002/Documents/GitHub/sunnyhills/personal_epochs/veronica/may/validation/tess nasa exo archive confirmed.csv')

    df = df.dropna(subset=['pl_orbper','sy_pnum','pl_controv_flag'])

    #adjusted = df.mask(df['pl_orbper']>.5) #and df['pl_orbper']<15)
    convert_to_num(df['pl_orbper'])


def convert_to_num(key):
    new_keys = []
    for i in key:
        i = i.strip('s') 
        new_keys.append(int(i))
    print(new_keys)
    
bls_check()