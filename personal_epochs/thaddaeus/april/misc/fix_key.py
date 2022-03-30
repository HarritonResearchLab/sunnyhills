import pandas as pd

df = pd.read_csv('./data/current/current_key_with_cadences.csv')

has_two_min = []
counter = 0 
for status in df['cadences']: 
    
    if '120' in status: 
        has_two_min.append(True)
        counter+=1
    else: 
        has_two_min.append(False)

print(counter)

df['has_two_min'] = has_two_min

df.to_csv('temp.csv', index=False)