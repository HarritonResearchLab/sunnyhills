def get_append_lines():
    header = []
    with open('/ar1/PROJ/fjuhsd/shared/github/sunnyhills/sunnyhills/pipeline/alpha.py', 'r') as f: 
        for line in f: 
            if '.append' in line: 
                col = line.replace('\n','').strip().split('.append(')[0]
                if col[-1].upper()=='S':
                    col = col[:-1]
                    header.append(col)
    print(len(header))
    print('['+','.join(header)+']')
        

get_append_lines()