def run_test(): 
    import time 
    baps = []
    for i in range(32): 
        baps.append(i)
        print(i)
        time.sleep(0.25)

    with open('personal_epochs/thaddaeus/june/topical/piping/out.txt', 'w') as f: 
        f.write(' '.join([str(i) for i in baps])) 

run_test()

