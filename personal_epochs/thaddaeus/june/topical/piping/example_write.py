def foo(): 
    for i in range(222):
        with open('personal_epochs/thaddaeus/june/topical/piping/out.txt', 'a') as f: 
            f.write(str(i)+'\n')

foo()