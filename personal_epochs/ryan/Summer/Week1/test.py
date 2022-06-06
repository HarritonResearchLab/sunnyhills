import random
with open('Compile.csv', 'w') as fp:
    for i in range(5):
        foo=random.randint(1,10)
        bar=str(foo)
        fp.write('\n')
        fp.write(bar)
fp.close()
