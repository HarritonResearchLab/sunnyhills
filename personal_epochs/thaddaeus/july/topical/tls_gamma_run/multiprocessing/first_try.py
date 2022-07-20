import numpy as np
import multiprocessing as mp 
import time 
import contextlib

from numpy.linalg import inv 

x = np.random.rand(72,72,72,72)

arr = [x for i in range(8)]

def first():

    a = None
    b = None 

    start = time.time()
    with mp.Pool(processes=32) as pool: 
        a = pool.map(inv, arr)
        print(time.time()-start)
    
first()

t = time.time()
for i in arr: 
    y = np.linalg.inv(i)

print(time.time()-t)