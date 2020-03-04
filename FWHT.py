# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 15:21:43 2015
Fast Walsh-Hadamard Transform with Sequency Order
Author: Ding Luo@Fraunhofer IOSB
"""
from math import log
import numpy as np
import GrayCode
from time import clock

def get_sequency_list(inputArray):
    """ Sort input 1D array into sequency order
    Utilizes gray code generation from a Python recipe from Internet.
    """
    length = inputArray.size
    bitlength = int(log(length,2))
    # Gray Code
    graycodes=GrayCode.GrayCode(bitlength)
    # Bitreverse of gray code
    bitreverse = [int(graycodes[i][::-1],2) for i in xrange(length)]
    
    outputArray = inputArray.copy()
    outputArray[bitreverse] = inputArray[:]

    return outputArray


def SFWHT(X):
    """ 'Slower' Fast Walsh-Hadamard Transform
    Step#1 Get sequency-ordered input
    Step#2 Perform Hadamard Transform
    """
    x=get_sequency_list(X)
    N = x.size
    M = int(log(N,2))
    out = x.copy()
    for m in xrange(M):
        outtemp = out.copy()         
        step = 2**m
        numCalc = 2**m
        for g in xrange(0,N,2*step): # number of groups
            for c in xrange(numCalc):
                index = g + c
                out[index] = outtemp[index] + outtemp[index+step]
                out[index+step] = outtemp[index] - outtemp[index+step]
    return out/float(N)


def FWHT(x):
    """ Fast Walsh-Hadamard Transform
    Based on mex function written by Chengbo Li@Rice Uni for his TVAL3 algorithm.
    His code is according to the K.G. Beauchamp's book -- Applications of Walsh and Related Functions.
    """
    x = x.squeeze()
    N = x.size
    G = N/2 # Number of Groups
    M = 2 # Number of Members in Each Group

    # First stage
    y = np.zeros((N/2,2))
    y[:,0] = x[0::2] + x[1::2]
    y[:,1] = x[0::2] - x[1::2]
    x = y.copy()
    # Second and further stage
    for nStage in xrange(2,int(log(N,2))+1):
        y = np.zeros((G/2,M*2))
        y[0:G/2,0:M*2:4] = x[0:G:2,0:M:2] + x[1:G:2,0:M:2]
        y[0:G/2,1:M*2:4] = x[0:G:2,1:M:2] - x[1:G:2,1:M:2]
        y[0:G/2,2:M*2:4] = x[0:G:2,0:M:2] - x[1:G:2,0:M:2]
        y[0:G/2,3:M*2:4] = x[0:G:2,1:M:2] + x[1:G:2,1:M:2]
        x = y.copy()
        G = G/2
        M = M*2
    x = y[0,:]
    # x = x.reshape((x.size,1))
    return x/float(N)

def myfwht(x):
    """A recursive implementation of the 1D Cooley-Tukey FFT"""
    # x = np.asarray(x, dtype=float)
    N = x.shape[0]
    
    if N == 1:
    	return x
    else:
        X_even = myfwht(x[0:(N/2)])
        X_odd = myfwht(x[(N/2):])
        return np.concatenate([(X_even + X_odd),
                               (X_even - X_odd)])

def myvfwht(x):
    """ Fast Walsh-Hadamard Transform
    Based on mex function written by Chengbo Li@Rice Uni for his TVAL3 algorithm.
    His code is according to the K.G. Beauchamp's book -- Applications of Walsh and Related Functions.
    """
    x = x.squeeze()
    N = x.size

    M = N/2
    G = 2
    # First stage
    y = np.zeros((M,G))

    y[:,0] = x[0:M] + x[M:]
    y[:,1] = x[0:M] - x[M:]

    x = y.copy()
    # Second and further stage
    for nStage in xrange(2,int(log(N,2))+1):

    	next_M = M/2
    	next_G = G*2

        y = np.zeros((next_M,next_G))
        
        y[0:next_M,0:next_G:2] = x[0:next_M,0:G] + x[next_M:,0:G]
        y[0:next_M,1:next_G:2] = x[0:next_M,0:G] - x[next_M:,0:G]

        x = y.copy()
        G = next_G
        M = next_M

    x = y[0,:]
    # x = x.reshape((x.size,1))
    return x/float(N)
   
if __name__ == "__main__":
    x = np.random.random(1024**2)
    
    t1 = clock()
    y1 = SFWHT(x)
    t1 = clock() - t1
    print(t1)
    
    t2 = clock()
    y2 = FWHT(x)
    t2 = clock() - t2
    print(t2)
