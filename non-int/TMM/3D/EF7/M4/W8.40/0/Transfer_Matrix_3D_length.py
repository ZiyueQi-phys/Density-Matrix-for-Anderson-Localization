


import math as math
import scipy as sp
import random
import numba as nb
from scipy.sparse.linalg import eigsh
import os
from numpy import conj
from math import *
import cmath
import matplotlib.cm as cm
from scipy.sparse import csr_matrix, identity
from scipy.sparse.linalg import eigs, eigsh
import time
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

def get_TMM(M,g,W1,W2,EF):
    TM = np.zeros((M*M*2, M*M*2),  dtype=np.complex128)

    for ii in range(M):
        for jj in range(M):
            TM[ii*M+jj, ii*M+jj] = (EF-random.uniform(W1,W2))/g

            TM[ii*M+jj, (ii+1)%M*M+jj] = -1
            TM[ii*M+jj, (ii-1)%M*M+jj] = -1

            TM[ii*M+jj, (ii)*M+(jj+1)%M] = -1
            TM[ii*M+jj, (ii)*M+(jj-1)%M] = -1
    
    for ii in range(M):
        for jj in range(M):
            TM[ii*M+jj, ii*M+jj+M*M] = -1
            TM[ii*M+jj+M*M, ii*M+jj] = 1

    return TM

@nb.jit()
def get_localization_length(M,g,W,EF,NL):

    TM = np.eye((2*M*M),dtype=complex)

    R_log_list = np.zeros((2*M*M),dtype=float)

    R_temp_list = np.zeros((2*M*M),dtype=float)
    R_accumulate = []
    R_average_list = [0]

    for ii in range(NL):
        #print("iter:",ii)
        TM_temp = get_TMM(M,g,-W,W,EF)

        TM = TM_temp.dot( ( TM ) )

        Q,R = np.linalg.qr(TM)
        #print(Q)
        #print(np.diag(R))
        #print(np.diag(np.dot(R.transpose(),R)))

        R_log_list = R_log_list + np.log(np.abs(np.diag(R)))/NL


        R_temp_list = R_temp_list + np.log(np.abs(np.diag(R)))
        temp = 1/min(np.abs(R_temp_list/(ii+1)))/M
        #print(  temp)
        
        R_accumulate.append(temp)
        #temp = np.average(R_accumulate)
        #if(abs(R_average_list[-1] - temp ) < 0.00001 ):
        #    break

        #R_average_list.append(temp)
        TM = Q
        #print(temp)
    #print("R_log_list")
    #print(np.sort(R_log_list))
    R_log_list = np.sort(np.abs(np.array(R_log_list)))
    #print("localization_length:",1/min(R_log_list)/M)
    #print(np.average(R_accumulate))

    return np.average(R_accumulate)
    '''
    R_log_list = np.sort(np.abs(np.array(R_log_list)))
    #print(np.diag(R))
    print(R_temp_list)
    
    print("lyapnov:",R_log_list)
    
    print("localization_length:",1/min(R_log_list)/M)
    '''

    #e,v = np.linalg.eigh(np.dot(TM.transpose(),TM))
    #print(e)



def main():

    g = 1
    W = 0
    EF = 7
    M = 4
    NL = 300000

    TM = get_TMM(M,g,-W,W,EF)
    e,v = np.linalg.eigh(TM)
    print(e)
    num_W = 1
    num_iter = 1

    #print(get_TMM(M,g,0,0,EF))

    #get_localization_length(M,g,W,EF,NL)
    
    color_list = ['r','k','b','m','pink','orange','g']
    
    W_list = [8.4]
    
    #plt.figure()
    
    for nn in range(1):
        M = 4+1*nn
        #W_list = []
        Length_list = []
        for ii in range(num_W):
            leng = 0
            W = 8.4

            print("W:",W)
            for ss in tqdm(range(num_iter)):
                print("num_iter:",ss)
                leng_temp = get_localization_length(M,g,W,EF,NL)
                leng = leng + leng_temp/num_iter*M
            #W_list.append(W)
            Length_list.append(leng)
        print("W_list:",W_list)
        print("L_list:",Length_list)
        plt.plot(W_list,Length_list,'-*',color=color_list[nn])
    #plt.show()
    #plt.savefig('TMM.pdf')
    #print(TM)
    
    #print("1")
    #print( (TM.transpose()).dot(TM) )
    #e,v = np.linalg.eigh((TM.transpose()).dot(TM) )
    
    #print(e)

if __name__ == '__main__':
    main()    
            




