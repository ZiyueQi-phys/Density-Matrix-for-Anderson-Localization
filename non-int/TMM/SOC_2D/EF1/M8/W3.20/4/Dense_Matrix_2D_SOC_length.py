

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
import scipy.stats as stats
from tqdm import tqdm

@nb.jit()
def get_TMM(M,g,W1,W2,EF, Vfront):

    TM = np.zeros((M*4, M*4),  dtype=np.complex128)

    Ham = np.zeros((M*2, M*2),  dtype=np.complex128)

    for ii in range(M):
        eps_i = random.uniform(W1,W2)
        Ham[2*ii, 2*ii] = eps_i
        Ham[2*ii+1, 2*ii+1] = eps_i

        #Ham[2*ii:2*ii+2, 2*ii:2*ii+2] = Ham[2*ii:2*ii+2, 2*ii:2*ii+2] + np.eye(2,dtype=complex)*eps_i

    for ii in range(M-1):
        alpha = random.uniform(0, 2*pi)
        beta  = random.uniform(0,1)
        beta = 1/2 * np.arccos(1 - 2 * beta)
        gamma = random.uniform(0, 2*pi)


        V_su2 = np.zeros((2,2), dtype=np.complex128)
        V_su2[0,0] = cmath.exp(1j*alpha)*cos(beta)
        V_su2[0,1] = cmath.exp(1j*gamma)*sin(beta)
        V_su2[1,0] = -cmath.exp(-1j*gamma)*sin(beta)
        V_su2[1,1] = cmath.exp(-1j*alpha)*cos(beta)

        Ham[2*ii:2*ii+2, 2*(ii+1):2*(ii+1)+2] = Ham[2*ii:2*ii+2, 2*(ii+1):2*(ii+1)+2] - V_su2
        Ham[2*(ii+1):2*(ii+1)+2, 2*ii:2*ii+2] = Ham[2*(ii+1):2*(ii+1)+2, 2*ii:2*ii+2] - V_su2.transpose().conj()
    
    
    alpha = random.uniform(0, 2*pi)
    beta  = random.uniform(0,1)
    beta = 1/2 * np.arccos(1 - 2 * beta)
    gamma = random.uniform(0, 2*pi)
    V_su2 = np.zeros((2,2), dtype=np.complex128)
    V_su2[0,0] = cmath.exp(1j*alpha)*cos(beta)
    V_su2[0,1] = cmath.exp(1j*gamma)*sin(beta)
    V_su2[1,0] = -cmath.exp(-1j*gamma)*sin(beta)
    V_su2[1,1] = cmath.exp(-1j*alpha)*cos(beta)

    Ham[2*(M-1):2*(M-1)+2, 0:2] = Ham[2*(M-1):2*(M-1)+2, 0:2] - V_su2
    Ham[0:2, 2*(M-1):2*(M-1)+2] = Ham[0:2, 2*(M-1):2*(M-1)+2] - V_su2.transpose().conj()


    V_front = np.zeros((M*2, M*2),  dtype=np.complex128)
    V_back = np.zeros((M*2, M*2),  dtype=np.complex128)

    for ii in range(M):
        alpha = random.uniform(0, 2*pi)
        beta  = random.uniform(0,1)
        beta = 1/2 * np.arccos(1 - 2 * beta)
        gamma = random.uniform(0, 2*pi)

        V_su2 = np.zeros((2,2), dtype=np.complex128)
        V_su2[0,0] = cmath.exp(1j*alpha)*cos(beta)
        V_su2[0,1] = cmath.exp(1j*gamma)*sin(beta)
        V_su2[1,0] = -cmath.exp(-1j*gamma)*sin(beta)
        V_su2[1,1] = cmath.exp(-1j*alpha)*cos(beta)

        V_front[2*ii:2*ii+2, 2*ii:2*ii+2] = V_front[2*ii:2*ii+2, 2*ii:2*ii+2] - V_su2

        alpha = random.uniform(0, 2*pi)
        beta  = random.uniform(0,1)
        beta = 1/2 * np.arccos(1 - 2 * beta)
        gamma = random.uniform(0, 2*pi)

        V_su2 = np.zeros((2,2), dtype=complex)
        V_su2[0,0] = cmath.exp(1j*alpha)*cos(beta)
        V_su2[0,1] = cmath.exp(1j*gamma)*sin(beta)
        V_su2[1,0] = -cmath.exp(-1j*gamma)*sin(beta)
        V_su2[1,1] = cmath.exp(-1j*alpha)*cos(beta)

        V_back[2*ii:2*ii+2, 2*ii:2*ii+2] = V_back[2*ii:2*ii+2, 2*ii:2*ii+2] - V_su2#.transpose().conj()

    #if(abs(np.trace(Vfront))>0.000001):
    #    V_back = Vback.transpose().conj()
    if(abs(np.trace(Vfront))>0.000001):
        #    print("wrong!!!")
        V_front = Vfront.transpose().conj()
    
    V_back_inv = np.linalg.inv(V_back)

    TM[0:2*M,0:2*M] = np.dot(V_back_inv, EF*np.eye(2*M,dtype=complex) - Ham )
    TM[0:2*M,2*M:4*M] = -np.dot(V_back_inv, V_front )
    TM[2*M:4*M, 0:2*M] = np.eye(2*M,dtype=complex)

    return TM, V_back


@nb.jit()
def get_localization_length(M,g,W,EF,NL):

    TM = np.eye((4*M),dtype=complex)

    R_log_list = np.zeros((4*M),dtype=float)

    R_temp_list = np.zeros((4*M),dtype=float)
    R_accumulate = []
    R_average_list = [0]
    V_front = np.zeros((M*2, M*2),  dtype=np.complex128)
    for ii in range(NL):
        #print("iter:",ii)
        TM_temp, V_front = get_TMM(M,g,-W,W,EF, V_front)

        TM = TM_temp.dot( ( TM ) )

        Q,R = np.linalg.qr(TM, mode='complete')

        #for ss in range(len(Q)):
        #    print("Q:",Q[ss])
        #    print("R:",R[ss])
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

        R_inv_abs = np.zeros((4*M,4*M),dtype=complex)
        for ss in range(4*M):
            if(abs(R[ss,ss])>0.000001):
                R_inv_abs[ss,ss] = 1/abs(R[ss,ss])
        R_diag = np.zeros((4*M,4*M),dtype=complex)
        for ss in range(4*M):
            R_diag[ss,ss] = R[ss,ss]
        TM = np.dot(Q, np.dot( R_diag,  R_inv_abs  )   )
        #print(temp)
    #print("R_log_list")
    print(np.sort(R_log_list))
    R_log_list = np.sort(np.abs(np.array(R_log_list)))
    print("localization_length:",1/min(R_log_list)/M)
    print(np.average(R_accumulate))

    return np.average(R_accumulate)



def main():

    g = 1
    W = 0
    EF = 1
    M = 4
    NL = 300000

    #TM = get_TMM(M,g,-W,W,EF)
    #e,v = np.linalg.eigh(TM)
    #print(e)
    num_W = 1
    num_iter = 1

    #print(get_TMM(M,g,0,0,EF))

    #get_localization_length(M,g,W,EF,NL)
    
    color_list = ['r','k','b','m','pink','orange','g']
    

    
    #plt.figure()
    

    for nn in range(1):
        M = 8+4*nn
        W_list = []
        Length_list = []
        for ii in range(num_W):
            leng = 0
            W = 3.2+ii*0.2

            print("W:",W)
            for ss in tqdm(range(num_iter)):
                print("num_iter:",ss)
                leng_temp = get_localization_length(M,g,W,EF,NL)
                leng = leng + leng_temp/num_iter*M
            W_list.append(W)
            Length_list.append(leng)
            
            print("W_list:",W_list)
            print("L_list:",[leng])
        
        #plt.plot(W_list,Length_list,'-*',color=color_list[nn])
    #plt.show()
    #plt.savefig('TMM.pdf')
    #print(TM)
    
    #print("1")
    #print( (TM.transpose()).dot(TM) )
    #e,v = np.linalg.eigh((TM.transpose()).dot(TM) )
    
    #print(e)

if __name__ == '__main__':
    main()    




