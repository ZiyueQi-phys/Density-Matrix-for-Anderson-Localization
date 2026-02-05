import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
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
import scipy.sparse
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import eigs, eigsh
import time
from scipy.optimize import curve_fit



#@nb.jit()
def get_Ham_sprase(g,W1,W2,ndim,ii):

    current_path = os.getcwd()
    #sample_ind = 66
    path = file_path = os.path.join(current_path, str(ii) ,'out')
    with open(path, 'r') as phase:
        orbital_value = phase.readlines()
    phase.close()

    L1 = 120
    L2 = 1
    W_list = []
    for ii in range(L1*L2):
        orbital_value[ii] = orbital_value[ii].split(':')
        
        string_temp = orbital_value[ii]
        #print(string_temp)
        #x = int(float(string_temp[0])) - 1
        #y = int(float(string_temp[1])) - 1
        mu = (float(string_temp[1][:-1])) 
        W_list.append(mu)


    #print("W_list:",W_list)
    h0 = np.mat(np.zeros((ndim,ndim),dtype=complex))
    h1 = np.mat(np.zeros((ndim,ndim),dtype=complex))


    for ii in range(ndim-1):
        h0[ii,(ii+1)%ndim] = -g
        h0[(ii+1)%ndim,ii] = -g
    #h0[ndim-1,0] = -1*g*cmath.exp(1j*pi*5/4)
    #h0[0,ndim-1] = -1*g*cmath.exp(-1j*pi*5/4)
    
    for ii in range(ndim):
        h1[ii,ii] = W_list[ii]
    return h0+h1


#@nb.jit()
def Diagonal_Ham(Ham_dense,ndim3):

    #e,v = scipy.sparse.linalg.eigsh(Ham_dense, k=ndim3-8, which='SA')
    e,v = np.linalg.eigh(Ham_dense)
    #print("e:",e)
    #print("v:",v)

    return e,v

def custom_func(x, a, b):
    return a*np.exp(-b*x)

#@nb.jit()
def Calc_corr(e,v,ndim3,num_iter,num_L):

    
    C_list = np.zeros((num_L),dtype=float)
    for ss in range(len(e)):
        if(ss==ndim3//2-1):
            v0 = v[:,ss]
            for nn in range(num_L):
                for mm in range(19,ndim3-40):
                    C_list[nn] = C_list[nn] + np.abs( v0[mm]*(v0[(nn+mm+1)%ndim3]) ) / (ndim3-60) / num_iter
    return C_list

def Correlation_length(g,W1,W2,ndim3):

    num_iter = 500

    num_L =20
    R_list = np.array([ii+1 for ii in range(num_L)])
    C_list = np.zeros((num_L),dtype=float)
    for ii in range(0,0+num_iter):
        print("num_iter:",ii)
        #if(ii==8 or ii==9 or ii==19 or ii==2 or ii==3 or ii==5 or ii==19 or ii==42 or ii==43 or ii==61 or ii==65 or ii==69 or ii==79 or ii==80 or ii==90 or ii==92 or ii==93 or ii==97 or ii==98 or ii==99):
        #if(ii==8 or ii==9 or  ii==10 or ii==11 or  ii==19 or (ii>=70 and ii<80)  or ii==2 or ii==5 or ii==3  or ii==42 or ii==43 or ii==61 or ii==65 or ii==69  or ii==80 or ii==90 or ii==92 or ii==93 or ii==97 or ii==98 or ii==99):
        #    continue
        #if(ii==2 or ii== 5 or ii== 8 or ii== 9 or ii== 10 or ii== 11 or ii== 17 or ii== 19 or ii== 21 or ii== 23 or ii== 28 or ii== 36 or ii== 40 or ii== 45 or ii== 50 or ii== 52 or ii== 57 or ii== 61 or ii== 63 or ii== 65 or ii== 68 or ii== 69 or ii== 70 or ii== 79 or ii== 88 or ii== 90 or ii== 91 or ii== 92 or ii== 93 or ii== 98 or ii== 99):
        #    continue
        #if(ii==42 or ii==43 or ii==97):
        #    continue
        #if(ii==8):
        #    continue
        Ham_sparse = get_Ham_sprase(g,W1,W2,ndim3,ii)

        e,v = Diagonal_Ham(Ham_sparse,ndim3)
        
        C_temp = Calc_corr(e,v,ndim3,num_iter,num_L)
        for ss in range(len(C_temp)):
            C_temp[ss] = C_temp[ss]**1  #C_temp[ss].conj()
        C_list = C_list + C_temp
        #for ss in range(len(e)):
        #    v0 = v[:,ss]
        #    for nn in range(num_L):
        #        for mm in range(ndim3-num_L-10):
        #            C_list[nn] = C_list[nn] + np.abs( v0[mm]*v0[nn+1+mm] )**2 / (ndim3-num_L-10) / num_iter
    
    C_list = C_list**2
    a4, b4= curve_fit(custom_func, R_list, C_list,maxfev=20000)[0]
    x4 = np.arange(1,num_L, 0.01)
    y4 = a4*np.exp(-b4*x4)
    print(a4,b4)


    print("C_final_list:", list(C_list))
    # 绘制自相关函数和拟合结果
    plt.figure(figsize=(8, 6))
    plt.scatter(R_list, C_list, label='C(r)')
    plt.plot(x4, y4, '--', label=f'Fit: ξ = {1/b4:.2f}')
    plt.xlabel('Distance r')
    plt.ylabel('C(r)')
    plt.legend()
    plt.title('Localization Length from Green\'s Function')
    plt.savefig('bench.pdf')
    

    

    return 1/b4
            

def Length_W(g,ndim3,color):
    num_W =1

    W_list = []
    L_list = []
    for ii in range(num_W):
        #print("num_W:",ii)
        W = 2
        W_list.append(W)
        L_temp = Correlation_length(g,-W,W,ndim3)

        L_list.append(L_temp)
    plt.figure()
    plt.plot(W_list,L_list,'-*',color = color, label='GS_'+str(ndim3))
    #plt.ylim(bottom=0)
    plt.show()

    print("W_list=", W_list)
    print("L_list=", L_list)

def calc_total_energy(g,W1,W2,ndim3,sss):
    Ham_sparse = get_Ham_sprase(g,W1,W2,ndim3,sss)

    e,v = Diagonal_Ham(Ham_sparse,ndim3)

    Ntot = ndim3
    Etot = 0
    #print("energy_list:", list(e))
    for ii in range(59-0):
        Etot = Etot + e[ii]
    #for ii in range(59-0):
    #    Etot = Etot + e[ii]
    #print(e[Ntot//2])
    #print(e[Ntot//2])
    print(sss, ",Etot:",Etot)


def main(ndim3,color):

    g = 1
    W1 = -2
    W2 = 2
    ndim1 = 1
    ndim2 = 1
    #ndim3 = 500
    #for ii in range(10):
    #    calc_total_energy(g,W1,W2,ndim3,ii)
    Length_W(g,ndim3,color)



    #Correlation_length(g,W1,W2,ndim1,ndim2,ndim3)


    #Ham_sparse = get_Ham_sprase(g,W1,W2,ndim1,ndim2,ndim3)

    #Diagonal_Ham(Ham_sparse,ndim1,ndim2,ndim3)


if __name__ == '__main__':
    main(120,'k')
