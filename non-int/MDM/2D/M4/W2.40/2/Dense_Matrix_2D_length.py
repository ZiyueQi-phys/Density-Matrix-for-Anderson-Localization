

import matplotlib.pyplot as plt
import numpy as np
import pylab as pl
import math as math
import scipy as sp
import random
#from qutip import Qobj, basis, SESolver, SolverOptions, MESolver, sesolve
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
from tqdm import tqdm

@nb.jit()
def get_Ham_sprase(g,W1,W2,ndim1,ndim2):
    row_indices = []
    col_indices = []
    data = []

    ### disorder term ###
    for ii in range(ndim1):
        for jj in range(ndim2):
            row_indices.append(ii*ndim2+jj)
            col_indices.append(ii*ndim2+jj)
            data.append(random.uniform(W1,W2))
    
    ### hopping term ###
    for ii in range(ndim1):
        for jj in range(ndim2):
            if(jj==ndim2*10):
                    continue
            row_indices.append(ii*ndim2+jj)
            col_indices.append(ii*ndim2+(jj+1)%ndim2)
            data.append(-g)
            
            col_indices.append(ii*ndim2+jj)
            row_indices.append(ii*ndim2+(jj+1)%ndim2)
            data.append(-g)

    for ii in range(ndim1):
        for jj in range(ndim2):

            if(ndim1 > 1.01):
                if(ii==ndim1*10 ):
                    1#print(ii)
                    continue
                #if(ndim1 == 2 and ii==1):
                #    continue
                row_indices.append(ii*ndim2+jj)
                col_indices.append((ii+1)%ndim1*ndim2+jj)
                data.append(-g)

                col_indices.append(ii*ndim2+jj)
                row_indices.append((ii+1)%ndim1*ndim2+jj)
                data.append(-g)
            '''
            elif(ndim1 < 2.01 and ndim1 > 1.01):
                row_indices.append(0*ndim2+jj)
                col_indices.append(1*ndim2+jj)
                data.append(-g)

                col_indices.append(0*ndim2+jj)
                row_indices.append(1*ndim2+jj)
                data.append(-g)
            '''
    row_indices = np.array(row_indices)
    col_indices = np.array(col_indices)
    data = np.array(data)

    Ham_sparse = csr_matrix((data, (row_indices, col_indices)), shape=(ndim1*ndim2, ndim1*ndim2))
    return Ham_sparse


@nb.jit()
def Diagonal_Ham(Ham_sparse,ndim1,ndim2):

    #e,v = scipy.sparse.linalg.eigsh(Ham_sparse, k=ndim1*ndim2//2, which='SA')
    #e,v = np.linalg.eigh(Ham_sparse.todense())
    e,v = eigs(
            Ham_sparse,
            k=1,
            sigma=0,
            which='LM',
            tol=1e-6
        )
    return e,v
    #print("e:",e)
    #print("v:",v)

    return e,v

def custom_func(x, a, b):
    return a*np.exp(-b*x)



@nb.jit()
def Calc_corr(e,v,ndim1,ndim2,num_iter,num_L):

    C_list = np.zeros((ndim1, ndim1, num_L),dtype=float)
    for ss in range(len(e)):
            v0 = v[:,ss]
            for nn in range(num_L):
                for aa in range(ndim1):
                    for cc in range(ndim1):
                        
                        for mm in range(9,ndim2-num_L-10):
                            #C_list[aa, cc, nn] = C_list[aa, cc, nn] + ( v0[aa*ndim2 + mm]*v0[cc*ndim2 + nn+1+mm] + v0[cc*ndim2 + mm]*v0[aa*ndim2 + nn+1+mm] )**1 / (ndim2-num_L-10) / num_iter
                            C_list[aa, cc, nn] = C_list[aa, cc, nn] + abs( v0[aa*ndim2 + mm]*v0[cc*ndim2 + nn+1+mm] )**1 / (ndim2-30-20) / num_iter

    return C_list


def Correlation_length(g,W1,W2,ndim1,ndim2):

    num_iter = 3000

    num_L = 30
    R_list = np.array([ii+1 for ii in range(num_L)])
    C_list = np.zeros((ndim1, ndim1, num_L),dtype=float)
    for ii in tqdm(range(num_iter)):
        #print("num_iter:",ii)
        Ham_sparse = get_Ham_sprase(g,W1,W2,ndim1,ndim2)

        e,v = Diagonal_Ham(Ham_sparse,ndim1,ndim2)
        
        C_temp = Calc_corr(e,v,ndim1,ndim2,num_iter,num_L)
        C_list = C_list + C_temp
        #for ss in range(len(e)):
        #    v0 = v[:,ss]
        #    for nn in range(num_L):
        #        for mm in range(ndim3-num_L-10):
        #            C_list[nn] = C_list[nn] + np.abs( v0[mm]*v0[nn+1+mm] )**2 / (ndim3-num_L-10) / num_iter
    
    C_final_list = []
    C_final_list_svd = []
    for ii in range(num_L):
        C_temp = C_list[:,:,ii]
        e,v = np.linalg.eigh(np.dot(C_temp,  C_temp.transpose()))
        #u,e,v = np.linalg.svd(np.dot(C_temp,  1))
        #e,v = np.linalg.eigh(np.dot(C_temp,  C_temp.transpose()))
        #print(ii,e)
        u, sigma, vt = np.linalg.svd(C_temp)
        #e,v = np.linalg.eigh(C_temp+ C_temp.transpose())
        #e = sigma
        #print("e:",e)
        #print("sigma:",sigma)
        C_final_list.append(((np.max(e))))
        #C_final_list_svd.append(np.sqrt((np.max(sigma)))/2)


    R_list_eff = []
    C_list_eff = []
    for ii in range(num_L):
        #if(ii%(ndim1//2) == 0):
        R_list_eff.append(R_list[ii])
        C_list_eff.append(C_final_list[ii])
    print("C_list:")
    for kk in range(num_L):
        for ii in range(ndim1):
            for jj in range(ndim1):
                1#print(ii,jj,kk,C_list[ii,jj,kk])
    #print(e)
    a4, b4= curve_fit(custom_func, R_list_eff, C_list_eff)[0]
    x4 = np.arange(1,num_L, 0.01)
    y4 = a4*np.exp(-b4*x4)
    print(a4, b4)

    
    '''
    # 绘制自相关函数和拟合结果
    plt.figure(figsize=(8, 6))
    print("C_final_list:",C_final_list)
    plt.scatter(R_list, C_final_list, color='r',label='positive')
    
    #plt.scatter(R_list, C_final_list_svd,color='b',label='svd')
    plt.plot(x4, y4, '--')
    plt.xlabel('Distance r')
    plt.ylabel('C(r)')
    plt.legend()
    plt.title('Localization Length from Green's Function')
    plt.show()
    
    '''
    
    return 1/b4


def Length_W(g,ndim1,ndim2,color):
    num_W = 1


    W_list = []
    L_list = []
    for ii in range(num_W):
        W = 2.4
        W_list.append(W)
        L_temp = Correlation_length(g,-W,W,ndim1,ndim2)

        L_list.append(L_temp)

    print("ndim2:",ndim2)
    print("W_list:",W_list)
    print("L_list:",L_list)
    #plt.figure()
    #plt.plot(W_list,np.array(L_list)/ndim1,'-*',color=color,label='M='+str(ndim1))
    #plt.show()
    #plt.savefig('Length_W.pdf')


def main():
    g = 1
    W1 = -8.5
    W2 = 8.5
    #plt.figure()
    dim_list = [4,8,16,32]
    color_list = ['r','k','b','m','pink','orange']
    for ii in range(1):
        ndim1 = 4
        ndim2 = 500
        Length_W(g,ndim1,ndim2,color_list[ii])



if __name__ == '__main__':
    main()


