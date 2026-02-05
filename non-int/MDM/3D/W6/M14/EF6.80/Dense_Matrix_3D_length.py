


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
from tqdm import tqdm
from joblib import Parallel, delayed

@nb.jit()
def get_Ham_sprase(g,W1,W2,ndim1,ndim2,ndim3):
    row_indices = []
    col_indices = []
    data = []

    ### disorder term ###
    for ii in range(ndim1):
        for jj in range(ndim2):
            for kk in range(ndim3):
                row_indices.append(ii*ndim2*ndim3+jj*ndim3+kk)
                col_indices.append(ii*ndim2*ndim3+jj*ndim3+kk)
                data.append(random.uniform(W1,W2))
    
    ### hopping term ###
    for ii in range(ndim1):
        for jj in range(ndim2):
            for kk in range(ndim3):
                row_indices.append(ii*ndim2*ndim3+jj*ndim3+kk)
                col_indices.append(ii*ndim2*ndim3+jj*ndim3+(kk+1)%ndim3)
                data.append(-g)

                col_indices.append(ii*ndim2*ndim3+jj*ndim3+kk)
                row_indices.append(ii*ndim2*ndim3+jj*ndim3+(kk+1)%ndim3)
                data.append(-g)

                if(ndim1 > 1.01 and ndim2 > 1.01):
                    row_indices.append(ii*ndim2*ndim3+jj*ndim3+kk)
                    col_indices.append((ii+1)%ndim1*ndim2*ndim3+jj*ndim3+kk)
                    data.append(-g)

                    col_indices.append(ii*ndim2*ndim3+jj*ndim3+kk)
                    row_indices.append((ii+1)%ndim1*ndim2*ndim3+jj*ndim3+kk)
                    data.append(-g)

                    row_indices.append(ii*ndim2*ndim3+jj*ndim3+kk)
                    col_indices.append(ii*ndim2*ndim3+(jj+1)%ndim2*ndim3+kk)
                    data.append(-g)

                    col_indices.append(ii*ndim2*ndim3+jj*ndim3+kk)
                    row_indices.append(ii*ndim2*ndim3+(jj+1)%ndim2*ndim3+kk)
                    data.append(-g)
    row_indices = np.array(row_indices)
    col_indices = np.array(col_indices)
    data = np.array(data)

    Ham_sparse = csr_matrix((data, (row_indices, col_indices)), shape=(ndim1*ndim2*ndim3, ndim1*ndim2*ndim3))
    return Ham_sparse


#@nb.jit()
def Diagonal_Ham(Ham_sparse,ndim1,ndim2,ndim3):

    #e,v = scipy.sparse.linalg.eigsh(Ham_sparse, k=ndim1*ndim2*ndim3//2, which='SA')
    
    e,v = eigs(
            Ham_sparse,
            k=1,
            sigma=6.800000000000001,
            which='LM',
            tol=1e-6
        )
    #print("e:",e)
    #print("v:",v)

    return e,v

def custom_func(x, a, b):
    return a*np.exp(-b*x)



@nb.jit()
def Calc_corr(e,v,ndim1,ndim2,ndim3,num_iter,num_L):

    C_list = np.zeros((ndim1*ndim2, ndim1*ndim2, num_L),dtype=float)
    for ss in range(len(e)):
            v0 = v[:,ss]
            for nn in range(num_L):
                for aa in range(ndim1):
                    for bb in range(ndim2):
                        for cc in range(ndim1):
                            for dd in range(ndim2):
                                if(1):
                                    for mm in range(10, ndim3-num_L-10):
                                        C_list[aa*ndim2+bb, cc*ndim2+dd, nn] = C_list[aa*ndim2+bb, cc*ndim2+dd, nn] + np.abs( v0[aa*ndim2*ndim3 + bb*ndim3 + mm]*v0[cc*ndim2*ndim3 + dd*ndim3 + nn+1+mm] )**1 / (ndim3-num_L-20) / num_iter
    return C_list





#@nb.jit()
def Correlation_length(g,W1,W2,ndim1,ndim2,ndim3):

    num_iter = 2500

    num_L = 30
    R_list = np.array([ii+1 for ii in range(num_L)])
    C_list = np.zeros((ndim1*ndim2, ndim1*ndim2, num_L),dtype=float)

    
    for ii in tqdm(range(num_iter)):
        #print("num_iter:",ii)
        Ham_sparse = get_Ham_sprase(g,W1,W2,ndim1,ndim2,ndim3)

        e,v = Diagonal_Ham(Ham_sparse,ndim1,ndim2,ndim3)
        
        C_temp = Calc_corr(e,v,ndim1,ndim2,ndim3,num_iter,num_L)
        C_list = C_list + C_temp
        #for ss in range(len(e)):
        #    v0 = v[:,ss]
        #    for nn in range(num_L):
        #        for mm in range(ndim3-num_L-10):
        #            C_list[nn] = C_list[nn] + np.abs( v0[mm]*v0[nn+1+mm] )**2 / (ndim3-num_L-10) / num_iter
    '''

    with Parallel(n_jobs=30, prefer="processes") as parallel:
        # 生成并行任务
        tasks = (delayed(process_iteration)(ii,g, W1, W2, ndim1, ndim2, ndim3,  num_iter, num_L) for ii in range(num_iter))
        # 使用tqdm包装任务迭代器以显示进度
        results = parallel(tqdm(tasks, total=num_iter, desc="Processing iterations"))

        # 累加结果
        for result in results:
            C_list += result
    '''

    C_final_list = []
    for ii in range(num_L):
        C_temp = C_list[:,:,ii]
        e,v = np.linalg.eigh(np.dot(C_temp, C_temp.transpose()))
        #e,v = np.linalg.eigh(C_temp+ C_temp.transpose())
        #print("e:",e)
        C_final_list.append(np.max(e))

    
    #print(e)
    a4, b4= curve_fit(custom_func, R_list, C_final_list)[0]
    x4 = np.arange(0,num_L, 0.01)
    y4 = a4*np.exp(-b4*x4)
    #print(a4,b4)


    '''
    # 绘制自相关函数和拟合结果
    plt.figure(figsize=(8, 6))
    plt.plot(R_list, C_final_list, label='C(r)')
    plt.plot(x4, y4, '--', label=str(1/b4))
    plt.xlabel('Distance r')
    plt.ylabel('C(r)')
    plt.legend()
    plt.title('Localization Length from Green's Function')
    plt.show()
    '''
    print("R_list=",R_list)
    print("C_final_list=",C_final_list)

    return 1/b4


def Length_W(g,ndim1,ndim2,ndim3):
    num_W = 1

    W_list = []
    L_list = []
    for ii in range(num_W):
        W = 6
        W_list.append(W)
        L_temp = Correlation_length(g,-W,W,ndim1,ndim2,ndim3)

        L_list.append(L_temp)

    print("W_list:",W_list)
    print("L_list:",L_list)
    '''
    plt.figure()
    plt.plot(W_list,L_list,'b-*')
    #plt.show()
    plt.savefig('Length_W.pdf')
    '''

def main():
    g = 1
    W1 = -8.5
    W2 = 8.5
    ndim1 = 14
    ndim2 = 14
    ndim3 = 300
    Length_W(g,ndim1,ndim2,ndim3)
    #Correlation_length(g,W1,W2,ndim1,ndim2,ndim3)


if __name__ == '__main__':
    main()
    


