

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

#@nb.jit()
def get_Ham_dense(g,W1,W2,ndim1,ndim2):
    
    Ham_dense = np.zeros((2*ndim1*ndim2, 2*ndim1*ndim2),  dtype=complex)

    for ii in range(ndim1):
        for jj in range(ndim2):
            eps_ij = random.uniform(W1,W2)
            Ham_dense[2*(ii*ndim2+jj), 2*(ii*ndim2+jj)] = eps_ij
            Ham_dense[2*(ii*ndim2+jj)+1, 2*(ii*ndim2+jj)+1] = eps_ij
    
    
    for ii in range(ndim1):
        for jj in range(ndim2):
            
            alpha = random.uniform(0, 2*pi)
            beta  = random.uniform(0,1)
            beta = 1/2 * np.arccos(1 - 2 * beta)
            gamma = random.uniform(0, 2*pi)

            V_su2 = np.zeros((2,2), dtype=np.complex128)
            V_su2[0,0] = cmath.exp(1j*alpha)*cos(beta)
            V_su2[0,1] = cmath.exp(1j*gamma)*sin(beta)
            V_su2[1,0] = -cmath.exp(-1j*gamma)*sin(beta)
            V_su2[1,1] = cmath.exp(-1j*alpha)*cos(beta)

            Ham_dense[2*(ii*ndim2+jj):2*(ii*ndim2+jj)+2, 2*(ii*ndim2+(jj+1)%ndim2):2*(ii*ndim2+(jj+1)%ndim2)+2] = -g*V_su2
            Ham_dense[2*(ii*ndim2+(jj+1)%ndim2):2*(ii*ndim2+(jj+1)%ndim2)+2, 2*(ii*ndim2+jj):2*(ii*ndim2+jj)+2] = -g*V_su2.transpose().conj()
    for ii in range(ndim1):
        for jj in range(ndim2):
            #if(ii==ndim1-1):
            #    continue

            alpha = random.uniform(0, 2*pi)
            beta  = random.uniform(0,1)
            beta = 1/2 * np.arccos(1 - 2 * beta)
            gamma = random.uniform(0, 2*pi)

            V_su2 = np.zeros((2,2), dtype=np.complex128)
            V_su2[0,0] = cmath.exp(1j*alpha)*cos(beta)
            V_su2[0,1] = cmath.exp(1j*gamma)*sin(beta)
            V_su2[1,0] = -cmath.exp(-1j*gamma)*sin(beta)
            V_su2[1,1] = cmath.exp(-1j*alpha)*cos(beta)

            Ham_dense[2*((ii)*ndim2+jj):2*((ii)*ndim2+jj)+2, 2*((ii+1)%ndim1*ndim2+jj):2*((ii+1)%ndim1*ndim2+jj)+2] = -g*V_su2
            Ham_dense[2*((ii+1)%ndim1*ndim2+jj):2*((ii+1)%ndim1*ndim2+jj)+2, 2*((ii)*ndim2+jj):2*((ii)*ndim2+jj)+2] = -g*V_su2.transpose().conj()
            

    return Ham_dense



#@nb.jit()
def Diagonal_Ham(Ham_dense,ndim1,ndim2):

    #e,v = scipy.sparse.linalg.eigsh(Ham_sparse, k=ndim1*ndim2//2, which='SA')
    #e,v = np.linalg.eigh(Ham_dense)
    #print("e:",e)
    #print("v:",v)
    Ham_sparse = csr_matrix(Ham_dense)




    e,v = eigs(
            Ham_sparse,
            k=1,
            sigma=0,
            which='LM',
            tol=1e-6
        )
    return e,v

def custom_func(x, a, b):
    return a*np.exp(-b*x)



#@nb.jit()
def Calc_corr(e,v,ndim1,ndim2,num_iter,num_L):

    C_list = np.zeros((ndim1*2, ndim1*2, num_L),dtype=float)
    for ss in range(len(e)):
        if(1):
            #print(e[ss])
            v0 = v[:,ss]
            for nn in range(num_L):
                for aa in range(ndim1):
                    for cc in range(ndim1):
                        for x in range(2):
                            for y in range(2):
                                for mm in range(30, ndim2-num_L-30):
                                    C_list[aa*2+x, cc*2+y, nn] = C_list[aa*2+x, cc*2+y, nn] + abs(v0[2*(aa*ndim2 + mm)+x].conj() * v0[2*(cc*ndim2 + nn+1+mm)+y] **1) / (ndim2-num_L-60) / num_iter

    return C_list



def process_iteration(ii,g, W1, W2, ndim1, ndim2, num_iter, num_L):
    
    # 生成哈密顿量矩阵
    Ham_sparse = get_Ham_dense(g, W1, W2, ndim1, ndim2)
    # 对角化哈密顿量
    e, v = Diagonal_Ham(Ham_sparse, ndim1, ndim2)
    # 计算关联函数
    return Calc_corr(e, v, ndim1, ndim2, num_iter, num_L)

def Correlation_length(g,W1,W2,ndim1,ndim2):

    num_iter = 2000

    num_L = 50
    R_list = np.array([ii+1 for ii in range(num_L)])
    C_list = np.zeros((ndim1*2, ndim1*2, num_L),dtype=float)
    
    '''
    for ii in tqdm(range(num_iter)):
        #print("num_iter:",ii)
        Ham_sparse = get_Ham_dense(g,W1,W2,ndim1,ndim2)

        e,v = Diagonal_Ham(Ham_sparse,ndim1,ndim2)
        #print(e)
        C_temp = Calc_corr(e,v,ndim1,ndim2,num_iter,num_L)
        C_list = C_list + C_temp
        #for ss in range(len(e)):
        #    v0 = v[:,ss]
        #    for nn in range(num_L):
        #        for mm in range(ndim3-num_L-10):
        #            C_list[nn] = C_list[nn] + np.abs( v0[mm]*v0[nn+1+mm] )**2 / (ndim3-num_L-10) / num_iter
    '''
    with Parallel(n_jobs=50, prefer="processes") as parallel:
        # 生成并行任务
        tasks = (delayed(process_iteration)(ii,g, W1, W2, ndim1, ndim2, num_iter, num_L) for ii in range(num_iter))
        # 使用tqdm包装任务迭代器以显示进度
        results = parallel(tqdm(tasks, total=num_iter, desc="Processing iterations"))

        # 累加结果
        for result in results:
            C_list += result
    
    C_final_list = []
    for ii in range(num_L):
        C_temp = C_list[:,:,ii]
        e,v = np.linalg.eigh(np.dot(C_temp, C_temp.transpose()))
        #e,v = np.linalg.eigh(C_temp+ C_temp.transpose())
        #print("e:",e)
        C_final_list.append((np.max(e)))
        

    
    #print(e)
    a4, b4= curve_fit(custom_func, R_list, C_final_list)[0]
    x4 = np.arange(0,num_L, 0.01)
    y4 = a4*np.exp(-b4*x4)
    #print(a4, b4)

    print("C_final_list:", C_final_list)
    '''    
    # 绘制自相关函数和拟合结果
    plt.figure(figsize=(8, 6))
    plt.scatter(R_list, C_final_list, label='C(r)')
 
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
        W = 4.2+ii/num_W*2.5
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
    color_list = ['r','k','b','m','pink','orange']
    for ii in range(1):
        ndim1 = 8+2*ii
        ndim2 = 500
        Length_W(g,ndim1,ndim2,color_list[ii])
    #plt.legend()
    #plt.show()
    #plt.savefig('Length_W_2D.pdf')
    #Correlation_length(g,W1,W2,ndim1,ndim2,ndim3)


if __name__ == '__main__':
    main()



