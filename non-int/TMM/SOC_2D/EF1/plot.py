import os
from re import L
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


current_path = os.getcwd()


M_list = [4,8,16,32]

W_list = [2.0+ii*0.2 for ii in range(12)]
#W_list = [7.0,7.4,7.8,8.0,8.2,8.4,8.6,9.0,9.4]
color_list = ['r','k','b','m','pink','orange']

fig = plt.figure(figsize=(4.5,4.5))
sub = fig.add_subplot(111)
plt.tick_params(size=8,width=1.5,labelsize = 12)
sub.spines['bottom'].set_linewidth(1.6)
sub.spines['left'].set_linewidth(1.6)
sub.spines['right'].set_linewidth(1.6)
sub.spines['top'].set_linewidth(1.6)

sub.spines['bottom'].set_color('black')
sub.spines['top'].set_color('black')
sub.spines['right'].set_color('black')
sub.spines['left'].set_color('black')


for ii in range(len(M_list)):
    M = M_list[ii]
    Length_list = []
    Err_list = []
    for jj in range(len(W_list)):
        W = W_list[jj]

        L_list = []
        for kk in range(5):
            path = os.path.join(current_path, 'M'+str(M), 'W'+"{:.2f}".format(W),str(kk), 'out')
            if os.path.exists(path):
                with open(path, 'r') as file:
                    orbital_value = file.readlines()
                    file.close()
                for ss in range(len(orbital_value)):
                    orbital_value[ss] = orbital_value[ss].split()
            
                for ss in range(len(orbital_value)):
                    if(len(orbital_value[ss]) > 0): 
                        if(orbital_value[ss][0] == 'L_list:'):
                            print(orbital_value[ss][1][1:-1])

                            L_list.append(float(orbital_value[ss][1][1:-1]))
        L_list = np.array(L_list)
        mean_value = np.mean(L_list)
        std_value = np.std(L_list)

        Length_list.append(mean_value)
        Err_list.append(std_value)

    plt.plot(1*np.array(W_list) ,np.array(Length_list)/M,'-',color=color_list[ii],linewidth = 1.5)

    plt.errorbar(1*np.array(W_list), np.array(Length_list)/M, 
                 yerr=np.array(Err_list)/M/sqrt(5),  # y方向的误差\
                 fmt='-',       # 线条和标记样式
                 color=color_list[ii], 
                 linewidth=1.5, 
                 capsize=3,      # 误差线端帽长度
                 capthick=2,     # 误差线端帽粗细
                 elinewidth=2,label = 'M='+str(M))   # 误差线粗细
#plt.legend(prop={'size':12})
leg = plt.legend(frameon=True, prop={'size':12}, loc='upper right')
leg.get_frame().set_edgecolor('none')
leg.get_frame().set_linewidth(0.0)
#plt.ylim([0.0,1.3])
plt.yticks([1,2,3,4])
plt.savefig('2D-SOC-EF1-TMM.svg')


