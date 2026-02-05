import os
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




num_iter = 10
num_W = 6
num_M = 3


color_list = ['r','b','m','m','pink','orange']

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


Length_list = np.zeros((num_W,num_M,num_iter), dtype=float)


W_list = np.array([2.3+0.1*ii for ii in range(6)])


Mean_value_list = np.zeros((num_W,num_M), dtype=float)
Err_value_list = np.zeros((num_W,num_M), dtype=float)



for ii in range(0,0+num_iter):
    
    #if(ii==30):
    #    continue
    path = os.path.join(current_path, str(ii+0), 'out')
        
    with open(path, 'r') as file:
        orbital_value = file.readlines()
        file.close()
    for ss in range(len(orbital_value)):
        orbital_value[ss] = orbital_value[ss].split()
    
    
    print(ii) 
    
    M_index_temp = 0

    for ss in range(len(orbital_value)):
        if(len(orbital_value[ss]) > 0): 
            if(orbital_value[ss][0] == 'L_list:'):
                #print(ss)
                Length_list[0,M_index_temp,ii] = float(orbital_value[ss][1][1:-1])
                Length_list[1,M_index_temp,ii] = float(orbital_value[ss][2][0:-1])
                Length_list[2,M_index_temp,ii] = float(orbital_value[ss][3][0:-1])
                Length_list[3,M_index_temp,ii] = float(orbital_value[ss][4][0:-1])
                Length_list[4,M_index_temp,ii] = float(orbital_value[ss][5][0:-1])
                Length_list[5,M_index_temp,ii] = float(orbital_value[ss][6][0:-1])
                M_index_temp = M_index_temp + 1
    if(M_index_temp!=3):
       # print(ii)
        print(M_index_temp)
        print("Wrong!!!!!!!!!!")
        exit()

for ii in range(num_W):
    for jj in range(num_M):
        M = jj+2
        mean_value = np.mean(Length_list[ii,jj,:]/M)
        err_value = np.std(Length_list[ii,jj,:]/M)
        Mean_value_list[ii,jj] = mean_value
        Err_value_list[ii,jj] = err_value

for ii in range(num_M):
    plt.plot(W_list, Mean_value_list[:,ii],color=color_list[ii], linewidth = 2)
    plt.errorbar( W_list, Mean_value_list[:,ii], yerr=Err_value_list[:,ii]/sqrt(num_iter), fmt='o', color=color_list[ii], ecolor=color_list[ii], capsize=5, label='M='+str(ii+2) )



                    
leg = plt.legend(frameon=True, prop={'size':12.6}, loc='upper right')
leg.get_frame().set_edgecolor('none')
leg.get_frame().set_linewidth(0.0)
#plt.legend(prop={'size':12})

plt.yticks([1.2,1.6,2.0,2.4])
#plt.ylim(top=4.5)
plt.savefig('length_scaling.svg')

