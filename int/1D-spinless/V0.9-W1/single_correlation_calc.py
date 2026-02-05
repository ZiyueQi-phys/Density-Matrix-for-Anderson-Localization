import matplotlib.ticker as mticker
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


def custom_func(x, a, b, c, xi):
    #return a*np.exp(-b*x)
    return a*np.sin((b)*(x))/(x**1)*np.exp(-x/xi)


def custom_func_single(x, a, b, c, d, f, xi):
    return a*np.exp(-b*x) + d*np.exp(-f*x)*np.cos(pi*x)
    #return a*np.sin((b)*(x))/(x**1)*np.exp(-x/xi)



len_list = 20

L_list = np.array([ii+1 for ii in range(len_list)], dtype=float)
C_list_single = np.zeros((len_list),dtype=float)
C_list = np.zeros((len_list),dtype=float)
num_iter = 500
for ii in range(0,0+num_iter):

    path = file_path = os.path.join(current_path, str(ii),'SDM', 'out')
    with open(path, 'r') as phase: 
    
        orbital_value = phase.readlines()
    phase.close()
    for ii in range(len(orbital_value)):
        orbital_value[ii] = orbital_value[ii].split()
    #for ii in range(len(orbital_value)):
    #    print(orbital_value[ii])
    #print(orbital_value[3])

    C_temp = orbital_value[-len_list*2:]
    #print(C_temp)

    for jj in range(len(C_temp)):
        
            
        C_temp[jj] = C_temp[jj][-1]        
    

    


    for jj in range(len_list*2):
        if(jj%2==0):
            C_list_single[jj//2] = C_list_single[jj//2] + abs(float(C_temp[jj]))**1/num_iter
        else:
            C_list[jj//2] = C_list[jj//2] + float(C_temp[jj])/60/2

C_list_single = C_list_single**2

bounds = ([0, 0, 0, 0],[np.inf, np.pi, np.inf, np.inf])
a4, b4, c4, xi4= curve_fit(custom_func, L_list, C_list,maxfev=20000)[0]
x4 = np.arange(1,len_list, 0.01)
y4 = a4*np.sin(b4*(x4))/(x4**1)*np.exp(-x4/xi4)
print(a4,b4,c4,xi4)
print("local_length:",xi4)
print("oscilation_frequency:", b4)




bounds = ([0, 0, 0, 0],[np.inf, np.inf, np.inf, np.inf])
a1, b1, c1, d1, f1, xi1= curve_fit(custom_func_single, L_list[1:], C_list_single[1:],maxfev=20000)[0]
x1 = np.arange(2.0,len_list, 0.01)
y1 = a1*np.exp(-b1*x1)  + d1*np.exp(-f1*x1)*np.cos(pi*x1)


print(a1,1/b1,c1,d1,1/f1,xi1)
print("single_local_length:",1/b1)


    
#plt.figure(figsize=(8, 6))
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
plt.scatter(L_list[0:], C_list_single[0:]*1, label='C(r)',color = 'red')



plt.plot(x1, y1*1, 'k--', label=f'Fit: ξ = {(1/b1):.2f}')


plt.yticks([2*1E-5,4*1E-5,6*1E-5,8*1E-5,10*1E-5])
#plt.yticks([0.2*1E-4*100000,0.6*1E-4*100000,1.0*1E-4*100000])
ax = plt.gca()
ax.set_yticks([2*1E-5,4*1E-5,6*1E-5,8*1E-5,10*1E-5])

ax.set_yticklabels([f'{v/1e-5:.0f}' for v in [2*1E-5,4*1E-5,6*1E-5,8*1E-5,10*1E-5]])

ax.text(0.0, 1.02, '1e-5', transform=ax.transAxes,
        ha='left', va='bottom')

leg = ax.legend(frameon=True, prop={'size':12.6}, loc='upper right')
leg.get_frame().set_edgecolor('none')
leg.get_frame().set_linewidth(0.0)



plt.savefig('C_list_single_V0.9_W1.svg')




print("C_list_single:",list(C_list_single))
print("C_list",list(C_list))

#
