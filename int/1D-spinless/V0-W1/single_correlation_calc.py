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
    #if(ii==8 or ii==9 or ii==10 or ii==11 or ii==19 or (ii>=70 and ii<80)  or ii==2 or ii==5 or ii==3  or ii==42 or ii==43 or ii==61 or ii==65 or ii==69  or ii==80 or ii==90 or ii==92 or ii==93 or ii==97 or ii==98 or ii==99):
    #    continue
    #if(ii==2 or ii== 5 or ii== 8 or ii== 9 or ii== 10 or ii== 11 or ii== 17 or ii== 19 or ii== 21 or ii== 23 or ii== 28 or ii== 36 or ii== 40 or ii== 45 or ii== 50 or ii== 52 or ii== 57 or ii== 61 or ii== 63 or ii== 65 or ii== 68 or ii== 69 or ii== 70 or ii== 79 or ii== 88 or ii== 90 or ii== 91 or ii== 92 or ii== 93 or ii== 98 or ii== 99):
    #    continue
    #if( ii==42 or ii==43 or ii==97):
    #    continue
    #if(ii==8):
    #    continue
    path = file_path = os.path.join(current_path, str(ii),'SDM', 'out')
    with open(path, 'r') as phase: 
    #with open('E:\\DFT_calculation\\wannier_tools\\KIn\\sigma_ahc_eta10.00meV.txt', 'r') as phase: 
    
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
a1, b1, c1, d1, f1, xi1= curve_fit(custom_func_single, L_list, C_list_single,maxfev=20000)[0]
x1 = np.arange(1,len_list, 0.01)
y1 = a1*np.exp(-b1*x1) + + d1*np.exp(-f1*x1)*np.cos(pi*x1)
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
plt.scatter(L_list, C_list_single, label='C(r)',color = 'red')


C_list_ED =  [3.6267705978406536e-05, 5.414946631085121e-05, 3.692257363121452e-05, 4.64725914515008e-05, 3.6049832571726507e-05, 4.148000697922487e-05, 3.439432944655676e-05, 3.7486704321821114e-05, 3.259554572132218e-05, 3.403818417664318e-05, 3.0796974430167626e-05, 3.131196955586237e-05, 2.883202681214743e-05, 2.8743514729676216e-05, 2.6754234777816656e-05, 2.630365148881059e-05, 2.4935629132694357e-05, 2.4218540260528637e-05, 2.3246175559230494e-05, 2.2401162578414162e-05]
plt.scatter(L_list, C_list_ED, label = 'ED',color = 'blue',marker = '*')

plt.plot(x1, y1, 'k--', label=f'Fit: ξ = {(1/b1):.2f}')
#plt.xlabel('Distance r')
#plt.ylabel('C(r)')
#plt.legend()
#plt.title('Localization Length from Green\'s Function')
#plt.legend(prop={'size':12})

leg = plt.legend(frameon=True, prop={'size':12.0}, loc='upper right')
leg.get_frame().set_edgecolor('none')
leg.get_frame().set_linewidth(0.0)

plt.yticks([2*1E-5,3*1E-5,4*1E-5,5*1E-5])

plt.savefig('C_list_single_V0_W1.svg')



print("C_list_single:",list(C_list_single))
print("C_list",list(C_list))

#
