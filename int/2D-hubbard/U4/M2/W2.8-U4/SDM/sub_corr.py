import os

current_path = os.getcwd()

path =  os.path.join(current_path,  'command_sub.dat')

#W = 2
num_iter = 100
with open(path, 'w') as file:
    for ii in range(0, num_iter+0):
        #if(ii!=47 and ii!=48 and ii!=52 and ii!=53 and ii!=56 and ii!=57 and ii!=58 and ii!=59):
        #    continue
        file.write('cd ' + str(ii) + '\n')
        #file.write('rm -r corr_series' + '\n')
        #file.write('mkdir SDM' + '\n')
        #file.write('cd SDM' + '\n')
        file.write('cp ../calc_correlation.jl .' + '\n')
        file.write('sleep 0.01' + '\n')
        file.write('cp ../anderson_corr.slum .' + '\n')
        file.write('sleep 0.01' + '\n')
        file.write('cp ../../GS_N_N-1/WF_half2500.h5 .' + '\n')
        file.write('sleep 0.01' + '\n')
        file.write('cp ../../GS_N_N-1/WF_half_m12500.h5 .' + '\n')
        file.write('sleep 0.01' + '\n')
        file.write('sbatch anderson_corr.slum' + '\n')
        file.write('cd ..' + '\n')
        file.write('cd ..' + '\n')
        print(ii)
file.close()

