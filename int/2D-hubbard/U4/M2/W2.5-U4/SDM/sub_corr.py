import os

current_path = os.getcwd()

path =  os.path.join(current_path,  'command_sub.dat')

#W = 2
num_iter = 50
with open(path, 'w') as file:
    for ii in range(80, num_iter+80):
        #if(ii==69):
        #    continue
        file.write('cd ' + str(ii) + '\n')
        #file.write('rm -r corr_series' + '\n')
        file.write('mkdir corr_series' + '\n')
        file.write('cd corr_series' + '\n')
        file.write('cp ../../calc_correlation.jl .' + '\n')
        file.write('sleep 0.1' + '\n')
        file.write('cp ../../anderson_corr.slum .' + '\n')
        file.write('sleep 0.1' + '\n')
        file.write('cp ../WF_half2000.h5 .' + '\n')
        file.write('sleep 1' + '\n')
        file.write('cp ../WF_half_m12000.h5 .' + '\n')
        file.write('sleep 1' + '\n')
        file.write('sbatch anderson_corr.slum' + '\n')
        file.write('cd ..' + '\n')
        file.write('cd ..' + '\n')
        
file.close()
