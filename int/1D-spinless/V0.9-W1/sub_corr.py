import os

current_path = os.getcwd()

path =  os.path.join(current_path,  'command_sub.dat')

#W = 2
num_iter = 300
with open(path, 'w') as file:
    for ii in range(500, num_iter+500):
        #if(ii==5):
        #    continue
        file.write('cd ' + str(ii) + '\n')
        #file.write('rm -r corr_series' + '\n')
        file.write('mkdir corr_series' + '\n')
        file.write('cd corr_series' + '\n')
        file.write('cp ../../calc_correlation.jl .' + '\n')
        file.write('sleep 0.01' + '\n')
        file.write('cp ../../anderson_corr.slum .' + '\n')
        file.write('sleep 0.01' + '\n')
        file.write('cp ../WF_half200.h5 .' + '\n')
        file.write('sleep 0.01' + '\n')
        file.write('cp ../WF_half_m1200.h5 .' + '\n')
        file.write('sleep 0.01' + '\n')
        file.write('sbatch anderson_corr.slum' + '\n')
        file.write('cd ..' + '\n')
        file.write('cd ..' + '\n')
        print(ii)
file.close()


