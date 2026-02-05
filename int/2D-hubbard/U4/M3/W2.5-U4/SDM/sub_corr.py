import os

current_path = os.getcwd()

path =  os.path.join(current_path,  'command_sub.dat')

#W = 2
num_iter = 20
with open(path, 'w') as file:
    for ii in range(80, num_iter+80):
        if(ii==80 or ii==87 ):
            continue
        
        #file.write('rm -r corr_series' + '\n')
        file.write('mkdir ' + str(ii) + '\n')
        file.write('cd ' + str(ii) + '\n')

        #file.write('cd corr_series' + '\n')
        file.write('cp ../calc_correlation.jl .' + '\n')
        file.write('sleep 0.1' + '\n')
        file.write('cp ../anderson_corr.slum .' + '\n')
        file.write('sleep 0.1' + '\n')
        file.write('cp ../../p4/' + str(ii) + '/WF_half4000.h5 .' + '\n')
        file.write('sleep 1' + '\n')
        file.write('cp ../../p3/' + str(ii) + '/WF_half_m14000.h5 WF_half_m14000.h5' + '\n')
        file.write('sleep 1' + '\n')
        file.write('sbatch anderson_corr.slum' + '\n')
        file.write('cd ..' + '\n')
        #file.write('cd ..' + '\n')
        
file.close()
