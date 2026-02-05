import os
current_path = os.getcwd()

directory_list = [f for f in os.listdir('.') if os.path.isdir(f)]
path = os.path.join(current_path, 'command_sub.dat')
with open(path, 'w') as file:
    
    for ii in range(len(directory_list)):
        file.write('cd '+ directory_list[ii] + '\n')
        for jj in range(10):
            file.write('mkdir ' + str(jj) + '\n')
            file.write('cd ' + str(jj) + '\n')
            file.write('cp ../../anderson.slum .' + '\n')
            file.write('dos2unix anderson.slum' + '\n')
            #file.write('cp ../WF_half5500.h5 .' + '\n')
            #file.write('cp ../WF_half_m15500.h5 .' + '\n')
            file.write('cp ../Transfer_Matrix_3D_length.py .' + '\n')

            file.write('sbatch anderson.slum' + '\n')
            file.write('cd ..' + '\n')
        file.write('cd ..' + '\n')
    
file.close()

#for ii in range(len(directory_list)):

import os

