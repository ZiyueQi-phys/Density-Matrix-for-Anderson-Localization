import os



current_path = os.getcwd()

path =  os.path.join(current_path,  'command_sub.dat')

#W = 2
num_iter = 300
with open(path, 'w') as file:
    for ii in range(500, num_iter+500):
        file.write('mkdir '+ str(ii) + '\n')
        file.write('sleep 0.01' + '\n')
        file.write('cp anderson.slum '+'./'+str(ii)+'/.' + '\n')
        file.write('sleep 0.01' + '\n')
        file.write('cp 1d_hubbard_ninj.jl '+'./'+str(ii)+'/.' + '\n')
        file.write('sleep 0.01' + '\n')
        #file.write('cp ../../U0/W1/'+ str(ii) +'/WF_half200.h5 '+'./'+str(ii)+'/WF_init0.h5' + '\n')
        #file.write('cp ../../U0/W1/'+ str(ii) +'/WF_half_m1200.h5 '+'./'+str(ii)+'/WF_init1.h5' + '\n')
        #file.write('cp ../../U0/W1/'+ str(ii) +'/out '+'./'+str(ii)+'/out1' + '\n')
        file.write('cd '+ str(ii) + '\n')
        file.write('sleep 0.01' + '\n')
        file.write('sbatch anderson.slum' + '\n')
        file.write('sleep 0.01' + '\n')
        file.write('cd ..'+'\n')
        file.write('sleep 0.01' + '\n')
file.close()

