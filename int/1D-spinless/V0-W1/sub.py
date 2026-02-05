import os



current_path = os.getcwd()

path =  os.path.join(current_path,  'command_sub.dat')

#W = 2
num_iter = 200
with open(path, 'w') as file:
    for ii in range(500, num_iter+500):
        file.write('mkdir '+ str(ii) + '\n')
        file.write('sleep 0.1' + '\n')
        file.write('cp anderson.slum '+'./'+str(ii)+'/.' + '\n')
        file.write('sleep 0.1' + '\n')
        file.write('cp 1d_hubbard_ninj.jl '+'./'+str(ii)+'/.' + '\n')
        file.write('sleep 0.1' + '\n')
        file.write('cd '+ str(ii) + '\n')
        file.write('sleep 0.1' + '\n')
        file.write('sbatch anderson.slum' + '\n')
        file.write('sleep 0.1' + '\n')
        file.write('cd ..'+'\n')
        file.write('sleep 0.1' + '\n')
file.close()
