import os



current_path = os.getcwd()

path =  os.path.join(current_path,  'command_sub.dat')

#W = 2
num_iter = 100
with open(path, 'w') as file:
    for ii in range(00, num_iter+00):
        #file.write('mkdir '+ str(ii) + '\n')
        #file.write('sleep 0.1' + '\n')
        #file.write('cp anderson.slum '+'./'+str(ii)+'/.' + '\n')
        #file.write('sleep 0.1' + '\n')
        #file.write('cp 2d_hubbard.jl '+'./'+str(ii)+'/.' + '\n')
        #file.write('sleep 0.1' + '\n')
        file.write('cd '+ str(ii) + '\n')
        file.write('sleep 0.1' + '\n')
        #file.write('cp ../../p4/'+str(ii)+'/out out1' + '\n')
        file.write('sbatch anderson.slum' + '\n')
        file.write('sleep 0.1' + '\n')
        file.write('cd ..'+'\n')
        file.write('sleep 0.1' + '\n')
file.close()

