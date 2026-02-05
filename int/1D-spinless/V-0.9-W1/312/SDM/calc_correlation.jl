

using ITensors, ITensorMPS
using Random
using HDF5

let

    t1 = 1.0
    U = 0.0
    W = 1.0
    N1 = 120
    N = N1
    L = N
    num_L = 20

    sites = siteinds("Electron", N; conserve_qns=true)

    

    f0 = h5open(string("WF_half",string(200),".h5") ,"r")
    psi0 = read(f0,"T",MPS)
    close(f0)

    f1 = h5open(string("WF_half_m1",string(200),".h5") ,"r")
    psi1 = read(f1,"T",MPS)
    close(f1)

    sites0 = siteinds(psi0)
    sites1 = siteinds(psi1)

    println("Norm difference: ", norm(psi0))
    println("Norm difference: ", norm(psi1))
    


    







    for x in 1:20
        expectation_value = 0
        expectation_value_single = 0
        for i in 20:(N-40)

        
            c_i_dagger_c_j = OpSum()
            c_i_dagger_c_j += 1, "Cdagup", i, "Cup", (i+x)%N
            c_i_dagger_c_j += 1, "Cdagdn", i, "Cdn", (i+x)%N
            zz0=MPO(c_i_dagger_c_j, sites0)
            zz1=MPO(c_i_dagger_c_j, sites1)
            expectation_value_single = expectation_value_single + abs( (inner(psi0',zz0,psi0)) - inner(psi1',zz1,psi1) )  / (N - 20 -40)
            expectation_value = expectation_value + (inner(psi0',zz0,psi0)) / (N - 20 -40)
        end
        println("<c_i" , "^dagger c_i+", x, "_single> = ", expectation_value_single)
        println("<c_i" , "^dagger c_i+", x, "> = ", expectation_value)
    end

end


