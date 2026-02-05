

using ITensors, ITensorMPS
using Random
using HDF5

let

    t1 = 1.0
    U = 4.0
    W = 2.3
    N1 = 60
    N2 = 3
    N = N1 * N2
    num_L = 20

    sites = siteinds("Electron", N; conserve_qns=true)

    
    f0 = h5open(string("WF_half",string(4000),".h5") ,"r")
    psi0 = read(f0,"T",MPS)
    close(f0)

    f1 = h5open(string("WF_half_m1",string(4000),".h5") ,"r")
    psi1 = read(f1,"T",MPS)
    close(f1)

    sites0 = siteinds(psi0)
    sites1 = siteinds(psi1)

    println("Norm difference: ", norm(psi0))
    println("Norm difference: ", norm(psi1))
    
    C_list = zeros(N2, N2, num_L)
    for x in 1:20
        for i1 in 1:N2, i2 in 1:N2
            println(i1," ",i2)
            sum_val = 0.0
            for jj in 10:(N1 - 20 - 10)
                site_i = (jj-1)*N2 + i1
                site_j = (jj+x-1)*N2 + i2
                

                c_i_dagger_c_j = OpSum()
                
                c_i_dagger_c_j += 1, "Cdagdn", (jj-1)*N2 + i1, "Cdn", (jj+x-1)*N2 + i2

                c_i_dagger_c_j += 1, "Cdagup", (jj-1)*N2 + i1, "Cup", (jj+x-1)*N2 + i2
                zz0=MPO(c_i_dagger_c_j, sites0)
                zz1=MPO(c_i_dagger_c_j, sites1)
                # 计算期望值
                val0 = inner(psi0', zz0, psi0)
                val1 = inner(psi1', zz1, psi1)
                sum_val += abs(val0-val1)

                
            end
            C_list[i1, i2, x] = sum_val / (N1 - 30 - 20 )  # 修正除数
        end
    end
    # 输出结果
    for x in 1:num_L, i1 in 1:N2, i2 in 1:N2
        println(i1, " ", i2, " ", x, " ", C_list[i1, i2, x])
    end

end


