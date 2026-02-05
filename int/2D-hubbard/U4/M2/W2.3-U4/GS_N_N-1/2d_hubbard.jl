using ITensors, ITensorMPS
using Random
using HDF5

let
    t1 = 1.0
    U = 4.0
    W = 2.3
    N1 = 60
    N2 = 2
    N = N1 * N2
    num_L = 30

    # 生成随机化学势并构建Hamiltonian
    sites = siteinds("Electron", N; conserve_qns=true)
    os = OpSum()
    # 跃迁项
    for a in 1:N1
        for b in 1:(N2-1)
            i = (a-1)*N2 + b
            j = i + 1
            os -= t1, "Cdagup", i, "Cup", j
            os -= t1, "Cdagup", j, "Cup", i
            os -= t1, "Cdagdn", i, "Cdn", j
            os -= t1, "Cdagdn", j, "Cdn", i
        end
    end
    for a in 1:(N1-1)
        for b in 1:N2
            i = (a-1)*N2 + b
            j = a*N2 + b
            os -= t1, "Cdagup", i, "Cup", j
            os -= t1, "Cdagup", j, "Cup", i
            os -= t1, "Cdagdn", i, "Cdn", j
            os -= t1, "Cdagdn", j, "Cdn", i
        end
    end
    # 随机化学势和U项
    #mus = [-W + 2 * W * rand() for _ in 1:N]
    
    for i in 1:N
        os += U, "Nupdn", i
    end

    for ii in 1:N1
        for jj in 1:N2
            #mu = cos( (1/(sqrt(5)-1)) * ((ii-1)*N2 + jj-1)  )
            
            nn = (ii-1)*N2+jj
            mu = -W+2*W*rand()
            println(ii," ",jj," ",mu)
            os += mu, "Cdagup", nn, "Cup", nn
            os += mu, "Cdagdn", nn, "Cdn", nn
        end
    end


    H = MPO(os, sites)

    # DMRG参数
    nsweeps = 25
    maxdim1 = [50, 100, 200, 400, 800, 1000, 1500,2000,2500]
    cutoff = [1E-12]
    noise = [1E-5]

    # 基态1: Npart = N
    state0 = fill("Emp", N)
    p = 64
    for i in N:-1:1
        if p > i
            state0[i] = "UpDn"
            p -= 2
        elseif p > 0
            state0[i] = isodd(i) ? "Up" : "Dn"
            p -= 1
        end
    end
    psi0_init = random_mps(sites, state0; linkdims=10)

    for ii in 1:size(maxdim1)[1]
        if ii > 6
            nsweeps = 5
        end
        maxdim = [maxdim1[ii], maxdim1[ii], maxdim1[ii]]
        print("maxdim:",maxdim)
        energy0, psi0 = dmrg(H, psi0_init; nsweeps, maxdim, cutoff, noise)
        
        f1 = h5open(string("WF_half",string(maxdim1[ii]),".h5") ,"w")
        write(f1,"T",psi0)
        close(f1)

        psi0_init = psi0
    end


    # 基态1: Npart = N-1
    state0 = fill("Emp", N)
    p = 63
    nsweeps = 25
    for i in N:-1:1
        if p > i
            state0[i] = "UpDn"
            p -= 2
        elseif p > 0
            state0[i] = isodd(i) ? "Up" : "Dn"
            p -= 1
        end
    end
    psi0_init = random_mps(sites, state0; linkdims=10)

    for ii in 1:size(maxdim1)[1]
        if ii > 6
            nsweeps = 5
        end
        maxdim = [maxdim1[ii], maxdim1[ii], maxdim1[ii]]
        print("maxdim:",maxdim)
        energy0, psi0 = dmrg(H, psi0_init; nsweeps, maxdim, cutoff, noise)
        f1 = h5open(string("WF_half_m1",string(maxdim1[ii]),".h5") ,"w")
        write(f1,"T",psi0)
        close(f1)

        psi0_init = psi0
    end


    

    
end
