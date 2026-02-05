using ITensors, ITensorMPS
using Random
using HDF5

let
    t1 = 1.0
    V = 0.9
    W = 1


    num_L = 30

    N = 120


    filepath = "out1"
    # 初始化数组和计数器
    W_list = Float64[]
    line_count = 0
    max_lines = N

    open(filepath, "r") do file
        for line in eachline(file)
            line_count += 1
            if line_count > max_lines
                break
            end
            # 分割行并提取最后一列
            fields = split(line,":")
            push!(W_list, parse(Float64, fields[end]))
        end
    end


    expectation_value = 0
    sites = siteinds("Electron", N; conserve_qns=true)
    #sites = siteinds("Electron", N)

    os = OpSum()
    for b in 1:(N - 1)
        os -= t1, "Cdagup", b, "Cup", b + 1
        os -= t1, "Cdagup", b + 1, "Cup", b
        os -= t1, "Cdagdn", b, "Cdn", b + 1
        os -= t1, "Cdagdn", b + 1, "Cdn", b

    end

    for i in 1:N
        #mu = -W+2*W*rand()
        mu = W_list[i]
        println("mu",i,":",mu)
        os += mu, "Cdagup", i, "Cup", i
        os += mu, "Cdagdn", i, "Cdn", i
    end

    for i in 1:(N - 1)
        os -= V, "Nup", i, "Nup" ,i+1
    end


    H = MPO(os, sites)

    # DMRG参数
    nsweeps = 100
    maxdim1 = [ 30, 50, 100, 200, 400 ]
    cutoff = [1E-12]
    noise = [1E-5]

    # 基态1: Npart = N
    state0 = fill("Emp", N)
    
    p = 60
    for i in N:-1:1
        if p > i
            state0[i] = "UpDn"
            p -= 2
        elseif p > 0
            state0[i] = isodd(i) ? "Up" : "Up"
            p -= 1
        end
    end

    psi0_init = random_mps(sites, state0; linkdims=10)
    println("sites:", typeof(sites))
    @show flux(psi0_init)

    for ii in 1:size(maxdim1)[1]
        if ii > 4
            nsweeps = 10
            cutoff = [1E-15]
        end
        if ii > 6
            nsweeps = 30
            cutoff = [1E-18]
        end
        maxdim = [maxdim1[ii], maxdim1[ii], maxdim1[ii]]
        print("maxdim:",maxdim)
        energy0, psi0 = dmrg(H, psi0_init; nsweeps, maxdim, cutoff, noise)
        if ii > 1
            f1 = h5open(string("WF_half",string(maxdim1[ii]),".h5") ,"w")
            write(f1,"T",psi0)
            close(f1)
        end
        psi0_init = psi0
    end

    


    # 基态1: Npart = N-1
    state0 = fill("Emp", N)
    p = 59
    nsweeps = 100
    cutoff = [1E-12]
    for i in N:-1:1
        if p > i
            state0[i] = "UpDn"
            p -= 2
        elseif p > 0
            state0[i] = isodd(i) ? "Up" : "Up"
            p -= 1
        end
    end
    
    psi0_init = random_mps(sites, state0; linkdims=10)
    println("sites:", typeof(sites))
    @show flux(psi0_init)

    for ii in 1:size(maxdim1)[1]
        if ii > 4
            cutoff = [1E-15]
            nsweeps = 10
        end

        if ii > 6
            nsweeps = 30
            cutoff = [1E-18]
        end
        maxdim = [maxdim1[ii], maxdim1[ii], maxdim1[ii]]
        print("maxdim:",maxdim)
        energy0, psi0 = dmrg(H, psi0_init; nsweeps, maxdim, cutoff, noise)
        if ii > 1
            f1 = h5open(string("WF_half_m1",string(maxdim1[ii]),".h5") ,"w")
            write(f1,"T",psi0)
            close(f1)
        end
        psi0_init = psi0
    end


    

    
end
