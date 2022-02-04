using Test
using AeroAcoustics, HDF5
import DSP

@testset "Test 2: Monopole source without flow:" begin
    csm_ref = h5open("data/test1_csm.h5", "r") do file
        read(file, "CsmData/csmReal")+im*read(file, "CsmData/csmImag")
    end

    # Test csm
    fs = h5readattr("data/test1_csm.h5", "CsmData")["fs"]
    n = h5readattr("data/test1_csm.h5", "CsmData")["n"]

    # Setup beamforming
    fc = h5read("data/test1_csm.h5", "CsmData")["binCenterFrequenciesHz"]
    z0 = h5readattr("data/test1_csm.h5", "CsmData")["z0"]
    micgeom = h5read("data/test1_csm.h5", "CsmData")["arrayGeom"]
    @test size(micgeom) == (3,84)
    weights_vector = ones(84)
    weights_vector[4:34] .= 0

    env = Environment(z0=z0,
                      micgeom=micgeom,
                      flim=(3000,4000),
                      wv=weights_vector,
                      Nx = 21,
                      Ny = 21,
                      xlim=(-0.5,0.5),
                      ylim=(-0.5,0.5),
                      CSM=FreqArray(csm_ref,fc))

    steeringvectors!(env)
    b = beamforming(env)
    bd = beamforming(env.CSM_s.arr,env.steeringvec.arr)
    @test b.arr ≈ bd
    idx = 1 # Frequency index
    s1,p1 = findmax(reshape(b[:,idx],env.Nx,env.Ny))
    bmax = ceil(SPL(s1./sqrt(2)))
    @test bmax == 49.0
    @test p1.I == (10,13) # (19,24) for n = 41
    p_1 = psf(env)[:,idx]
    s2,p2 = findmax(reshape(p_1,env.Nx,env.Ny))
    @test ceil(SPL(s2./sqrt(2))) == 94
    @test p2.I == (floor(Int,env.Nx/2)+1,floor(Int,env.Ny/2)+1)
    pcol_1 = zeros(env.N)
    AeroAcoustics.psf_col!(pcol_1,env.steeringvec.arr[:,:,idx],floor(Int,env.N/2)+1)
    @test pcol_1 ≈ p_1
    # DAMAS
    x = damas(env, b; maxiter = 10)
    x2 = damas(env, b,[b.fc[idx]]; maxiter = 10)
    @test x[:,idx] ≈ x2
    id1,id2 = UnitRange.(p1.I.-2,p1.I.+2)
    limits = [env.rx[id1][1],env.rx[id1][end],env.ry[id2][1],env.ry[id2][end]]
    @test abs.(bmax-SPL.(sourceintegration(x[:,idx],env,limits))) <= 1 # Within 1dB of beamforming is OK (increase number of iterations to get better estimate)
    #CLEAN-SC
    x_clean = cleanSC(env;maxiter=10,ϕ=0.5)
    @test abs.(bmax-SPL.(sourceintegration(x_clean[:,idx],env,limits)/sqrt(2))) <= 1
end
