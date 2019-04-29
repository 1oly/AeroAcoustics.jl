using Test
using AeroAcoustics, HDF5
import DSP

@testset "Test 1: Monopole source without flow:" begin
    csm_ref = h5open("data/test1_csm.h5", "r") do file
        read(file, "CsmData/csmReal")+im*read(file, "CsmData/csmImag")
    end
    t = h5open("data/test1_timeseries.h5", "r") do file
        read(file, "MicrophoneData/microphoneDataPa")
    end

    # Test csm
    fs = h5readattr("data/test1_csm.h5", "CsmData")["fs"]
    n = h5readattr("data/test1_csm.h5", "CsmData")["n"]
    @timeit "compute csm" csm_test = csm(t;n=n,noverlap=div(n,2),fs=fs,win=DSP.hanning(n))
    @test csm_test.arr ≈ csm_ref

    # Setup beamforming
    fc2 = h5read("data/test1_csm.h5", "CsmData")["binCenterFrequenciesHz"]
    @test csm_test.fc==fc2
    z0 = h5readattr("data/test1_csm.h5", "CsmData")["z0"]
    micgeom = h5read("data/test1_csm.h5", "CsmData")["arrayGeom"]
    @test size(micgeom) == (3,84)

    @timeit "compute env" env = Environment(z0=z0,
                      micgeom=micgeom,
                      flim=(3000,4000),
                      Nx = 21,
                      Ny = 21,
                      xlim=(-0.5,0.5),
                      ylim=(-0.5,0.5),
                      CSM=csm_test)

    @timeit "compute steeringvec" steeringvectors!(env)
    @timeit "compute beamforming" b = beamforming(env)
    idx = 1 # Frequency index
    s1,p1 = findmax(reshape(b[:,idx],env.Nx,env.Ny))
    bmax = ceil(SPL(sqrt(2).*s1))
    @test bmax == 52.0
    @test p1.I == (10,13) # (19,24) for n = 41
    @timeit "compute psf" p_1 = psf(env)[:,idx]
    s2,p2 = findmax(reshape(p_1,env.Nx,env.Ny))
    @test ceil(SPL(sqrt(2).*s2)) == 94
    @test p2.I == (floor(Int,env.Nx/2)+1,floor(Int,env.Ny/2)+1)
    pcol_1 = zeros(env.N)
    @timeit "compute psf_col" AeroAcoustics.psf_col!(pcol_1,env.steeringvec.arr[:,:,idx],floor(Int,env.N/2)+1)
    @test pcol_1 ≈ p_1
    # DAMAS
    x = zeros(size(b))
    @timeit "compute DAMAS" damas!(x, env, b; maxiter = 10)
    id1,id2 = UnitRange.(p1.I.-2,p1.I.+2)
    limits = [env.rx[id1][1],env.rx[id1][end],env.ry[id2][1],env.ry[id2][end]]
    @test abs.(bmax-SPL.(sourceintegration(x[:,idx],env,limits))) <= 1 # Within 1dB of beamforming is OK (increase number of iterations to get better estimate)
end
