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
                      flim=(100,10000),
                      Nx = 21,
                      Ny = 21,
                      xlim=(-0.5,0.5),
                      ylim=(-0.5,0.5),
                      CSM=csm_test)

    @timeit "compute steeringvec" steeringvectors!(env)
    @timeit "compute beamforming" s,p = findmax(reshape(beamforming(env)[:,10],21,21))
    @test ceil(SPL(s)) == 47
    @test p.I == (10,13)
    @timeit "compute psf" s,p = findmax(reshape(psf(env)[:,10],21,21))
    @test ceil(SPL(sqrt(2).*s)) == 94
    @test p.I == (11,11)
end
