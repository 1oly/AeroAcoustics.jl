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
    csm_test,fc1 = csm(t;n=n,noverlap=div(n,2),fs=fs,win=DSP.hanning(n))
    @test csm_test ≈ csm_ref

    # Setup beamforming
    fc2 = h5read("data/test1_csm.h5", "CsmData")["binCenterFrequenciesHz"]
    @test fc1==fc2
    z0 = h5readattr("data/test1_csm.h5", "CsmData")["z0"]
    micgeom = h5read("data/test1_csm.h5", "CsmData")["arrayGeom"]
    @test size(micgeom) == (3,84)

    env = Environment(z0=z0,
                      micgeom=micgeom,
                      Nx = 21,
                      Ny = 21,
                      flim=(100,10000),
                      f=fc1,
                      xlim=(-0.5,0.5),
                      ylim=(-0.5,0.5))

    v1 = steeringvectors(env)
    v1c = convert(Array{Complex{Float32},3},v1)
    b = beamforming(csm_ref[:,:,1:39],v1c)
    s,p = findmax(reshape(b[:,10],21,21))
    @test ceil(SPL(s)) == 48
    @test p.I == (10,12)
    println("Hello from branch MASTER")
end
