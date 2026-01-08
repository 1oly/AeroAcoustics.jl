using Test
using AeroAcoustics
using HDF5
using DelimitedFiles

const DATA_DIR = joinpath(@__DIR__, "data")
const CSM_FILE = joinpath(
    DATA_DIR,
    "DTU_PLCT_NACA63018_trip_5PS_5SS_U0_50_AoA_0_octave-12_CsmEss.h5",
)
const SRCINT_REF_FILE = joinpath(DATA_DIR, "task6_DTU_AeroAcousticsjl_CBF_srcint.csv")

function read_dtu_csm(path)
    h5open(path, "r") do file
        csm_data = read(file, "CsmData")
        fc = csm_data["binCenterFrequenciesHz"]
        csm = csm_data["csmReal"] + im * csm_data["csmImaginary"]
        # HDF5 can flip dimensions depending on storage order.
        if size(csm, 1) != size(csm, 2)
            csm = permutedims(csm, reverse(1:ndims(csm)))
        end
        meta = read(file, "MetaData")
        micgeom = meta["ArrayAttributes"]["microphonePositionsM"]
        if size(micgeom, 1) != 3
            micgeom = permutedims(micgeom)
        end
        return FreqArray(csm, fc), micgeom
    end
end

function read_z0(test_attrs)
    bounds = test_attrs["domainBoundsM"]
    if bounds isa AbstractMatrix
        return bounds[1, 3]
    end
    if bounds isa AbstractVector
        first_entry = bounds[1]
        if first_entry isa NamedTuple
            return first_entry[Symbol("3")]
        end
        return first_entry[3]
    end
    error("Unsupported domainBoundsM format: $(typeof(bounds))")
end

function read_mach(test_attrs)
    ma = test_attrs["machNumber"]
    if ma isa AbstractVector
        return ma[1]
    end
    return ma
end

@testset "Benchmark: DTU airfoil CBF" begin
    @test isfile(CSM_FILE)
    @test isfile(SRCINT_REF_FILE)

    csm, micgeom = read_dtu_csm(CSM_FILE)
    test_attrs = h5readattr(CSM_FILE, "MetaData/TestAttributes")
    measurement_data = h5readattr(CSM_FILE, "MeasurementData")
    z0 = read_z0(test_attrs)
    ma = read_mach(test_attrs)
    c = measurement_data["speedOfSoundMPerS"]
    if c isa AbstractVector
        c = c[1]
    end

    freq_targets = [500.0, 1000.0, 2000.0, 4000.0]
    freq_idx = [argmin(abs.(csm.fc .- f)) for f in freq_targets]
    csm_sel = FreqArray(csm.arr[:, :, freq_idx], csm.fc[freq_idx])

    env_flow = Environment(
        CSM = csm_sel,
        micgeom = micgeom,
        z0 = z0,
        dr = true,
        shear = true,
        Ma = -ma,
        c = c,
        h = 1.5,
        Nx = 81,
        Ny = 41,
        xlim = (-2, 2),
        ylim = (-1, 1),
        flim = (400, 4000),
    )
    steeringvectors!(env_flow)
    b_flow = beamforming(env_flow)

    srcint_ref = readdlm(SRCINT_REF_FILE, ';', Float64)
    srcint_freqs = srcint_ref[:, 1]
    srcint_vals = srcint_ref[:, 2]

    limits = [-0.5, 0.5, -0.4, 0.4]
    srcint_cbf = sourceintegration(b_flow, env_flow, limits)
    for (f, val) in zip(b_flow.fc, SPL.(srcint_cbf.arr))
        ref_idx = argmin(abs.(srcint_freqs .- f))
        @test abs.(val - srcint_vals[ref_idx]) <= 1.0
    end
end
