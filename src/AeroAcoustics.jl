#__precompile__()
module AeroAcoustics
using NLsolve,
using Distances
using DSP
using JuMP
using SCS
using HDF5
using ProximalOperators

import Base.length,
       Base.push!,
       Base.Threads

export Constants,
       CrossSpectralMatrix,
       Environment,
       SteeringMatrix,
       cmf,
       SPL,
       octavebandlimits,
       beamformer,
       beamformer_corr,
       beamformersetup,
       parseHDF5data,
       steeringvectors,
       pointspreadfunction,
       sourceintegration,
       fista,
       fistaprox!,
       nonneg!,
       csm,
#       diagrm!,
       fftnnls,
       shear

include("types.jl")
include("utils.jl")
include("steeringvectors.jl")
include("csm.jl")
include("beamformer.jl")
include("pointspreadfunction.jl")
#include("beamformer_corr.jl")
include("cmf.jl")
include("fista.jl")
include("fistaprox.jl")
include("fftnnls.jl")

end
