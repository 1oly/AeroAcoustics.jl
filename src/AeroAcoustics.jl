#__precompile__()
module AeroAcoustics
using NLsolve, Distances, DSP, JuMP, SCS, HDF5, ProximalOperators

import Base.length,
       Base.push!,
       Base.Threads
       #ImageFiltering.imfilter,
       #ImageFiltering.Fill

export Constants,
       CrossSpectralMatrix,
       Environment,
       SteeringMatrix,
       cmf,
       SPL,
       octavebandlimits,
       beamformer,
       beamformer_corr,
       #beamformer_old,
       beamformersetup,
       parseHDF5data,
       steeringvectors,
       pointspreadfunction,
       sourceintegration,
       fista,
       fistaprox!,
       csm,
       diagrm!,
       fftnnls,
       shear,
       fistalasso

include("types.jl")
include("utils.jl")
include("steeringvectors.jl")
include("csm.jl")
include("beamformer.jl")
include("pointspreadfunction.jl")
#include("beamformer_old.jl")
include("beamformer_corr.jl")
include("cmf.jl")
include("fista.jl")
include("fistaprox.jl")
include("fftnnls.jl")

end
