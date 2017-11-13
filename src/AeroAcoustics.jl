#__precompile__()
module AeroAcoustics
using NLsolve, Distances, DSP, JuMP, SCS, HDF5

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
       csm,
       diagrm!,
       fftnnls,
       shear,
       fistalasso

include("types.jl")
include("utils.jl")
include("csm.jl")
include("beamformer.jl")
include("pointspreadfunction.jl")
#include("beamformer_old.jl")
include("beamformer_corr.jl")
include("cmf.jl")
include("fista.jl")
include("fftnnls.jl")

end
