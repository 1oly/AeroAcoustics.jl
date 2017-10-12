__precompile__()
module AeroAcoustics
using NLsolve, Distances

import Base.length,
       Base.push!,
       Base.Threads
       #ImageFiltering.imfilter,
       #ImageFiltering.Fill

export Constants,
       Environment,
       SteeringMatrix,
       cmf,
       SPL,
       octavebandlimits,
       beamformer,
       beamformer_corr,
       beamformer2,
       beamformersetup,
       steeringvectors,
       pointspreadfunction,
       fista,
       fftnnls,
       shear,
       fistalasso

include("types.jl")
include("utils.jl")
include("beamformer.jl")
include("beamformer2.jl")
include("beamformer_corr.jl")
include("cmf.jl")
include("fista.jl")
include("fftnnls.jl")

# precompile beamformer
#precompile(beamformer2, (Environment{Float64},Constants{Float64},SteeringMatrix{Float64}));
end
