#__precompile__()
module AeroAcoustics
using ImageFiltering, NLsolve

import Base.length,
       Base.push!,
       Base.Threads,
       ImageFiltering.imfilter,
       ImageFiltering.Fill

export cmf,
       SPL,
       octavebandlimits,
       beamformer,
       beamformer_corr,
       fista,
       shear

include("beamformer.jl")
include("beamformer_corr.jl")
include("cmf.jl")
include("utils.jl")
include("fista.jl")

end
