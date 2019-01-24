#__precompile__()
module AeroAcoustics
using Distances
using LinearAlgebra # LinAlg in julia 1.0
using Statistics
using Parameters
import DSP

import Base.length,
       Base.push!,
       Base.Threads

export Environment,
       CrossSpectralMatrix,
       SPL,
       steeringvectors,
       csm,
       beamforming



include("types.jl")
include("steeringvectors.jl")
include("csm.jl")
include("beamforming.jl")

function SPL(p::Array{T}) where T <: Real
    s = similar(p)
    s[p.>0] = 10*log10.(p[p.>0]/4e-10)
    s[p.<=0] = -350
    return s
end
SPL(p::Number) = p > 0.0 ? 10*log10(p/4e-10) : -350.0

end
