module AeroAcoustics
using Distances
using LinearAlgebra
using Statistics
using NLsolve
using Parameters
using LazyArrays
import DSP

import Base.length,
       Base.push!,
       Base.reshape

export Environment,
       FreqArray,
       SPL,
       steeringvectors,
       steeringvectors!,
       csm,
       beamforming,
       psf,
       damas,
       sourceintegration,
       octavebandlimits,
       octavebands,
       narrow2oct,
       enbw,
       coherence_weights


include("utils.jl")
include("types.jl")
include("steeringvectors.jl")
include("csm.jl")
include("beamforming.jl")
include("psf.jl")
include("damas.jl")
include("sourceintegration.jl")


reshape(x::FreqArray,dims::Vararg{Int64,N}) where N = FreqArray(reshape(x.arr,dims...),x.fc)

"""
    narrow2oct(x::FreqArray,n,nomial::Bool=true;psd=false)

Sum narrow band spectrum to 1/n octave bands given narrow band frequencies `f` in x.
"""
function narrow2oct(x::FreqArray,n;nomial::Bool=true,psd::Bool=false)
    fc = octavebands(n,extrema(x.fc);nomial=nomial)
    fl, fu = octavebandlimits(fc,n)
    out = Array{Union{Missing, Float64},2}(undef,size(x,1),length(fc))
    for i in 1:length(fc)
        inds = selectfreqs(x.fc,(fl[i],fu[i]))
        if any(inds .== true)
            if psd
                out[:,i] = sum(x[:,inds],dims=2)/(x.fc[2]-x.fc[1])
            else
                out[:,i] = sum(x[:,inds],dims=2)
            end
        else
            @warn "No data in frequency band $(fc[i]) replacing with missing"
            out[:,i] = Array{Missing,1}(undef,size(x,1))
        end
    end
    return FreqArray(out,fc)
end

end