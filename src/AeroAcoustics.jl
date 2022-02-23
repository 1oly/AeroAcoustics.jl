module AeroAcoustics
using Distances
using LinearAlgebra
using Statistics
using NLsolve
using Parameters
using LazyArrays
using FFTW
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
       cleanSC,
       sourceintegration,
       octavebandlimits,
       octavebands,
       narrow2oct,
       enbw,
       coherence_weights,
       DR!,
       DR,
       SPI,
       fista


include("utils.jl")
include("types.jl")
include("steeringvectors.jl")
include("csm.jl")
include("beamforming.jl")
include("psf.jl")
include("damas.jl")
include("cleansc.jl")
include("sourceintegration.jl")
include("fista.jl")

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

"""
    DR!(CSM::FreqArray)
    
Apply diagonal removal to CSM.
"""
function DR!(CSM::FreqArray)
    T = eltype(CSM)
    M,M,Nf = size(CSM)
    for j = 1:Nf
        for i = 1:M
            CSM.arr[i,i,j] = zero(T)
        end
    end
end

"""
    DR(CSM::FreqArray)
    
Apply diagonal removal to CSM and new array.
"""
function DR(CSM::FreqArray)
    T = eltype(CSM)
    M,M,Nf = size(CSM)
    CSMd = copy(CSM)
    for j = 1:Nf
        for i = 1:M
            CSMd[i,i,j] = zero(T)
        end
    end
    return FreqArray(CSMd,CSM.fc)
end

end