module AeroAcoustics
using Distances
using LinearAlgebra
using Statistics
using NLsolve
using Parameters
import DSP

import Base.length,
       Base.push!,
       Base.Threads

export Environment,
       FreqArray,
       SPL,
       steeringvectors,
       steeringvectors!,
       csm,
       csm!,
       beamforming,
       psf,
       damas!,
       sourceintegration,
       octavebandlimits,
       octavebands,
       narrow2oct


include("types.jl")
include("steeringvectors.jl")
include("csm.jl")
include("beamforming.jl")
include("psf.jl")
include("damas.jl")
include("sourceintegration.jl")

SPL(p::Missing) = NaN
SPL(p::Number) = p > 0.0 ? 10*log10(p/4e-10) : NaN

"""
    octavebandlimits(fc,n)

Compute 1/n octave band limits given center frequencies fc.
"""
function octavebandlimits(fc,n)
    fl,fu = similar(fc)
    C = 10^(3. /10.)
    fl = fc*C^(-1/(2*n))
    fu = fc*C^(1/(2*n))
    return fl, fu
end

selectfreqs(f,flim) = (f.>=flim[1]) .& (f.<=flim[2])

function octavebands_nomial(n)
    if n == 1
        return [31.5, 63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0, 16000.0]
    elseif n == 3
        return [25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0,
            250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0,1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0,
            5000.0, 6300.0, 8000.0, 10000.0, 12500.0, 16000.0, 20000.0]
    else
        error("Value of `n` not supported. Try `n = 1` or `n = 3`")
    end
end

"""
    octavebands(n,flim=(25.,20000.),nomial::Bool=false)

Compute 1/n octave band center frequencies.
"""
function octavebands(n,flim=(25.,20000.),nomial::Bool=false)
    if nomial && (n==1 || n==3)
        fc = octavebands_nomial(n)
        f_inds = selectfreqs(fc,flim)
        return fc[f_inds]
    else
        fc3 = octavebands_nomial(3)
        smin,imin = findmin(abs.(fc3.-flim[1]))
        smax,imax = findmin(abs.(fc3.-flim[2]))
        (smin != 0) ? @warn("low frequency limit is not supported and has been replaced by $(fc3[imin])") : nothing
        (smax != 0) ? @warn("high frequency limit is not supported and has been replaced by $(fc3[imax])") : nothing
        flim = (fc3[imin],fc3[imax])
        G = 10^(3. /10.)
        fc = exp.(log(flim[1]):(log(G)/n):log(flim[2]))
        return fc
    end
end

"""
    narrow2oct(x::FreqArray,n,nomial::Bool=true)

Sum narrow band spectrum to 1/n octave bands given narrow band frequencies `f` in x.
"""
function narrow2oct(x::FreqArray,n,nomial::Bool=true)
    fc = octavebands(n,extrema(x.fc),nomial)
    fl, fu = octavebandlimits(fc,n)
    out = Array{Union{Missing, Float64},2}(undef,size(x,1),length(fc))
    for i in 1:length(fc)
        inds = selectfreqs(x.fc,(fl[i],fu[i]))
        if any(inds .== true)
            out[:,i] = sum(x[:,inds],dims=2)/(fu[i]-fl[i])
        else
            @warn "No data in frequency band $(fc[i]) replacing with missing"
            out[:,i] = Array{Missing,1}(undef,size(x,1))
        end
    end
    return FreqArray(out,fc)
end

end
