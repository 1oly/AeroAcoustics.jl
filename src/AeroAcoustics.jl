module AeroAcoustics
using Distances
using LinearAlgebra # LinAlg in julia 1.0
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
       narrow2oct


include("types.jl")
include("steeringvectors.jl")
include("csm.jl")
include("beamforming.jl")
include("psf.jl")
include("damas.jl")
include("sourceintegration.jl")

function SPL(p::Array{T}) where T <: Real
    s = similar(p)
    s[p.>0] = 10*log10.(p[p.>0]/4e-10)
    s[p.<=0] = -350
    return s
end
SPL(p::Number) = p > 0.0 ? 10*log10(p/4e-10) : -350.0

"""
    octavebandlimits(fc,n)

Compute `1/n` octave band limits given center frequencies fc.
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

function octavebands(n,flim=extrema(octavebands_nomial(n)))
    fc = octavebands_nomial(n)
    f_inds = selectfreqs(fc,flim)
    return fc[f_inds]
end

"""
    narrow2oct(x::FreqArray,fc,n)

Sum narrow band spectrum to 1/n octave bands given narrow band frequencies `f` in x
and octave band center frequencies `fc`.
"""
function narrow2oct(x::FreqArray,n)
    fc = octavebands(n,extrema(x.fc))
    fl, fu = octavebandlimits(fc[f_inds],n)
    out = similar(fc[f_inds])
    for i in 1:length(out)
        inds = selectfreqs(x.fc,(fl[i],fu[i]))
        out[i] = sum(x[:,inds])/(fu[i]-fl[i])
    end
    return FreqArray(out,fc[f_inds])
end

end
