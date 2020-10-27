SPL(p::Missing) = NaN
SPL(p::Number) = p > 0.0 ? 10*log10(p/4e-10) : NaN

"""
    octavebandlimits(fc,n)

Compute 1/n octave band limits given center frequencies fc.
"""
function octavebandlimits(fc::Vector{T},n) where T <: Real
    fl = fu = similar(fc)
    C = 10^(3. /10.)
    fl = fc*C^(-1/(2*n))
    fu = fc*C^(1/(2*n))
    return fl, fu
end

function octavebandlimits(fc::T,n) where T <: Real
    C = 10^(3. /10.)
    fl = fc*C^(-1/(2*n))
    fu = fc*C^(1/(2*n))
    return fl, fu
end

selectfreqs(f,flim) = (f.>=flim[1]) .& (f.<=flim[2])

function octavebands_nomial(n)
    if n == 1
        return [16, 31.5, 63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0, 16000.0, 32000.0]
    elseif n == 3
        return [12.5, 16, 20, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0,
            250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0,1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0,
            5000.0, 6300.0, 8000.0, 10000.0, 12500.0, 16000.0, 20000.0, 25000.0, 31500.0, 40000.0]
    else
        error("Value of `n` not supported. Try `n = 1` or `n = 3`")
    end
end

"""
    octavebands(n,flim=(16.,40000.),nomial::Bool=false)

Compute 1/n octave band center frequencies.
"""
function octavebands(n,flim=(16.,40000.);nomial::Bool=false)
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
    coherence_weights(CSM)

Calculate microphone weights according to mean mututal-coherence.

Amaral et al. “Slat Noise from an MD30P30N Airfoil at Extreme Angles of Attack,” AIAA Journal, vol. 56, no. 3, pp. 964–978, Mar. 2018.
"""
function coherence_weights(CSM)
    M,M,Nf = size(CSM)
    weight = zeros(M,Nf)
    for ω = 1:Nf
        for m = 1:M
            for n = 1:M
                if m != n
                    weight[m,ω] += abs(CSM[m,n,ω])^2 / abs(CSM[m,m,ω]*CSM[n,n,ω])
                end
            end
            weight[m,ω] /= M-1
        end
    end
    return weight
end

"""
    enbw(fs,win)

Equivalent Noise Band Width. Use to convert power spectrum to power spectral density and reverse. PS = PSD*ENBW
"""
enbw(fs,win) = fs*sum(abs2,win)/sum(win)^2
