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
       damas!


include("types.jl")
include("steeringvectors.jl")
include("csm.jl")
include("beamforming.jl")
include("psf.jl")
include("damas.jl")

function SPL(p::Array{T}) where T <: Real
    s = similar(p)
    s[p.>0] = 10*log10.(p[p.>0]/4e-10)
    s[p.<=0] = -350
    return s
end
SPL(p::Number) = p > 0.0 ? 10*log10(p/4e-10) : -350.0

function shear!(fvec,xi,xm,Ma,h)
    a = sqrt(xi[1]^2+(1-Ma^2)*(xi[2]^2+xi[3]^2))
    b = norm(xm-xi)
    fvec[1] = (1/a)*xi[1]-(1/b)*(1-Ma^2)*(xm[1]-xi[1]) - Ma
    fvec[2] = (1/a)*xi[2]-(1/b)*(xm[2]-xi[2])
    fvec[3] = xi[3] - h
end

function propagation_time(E::Environment)
    @unpack N,M,micgeom,Rxy,Ma,h,c = E
    ta = Array{Float64,2}(undef,N,M)
    xm = zeros(3,M)
    for n = 1:N
        xm[1:2,:] .= micgeom[1:2,:] .- Rxy[1:2,n]
        for m = 1:M
            res = nlsolve((fvec,x)->shear!(fvec,x,xm[:,m],Ma,h),[.1,.1,h])
            xi = res.zero
            r1 = norm(xi)
            d = xi[1]*Ma*c/r1
            c1 = d + sqrt(d^2+c^2-(Ma*c)^2)
            r2 = norm(xm[:,m]-xi)
            ta[n,m] = r1/c1 + r2/c
        end
    end
    return ta
end

end
