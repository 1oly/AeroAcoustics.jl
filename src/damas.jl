import Base: getindex, iterate

"""
    damas(env::Environment, b [,f<:AbstractArray, ω::Real]; maxiter=10)
Performs exactly `maxiter` DAMAS iterations for all frequency bins in `env.fn` or frequencies in `f`. `f` must contain values that match exact a (sub)set of the values in `env.fn`. Successive Over Relaxation (SOR) can be set with relaxation parameter `ω`. Default is `ω=1` corresponding to no relaxation.
"""
function damas(env::Environment, b::FA, f::T=env.fn, ω::Real=1.0; maxiter::Int=10) where {T <: AbstractArray, FA <: FreqArray}
    env.fn == f ? (f_inds=1:env.Nf) : (f_inds = findall(x->x in f, env.fn))
    x = zeros(size(b,1),length(f_inds))
    for (i,fn) in enumerate(f_inds)
        @views damas!(x[:,i], env.steeringvec.arr[:,:,fn], b[:,fn], ω; maxiter=maxiter)
    end
    return FreqArray(x,env.fn[f_inds])
end

function damas!(x, steer, b, ω::Real=1.0; maxiter::Int=10)
    iterable = DAMASSORIterable(x, steer, similar(x), similar(x), b, ω, maxiter)
    for _ = iterable end
    x
end

mutable struct DAMASSORIterable{solT,steerT,vecT,rhsT,numT}
    x::solT
    steer::steerT
    tmp1::vecT
    tmp2::vecT
    b::rhsT
    ω::numT
    maxiter::Int
end

start(::DAMASSORIterable) = 1
done(it::DAMASSORIterable, iteration::Int) = iteration > it.maxiter
function iterate(s::DAMASSORIterable, iteration::Int=start(s))
    if done(s, iteration) return nothing end

    n = size(s.b, 1)

    for col = 1 : n
        AeroAcoustics.psf_col!(s.tmp2,s.steer,col)
        @simd for row = 1 : col - 1
            @inbounds s.tmp1[row] = max(0.0,s.tmp1[row] - s.tmp2[row] * s.x[col])
        end

        s.tmp1[col] = s.b[col]
    end

    for col = 1 : n
        AeroAcoustics.psf_col!(s.tmp2,s.steer,col)
        @inbounds s.x[col] += s.ω * (s.tmp1[col] / s.tmp2[col] - s.x[col])
        @simd for row = col + 1 : n
            @inbounds s.tmp1[row] = max(0.0,s.tmp1[row] - s.tmp2[row] * s.x[col])
        end
    end

    nothing, iteration + 1
end
