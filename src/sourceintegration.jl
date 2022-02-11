"""
    sourceintegration(x::Vector,env::Environment,limits::Vector{T}) where T
    sourceintegration(x::Vector,env::Environment,limits::Vector{Vector{T}}) where T
    sourceintegration(x::FreqArray,env::Environment,limits::Vector{T}) where T
    sourceintegration(x::FreqArray,env::Environment,limits::Vector{Vector{T}}) where T

Source integration of `x` from `limits` given in cartesian coordinates `[xmin,xmax,ymin,ymax]`.
```
src1 = [-0.1,0.1,-0.1,0.1]
src2 = [-0.2,0.0,-0.1,0.1]
limits = [src1,src2]
srcint = sourceintegration(x,env,limits)
```

Optionally supply keyword `int_thres` to set a threshold on source integration. Default is set to `Inf` which sum all values within integration region.
"""
function sourceintegration(x::Vector{T1},env::Environment,limits::Vector{T2};int_thres=Inf) where {T1,T2}
    xi = findall((env.rx.>=limits[1]) .& (env.rx.<=limits[2]))
    yi = findall((env.ry.>=limits[3]) .& (env.ry.<=limits[4]))
    I = LinearIndices((env.Nx,env.Ny))[xi,yi]
    pos_idx = x[I] .> 0
    if iszero(pos_idx)
        return NaN
    else
        fmax,g_idx = findmax(x[I][pos_idx])
        thres_idx = SPL.(x[I][pos_idx]) .- SPL.(fmax) .<= int_thres
        return sum(x[I][pos_idx][thres_idx])
    end
end

function sourceintegration(x::Vector{T1},env::Environment,limits::Vector{Vector{T2}};int_thres=Inf) where {T1,T2}
    out = Float64[]
    for src in limits
        push!(out,sourceintegration(x,env,src;int_thres=int_thres))
    end
    return out
end

function sourceintegration(x::FreqArray,env::Environment,limits::Vector{Vector{T}};int_thres=Inf) where T
    out = Array{Float64}(undef,length(limits),length(x.fc))
    for i in 1:length(x.fc)
        tmp = Float64[]
        for src in limits
            push!(tmp,sourceintegration(x.arr[:,i],env,src;int_thres=int_thres))
        end
        out[:,i] = tmp
    end
    return FreqArray(out,x.fc)
end

function sourceintegration(x::FreqArray,env::Environment,limits::Vector{T};int_thres=Inf) where T
    out = Array{Union{Missing, Float64}}(undef,length(x.fc))
    for i in 1:length(x.fc)
        out[i] = sourceintegration(x.arr[:,i],env,limits;int_thres=int_thres)
    end
    return FreqArray(out,x.fc)
end


"""
    point_to_region(src,dxdy)

Helper function to get integration limits from `src` coordinates and size of area `dxdy`
```
sources = [(0.1,0.1),(-0.1,0.1),(-0.1,-0.1),(0.1,-0.1)]
dxdy = (0.1,0.1)
limits = point_to_region(sources,dxdy)
srcint = sourceintegration(x,env,limits)
```
"""
function point_to_region(src::NTuple{2},dxdy::Tuple{Vararg{T,N}}) where {T,N}
    N == 1 && (dxdy = (dxdy[1],dxdy[1]))
    xmin = src[1]-dxdy[1]/2
    xmax = src[1]+dxdy[1]/2
    ymin = src[2]-dxdy[2]/2
    ymax = src[2]+dxdy[2]/2
    return [xmin,xmax,ymin,ymax]
end
function point_to_region(src::NTuple{2},dxdy::T) where T <: Real
    point_to_region(src::NTuple{2},(dxdy[1],dxdy[1]))
end
function point_to_region(src::T,dxdy) where T <: AbstractArray
    eltype(src) <: NTuple || (src = Tuple.(src)) # Convertion to array of tuples
    out = Array{Array{Float64,1},1}(undef,0)
    for i = 1:length(src)
        push!(out,point_to_region(src[i],dxdy))
    end
    return out
end

"""
    SPI(b::FreqArray,E::Environment,dxdy,limits)

Convenience function to normalize source integration of a beamforming result.
"""
function SPI(b::FreqArray,E::Environment,dxdy,limits)
    midpoint = (sum(extrema(E.rx))/2,sum(extrema(E.ry))/2)
    limits_psf  = AeroAcoustics.point_to_region(midpoint,dxdy)
    srcint_psf = sourceintegration(psf(E),E,limits_psf)
    srcint_b = sourceintegration(b,E,limits)
    return FreqArray(srcint_b./srcint_psf,srcint_b.fc)
end