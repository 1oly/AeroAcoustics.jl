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
"""
function sourceintegration(x::Vector{T1},env::Environment,limits::Vector{T2}) where {T1,T2}
    xi = findall((env.rx.>=limits[1]) .& (env.rx.<=limits[2]))
    yi = findall((env.ry.>=limits[3]) .& (env.ry.<=limits[4]))
    I = LinearIndices((env.Nx,env.Ny))[CartesianIndex.(xi,yi)]
    return sum(x[I])
end

function sourceintegration(x::Vector{T1},env::Environment,limits::Vector{Vector{T2}}) where {T1,T2}
    out = Float64[]
    for src in limits
        push!(out,sourceintegration(x,env,src))
    end
    return out
end

function sourceintegration(x::FreqArray,env::Environment,limits::Vector{Vector{T}}) where T
    out = Array{Float64}(undef,length(x.fc),length(limits))
    for i in 1:length(x.fc)
        tmp = Float64[]
        for src in limits
            push!(tmp,sourceintegration(x.arr[:,i],env,src))
        end
        out[i,:] = tmp
    end
    return FreqArray(out,x.fc)
end

function sourceintegration(x::FreqArray,env::Environment,limits::Vector{T}) where T
    out = Array{Union{Missing, Float64}}(undef,length(x.fc))
    for i in 1:length(x.fc)
        out[i] = sourceintegration(x.arr[:,i],env,limits)
    end
    return FreqArray(out,x.fc)
end
