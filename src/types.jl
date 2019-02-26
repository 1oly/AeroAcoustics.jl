abstract type AeroAcousticType end

struct FreqArray{T, N, AA<:AbstractVector} <: AbstractArray{T, N}
    arr::Array{T, N}
    fc::AA
end

# Standard array functions and indexing
Base.size(A::FreqArray) = size(A.arr)
Base.axes(A::FreqArray) = axes(A.arr)
Base.getindex(A::FreqArray, i::Int) = A.arr[i]
Base.getindex(A::FreqArray{T, N}, I::Vararg{Int, N}) where {T, N} = A.arr[I...]

@with_kw mutable struct Environment <: AeroAcousticType
    micgeom::Matrix{<:AbstractFloat}
    z0::Float64
    CSM::FreqArray
    flim::NTuple{2,Real} = extrema(CSM.fc)
    Cinds = (CSM.fc.>=flim[1]) .& (CSM.fc.<=flim[2])
    fn = CSM.fc[Cinds]
    Nf::Int = length(fn)
    Nx::Int = 21
    Ny::Int = 21
    @assert Nx>0 && Ny>0 && z0>0 "Nx, Ny, z0 must be positive"
    xlim::NTuple{2,Real} = (-1.,1.)
    ylim::NTuple{2,Real} = (-1.,1.)
    c::Real = 343.0 # Speed of sound
    Ma::Real = 0.0 # Flow Mach speed
    ### Compute extra parameters
    M::Int = size(micgeom,2)
    N = Nx*Ny
    rx = range(xlim[1],stop = xlim[2], length = Nx)
    ry = range(ylim[1],stop = ylim[2], length = Ny)
    Rxy = hcat([[x, y, z] for x in rx, y in ry, z in z0]...)
    D0 = colwise(Euclidean(), Rxy, [0,0,0]) # Distance from center of array to grid points
    D = pairwise(Euclidean(), Rxy, micgeom; dims=2) # Distance from each mic to grid points
    steeringvec = nothing
end
