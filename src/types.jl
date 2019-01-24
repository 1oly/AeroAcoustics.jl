abstract type AeroAcousticType end

@with_kw struct Environment{R<:Real,T<:AbstractFloat,T1<:AbstractArray} <: AeroAcousticType
    micgeom::Matrix{T}
    z0::R
    f::T1
    Nx::Int = 21
    Ny::Int = 21
    @assert Nx>0 && Ny>0 && z0>0 "Nx, Ny, z0 must be positive"
    xlim::Tuple = (-1.,1.)
    ylim::Tuple = (-1.,1.)
    flim::Tuple = extrema(f)
    c::R = 343.0 # Speed of sound
    Ma::R = 0.0 # Flow Mach speed
    ### Compute extra parameters
    fn = f[(f.>=flim[1]) .& (f.<=flim[2])]
    M::Int = size(micgeom,2)
    N = Nx*Ny
    Nf::Int = length(fn)
    rx = range(xlim[1],stop = xlim[2], length = Nx)
    ry = range(ylim[1],stop = ylim[2], length = Ny)
    Rxy = hcat([[x, y, z] for x in rx, y in ry, z in z0]...)
    D0 = colwise(Euclidean(), Rxy, [0,0,0]) # Distance from center of array to grid points
    D = pairwise(Euclidean(), Rxy, micgeom) # Distance from each mic to grid points
end

struct CrossSpectralMatrix{T1<:AbstractFloat,T2<:AbstractArray} <: AeroAcousticType
    csm::Array{Complex{T1},3}
    fc::T2
    diagrm::Bool
end
