abstract type WindTunnelType end
abstract type SteeringVectorType end

struct Type2 <: SteeringVectorType end
struct Type3 <: SteeringVectorType end
struct Shear <: SteeringVectorType end
struct Uniform <: SteeringVectorType end


# TODO: From imfilter.jl, dispatch on string arguments...
#const valid_types = ("type2", "type3", "shear")

#function typeinstance(types::AbstractString)
#    if type ∈ valid_types
#        return Kind(Symbol(type))
#    else
#        throw(ArgumentError("$type not a recognized steering vector type"))
#    end
#end
#typeinstance(b::SteeringVectorType) = b

struct Constants{T<:Real} <: WindTunnelType
    # TODO: promote
    # TODO: direction of flow
    Ma::T   # Mach number
    c::T    # Speed of sound
    h::T    # Shear layer distance from mic array
end

# If no distance or Ma is set, default is zero
Constants(a::T) where T<:Real = Constants(zero(T),a::T,zero(T))
# If no distance is set, default is zero
Constants(a::T,b::T) where T<:Real = Constants(a::T,b::T,zero(T))

struct CrossSpectralMatrix{T<:AbstractFloat} <: WindTunnelType
    csmReal::Array{T,3}
    csmImag::Array{T,3}
    binCenterFrequenciesHz::Array{T,1}
    diagrm::Bool
end

# TODO: Should Environment have a dx,dy field?
struct Environment{T<:AbstractFloat} <: WindTunnelType
    N::Int64
    M::Int64
    Nx::Int64
    Ny::Int64
    Nz::Int64
    Nf::Int64
    f::Array{T,1}          # frequency vector
    micgeom::Array{T,2}    # microphone coordinates
    # TODO: Or just a vector e.g. linspace
    rx::StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T}}              # x range
    ry::StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T}}              # y range
    rz::T              # z distance
    Rxy::Array{T,2}
    D0::Array{T,1}
    D::Array{T,2}
    CSM::CrossSpectralMatrix{T}
end

struct SteeringMatrix{T<:AbstractFloat,S<:SteeringVectorType} <: WindTunnelType
    v::Array{Complex{T},3}
    kind::S
end

# input data and results from deconvolution:
# TODO: make fistaprox! dispatch on data::DeconvData and return res::DeconvResults

struct DeconvData{T<:AbstractFloat, R<:Real, I<:Int} <: WindTunnelType
    x::Array{T,3}
    psf::Array{T,3}
    b::Array{T,3}
    tol::R
    maxit::I
    lambda::R
    function DeconvData{T,R,I}(
        x::Array{T,3},
        psf::Array{T,3},
        b::Array{T,3},
        tol::R,
        maxit::I,
        lambda::R) where {T<:AbstractFloat, R<:Real, I<:Int}
        new(x,psf,b,tol,maxit,lambda)
    end
end

DeconvData{T<:AbstractFloat, R<:Real, I<:Int}(
    x::Array{T,3} = zeros(Float64,1,1,1),
    psf::Array{T,3} = zeros(Float64,1,1,1),
    b::Array{T,3} = zeros(Float64,1,1,1),
    tol::R = 1e-3,
    maxit::I = 5000,
    lambda::R = 0.01) = DeconvData{T,R,I}(x,psf,b,tol,maxit,lambda)

mutable struct DeconvResults{T<:AbstractFloat, R<:Real, I<:Int} <: WindTunnelType
    xopt::Array{T,3}
    tol::R
    maxit::Vector{I}
    lambda::R
    prox::Symbol
    function DeconvResults{T,R,I}(
        xopt::Array{T,3},
        tol::R,
        maxit::Vector{I},
        lambda::R,
        prox::Symbol) where {T<:AbstractFloat, R<:Real, I<:Int}
        new(xopt,tol,maxit,lambda,prox)
    end
end

DeconvResults{T<:AbstractFloat, R<:Real, I<:Int}(
    xopt::Array{T,3} = zeros(1,1,1),
    tol::R = zero(Float64),
    maxit::Vector{I} = Int64[],
    lambda::R = zero(Float64),
    prox::Symbol = Symbol()) = DeconvResults{T,R,I}(xopt,tol,maxit,lambda,prox)

DeconvResults(prox::Symbol) = DeconvResults(zeros(1,1,1),zero(Float64),Int64[],zero(Float64),prox)
