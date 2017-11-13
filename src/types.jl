abstract type WindTunnelType end
abstract type SteeringVectorType end

struct Kind <: SteeringVectorType
    style::Symbol
end

Kind() = Kind(:type3)

const valid_types = ("type2", "type3", "shear")

function typeinstance(types::AbstractString)
    if type ∈ valid_types
        return Kind(Symbol(type))
    else
        throw(ArgumentError("$type not a recognized steering vector type"))
    end
end
typeinstance(b::SteeringVectorType) = b

struct Constants{T<:Real} <: WindTunnelType
    Ma::T   # Mach number
    c::T    # Speed of sound
    h::T    # Shear layer distance from mic array
end

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
    rx::StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T}}              # x range
    ry::StepRangeLen{T,Base.TwicePrecision{T},Base.TwicePrecision{T}}              # y range
    rz::T              # z distance
    Rxy::Array{T,2}
    D0::Array{T,1}
    D::Array{T,2}
    #CSM::Array{Complex{T},3}
    CSM::CrossSpectralMatrix{T}
end

struct SteeringMatrix{T<:AbstractFloat} <: WindTunnelType
    v::Array{Complex{T},3}
    kind::Kind
end
