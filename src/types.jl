abstract type WindTunnelType end

struct Constants{T<:Real} <: WindTunnelType
    Ma::T   # Mach number
    c::T    # Speed of sound
    h::T    # Shear layer distance from mic array
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
    CSM::Array{Complex{T},3}
end

struct SteeringMatrix{T<:AbstractFloat} <: WindTunnelType
    v::Array{Complex{T},3}
    kind::AbstractString
end
