abstract type WindTunnelType end

struct Constants{T<:Real} <: WindTunnelType
    Ma::T   # Mach number
    c::T    # Speed of sound
end

mutable struct Environment{N<:Real,T<:AbstractFloat} <: WindTunnelType
    f::AbstractVector{T}    # frequency vector
    mic::AbstractMatrix{T}  # microphone coordinates
    rx::AbstractVector{N}   # x range
    ry::AbstractVector{N}   # y range
    z0::N                   # z distance
end
