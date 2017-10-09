abstract type WindTunnelType end

struct Constants{T<:Real} <: WindTunnelType
    Ma::T   # Mach number
    c::T    # Speed of sound
end

mutable struct Environment <: WindTunnelType
    N::Int64
    M::Int64
    Nx::Int64
    Ny::Int64
    Nz::Int64
    Nf::Int64
    f::Array{Float64,1}          # frequency vector
    micgeom::Array{Float64,2}    # microphone coordinates
    rx::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}              # x range
    ry::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}              # y range
    rz::Float64              # z distance
    Rxy::Array{Float64,2}
    D0::Array{Float64,1}
    D::Array{Float64,2}
    CSM::Array{Complex{Float64},3}
end

mutable struct SteeringMatrix <: WindTunnelType
    v::Array{Complex{Float64},3}
end
