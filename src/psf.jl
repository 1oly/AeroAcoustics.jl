"""
    psf(Environment[,cent])

Calculate frequency-domain point spread function using the Environment struct to access
steeringvectors. Optionally, supply the index where the psf is centered, default is (N/2)+1.
"""
function psf(E::Environment,cent::Int64=floor(Int,E.N/2)+1)
    @unpack steeringvec,M,N,Nf,fn = E
    p = Array{Float64, 2}(undef, N, Nf)
    @views @inbounds for j in 1:Nf
        AeroAcoustics.psf_col!(p[:,j],steeringvec.arr[:,:,j],cent)
    end
    return FreqArray(p,fn)
end

"""
    psf_col!(p,steeringvec,cent)

Calculate single frequency point spread function as a column vector.
"""
function psf_col!(p::A,steeringvec::B,cent::Int64) where {T<:AbstractFloat, N, A<:AbstractArray{T,1}, B<:AbstractArray{Complex{T},N}}
    M = size(steeringvec,1)
    @views p .= 0.5.*M^2 .*abs2.(steeringvec'*steeringvec[:,cent])
end
