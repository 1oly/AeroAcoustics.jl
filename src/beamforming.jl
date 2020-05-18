"""
    beamforming(csm,v)

Calculate frequency-domain beamforming using cross-spectral matrix `csm` and steering vector `v`.
`csm` must be a square (Hermitian) matrix optionally with a third (frequency) dimension.
First dimension of `v` and `csm` must be equal.
"""
function beamforming(csm::A,v::B) where {T1<:AbstractFloat,T2<:AbstractFloat, N,  A<:AbstractArray{Complex{T1},N},B<:AbstractArray{Complex{T2},N}}
    #@assert size(csm,1) == size(csm,2) "csm must be square with dimensions M x M (x Nf)!"
    #@assert size(v,1) == size(csm,1) "First dimension of v and csm must be equal!"
    b = Array{T2, 2}(undef, size(v,2), size(csm,3))
    for j in axes(csm,3)
        csmd = @view csm[:,:,j]
        for i in axes(v,2)
            vd = @view v[:,i,j]
            b[i,j] = real(vd'*csmd*vd)
        end
    end
    return b
end

"""
    beamforming(Environment)

Calculate frequency-domain beamforming using cross-spectral matrix `csm` and steering vector `v`
stored in Environment struct
"""
function beamforming(E::Environment)
    @unpack CSM_s,steeringvec,M,N,Nf,fn = E
    b = Array{Float64}(undef, N, Nf)
    @views @inbounds for j in 1:Nf
        AeroAcoustics.bf_col!(b[:,j],steeringvec.arr[:,:,j],selectdim(CSM_s.arr, 3, j))
    end
    return FreqArray(b,fn)
end

function bf_col!(b::A,st::B,csm::C) where {T1<:AbstractFloat,T2<:AbstractFloat, N, A<:AbstractArray{T1,1}, B<:AbstractArray{Complex{T1},N},C<:AbstractArray{Complex{T2},N}}
    @views @inbounds for i in 1:length(b)
        b[i] = real(st[:,i]'*csm*st[:,i])
    end
end
