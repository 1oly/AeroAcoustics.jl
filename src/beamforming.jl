"""
    beamforming(csm,v)

Calculate frequency-domain beamforming using cross-spectral matrix `csm` and steering vector `v`.
`csm` must be a square (Hermitian) matrix optionally with a third (frequency) dimension.
First dimension of `v` and `csm` must be equal.
"""
function beamforming(csm::AbstractArray{Complex{T},N},v::AbstractArray{Complex{T},N}) where {T <: AbstractFloat,N}
    @assert size(csm,1) == size(csm,2) "csm must be square with dimensions M x M (x Nf)!"
    @assert size(v,1) == size(csm,1) "First dimension of v and csm must be equal!"
    b = Array{T, 2}(undef, size(v,2), size(csm,3))
    for j in axes(csm,3)
        csmd = @view csm[:,:,j]
        for i in axes(v,2)
            vd = @view v[:,i,j]
            b[i,j] = real(vd'*csmd*vd)
        end
    end
    return b
end
