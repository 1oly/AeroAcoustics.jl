function beamformer(E::Environment{T},C::Constants,V::SteeringMatrix{T,S}) where {T <: AbstractFloat, S <: SteeringVectorType}
    b = Array{T}(E.N,E.Nf)
    vd = Array{Complex{T}}(E.M)
    C = Array{Complex{T}}(E.M,E.M)
    # Compute beamforming
    for j in 1:length(E.f)
        C = E.CSM.csmReal[j,:,:] + im*E.CSM.csmImag[j,:,:]
        for i in eachindex(E.D0)
            vd = V.v[:,i,j]
            b[i,j] = real(vd'*C*vd)
        end
    end
    return reshape(b,E.Nx,E.Ny,E.Nf)
end
