function beamformer(E::Environment{T},C::Constants,V::SteeringMatrix{T}) where T <: AbstractFloat
    b = Array{T}(E.N,E.Nf)
    vd = Array{Complex{T}}(E.M)
    CSM = Array{Complex{T}}(E.M,E.M)
    # Compute beamforming
    for j in 1:length(E.f)
        CSM = E.CSM[j,:,:]
        for i in eachindex(E.D0)
            vd = V.v[:,i,j]
            b[i,j] = real(vd'*CSM*vd)
        end
    end
    return reshape(b,E.Nx,E.Ny,E.Nf)
end
