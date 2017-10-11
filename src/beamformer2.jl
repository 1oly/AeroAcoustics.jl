function beamformer2(E::Environment{T},C::Constants,V::SteeringMatrix{T}) where T <: AbstractFloat
    b = Array{T}(E.N,E.Nf)

    # Compute beamforming
    for j in 1:length(E.f)
        for i in eachindex(E.D0)
            vd = @view V.v[:,i,j]
            b[i,j] = real(vd'*E.CSM[j,:,:]*vd)
        end
    end
    return reshape(b,E.Nx,E.Ny,E.Nf)
end

function pointspreadfunction(E::Environment{T},C::Constants,V::SteeringMatrix{T}) where T <: AbstractFloat
    PSF = Array{T}(E.N,E.Nf)
    gm = Array{Complex{T}}(E.M)
    cent = round(Int64,E.N/2)+1

    # Compute PSF
    for j in 1:length(E.f)
        gm = E.M*V.v[:,cent,j].*(E.D0[cent]./E.D[cent,:]) # Compensate to get correct level.
        for i in eachindex(E.D0)
            PSF[i,j] = abs2(V.v[:,i,j]'*gm)
        end
    end
    return reshape(PSF,E.Nx,E.Ny,E.Nf)
end
