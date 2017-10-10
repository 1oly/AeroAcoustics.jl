function beamformer2(E::Environment,C::Constants,V::SteeringMatrix;psf::Bool=false)
    b = Array{Float64}(E.N,E.Nf)

    # Compute beamforming
    for j in 1:length(E.f)
        for i in eachindex(E.D0)
            b[i,j] = real(V.v[:,i,j]'*E.CSM[j,:,:]*V.v[:,i,j])
        end
    end
    return reshape(b,E.Nx,E.Ny,E.Nf)
end

function beamformer2(E::Environment,C::Constants,V::SteeringMatrix;psf::Bool=true)
    b = Array{Float64}(E.N,E.Nf)
    PSF = similar(b)
    gm = Array{Complex{Float64}}(E.M)

    cent = round(Int64,E.N/2)+1

    # Compute beamforming
    for j in 1:length(E.f)
        for i in eachindex(E.D0)
            b[i,j] = real(V.v[:,i,j]'*E.CSM[j,:,:]*V.v[:,i,j])
        end
    end

    # Compute PSF
    for j in 1:length(E.f)
        gm = E.M*V.v[:,cent,j].*(E.D0[cent]./E.D[cent,:]) # Compensate to get correct level.
        for i in eachindex(E.D0)
            PSF[i,j] = abs(V.v[:,i,j]'*gm)^2
            #for m in 1:E.M
                #PSF[i,j] = (E.M*E.D0[i]/E.D[i,m])^2*abs(V.v[:,i,j]'*gm)^2
            #end
        end
    end
    return reshape(b,E.Nx,E.Ny,E.Nf),reshape(PSF,E.Nx,E.Ny,E.Nf)
end
