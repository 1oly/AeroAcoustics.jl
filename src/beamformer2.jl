@traceable function beamformer2(E::Environment,C::Constants,V::SteeringMatrix;psf::Bool=false)
    b = Array{Float64}(E.N,E.Nf)

    # Compute beamforming
    for j in 1:length(E.f)
        for i in eachindex(E.D0)
            b[i,j] = abs(V.v[:,i,j]'*E.CSM[j,:,:]*V.v[:,i,j])
        end
    end

    if psf
        PSF = similar(b)
        for j in 1:length(E.f)
            gm = vec(V.v[:,round(Int64,E.N/2)+1,j])
            for i in eachindex(E.D0)
                for m in 1:E.M
                    PSF[i,j] = (E.M*E.D0[i]/E.D[i,m])^2*abs(V.v[:,i,j]'*gm)^2
                end
            end
        end

        return reshape(b,E.Nx,E.Ny,E.Nf),reshape(PSF,E.Nx,E.Ny,E.Nf)
    else
        return reshape(b,E.Nx,E.Ny,E.Nf)
    end
end
