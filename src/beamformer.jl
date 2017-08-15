function beamformer(Nx::Int64,Ny::Int64,X,Y,z0,f,rn,CSM;psf=False)
    const M = size(rn,1)    # Number of microphones
    const omega = 2pi*f     # Angular frequency
    const c = 343           # Speed of sound

    #r0 = sqrt.(X.^2 + Y.^2 .+ z0^2)

    # Allocation of arrays
    gj = Array{Complex128}(M)
    gjs = Array{Complex128}(Nx,Ny,M)
    b = Array{Float64}(Nx,Ny)

    # Compute transfer functions
    #
    for i in 1:Nx
        for j in 1:Ny
            r0 = sqrt(X[i,j]^2 + Y[i,j]^2 + z0^2)
            Threads.@threads for m in 1:M
                #rmn[i,j,m] = sqrt.((X[i,j]-rn[m,1])^2+(Y[i,j]-rn[m,2])^2 + z0^2)
                #gj[i,j,m] = (r0/rmn[i,j,m])*exp(im*omega*(rmn[i,j,m]-r0)/c)
                rm = sqrt((X[i,j]-rn[m,1])^2+(Y[i,j]-rn[m,2])^2 + z0^2)
                gj[m] = (1/(Nx*Ny))*(rm/r0)*exp(-im*omega*(rm-r0)/c)
            end
            gjs[i,j,:] = gj
            b[i,j] = real(gj'*CSM*gj)
        end
    end

    # CSM[eye(Bool,M)] = 0;    # Diagonal removal
    if psf
        PSF = Array{Float64}(Nx,Ny)
        if iseven(Nx)
            midx = Nx/2
        elseif isodd(Nx)
            midx = round(Int64,Nx/2)+1
        end
        if iseven(Ny)
            midy = Ny/2
        elseif isodd(Ny)
            midy = round(Int64,Ny/2)+1
        end
        g1 = vec(gjs[midx,midy,:])
        Threads.@threads for i in 1:Nx
            for j in 1:Ny
                PSF[i,j] = real(vec(gjs[i,j,:])'*g1)
            end
        end
        return b,gjs,PSF
    else
        return b,gjs
    end
end
