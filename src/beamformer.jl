function beamformer{T,C}(Nx::Int64,Ny::Int64,X,Y,z0,f::T,rn,CSM::Array{C,2};psf::Bool=false)
    const M = size(rn,1)    # Number of microphones
    const omega = 2pi*f     # Angular frequency
    const c = 343           # Speed of sound
    # CSM[eye(Bool,M)] = 0;    # Naive diagonal removal

    # Allocation of arrays
    gj = Array{C}(M)
    gjs = Array{C}(Nx,Ny,M)
    b = Array{T}(Nx,Ny)
    # Compute transfer functions
    for i in 1:Nx
        for j in 1:Ny
            r0 = sqrt(X[i,j]^2 + Y[i,j]^2 + z0^2)
            Threads.@threads for m in 1:M
                rm = sqrt((X[i,j]-rn[m,1])^2+(Y[i,j]-rn[m,2])^2 + z0^2)
                # TYPE II? Steering vector:
                gj[m] = (1/M)*(rm/r0)*exp(-im*omega*(rm-r0)/c)
            end
            gjs[i,j,:] = gj
            b[i,j] = real(gj'*CSM*gj)
        end
    end


    if psf
        PSF = Array{T}(Nx,Ny)
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

function beamformer{T,C}(Nx::Int64,Ny::Int64,X,Y,z0,f::Array{T,1},rn,CSM::Array{C,3};psf::Bool=false)
    const M = size(rn,1)    # Number of microphones
    Nf = length(f)
    b = Array{T}(Nx,Ny,Nf)
    gjs = Array{C}(Nx,Ny,M,Nf)
    if psf
        PSF = Array{T}(Nx,Ny,Nf)
        Threads.@threads for i in 1:Nf
            b[:,:,i],gjs[:,:,:,i],PSF[:,:,i] = beamformer(Nx,Ny,X,Y,z0,f[i],rn,CSM[i,:,:];psf=true)
        end
        return b,gjs,PSF
    else
        Threads.@threads for i in 1:Nf
            b[:,:,i],gjs[:,:,:,i] = beamformer(Nx,Ny,X,Y,z0,f[i],rn,CSM[i,:,:];psf=false)
        end
        return b,gjs
    end

end
