function beamformer_corr{T,C}(
        Nx::Int64,
        Ny::Int64,
        X::Array{T,2},
        Y::Array{T,2},
        z0::T,
        f::T,
        rn::Array{T,2},
        CSM::Array{C,2};
        psf::Bool=false)
    const M::Int64 = size(rn,1)    # Number of microphones
    const omega::T = 2pi*f         # Angular frequency
    const c::Float64 = 343.0       # Speed of sound
    const Ma::Float64 = 0.2
    const h::Float64 = 0.25
    # CSM[eye(Bool,M)] = 0;        # Naive diagonal removal

    # Allocation of arrays
    gj = Array{C}(M)
    gjs = Array{C}(Nx,Ny,M)
    b = Array{T}(Nx,Ny)
    xm = zeros(3)
    # Compute transfer functions
    for i in 1:Nx
        for j in 1:Ny
            r0::Float64 = sqrt(X[i,j]^2 + Y[i,j]^2 + z0^2)
            for m in 1:M
                #rm::Float64 = sqrt((X[i,j]-rn[m,1])^2+(Y[i,j]-rn[m,2])^2 + z0^2)
                # Shear layer correction:
                xm = [rn[m,1]-X[i,j], rn[m,2]-Y[i,j], 0]
                res = nlsolve((x,fvec)->shear(x,fvec,xm,Ma,h), [ 1.; 1.; 1.])
                xi = res.zero
                r1 = sqrt(xi[1]^2+xi[2]^2+xi[3]^2)
                d = xi[1]*Ma*c/r1
                c1 = d + sqrt(d^2+c^2-(Ma*c)^2)
                r2 = sqrt((xm[1]-xi[1])^2+(xm[2]-xi[2])^2+(xm[3]-xi[3])^2)
                ta = r1/c1 + r2/c
                gj[m] = (1/M)*(ta*c/r0)*exp(-im*omega*ta) # TYPE II? Steering vector:
            end
            gjs[i,j,:] = gj
            b[i,j] = real(gj'*CSM*gj)
        end
    end


    if psf
        PSF = Array{T}(Nx,Ny)
        g1 = Array{C}(M)
        midx::Int64 = 0
        midy::Int64 = 0
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

function beamformer_corr{T,C}(Nx::Int64,Ny::Int64,X,Y,z0,f::Array{T,1},rn,CSM::Array{C,3};psf::Bool=false)
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
