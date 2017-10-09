function SPL{T}(p::Array{T})
    s = similar(p)
    s[p.>0] = 10*log10.(p[p.>0]/4e-10)
    s[p.<=0] = -350
    return s
end
SPL(p::Number) = 10*log10(p/4e-10)

function shear(xi,fvec,xm,M,h)
    a = sqrt(xi[1]^2+(1-M^2)*(xi[2]^2+xi[3]^2))
    b = sqrt((xm[1]-xi[1])^2+(xm[2]-xi[2])^2+(xm[3]-xi[3])^2)
    fvec[1] = (1/a)*xi[1]-(1/b)*(1-M^2)*(xm[1]-xi[1]) - M
    fvec[2] = (1/a)*xi[2]-(1/b)*(xm[2]-xi[2])
    fvec[3] = xi[3] - h
end

function octavebandlimits(fc,kind)
    fl,fu = similar(fc)
    const C = 10^(3./10.)
    fl = fc*C^(-1/(2*kind))
    fu = fc*C^(1/(2*kind))
    return fl, fu
end

function beamformersetup(dx,dy,x,y,z,f,micgeom,csmdata)
    const Nm, M = size(micgeom)
    # Region of interest
    rx = x[1]:dx:x[2]
    ry = y[1]:dy:y[2]
    rz = z  # TODO: Make 3D capabilities?

    Nx = length(rx)
    Ny = length(ry)
    Nz = length(rz)

    fl = f[1]
    fu = f[end]

    Rxy = hcat([[x, y, z] for x in rx, y in ry, z in rz]...)
    D0 = colwise(Euclidean(), Rxy, [0,0,0]) # Distance from center of array to grid points
    if Nm == 2
        micgeom = [micgeom; zeros(M)']
    end
    D = pairwise(Euclidean(), Rxy, micgeom) # Distance from each mic to grid points
    const N = size(D,1)
    fc = csmdata["binCenterFrequenciesHz"]
    fn = fc[(fc.>=fl) .& (fc.<=fu)]
    ind = findin(fc,fn)
    csm = convert(Array{Complex{Float64},3},csmdata["csmReal"][ind,:,:]+im*csmdata["csmImaginary"][ind,:,:])
    Nf = length(ind)
    return Environment(N,M,Nx,Ny,Nz,Nf,fn,micgeom,rx,ry,rz,Rxy,D0,D,csm)
end

function steeringvectors(E::Environment,C::Constants,kind::String="II")
    kw = 2pi*E.f/C.c
    if kind == "II"
        vi = Array{Complex{Float64}}(E.M,E.N,length(E.f))
        for j in 1:length(E.f)
            for i in eachindex(E.D0)
                for m in 1:E.M
                    vi[m,i,j] = (1/E.M)*(E.D[i,m]/E.D0[i])*exp(im*kw[j]*(E.D0[i]-E.D[i,m]))
                end
            end
        end
    end
    return SteeringMatrix(vi) # [mics,gridpoints,freqs]
end
