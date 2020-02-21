function shear!(fvec,xi,xm,Ma,h)
    a = sqrt(xi[1]^2+(1-Ma^2)*(xi[2]^2+xi[3]^2))
    b = norm(xm-xi)
    fvec[1] = (1/a)*xi[1]-(1/b)*(1-Ma^2)*(xm[1]-xi[1]) - Ma
    fvec[2] = (1/a)*xi[2]-(1/b)*(xm[2]-xi[2])
    fvec[3] = h - xi[3]
end

"""
    shearlayercorrection(E::Environment)

Compute shear layer correction using Amiet's derivation for a zero-thickness shear layer. Returns time-delay and amplitude correction. Relies on `shear=true`, `Ma` and `h` set in `Environment()`. The coordinate system is assumed to be at the microphone array center and `h` is the distance from the array center to the shear layer (this is contrary to the references).

References:
-	R. K. Amiet, “Refraction of sound by a shear layer,” Journal of Sound and Vibration, vol. 58, no. 4, pp. 467–482, 1978.
-	C. Bahr, et. al “Shear Layer Correction Validation Using A Non-Intrusive Acoustic Point Source,” presented at the 16th AIAA/CEAS Aeroacoustics Conference, Reston, Virigina, 2012.
"""
function shearlayercorrection(E::Environment)
    @unpack N,M,micgeom,Rxy,Ma,h,c,z0,ampcorr = E
    ta = Array{Float64,2}(undef,N,M)
    pcpm2 = Array{Float64,2}(undef,N,M)
    xn = zeros(3,M)
    for n = 1:N
        xn .= Rxy[:,n] .- micgeom
        for m = 1:M
            res = nlsolve((fvec,x)->shear!(fvec,x,xn[:,m],-Ma,h),[.1,.1,h])
            xi = res.zero
            r1 = norm(xi)
            d = -xi[1]*Ma*c/r1
            c1 = d + sqrt(d^2+c^2-(Ma*c)^2)
            r2v = xn[:,m]-xi
            r2 = norm(r2v)
            ta[n,m] = r1/c + r2/c1
            if ampcorr
                theta = acos(-sign(Ma)*xi[1])
                hH = (z0-h)/z0
                zeta = sqrt((1 - abs(Ma)*cos(theta))^2 - cos(theta)^2)
                pcpm2[n,m]=1/4/zeta^2*hH^2*(1+(1/hH-1)*zeta.^3/sin(theta)^3)*(1+(1/hH-1)*zeta/sin(theta))*(zeta+sin(theta)*(1-abs(Ma)*cos(theta))^2)^2
            end
        end
    end

    return ta,pcpm2
end

"""
    steeringvectors(E::Environment)

Pre-compute steeringvectors for beamforming using an `Environment` with the needed parameters.
"""
function steeringvectors(E::Environment)
    @unpack fn,c,M,N,Nf,D,D0 = E
    kw = 2pi*fn/c
    vi = Array{ComplexF64,3}(undef,N,M,Nf)
    Threads.@threads for j in 1:Nf
        vi[:,:,j] .= 1 ./(D.*D0.*sum(1 ./D.^2,dims=2)).*exp.(-im.*kw[j].*(D.-D0))
    end
    steeringvec = FreqArray(permutedims(vi,(2,1,3)),fn)
end

function steeringvectors!(E::Environment)
    @unpack fn,c,M,N,Nf,D,D0,Ma,h = E
    vi = Array{ComplexF64,3}(undef,N,M,Nf)

    if E.shear
        w = 2pi*fn
        ta,pcpm2 = AeroAcoustics.shearlayercorrection(E)
        for j in 1:Nf
            vi[:,:,j] .= 1 ./(ta.*c.*D0.*sum(1 ./(ta.*c).^2,dims=2)).*exp.(-im.*w[j].*ta)
            #vi[:,:,j] .= (1 ./M).*(ta.*c./D0).*exp.(-im.*w[j].*ta)
        end
        if E.ampcorr
            for j in 1:Nf
                vi[:,:,j] .*= sqrt.(pcpm2)
            end
        end
    else
        kw = 2pi*fn/c
        for j in 1:Nf
            #vi[:,:,j] .= (1 ./M).*(D./D0).*exp.(-im.*kw[j].*(D.-D0))
            vi[:,:,j] .= 1 ./(D.*D0.*sum(1 ./D.^2,dims=2)).*exp.(-im.*kw[j].*(D.-D0))
        end
    end

    E.steeringvec = FreqArray(permutedims(vi,(2,1,3)),fn)
    return nothing
end
