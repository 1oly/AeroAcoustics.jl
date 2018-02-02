function steeringvectors(E::Environment{T},C::Constants{T},kind::Type3) where T <: AbstractFloat
    kw = 2pi*E.f/C.c
    vi = Array{Complex{T}}(E.M,E.N,length(E.f))
    Dsum = sum(1./E.D.^2,2)
    for j in 1:length(E.f)
        for i in eachindex(E.D0)
            for m in 1:E.M
                vi[m,i,j] = 1/(E.D[i,m]*E.D0[i]*Dsum[i])*exp(-im*kw[j]*(E.D[i,m]-E.D0[i]))
            end
        end
    end
    return SteeringMatrix(vi,kind)
end

function steeringvectors(E::Environment{T},C::Constants{T},kind::Type2) where T <: AbstractFloat
    kw = 2pi*E.f/C.c
    vi = Array{Complex{T}}(E.M,E.N,length(E.f))

    for j in 1:length(E.f)
        for i in eachindex(E.D0)
            for m in 1:E.M
                vi[m,i,j] = (1/E.M)*(E.D[i,m]/E.D0[i])*exp(-im*kw[j]*(E.D[i,m]-E.D0[i]))
            end
        end
    end
    return SteeringMatrix(vi,kind)
end

function steeringvectors(E::Environment{T},C::Constants{T},kind::Shear) where T <: AbstractFloat
    w = 2pi*E.f
    t = Array{T}(E.M,E.N)
    xm = Array{T}(3)
    xm[3] = 0.0
    vi = Array{Complex{T}}(E.M,E.N,length(E.f))
    for i in eachindex(E.D0)
        for m in 1:E.M
            xm[1:2] .= E.micgeom[1:2,m]-E.Rxy[1:2,i]
            #f_c(x,fvec) = shear!(x,fvec,xm,C) # TODO: Maybe slightly faster
            #res = nlsolve(f_c,[1.;1.;1.])
            res = nlsolve((x,fvec)->shear!(x,fvec,xm,C), [ 1.; 1.; 1.])
            xi = res.zero
            r1 = sqrt(xi[1]^2+xi[2]^2+xi[3]^2)  # From 0 to shear interface
            d = xi[1]*C.Ma*C.c/r1
            c1 = d + sqrt(d^2+C.c^2-(C.Ma*C.c)^2)
            r2 = sqrt((xm[1]-xi[1])^2+(xm[2]-xi[2])^2+(xm[3]-xi[3])^2)
            t[m,i] = r1/c1 + r2/C.c
        end
    end
    for j in 1:length(E.f)
        for i in eachindex(E.D0)
            for m in 1:E.M  # Type II with shear layer
                vi[m,i,j] = (1/E.M)*(t[m,i]*C.c/E.D0[i])*exp(-im*w[j]*t[m,i])
            end
        end
    end
    return SteeringMatrix(vi,kind)
end

# Default steeringvector
function steeringvectors(E::Environment{T},C::Constants{T}) where T <: AbstractFloat
    steeringvectors(E,C,Type3())
end

function steeringvectors(E::Environment{T},C::Constants{T},kind) where T <: AbstractFloat
    !isa(kind, SteeringVectorType) && throw(ArgumentError("Unknown steering vector type. Use e.g. call steeringvectors(E::Environment{T},C::Constants{T},Type3())"))
    return nothing
end

# Still WIP:
function steeringvectors(E::Environment{T},C::Constants{T},kind::Uniform) where T <: AbstractFloat
    kw = 2pi*E.f/((1-C.Ma^2)*C.c)
    vi = Array{Complex{T}}(E.M,E.N,length(E.f))
    # Rescale distances to take doppler shift into account
    Rxy = E.Rxy./sqrt.(1-C.Ma.^2)
    micgeom = E.micgeom./sqrt.(1-C.Ma.^2)
    D = pairwise(Euclidean(), Rxy, micgeom)
    D0 = colwise(Euclidean(), Rxy, [0,0,0])
    for j in 1:length(E.f)
        for i in eachindex(E.D0)
            for m in 1:E.M
                vi[m,i,j] = (1/E.M)*(D[i,m]/D0[i])*exp(-im*kw[j]*((1-C.Ma)*D[i,m]-D0[i]))
            end
        end
    end
    return SteeringMatrix(vi,kind)
end

# If dispatch of typeinstace
#function steeringvectors(E::Environment{T},C::Constants{T},kind::AbstractString) where T <: AbstractFloat
#    steeringvectors(E,C,typeinstance(kind))
#end
