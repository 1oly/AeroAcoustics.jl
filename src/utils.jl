function SPL{T}(p::Array{T})
    s = similar(p)
    s[p.>0] = 10*log10.(p[p.>0]/4e-10)
    s[p.<=0] = -350
    return s
end
SPL(p::Number) = 10*log10(p/4e-10)

function shear!(xi,fvec,xm,C::Constants{T}) where T <: AbstractFloat
    a = sqrt(xi[1]^2+(1-C.Ma^2)*(xi[2]^2+xi[3]^2))
    b = sqrt((xm[1]-xi[1])^2+(xm[2]-xi[2])^2+(xm[3]-xi[3])^2)
    fvec[1] = (1/a)*xi[1]-(1/b)*(1-C.Ma^2)*(xm[1]-xi[1]) - C.Ma
    fvec[2] = (1/a)*xi[2]-(1/b)*(xm[2]-xi[2])
    fvec[3] = xi[3] - C.h
end

function octavebandlimits(fc,kind)
    fl,fu = similar(fc)
    const C = 10^(3./10.)
    fl = fc*C^(-1/(2*kind))
    fu = fc*C^(1/(2*kind))
    return fl, fu
end

function promote_array(arrays...)
    # see: https://stackoverflow.com/a/31235569
    eltype = Base.promote_eltype(arrays...)
    tuple([convert(Array{eltype}, array) for array in arrays]...)
end

function promote_array(T,arrays...)
    # see: https://stackoverflow.com/a/31235569
    tuple([convert(Array{T}, array) for array in arrays]...)
end

function parseHDF5data(filename::AbstractString)
    parseHDF5data(Float64,filename::AbstractString)
end

function parseHDF5data(::Type{T},filename::AbstractString) where T<:AbstractFloat
    CSM = h5open(filename, "r") do file
        read(file, "CsmData")
    end
    CSMreal,CSMimag,binfc = promote_array(T,CSM["csmReal"],CSM["csmImaginary"],CSM["binCenterFrequenciesHz"])
    return CrossSpectralMatrix(CSMreal,CSMimag,binfc,false)
end

function beamformersetup(dx::T,dy::T,x::Vector{T},y::Vector{T},z::T,f::Vector{T},micgeom::Matrix{T},csmdata::CrossSpectralMatrix{T}) where T <: AbstractFloat
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
        micgeom = [micgeom; zeros(T,M)']
    end
    D = pairwise(Euclidean(), Rxy, micgeom) # Distance from each mic to grid points
    const N = size(D,1)
    fc = csmdata.binCenterFrequenciesHz
    fn = fc[(fc.>=fl) .& (fc.<=fu)]
    ind = findin(fc,fn)
    C = CrossSpectralMatrix(csmdata.csmReal[ind,:,:],csmdata.csmImag[ind,:,:],fn,false)
    Nf = length(ind)
    return Environment(N,M,Nx,Ny,Nz,Nf,fn,micgeom,rx,ry,rz,Rxy,D0,D,C)
end

function beamformersetup(dx::T,dy::T,x::Vector{T},y::Vector{T},z::T,f::Vector{T},micgeom::Matrix{T},csmdata::AbstractString) where T<:AbstractFloat
    time_data = parseHDF5data(T,csmdata)
    beamformersetup(dx,dy,x,y,z,f,micgeom,time_data)
end


function beamformersetup(dx::A,dy::B,x::Vector{C},y::Vector{D},z::E,f::Vector{F},micgeom::Matrix{G},csmdata::AbstractString) where {A,B,C,D,E,F,G}
    Tc = promote_type(A,B,C,D,E,F,G)
    time_data = parseHDF5data(Tc,csmdata)
    beamformersetup(Tc(dx),Tc(dy),Vector{Tc}(x),Vector{Tc}(y),Tc(z),Vector{Tc}(f),Matrix{Tc}(micgeom),time_data)
end


function sourceintegration(res::Array{T,3},SourcePositions::S,E::Environment{T},fco::Array{T,1},intarea::Tuple) where {T <: AbstractFloat, S<: Dict}
    SourceInt = Dict{String, Array{T,1}}()
    SourceInt["fco"] = fco
    fl,fu = octavebandlimits(fco,3) # Calculate third-octave band limits
    idx,idy = Int.(round.(1/E.rx.step.hi)),Int.(round.(1/E.ry.step.hi))
    dxint,dyint = Int.(round.(intarea[1]/2E.rx.step.hi)),Int.(round.(intarea[2]/2E.ry.step.hi))
    srcint = similar(fco)
    for (key,value) in SourcePositions
        indx,indy = Int.(round.(idx*abs(E.rx[1]-value["xy"][1])))+1,Int.(round.(idy*abs(E.ry[1]-value["xy"][2])))+1
        for i in 1:length(fco)
            fn = E.f[(E.f.>=fl[i]) .& (E.f.<=fu[i])]
            ind = findin(E.f,fn)
            copy!(srcint[i],SPL(sum(res[indx-dxint:indx+dxint,indy-dyint:indy+dyint,ind])))
        end
        SourceInt[key] = srcint
    end
    return SourceInt
end
