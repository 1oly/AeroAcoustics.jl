using HDF5, AeroAcoustics, BenchmarkTools

dir = "/Users/oliver/Documents/AeroAcoustic_benchmarks/simulated/b7"
data = read(h5open(joinpath(dir,"ab7aCsmEss.h5")))
dataref = read(h5open(joinpath(dir,"ab7aCsmOpt.h5")))

micgeom = data["MetaData"]["ArrayAttributes"]["microphonePositionsM"]
csmdata = data["CsmData"]
id = 70
CSM = csmdata["csmReal"][id,:,:]+im*csmdata["csmImaginary"][id,:,:]
CSM = convert(Array{Complex{Float64},2},CSM)
fc = csmdata["binCenterFrequenciesHz"][id]

rn = micgeom'
dx = 0.05 # Should be 0.025 according to benchmark
dy = 0.05
Nx = round(Int64,(1/dx) + 1)
Ny = Nx
M = size(rn,1)

# Region of interest
x0,y0 = 0.5,0.5
rx = -x0:dx:x0
ry = -y0:dy:y0
z0 = 0.75
f = fc

Rxy = hcat([[x, y, z] for x in rx, y in ry, z in z0]...)
micgeom
using Distances
D = pairwise(Euclidean(), Rxy, micgeom)
D0 = colwise(Euclidean(), Rxy, [0,0,0])

gj = Array{Complex128}(size(D))
gj = (1/M)*(D./D0).*exp(-im*1000.*(D.-D0)/343.)

D./D0
D.-D0

Xs,Ys = (Float64[i for i in rx, j in ry],Float64[j for i in rx, j in ry])

function foo(Nx,Ny,X,Y,z0,f,rn,CSM)
    const M::Int64 = size(rn,1)    # Number of microphones
    const omega::Float64 = 2pi*f         # Angular frequency
    const c::Float64 = 343.0       # Speed of sound
    # CSM[eye(Bool,M)] = 0;        # Naive diagonal removal

    # Allocation of arrays
    gj = Array{Complex64}(M)
    gjs = Array{Complex64}(Nx,Ny,M)
    b = Array{Float64}(Nx,Ny)

    for i in 1:Nx
        for j in 1:Ny
            r0::Float64 = sqrt(X[i,j]^2 + Y[i,j]^2 + z0^2)
            for m in 1:M
                rm::Float64 = sqrt((X[i,j]-rn[m,1])^2+(Y[i,j]-rn[m,2])^2 + z0^2)
                gj[m] = (1/M)*(rm/r0)*exp(-im*omega*(rm-r0)/c) # TYPE II? Steering vector:
            end
            gjs[i,j,:] = gj
            b[i,j] = real(gj'*CSM*gj)
        end
    end
    return b, gjs
end

function foo2(D,D0,f,CSM)
    const M::Int64 = size(D,2)    # Number of microphones
    const N::Int64 = size(D,1)
    const omega::Float64 = 2pi*f         # Angular frequency
    const c::Float64 = 343.0       # Speed of sound
    # CSM[eye(Bool,M)] = 0;        # Naive diagonal removal

    # Allocation of arrays
    gj = Array{Complex128}(M,N)
    b = similar(D0)

    for i in eachindex(D0)
        for m in 1:M
            gj[m,i] = (1/M)*(D[i,m]/D0[i])*exp(-im*omega*(D[i,m]-D0[i])/c)
        end
            b[i] = real(gj[:,i]'*CSM*gj[:,i])
    end

    return b
end

function foo3(Nx,Ny,X,Y,z0,f,rn,CSM)
    const M::Int64 = size(rn,1)    # Number of microphones
    const omega::Float64 = 2pi*f         # Angular frequency
    const c::Float64 = 343.0       # Speed of sound
    # CSM[eye(Bool,M)] = 0;        # Naive diagonal removal

    # Allocation of arrays
    gj = Array{Complex64}(M)
    gjs = Array{Complex64}(Nx,Ny,M)
    b = Array{Float64}(Nx,Ny)

    # Uniform flow:
    e = [1,0,0] # Gives direction of flow
    ma = 0.0
    macostheta = ma*()

    for i in 1:Nx
        for j in 1:Ny
            r0::Float64 = sqrt(X[i,j]^2 + Y[i,j]^2 + z0^2)
            for m in 1:M
                rm::Float64 = sqrt((X[i,j]-rn[m,1])^2+(Y[i,j]-rn[m,2])^2 + z0^2)
                gj[m] = (1/M)*(rm/r0)*exp(-im*omega*(rm-r0)/c) # TYPE II? Steering vector:
            end
            gjs[i,j,:] = gj
            b[i,j] = real(gj'*CSM*gj)
        end
    end
    return b, gjs
end

@time foo(Nx,Ny,Xs,Ys,z0,f,rn,CSM)

@time foo2(D,D0,f,CSM)

foo3(Nx,Ny,Xs,Ys,z0,f,rn,CSM;)


using GR
inline("atom")
contourf(reshape(b,Nx,Ny))
