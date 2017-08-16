using HDF5, AeroAcoustics

dir = "/Users/oliver/Documents/AeroAcoustic_benchmarks/simulated/b7"
data = read(h5open(joinpath(dir,"ab7aCsmEss.h5")))
dataref = read(h5open(joinpath(dir,"ab7aCsmOpt.h5")))

micgeom = data["MetaData"]["ArrayAttributes"]["microphonePositionsM"][1:2,:]
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

Xs,Ys = (Float64[i for i in rx, j in ry],Float64[j for i in rx, j in ry])

#function beamformer{T}(Nx::Int64,Ny::Int64,X,Y,z0,f::Float64,rn,CSM::Array{T,2};psf=false)
M = size(rn,1)    # Number of microphones
const omega = 2pi*f     # Angular frequency
const c = 343           # Speed of sound

    #r0 = sqrt.(X.^2 + Y.^2 .+ z0^2)

    # Allocation of arrays
gj = Array{Complex128}(M)
gjs = Array{Complex128}(Nx,Ny,M)
b = Array{Float64}(Nx,Ny)
    # Compute transfer functions
    #
#for i in 1:Nx
#    for j in 1:Ny
i = 1
j = 1
r0 = sqrt(Xs[i,j]^2 + Ys[i,j]^2 + z0^2)
#        Threads.@threads for m in 1:M
            #rmn[i,j,m] = sqrt.((X[i,j]-rn[m,1])^2+(Y[i,j]-rn[m,2])^2 + z0^2)
            #gj[i,j,m] = (r0/rmn[i,j,m])*exp(im*omega*(rmn[i,j,m]-r0)/c)
rm = sqrt((Xs[i,j]-rn[1,1])^2+(Ys[i,j]-rn[1,2])^2 + z0^2)
gj = (1/(Nx*Ny))*(rm/r0)*exp(-im*omega*(rm-r0)/c)

gjs[i,j,:] = gj
b[i,j] = real(gj'*CSM*gj)
real(gj'*CSM*gj)
