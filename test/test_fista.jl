using HDF5, AeroAcoustics, GR
inline("atom")

dir = "/Users/oliver/Documents/AeroAcoustic_benchmarks/simulated/b7"
data = read(h5open(joinpath(dir,"ab7aCsmEss.h5")))
dataref = read(h5open(joinpath(dir,"ab7aCsmOpt.h5")))

micgeom = data["MetaData"]["ArrayAttributes"]["microphonePositionsM"][1:2,:]
csmdata = data["CsmData"]
CSM = csmdata["csmReal"]+im*csmdata["csmImaginary"]
fc = csmdata["binCenterFrequenciesHz"]

rn = micgeom'
dx = 0.05 # Should be 0.025 according to benchmark
dy = 0.05
M = size(rn,1)

# Region of interest
x0,y0 = 0.5,0.5
rx = -x0:dx:x0
ry = -y0:dy:y0
Nx = length(rx)
Ny = length(ry)
z0 = 0.75
fl = 2500
fu = 4000
f = fc[(fc.>fl) .& (fc.<fu)]
ind = findin(fc,f)
CSM = convert(Array{Complex{Float64},3},CSM[ind,:,:])

Xs,Ys = (Float64[i for i in rx, j in ry],Float64[j for i in rx, j in ry])

b,gjs,PSF = beamformer(Nx,Ny,Xs,Ys,z0,f,rn,CSM,psf=true)

X0 = zeros(Nx,Ny)
kmax = 100
x_fista = fista(PSF, b, X0,kmax)
x_fista
x_fista_sum = reshape(sum(x_fista,3),Nx,Ny)
