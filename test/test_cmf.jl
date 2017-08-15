using HDF5, AeroAcoustics

dir = "/Users/oliver/Documents/AeroAcoustic_benchmarks/simulated/b7"
data = read(h5open(joinpath(dir,"ab7aCsmEss.h5")))
dataref = read(h5open(joinpath(dir,"ab7aCsmOpt.h5")))

micgeom = data["MetaData"]["ArrayAttributes"]["microphonePositionsM"][1:2,:]
csmdata = data["CsmData"]
id = 70
CSM = csmdata["csmReal"][id,:,:]+im*csmdata["csmImaginary"][id,:,:]
fc = csmdata["binCenterFrequenciesHz"][id]

rn = micgeom'
dx = 0.05 # Should be 0.025 according to benchmark
dy = 0.05
Ns = round(Int64,(1/dx) + 1)
M = size(rn,1)

# Region of interest
x0,y0 = 0.5,0.5
rx = -x0:dx:x0
ry = -y0:dy:y0
z0 = 0.75
f = fc

Xs,Ys = (Float64[i for i in rx, j in ry],Float64[j for i in rx, j in ry])

b,gjs,PSF = beamformer(Ns,Xs,Ys,z0,f,rn,CSM,psf=true)

Gh = reshape(gjs,(Ns*Ns,M))
G = Gh'
N = Ns^2
X = eye(Float32,N)
L = norm(Gh*G,2)^2
t = 1/L
lambda = 1e-5
kmax = 1
GhG = Array{Complex{Float32},2}(Gh*G)
GhCSMG = Array{Complex{Float32},2}(Gh*CSM*G)
cmf(GhG,X,GhCSMG,t,lambda,1);
@time cmf(GhG,X,GhCSMG,t,lambda,kmax);

#x = (reshape(diag(real(Xhat)),(Ns,Ns)))
#heatmap(x,title="CMF")
