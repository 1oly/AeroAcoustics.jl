using HDF5, AeroAcoustics, TraceCalls

dir = "/Users/oliver/Documents/AeroAcoustic_benchmarks/simulated/b7"
data = read(h5open(joinpath(dir,"ab7aCsmEss.h5")))
dataref = read(h5open(joinpath(dir,"ab7aCsmOpt.h5")))

micgeom = data["MetaData"]["ArrayAttributes"]["microphonePositionsM"]
csmdata = data["CsmData"]

dx, dy = 0.025,0.025
x = [-0.5,0.5]
y = x
z = 0.75
f = [100.,300.]

Env = beamformersetup(dx,dy,x,y,z,f,micgeom,csmdata)
Const = Constants(0.0,343.0)
V = steeringvectors(Env,Const,"II")
@trace beamformer2(Env,Const,V;psf=true)
