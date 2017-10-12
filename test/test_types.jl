using HDF5, AeroAcoustics

dir = "/Users/oliver/Documents/AeroAcoustic_benchmarks/simulated/b7"
data = read(h5open(joinpath(dir,"ab7aCsmEss.h5")))
dataref = read(h5open(joinpath(dir,"ab7aCsmOpt.h5")))

micgeom = data["MetaData"]["ArrayAttributes"]["microphonePositionsM"]
csmdata = data["CsmData"]

dx, dy = 0.025,0.025
x = [-0.5,0.5]
y = x
z = 0.75
f = [100.,5000.]

Env = beamformersetup(dx,dy,x,y,z,f,micgeom,csmdata)
Const = Constants(0.0,343.0)
V = steeringvectors(Env,Const,"II")

@time beamformer2(Env,Const,V)

Profile.clear()
@profile beamformer2(Env,Const,V)
open("profile.txt", "w") do s
    Profile.print(IOContext(s, :displaysize => (300, 2000)))
end

#=
using GR
inline("atom")
figure(size=(800,300))
subplot(1,2,1)
heatmap(b[:,:,70])
subplot(1,2,2)
heatmap(PSF[:,:,70])

using ImageFiltering

Env.N
d = padarray(b[:,:,10], Fill(zero(eltype(b)),(2Env.N-1,2Env.N-1)))
figure(size=(300,300))
heatmap(d.parent)
size(d.parent)
Env.N

2Env.N-1-Env.N

d
Fill(zero(eltype(b)),(2Env.N-1,2Env.N-1))
Env.M

2pi*Env.f/Const.c
=#
