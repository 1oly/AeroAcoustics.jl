using HDF5, GR, AeroAcoustics
inline("atom")

#FFTW.set_num_threads(2)

# load data:
csm_path = "/Users/oliver/Documents/AeroAcoustic_benchmarks/Simulated/b7/ab7aCsmEss.h5"
csm_Ref = "/Users/oliver/Documents/AeroAcoustic_benchmarks/Simulated/b7/ab7aCsmEss.h5"
time_path = "/Users/oliver/Documents/AeroAcoustic_benchmarks/Simulated/b7/ab7aTimeSeries.h5"

micgeom = h5open(csm_path, "r") do file
    read(file, "MetaData/ArrayAttributes/microphonePositionsM")
end

SourcePos = Dict{String, Dict{String,Tuple}}()
SourcePos["Source3"] = Dict([("xy",(0.1,0.1))])
SourcePos["Source0"] = Dict([("xy",(0.1,-0.1))])
SourcePos["Source1"] = Dict([("xy",(-0.1,-0.1))])
SourcePos["Source2"] = Dict([("xy",(-0.1,0.1))])

# Set frequency range:
fco = 10.^(3.2:0.1:4.31)        # Center frequencies of interest
fl,fu = octavebandlimits(fco,3) # Calculate third-octave band limits

# Setup grid and limits:
dx, dy = 0.025,0.025    # Grid size
x = [-0.5,0.5]          # x limits
y = x                   # y limits
z = 0.75                # Measurement distance

# Create beamformer environment:
Env1 = beamformersetup(dx,dy,x,y,z,[fl[1];fu[end]],micgeom,csm_Ref)

# or alternatively use time_data:
#C = csm(time_path)

#Env2 = beamformersetup(dx,dy,x,y,z,[fl[1];fu[end]],micgeom,C)

# Declare constants:
Const = Constants(0.0,343.0)

# Calculate steering vectors:
V1 = steeringvectors(Env1,Const)
#V2 = steeringvectors(Env2,Const)

# Calculate beamformer output and point-spread function:
b1 = beamformer(Env1,Const,V1)
#b2 = beamformer(Env2,Const,V2)
PSF = pointspreadfunction(Env1,Const,V1)

indf = 40
figure(size=(600,600))
#subplot(1,2,1)
contourf(Env1.rx,Env1.ry,SPL(b1[:,:,indf]),title="Beamformer at $(Env1.f[indf]) Hz")
#subplot(1,2,2)
#contourf(Env1.rx,Env1.ry,SPL(b2[:,:,indf]))

# Deconvolution

# Deconvolution without boundary-conditions
X0 = zeros(Env1.Nx,Env1.Ny)
x_fistaprox = zeros(Env1.Nx,Env1.Ny,Env1.Nf)
kmax = 20
tol = -Inf
@time fistaprox!(x_fistaprox, PSF, b1, NormL1Pos(1e-2); tol=tol, maxit=kmax);
@time x_fista,objs = fista(PSF, b1, X0, kmax,tol);

figure(size=(700,300))
fin = 10
subplot(1,2,1)
contourf(Env1.rx,Env1.ry,x_fista[:,:,fin])
subplot(1,2,2)
contourf(Env1.rx,Env1.ry,x_fistaprox[:,:,fin])
