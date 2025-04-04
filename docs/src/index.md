# AeroAcoustics

AeroAcoustics.jl is a package for post-processing of microphone array measurements.
It aims to provide basic functionality for computing cross-spectral matrices and frequency-domain beamforming.

The package repository is here: <https://github.com/1oly/AeroAcoustics.jl>

This package has been developed over several years, but has reached a stable state. A brief introduction to the usage is given below.

## Installation
From the julia REPL:  
```
pkg> add AeroAcoustics
```   
and optionally running the tests:  
```
pkg> test AeroAcoustics
```  
This runs the tests and downloads a measurement file to the `/data` directory.

## Quick start
The package is structured around the struct `Environment`, that takes all neccessary variables associated with an acoustic image as input. We will use the measurement data downloaded for the unit tests, which is located in `/data`. This file is a cross-spectral matrix (csm), a frequency domain representation of the measurement. This example is also available as a Jupyter notebook in the `examples` directory. First we import packages for reading data and plotting:
```julia
using HDF5, AeroAcoustics, PyPlot
```
Open the hdf5 file from the test directory and assemble the data to a Complex array.
```julia
csm_file = joinpath(replace(dirname(pathof(AeroAcoustics)),"src"=>"test"),"data","test1_csm.h5")
csm_ref = h5open(csm_file, "r") do file
    read(file, "CsmData/csmReal")+im*read(file, "CsmData/csmImag")
end
```
and get the associated frequencies, microphone array and measurement distance `z0`:
```julia
fc = h5read(csm_file, "CsmData")["binCenterFrequenciesHz"]
micgeom = h5read(csm_file, "CsmData")["arrayGeom"]
z0 = h5readattr(csm_file, "CsmData")["z0"]
```
We now have all the data to populate the `Environment` struct. First, the csm is constructed as a `FreqArray`:
```julia
CSM = FreqArray(csm_ref,fc)
```
which holds the array and associated frequency bins. The `Environment` is defined:
```julia
E = Environment(
    z0=z0,
    micgeom=micgeom,
    CSM=CSM,
    flim=(3000,4000),
    Nx = 21,
    Ny = 21,
    xlim=(-0.5,0.5),
    ylim=(-0.5,0.5),
    multi_thread = true # multi-threading can be enabled globally like this
    )
```
Where the measurement distance `z0`, the microphone geometry `micgeom`, and the csm `CSM` are required variables. 
See the optional inputs for `Environment` by typing:
```julia
help?> Environment
```
Now, we need to assign steering vectors (transfer functions) between the grid points defined in the environment `E` and the microphone locations in `micgeom`, this is done in a simple manner:
```julia
steeringvectors!(E)
```
where the "!" mutates the environment `E` and stores the steering vectors associated with the Environment. If a flow field is defined in the environment, the correct steering vectors will automatically be calculated. Next, we calculate the beamforming image:
```julia
b = beamforming(E)
```
the output is a `FreqArray` of size `E.Nx*E.Ny` times the number of frequency bins within the limits defined in `E`. To plot the acoustic image, reshape the beamforming result and convert to dB:
```julia
bdB = SPL.(reshape(b[:,1],E.Nx,E.Ny))
pcolormesh(E.rx,E.ry,bdB)
colorbar()
```
Check out `examples/Quick_start.ipynb` to see the output image.

### Advanced methods
Two advanced and widely used methods for improving the beamforming image are
Clean-SC and DAMAS, which can be easily called by:   
```julia
xSC = cleanSC(E)
xD = damas(E,b)
```
Check out `examples/Quick_start.ipynb` to see examples of the acoustic images.

### Source integration
In a typical workflow, the acoustic images are used to focus on selected regions
of the spatial domain and extract a spectrum. The source integration is used on an acoustic map, e.g., `b`, `xSC`, or `xD`, computed above.
```julia
sourceintegration(b,E,int_region)
```
where `b` is the source map, `E` the environment struct and `int_region` describe the limits of
a square to integrate over. A utility function `AeroAcoustics.point_to_region` can help define limits 
by giving a point and extent as input, e.g.,   
```julia
src_pos = (0.0, 1.0) # Point to integrate
dxdy = (0.5,0.5) # Size of square extenting from src_pos 
int_region = AeroAcoustics.point_to_region(src_pos,dxdy)
4-element Vector{Float64}:
 -0.25
  0.25
  0.75
  1.25
```

```
