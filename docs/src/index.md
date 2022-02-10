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
```
using HDF5, AeroAcoustics, PyPlot
```
Open the hdf5 file from the test directory and assemble the data to a Complex array.
```
csm_file = joinpath(dirname(@__DIR__),"test","data","test1_csm.h5")
csm_ref = h5open(csm_file, "r") do file
    read(file, "CsmData/csmReal")+im*read(file, "CsmData/csmImag")
end
```
and get the associated frequencies, microphone array and measurement distance `z0`:
```
fc = h5read(csm_file, "CsmData")["binCenterFrequenciesHz"]
micgeom = h5read(csm_file, "CsmData")["arrayGeom"]
z0 = h5readattr(csm_file, "CsmData")["z0"]
```
We now have all the data to populate the `Environment` struct. First, the csm is constructed as a `FreqArray`:
```
CSM = FreqArray(csm_ref,fc)
```
which holds the array and associated frequency bins. The `Environment` is defined:
```
E = Environment(
    z0=z0,
    micgeom=micgeom,
    CSM=CSM,
    flim=(3000,4000),
    Nx = 21,
    Ny = 21,
    xlim=(-0.5,0.5),
    ylim=(-0.5,0.5)
    )
```
Where the measurement distance `z0`, the microphone geometry `micgeom`, and the csm `CSM` are required variables. 
See the optional inputs for `Environment` by typing:
```
help?> Environment
```
Now, we need to assign steering vectors (transfer functions) between the grid points defined in the environment `E` and the microphone locations in `micgeom`, this is done in a simple manner:
```
steeringvectors!(E)
```
where the "!" mutates the environment `E` and stores the steering vectors associated with the Environment. If a flow field is defined in the environment, the correct steering vectors will automatically be calculated. Next, we calculate the beamforming image:
```
b = beamforming(E)
```
the output is a `FreqArray` of size `E.Nx*E.Ny` times the number of frequency bins within the limits defined in `E`. To plot the acoustic image, reshape the beamforming result and convert to dB:
```
bdB = SPL.(reshape(b[:,1],E.Nx,E.Ny))
pcolormesh(E.rx,E.ry,bdB)
colorbar()
```
Check out `examples/Quick_start.ipynb` to see the output image.