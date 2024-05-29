# AeroAcoustics.jl
[![CI](https://github.com/1oly/AeroAcoustics.jl/workflows/CI/badge.svg)](https://github.com/1oly/AeroAcoustics.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/1oly/AeroAcoustics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/1oly/AeroAcoustics.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://1oly.github.io/AeroAcoustics.jl/dev)
[![status](https://joss.theoj.org/papers/9e20f1ec29f69e94bf0c9f1d2c22fa0d/status.svg)](https://joss.theoj.org/papers/9e20f1ec29f69e94bf0c9f1d2c22fa0d)

A [Julia](http://julialang.org) package for Aeroacoustics and acoustic imaging.

<p align="center">
  <a href="#installation">Installation</a> •
  <a href="#tests">Tests</a> •
  <a href="#documentation">Documentation</a> •
  <a href="#citation">Citation</a>
</p>

![Image](presentation.png?raw=true "Title")

## Overview

This package provide methods for working with microphone array measurements. Utilities for processing
beamforming and other acoustic imaging methods are collected in this package, and it is the intention, that 
a suite of both well-known and state of the art methods is continuously updated. 

The current set of methods include conventional frequency domain beamforming (CBF) with source power integration (SPI), and the following advanced methods:

| Method  | Reference  |
|---------|------------|
|DAMAS   | Brooks, T. F. et al. (2006). *A deconvolution approach for the mapping of acoustic sources (DAMAS) determined from phased microphone arrays*. J. Sound Vib. 294(4), 856–879. https://doi.org/10.1016/j.jsv.2005.12.046  |
|CLEAN-SC   |  Sijtsma, P. (2007). *CLEAN based on spatial source coherence*. Int.J.Aeroacou. 6(4), 357–374. https://doi.org/10.1260/147547207783359459 |
|FISTA   | - Beck, A., & Teboulle, M. (2009). *A fast iterative shrinkage-thresholding algorithm for linear inverse problems.* SIAM J. Imag. Sci., 2(1), 183-202. https://doi.org/10.1137/080716542 <br />- Lylloff, O., et al. (2015). *Improving the efficiency of deconvolution algorithms for sound source localization*. J. Acou. Soc. Am, 138(1), 172-180. https://doi.org/10.1121/1.4922516 |

On the roadmap is:   

| Method  | Reference  |
|---------|------------|
Functional beamforming | Dougherty, R.P. (2014). *Functional beamforming*. 5th Berlin Beamforming Conference, February 19–20 2014, Berlin, Germany, GFaI, e.V., Berlin. |
Orthogonal beamforming | Sarradj, E. (2010). *A fast signal subspace approach for the determination of absolute levels from phased microphone array measurements*. J. Sound Vib. 329(9), 1553–1569. https://doi.org/10.1016/j.jsv.2009.11.009|
Spectral Estimation Method (SEM) / <br> Covariance Matrix Fitting (CMF) | Blacodon, D. et al. (2004). *Level estimation of extended acoustic sources using a parametric method*. J. Airc. 41, 1360–1369. https://doi.org/10.2514/1.3053 <br> Yardibi, T. et al. (2010). *A covariance fitting approach for correlated acoustic source mapping*. J. Acoust. Soc. Am. 127(5), 2920–2931. https://doi.org/10.1121/1.3365260|
Generalized Inverse Beamforming (GIBF) | Suzuki, T. (2011). *L1 generalized inverse beamforming algorithm resolving coherent/incoherent, distributed and multipole sources*. J. Sound Vib. 330(24), 5835–5851. https://doi.org/10.1016/j.jsv.2011.05.021|

Additional methods can also be added by contributors to this repository. The code structure enables an easy and modular addition of new methods. 

Denoising of the cross-spectral matrix is another important part of succesful acoustic imaging. Several different methods are on the roadmap to be implemented.

Source integration of acoustic images is another important feature of this package. 
The output can be produced in narrow-band, 1/3rd or 1/12th octave bands. 

## Installation
First install [Julia](http://julialang.org) and start Julia in a terminal, [VS code](https://www.julia-vscode.org), [Jupyter](https://github.com/JuliaLang/IJulia.jl) or another application that can run Julia. This package is registered and can be installed with

```
pkg> add AeroAcoustics
```
the package manager `pkg>` can be accessed by typing `]`.
## Tests
This package constains tests used for CI, but can also be used to check if the package is working properly when installed. To run the tests, after adding the package, activate the package manager, by typing `]`, and write
```
pkg> test AeroAcoustics
```
The test suite will download an external file, that is stored in `test/data`.
## Documentation
The [documentation](https://1oly.github.io/AeroAcoustics.jl/dev/) gives an introduction to the packages and the API.

## Community and contributions
Any question regarding the installation, use or extension of the code can be posted in the [Github issue tracker](https://github.com/1oly/AeroAcoustics.jl/issues).
If you want to add an new algorithm, you can fork this package and start developing your code and test it. If you wish to contribute to the general development, you can look for open issues or reach out to discuss further. Contributions are very welcome!

## Related packages
Another noteworthy library for microphone array measurements is [Acoular](http://www.acoular.org), written in Python. AeroAcoustics.jl draws inspiration from Acoular but focusses only on the processing of *measurement data*, while Acoular also has functionality for generating simulated data.

## Citation
If you use this package in your work, please cite the following:
```  
@article{lylloff2024aeroacoustics, 
doi = {10.21105.joss.06390}, 
url = {https://doi.org/10.21105/joss.06390}, 
year = {2024}, 
publisher = {The Open Journal}, 
journal = {Journal of Open Source Software},
volume = {9}, 
number = {97}, 
pages = {6390}, 
author = {Oliver Lylloff}, 
title = {AeroAcoustics.jl: A Julia package for aeroacoustics} }
``` 