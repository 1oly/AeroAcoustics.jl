# AeroAcoustics.jl
[![CI](https://github.com/1oly/AeroAcoustics.jl/workflows/CI/badge.svg)](https://github.com/1oly/AeroAcoustics.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/1oly/AeroAcoustics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/1oly/AeroAcoustics.jl)
[![](https://img.shields.io/badge/docs-latest-blue.svg)](https://1oly.github.io/AeroAcoustics.jl/dev)

A [julia](http://julialang.org) package for Aeroacoustics and acoustic imaging.

![Image](presentation.png?raw=true "Title")

## Overview

This package provide methods for working with microphone array measurements. Utilities for processing
beamforming and other acoustic imaging methods are collected in this package, and it is the intention, that 
a suite of both well-known and state of the art methods is continuously updated. 

The current set of methods include conventional frequency domain beamforming (CBF) with source power integration (SPI), and the following advanced methods:

| Method  | Reference  |
|---------|------------|
|DAMAS   | Brooks, T. F. et al. (2006). *A deconvolution approach for the mapping of acoustic sources (DAMAS) determined from phased microphone arrays*. J. Sound. Vib. 294(4), 856–879. https://doi.org/10.1016/j.jsv.2005.12.046  |
|CLEAN-SC   |  Sijtsma, P. (2007). *CLEAN based on spatial source coherence*. Int.J.Aeroacou. 6(4), 357–374. https://doi.org/10.1260/147547207783359459 |
|FISTA   | - Beck, A., & Teboulle, M. (2009). *A fast iterative shrinkage-thresholding algorithm for linear inverse problems.* SIAM J. Imag. Sci., 2(1), 183-202. https://doi.org/10.1137/080716542 <br />- Lylloff, O., et al. (2015). *Improving the efficiency of deconvolution algorithms for sound source localization*. J. Acou. Soc. Am, 138(1), 172-180. https://doi.org/10.1121/1.4922516 |

On the roadmap is:   
- Functional/adaptive/orthogonal beamforming
- CMF

Additional methods can also be added by contributors to this repository. The code structure enables an easy and modular addition of new methods. 

Source integration of acoustic images is another important feature of this package. 
The output can be produced in narrow-band, 1/3rd or 1/12th octave bands. 

## Installation
First install [julia](http://julialang.org) and start julia in a terminal, [VS code](https://www.julia-vscode.org), [Jupyter](https://github.com/JuliaLang/IJulia.jl) or another application that can run julia. This package is registered and can be installed with

```
pkg> add AeroAcoustics
```
the package manager `pkg>` can be accessed by typing `]`.
## Contribution
Contributions are welcome! Issues are tracked on [Github issue tracker](https://github.com/1oly/AeroAcoustics.jl/issues). If you want to add an new algorithm, you can fork this package and start developing your code and test it.

## Related packages
Another noteworthy library for microphone array measurements is [Acoular](http://www.acoular.org), written in Python. AeroAcoustics.jl draws inspiration from Acoular but focusses only on the processing of *measurement data*, while Acoular also has functionality for generating simulated data.

