# AeroAcoustics
[![CI](https://github.com/1oly/AeroAcoustics.jl/workflows/CI/badge.svg)](https://github.com/1oly/AeroAcoustics.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/1oly/AeroAcoustics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/1oly/AeroAcoustics.jl)

A [julia](http://julialang.org) package for Aeroacoustics.

## Overview

This package provide methods for working with microphone array measurements.
Centered around frequency-domain beamforming, methods for source localization and
quantification have been collected here. The package is under active development.

Another noteworthy library for microphone array measurements is [Acoular](http://www.acoular.org), written in Python. AeroAcoustics.jl draws inspiration from Acoular but focusses only on the processing of *measurement data*, while Acoular also has functionality for generating simulations.

## Installation

This package is not yet registered but can be installed with

```
pkg> add https://gitlab.windenergy.dtu.dk/ollyl/AeroAcoustics.jl
```

## Contribution
Contributions are welcome!

