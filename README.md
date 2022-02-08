# AeroAcoustics.jl
[![CI](https://github.com/1oly/AeroAcoustics.jl/workflows/CI/badge.svg)](https://github.com/1oly/AeroAcoustics.jl/actions?query=workflow%3ACI)
[![Codecov](https://codecov.io/gh/1oly/AeroAcoustics.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/1oly/AeroAcoustics.jl)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://1oly.github.io/AeroAcoustics.jl/dev)

A [julia](http://julialang.org) package for Aeroacoustics.

## Overview

This package provide methods for working with microphone array measurements.
Centered around frequency-domain beamforming, methods for source localization and
quantification have been collected here. The package is under active development.

Another noteworthy library for microphone array measurements is [Acoular](http://www.acoular.org), written in Python. AeroAcoustics.jl draws inspiration from Acoular but focusses only on the processing of *measurement data*, while Acoular also has functionality for generating simulations.

## Installation
First install [julia](http://julialang.org) and start julia in a terminal, [VS code](https://www.julia-vscode.org), [Jupyter](https://github.com/JuliaLang/IJulia.jl) or another application that can run julia. This package is not yet registered yet but can be installed with

```
pkg> add https://gitlab.windenergy.dtu.dk/ollyl/AeroAcoustics.jl
```
the package manager `pkg>` can be accessed by typing `]`.
## Contribution
Contributions are welcome!

