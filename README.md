# AeroAcoustics

[![][pipeline-img]][pipeline-url] [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][coverage-img]][pipeline-url] [![][docs-img]][docs-url]

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


[pipeline-img]: https://gitlab.windenergy.dtu.dk/ollyl/AeroAcoustics.jl/badges/master/pipeline.svg
[pipeline-url]: https://gitlab.windenergy.dtu.dk/ollyl/AeroAcoustics.jl/commits/master

[coverage-img]: https://gitlab.windenergy.dtu.dk/ollyl/AeroAcoustics.jl/badges/master/coverage.svg

[docs-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-url]: https://ollyl.pages.windenergy.dtu.dk/AeroAcoustics.jl/
