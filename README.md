# AeroAcoustics

[![][pipeline-img]][pipeline-url] [![][travis-img]][travis-url] [![][appveyor-img]][appveyor-url] [![][coverage-img]][pipeline-url] [![][docs-img]][docs-url]

A [julia](http://julialang.org) package for Aeroacoustics.

## Overview

This package provide methods for working with microphone array measurements.
Centered around frequency-domain beamforming, methods for source localization and
quantification have been collected here. The package is under active development.

Here is a small [demonstration](https://nbviewer.jupyter.org/urls/gitlab.windenergy.dtu.dk/ollyl/AeroAcoustics.jl/raw/master/examples/Introduction.ipynb) of the current functionality.

Another noteworthy library for microphone array measurements is [Acoular](http://www.acoular.org), written in Python. AeroAcoustics.jl draws inspiration from Acoular but focusses only on the processing of *measurement data*, while Acoular also has functionality for generating simulations.

## Installation

This package is not yet registered but can be installed with

```
pkg> add https://gitlab.windenergy.dtu.dk/ollyl/AeroAcoustics.jl
```

## Contribution
Contributions are welcome, the roadmap and todos are tracked in the [roadmap meta-issue](https://gitlab.windenergy.dtu.dk/ollyl/AeroAcoustics.jl/issues/1).


[pipeline-img]: https://gitlab.windenergy.dtu.dk/ollyl/AeroAcoustics.jl/badges/master/pipeline.svg
[pipeline-url]: https://gitlab.windenergy.dtu.dk/ollyl/AeroAcoustics.jl/commits/master

[coverage-img]: https://gitlab.windenergy.dtu.dk/ollyl/AeroAcoustics.jl/badges/master/coverage.svg

[travis-img]: https://travis-ci.org/1oly/AeroAcoustics.jl.svg?branch=master
[travis-url]: https://travis-ci.org/1oly/AeroAcoustics.jl

[appveyor-img]: https://ci.appveyor.com/api/projects/status/l5igmy3p3q5f4n6d/branch/master?svg=true
[appveyor-url]: https://ci.appveyor.com/project/1oly/aeroacoustics-jl/branch/master

[docs-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-url]: https://ollyl.pages.windenergy.dtu.dk/AeroAcoustics.jl/
