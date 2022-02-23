---
title: 'AeroAcoustics.jl: A julia package for aeroacoustics'
tags:
  - julia
  - acoustics
  - aerodynamics
  - acoustic imaging
  - beamforming
authors:
  - name: Oliver Lylloff^[corresponding  author] 
    orcid: 0000-0002-1596-107X
    affiliation: 1
affiliations:
 - name: Technical University of Denmark (DTU), Dept. of Wind Energy
   index: 1
date: 10 February 2022
bibliography: paper.bib
---

# Summary

Flow-induced noise is an important research topic in automotive, aerospace and
wind energy sciences. Aeroacoustics is the interdiciplinary research field between
aerodynamics and acoustics, and can stretch from numerical simulations to 
experimental wind tunnel tests. In the experimental domain, the use of wind tunnels
is a common strategy for the assessment of aerodynamic performance and noise. 
The aeroacoustic noise, arising from the test item, must be isolated in order to obtain
a trustworthy evaluation due to high levels of background noise. A common
strategy is to use *acoustic images*, that are computed using microphone arrays,
to visualize and quantify the spatial distribution of noise. These methods for computing 
acoustic images often use the same code base and utility functions, and can
benefit from a modular software package to benchmark exisisting methods and develop new ones.

# Statement of need

`AeroAcoustics.jl` is a julia [@Bezanson2017] package for working with
microphone arrays and acoustic imaging (e.g., beamforming). 

The research community around acoustic imaging have developed benchmark test 
cases and compared results across different research groups [@Bahr2017;@Sarradj2017].
However, the software code was not compared or shared online.


# References