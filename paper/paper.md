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
 - name: Technical University of Denmark (DTU), Dept. of Wind and Energy Systems
   index: 1
date: 2 May 2023
bibliography: paper.bib
---

# Summary

Aeroacoustics is the interdiciplinary research field between
aerodynamics and acoustics. It is concerned with the noise created by the motion of air flowing around a solid body, and
is an important research topic in the automotive, aerospace and wind energy sciences []. 

The aeroacoustic noise, arising from the test item, is typically contaminated by the surrounding background noise, 
and must be isolated in order to obtain a trustworthy evaluation. A common
strategy is to use *acoustic images*, that are computed using microphone arrays,
to visualize and quantify the spatial distribution of noise. These methods for computing 
acoustic images often use the same code base and utility functions, and can
benefit from a modular software package to benchmark exisisting methods and develop new ones.

The AeroAcoustics.jl package provides a basis for working with microphone array measurements.

# Statement of need

`AeroAcoustics.jl` is a julia [@Bezanson2017] package for working with
microphone arrays and acoustic imaging (e.g., beamforming). 

The research community around acoustic imaging have developed benchmark test 
cases and compared results across different research groups [@Bahr2017;@Sarradj2017].
However, the software code was not compared or shared online. 

Microphone arrays and the algorithms used to create acoustic images, are also being used in domains
without flow. AeroAcoustics.jl can also be used here.

A noteworthy software packagage is Acoular [@acoular], that has similar features but written in Python.

# References