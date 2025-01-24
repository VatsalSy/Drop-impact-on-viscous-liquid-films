# Drop Impact on Viscous Liquid Films

## Overview
This repository contains the code for simulating drop impacts on viscous liquid films. The code used in the paper "Drop impact on viscous liquid films" by Sanjay et al. (2022) is available as version [v1.0.0](https://github.com/VatsalSy/Drop-impact-on-viscous-liquid-films/releases/tag/v1.0.0).

## Latest Version (v2.0)
Released on January 24, 2025, this version includes several improvements and changes from the published version:

### New Features
- Independent control of film density separate from drop density
- Independent control of surface tension at the film-air interface
- More flexible parameter configuration

### Repeating Variables
1. Density of the drop
2. Radius of the drop
3. Surface tension of the drop-air interface

### Technical Updates
In version 2.0:
- Improved default value configuration
- Removed `omega` adaptation (incompatible with newest Basilisk version for axi cases)
- Removed `adapt_wavelet_limited` due to incompatibility with newest Basilisk
  - For `adapt_wavelet_limited` implementation, please refer to: https://github.com/comphy-lab/adapt-wavelet-limited
- Introduced `dirichlet(...)` for setting boundary conditions of `f1` and `f2`
  - While both `f...[left] = ...;` and `f...[left] = dirichlet(...);` work similarly at high resolutions, `dirichlet(...)` performs better at limited resolutions

## License Information
The codes are developed using [Basilisk C](http://basilisk.fr), which is a [free software program](https://en.wikipedia.org/wiki/Free_software). In that spirit, this repository is also part of the free software program.
