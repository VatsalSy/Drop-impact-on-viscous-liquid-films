# Drop Impact on Viscous Liquid Films

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://zenodo.org/badge/DOI/10.1017/jfm.2023.13.svg)](https://doi.org/10.1017/jfm.2023.13)
[![GitHub release](https://img.shields.io/github/v/release/VatsalSy/Drop-impact-on-viscous-liquid-films)](https://github.com/VatsalSy/Drop-impact-on-viscous-liquid-films/releases)
[![Basilisk](https://img.shields.io/badge/Powered%20by-Basilisk-green)](http://basilisk.fr)
[![OpenMP](https://img.shields.io/badge/Parallel-OpenMP-orange)](http://openmp.org)
[![MPI](https://img.shields.io/badge/Parallel-MPI-orange)](https://www.open-mpi.org)
[![Issues](https://img.shields.io/github/issues/VatsalSy/Drop-impact-on-viscous-liquid-films)](https://github.com/VatsalSy/Drop-impact-on-viscous-liquid-films/issues)

## Overview
This repository contains the code for simulating drop impacts on viscous liquid films. The code used in the paper "Drop impact on viscous liquid films" by Sanjay et al. (2022) is available as version [v1.0.0](https://github.com/VatsalSy/Drop-impact-on-viscous-liquid-films/releases/tag/v1.0.0).

## Latest Version (v2.0)
Released on January 24, 2025, this version includes several improvements and changes from the published version:

### New Features
- Independent control of film density separate from drop density
- Independent control of surface tension at the film-air interface
- More flexible parameter configuration
- Added `reset_install_requirements.sh`, helping ensure codebase compatibility with the latest Basilisk version

### Repeating Variables
1. Density of the drop
2. Radius of the drop
3. Surface tension of the drop-air interface

### Technical Updates
- Improved default value configuration
- Removed `omega` adaptation (incompatible with the newest Basilisk version for axisymmetric cases)
- Removed `adapt_wavelet_limited` due to incompatibility with the newest Basilisk
  - For an `adapt_wavelet_limited` implementation, please refer to: [https://github.com/comphy-lab/adapt-wavelet-limited](https://github.com/comphy-lab/adapt-wavelet-limited)
- Introduced `dirichlet(...)` for boundary conditions of `f1` and `f2`
  - For sufficiently high resolution, `f...[left] = ...;` and `f...[left] = dirichlet(...);` behave similarly. At limited resolutions, however, `dirichlet(...)` is more robust.

## Running the Code

The simulation can be run in three different modes: serial, OpenMP (shared memory parallelism), or MPI (distributed memory parallelism). Example scripts for OpenMP and MPI are provided.

### Parameters
The code accepts the following command-line arguments in order:
```
./dropFilm [MAXlevel] [tmax] [We] [Ohd] [Bo] [Ohf] [hf] [rhof] [sigma23]
```

Default values are provided in the job scripts:
- MAXlevel = 12 (grid refinement level)
- tmax = 2.5 (simulation end time)
- We = 4.0 (Weber number)
- Ohd = 0.034 (Ohnesorge number for drop)
- Bo = 0.5 (Bond number)
- Ohf = 0.670 (Ohnesorge number for film)
- hf = 0.10 (film height)
- rhof = 1.0 (film density)
- sigma23 = 1.0 (surface tension at film-air interface)

### 1. Serial Run
```shell
qcc -O2 -Wall -disable-dimensions dropFilm.c -o dropFilm -lm
./dropFilm
```

### 2. OpenMP Run
```shell
export OMP_NUM_THREADS=8  # Set desired number of threads
qcc -O2 -Wall -disable-dimensions dropFilm.c -o dropFilm -lm
./dropFilm
```
Alternatively, use the provided script:
```shell
bash job-openMP.sh
```

### 3. MPI Run
```shell
CC99='mpicc -std=c99' qcc -Wall -O2 -D_MPI=1 -disable-dimensions dropFilm.c -o dropFilm -lm
mpirun -np 4 ./dropFilm  # Adjust number of processes as needed
```
Alternatively, use the provided script:
```shell
bash job-openMPI.sh
```

## License Information
The codes are developed using [Basilisk C](http://basilisk.fr), which is a [free software program](https://en.wikipedia.org/wiki/Free_software). In that spirit, this repository is also part of the free software program.

## Contributing and Feedback

We welcome various types of contributions and feedback. Please use our issue templates to ensure your feedback is addressed effectively:

### Issue Templates

1. **🐛 Bug Reports** - [`[New Bug Report]`](../../issues/new?template=bug_report.md)
   - For reporting code issues, compilation errors, or unexpected behavior
   - Include system information and steps to reproduce
   - Attach relevant logs and screenshots

2. **✨ Feature Requests** - [`[New Feature Request]`](../../issues/new?template=feature_request.md)
   - Suggest new features or improvements
   - Describe the problem you're trying to solve
   - Share implementation ideas if you have any

3. **🔬 Scientific Comments** - [`[New Scientific Comment]`](../../issues/new?template=scientific_comment.md)
   - Discuss scientific aspects of the simulation
   - Share observations about numerical methods
   - Suggest improvements to physical modeling
   - Provide references and supporting evidence

4. **💭 Open Discussions** - [`[New Discussion]`](../../issues/new?template=open_discussion.md)
   - Start general discussions about the project
   - Ask questions about implementation details
   - Share your experience using the code

Click on the links above to create a new issue using the corresponding template. This helps us better understand and address your feedback.
