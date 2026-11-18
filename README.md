# COMTAILS-POL

**COMTAILS-POL** generates synthetic images of the dust tails and comae of comets and other active small bodies (e.g. active asteroids): brightness, optical depth and the degree of linear polarization of the scattered light. It is the successor of COMTAILS.

Author: Fernando Moreno, Instituto de Astrofísica de Andalucía (IAA-CSIC), Granada, Spain.

## Download

The complete package is attached to the latest release as a tarball. It contains the tutorial, the FORTRAN and Python codes, and example input files.

**→ [Releases](https://github.com/FernandoMorenoDanvila/COMTAILS/releases)**

```bash
tar xzf COMTAILS-POL_v1.0.tar.gz
```

## What the code does

At each time step between the onset of activity and the observation date, dust is ejected from the nucleus. The ejection can be isotropic, from the sunward hemisphere, or from an active area on a rotating nucleus. Each grain follows its own Keplerian orbit under solar gravity and radiation pressure. Its projected sky-plane position at the observation date is added to the output images. The orbits of the target and of the observer are taken from JPL-Horizons.

New in COMTAILS-POL:

- **Size-dependent scattering properties** read from a database of geometric albedo and degree of linear polarization as a function of grain size and phase angle. Laboratory (CoDuLab) or numerical databases (e.g. TAMUdust2020 by Saito et al. 2021) data can be used.
- **Polarization maps**: total, polarized and brightness-limited, formed in a PSF-consistent way. They include the nucleus, with a trigonometric polarization phase law.
- **Aperture polarimetry**, together with Afρ and magnitude, and the **size distribution of the dust within the aperture**.
- **Random grain radius within each size bin**, which removes discrete syndyne artefacts.
- **Optical depth with the size-dependent extinction efficiency**, plus a corrected normalization.
- **Ephemeris cache**, a faster separable convolution, faster orbit and Kepler solvers, and near-linear MPI scaling.

## Versions

The two versions use the same input and output files, and their results agree within the Monte Carlo noise.

| | FORTRAN | Python |
|---|---|---|
| Files | `COMTAILS-POL.for`, `heliorbit.for` | `comtails.py`, `comtails_core.py` |
| Parallelization | MPI | mpi4py |
| Requirements | FORTRAN compiler (e.g. gfortran), MPI library (Open MPI or MPICH), CFITSIO; PGPLOT + X11 only for optional on-screen graphics; `curl` and `wget` | Python ≥ 3.10, numpy, numba, astropy, requests; mpi4py + MPI library for parallel runs |
| Run | `mpirun -np 8 ./COMTAILS-POL.exe` | `mpirun -np 8 python comtails.py` |

FORTRAN compilation, for example with gfortran and Open MPI:

```bash
mpif90 -O3 -march=native -std=legacy -fallow-argument-mismatch \
       -o COMTAILS-POL.exe COMTAILS-POL.for heliorbit.for \
       -lpgplot -lX11 -lcfitsio -lm
```

Python installation:

```bash
pip install numpy numba astropy requests mpi4py
```

See the tutorial (`Manual_COMTAILS-POL.pdf`) for a full description of the input files, outputs, compilation options and examples.

## Important

- Check that the target record number (line 1 of `TAIL_INPUTS.dat`) gives the desired object and orbit in JPL-Horizons, since these identifiers may change.
- Input files from COMTAILS must be updated: line 26 (scattering model) is new, and so are lines 27–30 when the scattering database is used. See the tutorial.

## Earlier versions

The original COMTAILS code, described in Moreno (2025), is still available:

- FORTRAN serial version: branch [`FORTRAN_SERIAL`](https://github.com/FernandoMorenoDanvila/COMTAILS/tree/FORTRAN_SERIAL)
- FORTRAN MPI version: branch [`FORTRAN_PARALLEL`](https://github.com/FernandoMorenoDanvila/COMTAILS/tree/FORTRAN_PARALLEL)
- Python version: repository [py_COMTAILS](https://github.com/FernandoMorenoDanvila/py_COMTAILS)

## IMPORTANT NOTE ON EARLIER VERSION: 
In the previous version of the code, the optical depth was overestimated by a
factor $N_\mathrm{ev}$, the number of directional Monte Carlo events. Each
directional sample represents $N_\mathrm{bin}/N_\mathrm{ev}$ real particles of
its size bin, a weight that was correctly applied to the scattered flux but
omitted in the computation of the optical depth. 

Also, in the previous parallel version of the code, the file
  $\mathrm{nm.dat}$, which gives the sky-plane coordinates when IPRN=1 (line 24 of file TAIL_INPUTS.dat), radius and ejection
  time of every Monte Carlo particle, contained only part of the particles,
  because it was written by a single process.

## Citation

If you use this code, please cite:

Moreno, F. 2025, *A&A*, 695, A263, COMetary dust TAIL Simulator (COMTAILS): A computer code to generate comet dust tail brightness images

## Contact

Please report bugs and questions through the [Issues](https://github.com/FernandoMorenoDanvila/COMTAILS/issues) page of this repository, or directly to fernando (at) iaa.es.
