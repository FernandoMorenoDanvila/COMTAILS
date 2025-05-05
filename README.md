# COMTAILS: COMet dust TAIL Simulator
### This FORTRAN parallel code is created by **Fernando Moreno**, (fernando at iaa.es) of the Instituto de Astrofísica de Andalucía, CSIC, Granada, Spain.
### It is desgined to generate synthetic dust tail images from comets and active asteroids.
### This is a Message Passing Interface (MPI) parallel version of COMTAILS that can be easily compiled with GNU free compiler mpif77 on Linux stations. FITSIO and PGPLOT libraries must be installed on the system.
### Example command to compile this parallel version:  
### mpif77 -fallow-argument-mismatch -o COMTAILS_PAR.exe COMTAILS_PAR.for heliorbit.for -L/usr/local/pgplot64 -lpgplot -/usr/X11R6/lib64 -LX11 -L. -lcfitsio -lm
### And to run the code: mpirun -np 8 ./COMTAILS_PAR.exe, where 8 is the number of available slots in your system.
### A manual describing the input files, compilation instructions, and example executions can be found in file Manual_COMTAILS.pdf (in the main files). 
## IMPORTANT:ALWAYS CHECK THE INPUT # Rec (with the JPL-Horizons system) IN FILE TAIL_INPUTS.dat, AS IT MAY CHANGE EVENTUALLY FOR THE DESIRED TARGET
### The code is described in Moreno, F., Astronomy and Astrophysics, 695, 263 (2025) 
