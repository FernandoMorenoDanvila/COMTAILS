mpif77 -fallow-argument-mismatch -o COMTAILS_PAR.exe COMTAILS_PAR.for heliorbit.for -L/usr/local/pgplot64 -lpgplot -L/usr/X11R6/lib64 -LX11 -L. -lcfitsio -lm
