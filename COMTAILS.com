gfortran COMTAILS.for heliorbit.for -o  COMTAILS.exe  libcfitsio.a  -L/usr/local/pgplot64 -lpgplot -L/usr/X11R6/lib64 -LX11 
