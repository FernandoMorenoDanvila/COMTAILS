C
C -----------------------------------------------------------------------
C CODE TO GENERATE SYNTHETIC DUST TAIL BRIGHTNESS FROM ACTIVE SMALL BODIES
C WRITTEN BY Fernando Moreno IAA-CSIC (fernando@iaa.es)      
C -----------------------------------------------------------------------

         IMPLICIT NONE
         INCLUDE 'image.inc'
         INCLUDE 'timemax.inc'
         INCLUDE 'ninputs.inc'
         INCLUDE 'starmax.inc'
C MISCELLANEOUS VARIABLES         
         CHARACTER*14 RECID
         REAL*8 PDEN,PV0,PV,PHASE_COEFF,RNUCLEUS,BRNUCLEUS
         REAL*8 PV0_NUC,PV_NUC,PHASE_COEFF_NUC
         REAL*8 POWER,START_JD,END_JD
         REAL*8 NMIN,NMAX,MMIN,MMAX,AUGX,AUGY
         REAL*8 PI,HALFPI,TWOPI,FOURPI,AUKM,CTEVEL,TORAD
         REAL*8 RADMIN,RADMAX,CP,QP,BETA,V0,GAMMA,KAPPA
         REAL*8 XCOBS,YCOBS,ZCOBS,XC_EJEC,YC_EJEC,ZC_EJEC,XE,YE,ZE
         REAL*8 DELTA,RCOBS,THETACOBS,XOUTE,YOUTE,ZOUTE
         REAL*8 PSANG,PA,PHAS_ANG,CHITAE,ETAE,GITAE,PER_JD,TC,DELTAT
         REAL*8 R1LOG,R2LOG,STEPRLOG,POWERINT,TAU
         REAL*8 XIT,YIT,ZIT,VEJ,VXIT,VYIT,VZIT
         REAL*8 XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
         REAL*8 XXP,XYP,XZP,YXP,YYP,YZP,ZXP,ZYP,ZZP
         REAL*8 INC,OM,WPER,MU,CTE_MAG,ORB_PER,HALF_PER,SEMIAXIS
         REAL*4 RANDOM1,RANDOM2,TIME_START,TIME_END,MAGSUN
         REAL*8 GRDSIZ,DMDT,XMASS,TIME,THETAC_EJEC,RC_EJEC
         REAL*8 RLOG,R1,R2,RAD,VXC_EJEC,VYC_EJEC,VZC_EJEC
         REAL*8 ENEINT,FLUX_OUT
         REAL*8 TOTAL_MASS,SFWHM,LAT,LON,COSZ,EXPOCOS
         REAL*8 RSUNKM,SOLID_ANGLE,ARCSEC_RAD,AREAPIX,CTE_PART,CTE_NUC
         REAL*8 RHO_AP,MAG,AFRHO,AFRHO_0,PCP
         INTEGER NX,NY,NSIZES,NTIMES,ITIME,IRLOG,IAP
         INTEGER K,II,JJ,NEVENT,NTOTMC,INUC,JNUC,IEJEC_MODE,ISUN,ICONV
         INTEGER KK,NSHADOW,IPRN,IGRAPHO
C VARIABLES FOR TRANSFORMATION TO SKY PLANE
         REAL*8 NMPAR,NMPAR1,NMPAR2,NMPAR3,NMPAR4,NMPAR5
         REAL*8 NMPAR6,NMPAR7,NMPAR8,MROT,NROT,MPAR,NPAR,LPAR
C VARIABLES TO BUILDUP DUST TAIL IMAGE/STAR FIELD IMAGE        
         REAL*8 NGRID(NXMAX),MGRID(NYMAX),SCALE,OPT_DEPTH_NUC
         REAL*4 FLUX(NXMAX,NYMAX),FLUX_STAR(NXMAX,NYMAX)
         REAL*4 CONVOLU(NXMAX,NYMAX),OPT_DEPTH(NXMAX,NYMAX)
         REAL*8 ED,QD,WD,IND,ND
         REAL*8 MAGLIM
         INTEGER ISTAR
         CHARACTER*80 OUTIMAGE
C JPL-EPHEMERIS DERIVED PARAMETERS
         REAL*8 XARR(3,NTIMEMAX),VARR(3,NTIMEMAX),XYZ_EARTH(3)
         REAL*8 TIMES(NTIMEMAX),EC,QR
         REAL*8 TRUEANARR(NTIMEMAX)
         REAL*8 RA,DEC,RAH,RAM,RAS,DED,DEM,DES,PLANG,XTEMP,YTEMP
         REAL*8 RA_STAR,DE_STAR,GMAG,RPMAG,BPMAG,RMAG
         REAL*8 START_D,START_M,START_Y,OBS_D,OBS_M,OBS_Y
C STAR FIELD PARAMETERS
         REAL*8 STARMAG(NSTARMAX),STARFLUX(NSTARMAX)
         INTEGER POSX(NSTARMAX),POSY(NSTARMAX)         
C PARAMETERS OF THE INPUT FILE, dM/dt, vel. factor, power index, and radii
         REAL*8 DTIME(NINPUTSMAX),DMDTLOG(NINPUTSMAX),VELFAC(NINPUTSMAX)
         REAL*8 POWERA(NINPUTSMAX)
         REAL*8 RADIOMIN(NINPUTSMAX),RADIOMAX(NINPUTSMAX)
         REAL*8 VFAC,XDTIME,DMDTL
         INTEGER NINPUTS,INDEX,NEVAL
C ROTATING NUCLEUS AND ACTIVE AREA PARAMETERS
         REAL*8 NUC_PHI,NUC_INC,CL1,CL2
         REAL*8 AREA_LONGITUDE,AREA_LATITUDE,PERIOD,UR,UTHETA,UZ
         REAL*8 VXOP,VYOP,VZOP,EJMOD,THETA0
         REAL*8 LAT_MIN_DEG,LAT_MAX_DEG
         REAL*8 LON_MIN_DEG,LON_MAX_DEG
         REAL*8 LAT_MIN,LAT_MAX,LON_MIN,LON_MAX,SUBSOLAT
         REAL*8 RED_FACTOR
         REAL*8 CSH(7)  ! Schleicher phase function polynomial fit
         
         CHARACTER*8 RA_CHAR,DE_CHAR
         CHARACTER LCOM1*36,LCOM2*56,LCOM3*51,LCOM4*50,LCOM*222,BLANK*33
         CHARACTER LCOM5*29
C KEPLER FUNCTIONS
         EXTERNAL FKEPLER,FKEPLERHMENOS,FKEPLERH,RTBIS,EKEPL2
         EXTERNAL HKEPLER
         
C COMMON PARAMETERS AND ARRAYS         
         COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
         COMMON /HELIOMATRIX/XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
         COMMON/HELIOMATRIX_PARTICLE/XXP,XYP,XZP,YXP,YYP,YZP,ZXP,ZYP,ZZP
         COMMON/ORBITAL_ELEMENTS_PARTICLE/QD,ED,IND,WD,ND
         COMMON/PARAM_NM/DELTA,NMPAR,NMPAR1,NMPAR2,NMPAR3,NMPAR4,
     &  NMPAR5,NMPAR6,NMPAR7,NMPAR8
         COMMON/INTERP/DTIME,DMDTLOG,VELFAC,POWERA,RADIOMIN,RADIOMAX
         COMMON /ORBITAL_PERIOD/ORB_PER
C COMMON FOR OUTPUT FITS IMAGE HEADER
         COMMON/FITSIMAGE/START_JD,END_JD,GRDSIZ,NTOTMC
         COMMON/OUTPUT_IMAGE/FLUX
         
C SIXTH ORDER POLYNOMIAL FITTING TO SCHLEICHER LOG10(PHASE FUNCTION)
        DATA CSH/-7.4004978D-03,-1.6492566D-02,1.0950353D-04,
     & 8.3640600D-07,1.0157539D-09,-9.6882641D-11,4.4184372D-13/

C -----------------------------------------------------------------
        
C CPU COMPUTATION TIME:
        CALL CPU_TIME (TIME_START)

        
C INPUT FILE: PARAMETERS THAT VARY WITH HELIOCENTRIC DISTANCE
C READ INPUT FILE: dm/dt(t), vfactor(t), power(t), rmin(t), rmax(t):
        OPEN (10,FILE='dmdt_vel_power_rmin_rmax.dat',STATUS='OLD')
        READ(10,*) ! Blank line
        DO II=1,NINPUTSMAX
        READ (10,*,END=50)DTIME(II),DMDTLOG(II),VELFAC(II),POWERA(II),
     &RADIOMIN(II),RADIOMAX(II)
        ENDDO
        STOP ' Too many input values:INCREASE NINPUTSMAX (ninputs.inc) '
 50     CLOSE(10)
        NINPUTS=II-1
C INITIALIZE ARRAY OF FLUXES AND OPTICAL DEPTH:
         DO II=1,NX
            DO JJ=1,NY
               FLUX(II,JJ)=0
               OPT_DEPTH(II,JJ)=0
               FLUX_STAR(II,JJ)=0
            ENDDO
         ENDDO
         OPT_DEPTH_NUC=0
            
        
C --- SOME USEFUL CONSTANTS:
        PI=2.D0*DACOS(0.0D0)
        HALFPI=PI/2.D0
        TWOPI=2.D0*PI
        FOURPI=4.D0*PI
        TORAD=PI/180.D0
        AUKM=1.4959787D8        ! Astronomical unit (km)  
        CTEVEL=5.775483D-4      ! km/s to AU/day conversion factor
        MU=2.959122082855911D-4 ! Standard gravitational parameter GM au³/day²
        CP=1.191D-3             ! Radiation pressure coefficient
        QP=1.D0                 ! Scattering efficiency for radiation pressure
        RSUNKM=695660.D0                 ! Solar radius, km
        SOLID_ANGLE=PI*(RSUNKM/AUKM)**2  ! Solid ang of sun disk from Earth, sr
        ARCSEC_RAD=180.D0*3600.D0/PI          ! arcsec per radian 
        SOLID_ANGLE=SOLID_ANGLE*ARCSEC_RAD**2 ! arcsec**2
        CTE_MAG=2.5D0*DLOG10(SOLID_ANGLE)
        
C -- INPUT FILE, GENERAL INPUTS:
        OPEN (1,FILE='TAIL_INPUTS.dat',STATUS='old')
C FIRST LINE IS THE OBJECT IDENTIFICATION CODE (REC # in JPL-HORIZONS):
         READ (1,20) RECID    ! Check this out, as sometimes changes with time 
	 READ (1,*) PDEN            ! Particle density (kg/m3)
         READ (1,*) PV0,PHASE_COEFF ! Particle geom. albedo & linear phase coeff
         READ (1,*) IEJEC_MODE  ! Ejection mode: 1-> isotropic; 2-> sunward
                                ! 3-> Actve areas on rotating spherical nuc
         READ (1,*) NUC_INC,NUC_PHI,PERIOD
         READ (1,*) LAT_MIN_DEG,LAT_MAX_DEG
         READ (1,*) LON_MIN_DEG,LON_MAX_DEG
         READ (1,*) ISUN
             ! Switch between solar(1)/no solar(0) dependence of ejection
         NUC_INC=NUC_INC*TORAD
         NUC_PHI=NUC_PHI*TORAD
         LON_MIN=LON_MIN_DEG*TORAD
         LON_MAX=LON_MAX_DEG*TORAD
         LAT_MIN=LAT_MIN_DEG*TORAD
         LAT_MAX=LAT_MAX_DEG*TORAD
         READ (1,*) V0,GAMMA,KAPPA  ! Particle speed parameters (v0 in km/s)
         READ (1,*) EXPOCOS         ! Velocity dependence of cos(zenith_angle)**EXPOCOS
         READ (1,*) MAGSUN          ! Magnitude of the Sun in the filter
         READ (1,*) START_JD        ! Start of integration (JD)
         READ (1,*) END_JD          ! End of integration (JD)
         READ (1,*) NX,NY           ! Dust tail image dimensions (px)
         READ (1,*) SCALE           ! Plate scale (arcsec/px)
         READ (1,*) NTIMES          ! Number of time bins
         READ (1,*) NSIZES          ! Number of size bins
         READ (1,*) NEVENT          ! Number of directional event samples
         READ (1,*) IAP,RHO_AP      ! Ap. radius (1>km, 2>arcsec)for mag & Afrho computation               
         READ (1,*) RNUCLEUS        ! Nucleus radius (m)
         READ (1,*) PV0_NUC,PHASE_COEFF_NUC ! Nuc. geom. albedo & phase coeff
         READ (1,*) ICONV,SFWHM     ! Convol 1=Yes/0=No ;'seeing' FWHM
         READ (1,*) ISTAR,MAGLIM    ! Background star field (1=Yes/0=No), limiting magnitude
         READ (1,*) IPRN            ! IPRN=1, creates file of particle positions in sky (nm.dat). Caution: it might be too large.
         READ (1,*) IGRAPHO,PCP      ! IGRAPHO=1, open interactive graphics output, PCP is the percent of particles to be plotted
         CLOSE(1)
         
         PCP=1.D0-PCP/100.D0
         
         IF (NX.GT.NXMAX.OR.NY.GT.NYMAX)
     &   STOP ' REDIM PARAMETERS IN image.inc'
         IF (START_JD.GE.END_JD)
     &   STOP ' Bad choice of starting/ending times'

         NTOTMC=NSIZES*NTIMES*NEVENT ! Total Monte Carlo events
         
         IF (ISTAR.EQ.1) OPEN (88,FILE='starpos.dat',STATUS='UNKNOWN')          
         IF (IPRN.EQ.1) OPEN (16,FILE='nm.dat',STATUS='UNKNOWN')
         IF (IEJEC_MODE.EQ.3) OPEN (29,FILE='subsolar.dat',
     &                            STATUS='UNKNOWN')         
         OPEN (77,FILE='dustlossrate.dat',STATUS='UNKNOWN')         
       
C TIME ARRAY         
         DELTAT=(END_JD-START_JD)/DFLOAT(NTIMES-1)
         TIMES(1)=START_JD
         DO II=2,NTIMES
            TIMES(II)=TIMES(II-1)+DELTAT
         ENDDO
         TIMES(NTIMES)=END_JD   ! OBSERVATION TIME
                  
          
        WRITE (*,*) ' Downloading ephemeris from JPL-Horizons ... ' 
        CALL GET_COORVEL (RECID,END_JD,NTIMES,TIMES,
     &    XARR,VARR,XYZ_EARTH,QR,EC,INC,OM,WPER,TRUEANARR,
     &    PER_JD,RA,DEC,DELTA,PSANG,PLANG,PHAS_ANG)

        CALL CALDATE (START_JD,START_D,START_M,START_Y)        
        CALL CALDATE (END_JD,OBS_D,OBS_M,OBS_Y)

        WRITE (*,8) START_D,INT(START_M),INT(START_Y),
     & OBS_D,INT(OBS_M),INT(OBS_Y)
        
        WRITE (*,9) START_JD,START_JD-PER_JD,END_JD,END_JD-PER_JD
        
         IF (EC.LT.1.D0) THEN                   ! elliptic orbit
            SEMIAXIS=QR/(1.D0-EC)               ! semi-major axis (au)
            ORB_PER=TWOPI*DSQRT(SEMIAXIS**3/MU) ! orbital period (days)
            HALF_PER=0.5D0*ORB_PER             
            WRITE (*,10) ORB_PER/365.25D0
         ENDIF

        
C PHASE FUNCTION CORRECTION IS APPROXIMATE: IN THE IDEAL CASE EACH PARTICLE
C WOULD HAVE ITS OWN PHASE FUNCTION, BUT HERE WE ASSUME ALL PARTICLE BEHAVE
C EQUALLY, AND ARE ASSUMED SPHERICAL. OTHERWISE, ONE MUST RUN ELECTROMAGNETIC
C SCATTERING CODES (e.g. DDA) TO OBTAIN PHASE FUCTIONS DEPENDING ON PARTICLE
C SIZE, SHAPE, AND REFRACTIVE INDEX, OR RELY ON LAB SCATTERING DATA ...
C IN THE FUTURE, THIS WOULD ALLOW THE GENERATION OF LINEAR POLARIZATION MAPS
        
C LINEAR PHASE COEFF CORRECTION TO GEOMETRIC ALBEDO (PARTICLES)
C APPLICABLE FOR PHASE ANGLE<50 DEG           
C           PV=PV0*10**(-PHAS_ANG*PHASE_COEFF/2.5D0)
C           WRITE (*,*) PHAS_ANG,PV
C CODE IS USING SCHLEICHER PHASE FUNCTION CORRECTION FOR ALL PARTICLES:
           PV=PV0*10**(CSH(1)+CSH(2)*PHAS_ANG+CSH(3)*PHAS_ANG**2+
     &        CSH(4)*PHAS_ANG**3+CSH(5)*PHAS_ANG**4+CSH(6)*PHAS_ANG**5+
     &        CSH(7)*PHAS_ANG**6)

C LINEAR PHASE COEFF CORRECTION TO GEOMETRIC ALBEDO (NUCLEUS)
C WE USE THE LINEAR APPROXIMATION FOR ALL NUCLEUS PHASE ANGLES
           PV_NUC=PV0_NUC*10**(-PHASE_COEFF_NUC*PHAS_ANG/2.5D0)
        
C --- SET UP DUST TAIL IMAGE GRID PARAMETERS ----------        
         CALL IMAGE_GRID (NX,NY,SCALE,DELTA,GRDSIZ,ARCSEC_RAD,
     &    NMIN,NMAX,MMIN,MMAX,NGRID,MGRID,INUC,JNUC)
          AUGX=DFLOAT(NX)/(NMAX-NMIN)
          AUGY=DFLOAT(NY)/(MMAX-MMIN)
          AREAPIX=(GRDSIZ*ARCSEC_RAD/(DELTA*AUKM))**2


          IF (IGRAPHO.EQ.1) THEN
	call pgbegin (0,'/xwindow',1,1)
	call pgenv (sngl(NMIN),sngl(NMAX),sngl(MMIN),sngl(MMAX),1,0)
	call pglabel ('Projected distance on RA axis [km]',
     & 'Projected distance on DEC axis [km]',' ')
        ENDIF
          
          

          
c Convert to radians some of the outputs:       
      INC=INC*TORAD
      WPER=WPER*TORAD
      OM=OM*TORAD
      PA=PSANG*TORAD
      
C EQUATORIAL COORDINATES       
      RAH=INT(RA/15.D0)
      RAM=(RA/15.D0-RAH)*60.D0
      RAS=(RAM-INT(RAM))*60.D0
      DED=INT(DEC)
      DEM=(DEC-DED)*60.D0
      DES=(DEM-INT(DEM))*60.D0

      CALL CALDATE (END_JD,OBS_D,OBS_M,OBS_Y)


      WRITE (*,11) XYZ_EARTH(1),XYZ_EARTH(2),XYZ_EARTH(3)
      WRITE (*,12)' TARGET RA (h,m,s): ',INT(RAH),INT(RAM),RAS
      WRITE (*,12)' TARGET DEC(d,m,s): ',INT(DED),ABS(INT(DEM)),ABS(DES)

      WRITE (*,13) PHAS_ANG,PSANG,PLANG

       IF (ISTAR.EQ.0) GOTO 881 ! to skip star field section
      
C NEW SECTION: STAR FIELD WITHIN IMAGE (LIMITED TO 1 DEGREE IN RADIUS)
      WRITE (RA_CHAR,500) RA
      WRITE (DE_CHAR,500) DEC
 500  FORMAT(F7.3)
      RA_CHAR=RA_CHAR(1:7)//'d'
      DE_CHAR=DE_CHAR(1:7)//'d'
      IF (DEC.LT.0.) DE_CHAR=' '//DE_CHAR
      WRITE (*,*) ' Downloading star field coordinate/mag table ...' 
      LCOM1='wget "https://irsa.ipac.caltech.edu/'
      LCOM2='cgi-bin/Gator/nph-query?outfmt=1&objstr='//RA_CHAR//DE_CHAR
      LCOM3='&spatial=cone&radius=3600&catalog=gaia_edr3_source'
      LCOM4='&selcols=ra,dec,phot_g_mean_mag,phot_rp_mean_mag,'
      LCOM5='phot_bp_mean_mag" -O star.dat'
      LCOM=LCOM1//LCOM2//LCOM3//LCOM4//LCOM5
      write (*,*) 'lcom=',lcom
      CALL EXECUTE_COMMAND_LINE(LCOM)
      
      OPEN(1,FILE='star.dat',STATUS='OLD')
      DO II=1,36
         READ (1,*) 
      ENDDO
      K=1

      
      DO II=1,NSTARMAX
         READ (1,22,END=903,ERR=902) RA_STAR,DE_STAR,BLANK,
     &                               GMAG,RPMAG,BPMAG
c TRANSFORMATION TO STANDARD COORDINATES --> PIXEL COORDINATES ON IMAGE 
         CALL STDCOOR (RA,DEC,RA_STAR,DE_STAR,XTEMP,YTEMP)
         POSX(K)=INT(XTEMP*ARCSEC_RAD/SCALE)+INUC
         POSY(K)=INT(YTEMP*ARCSEC_RAD/SCALE)+JNUC
         POSX(K)=NX-POSX(K)
         POSY(K)=NY-POSY(K) 
C     MAGNITUDES ARE CONVERTED TO R-COUSINS BY USING GBP-GRP VS. G-R FIT
C     Busso et al. 2022:Gaia DR3 documentation Chapter 5: Photometric data 
C     Section 5.5.1 by Carrasco & Bellazzini, Fig. 5.39
         XTEMP=BPMAG-RPMAG
         YTEMP=0.02275D0+0.3691D0*XTEMP-0.1243D0*XTEMP**2-
     &        0.01396D0*XTEMP**3+0.003775D0*XTEMP**4
         RMAG=GMAG-YTEMP

         STARMAG(K)=RMAG  ! Cousins R

         STARFLUX(K)=10**(-0.4D0*(STARMAG(K)+10.925)) !

         IF (POSX(K).GE.1.AND.POSX(K).LE.NX.AND.
     &        POSY(K).GE.1.AND.POSY(K).LE.NY) THEN

            IF (STARMAG(K).LT.MAGLIM) THEN
            WRITE (88,880) POSX(K),POSY(K),STARMAG(K),STARFLUX(K)
         
         
            FLUX_STAR(POSX(K),POSY(K))=FLUX_STAR(POSX(K),POSY(K))+
     &                                 SNGL(STARFLUX(K))
        ENDIF
       ENDIF

         K=K+1
 902  ENDDO
      STOP ' INCREASE NSTARMAX,DIMENSIONS IN POSX,POSY, AND STARMAG '
 903  CONTINUE  

      
  881   CONTINUE
      
      
C EARTH COORDINATES AT OBS_DATE:
      XE=XYZ_EARTH(1)
      YE=XYZ_EARTH(2)
      ZE=XYZ_EARTH(3)
C COMET COORDINATES AT OBS_DATE:     
      XCOBS=XARR(1,NTIMES)
      YCOBS=XARR(2,NTIMES)
      ZCOBS=XARR(3,NTIMES)
C COMET TRUE ANOMALY AT OBS_DATE:     
c True anomaly at observation date, radians:
      THETACOBS=TRUEANARR(NTIMES) 

      WRITE (*,14) EC,QR,PER_JD
      WRITE (*,15) OM/TORAD,WPER/TORAD,INC/TORAD

      WRITE (*,16) GRDSIZ,SCALE,NX*SCALE/60.D0,NY*SCALE/60.D0

C COMET HELIOCENTRIC DISTANCE IN DATE-OBS in au:         
        RCOBS=DSQRT(XCOBS**2+YCOBS**2+ZCOBS**2)
        WRITE (*,17) RCOBS,DELTA

        CTE_PART=PV*SOLID_ANGLE/
     &  ((AUKM*1.D3)**2*RCOBS**2*DELTA**2*AREAPIX)

        CTE_NUC=PV_NUC*SOLID_ANGLE/
     &  ((AUKM*1.D3)**2*RCOBS**2*DELTA**2*AREAPIX)
C NUCLEUS BRIGHTNESS
        BRNUCLEUS=CTE_NUC*RNUCLEUS**2

C CONVERT TO HELIOCENTRIC ORBIT PLANE COORDINATES:
        CALL HE_TO_HPO(XE,YE,ZE,XOUTe,YOUTe,ZOUTe)
C COMETOCENTRIC COORDINATES OF THE EARTH AT OBS. DATE
        CHITAE=XOUTe*DCOS(THETACOBS)+YOUTe*DSIN(THETACOBS)-RCOBS
        ETAE=  XOUTe*DSIN(THETACOBS)-YOUTe*DCOS(THETACOBS)
        GITAE=ZOUTe
        
C CONVERSION PARAMETERS TO SKY PLANE (N-M)
        NMPAR=DSQRT(ETAE**2+GITAE**2)
        NMPAR1=NMPAR/DELTA
        NMPAR2=CHITAE*ETAE/(NMPAR*DELTA)
        NMPAR3=CHITAE*GITAE/(NMPAR*DELTA)
        NMPAR4=GITAE/NMPAR
        NMPAR5=ETAE/NMPAR

        NMPAR6=CHITAE/DELTA
        NMPAR7=ETAE/DELTA
        NMPAR8=GITAE/DELTA
        
C OBSERVATION TIME SINCE PERIHELION PASSAGE IN DAYS
        TC=(END_JD-PER_JD)
C SIZE STEP
        R1LOG=LOG10(RADMIN)
        R2LOG=LOG10(RADMAX)	
        STEPRLOG=(R2LOG-R1LOG)/DFLOAT(NSIZES)

        WRITE (*,*) ' Building up dust tail ...'
        
C ----------------------------------------
C -------------- TIME LOOP ---------------
C ----------------------------------------

        TOTAL_MASS=0.0D0


        DO 7776 ITIME=1,NTIMES-1
           
         RED_FACTOR=1.D0  ! Reduction factor, see below when IEJEC_MODE=3
           
c$$$    IF (IGRAPHO.EQ.1) THEN
c$$$    call pgask (.false.)
c$$$	call pgenv (sngl(NMIN),sngl(NMAX),sngl(MMIN),sngl(MMAX),1,0)
c$$$    ENDIF
      
        
         TIME=TIMES(ITIME)
         TAU=END_JD-TIME
         XC_EJEC=XARR(1,ITIME)
         YC_EJEC=XARR(2,ITIME)
         ZC_EJEC=XARR(3,ITIME)

         VXC_EJEC=VARR(1,ITIME)
         VYC_EJEC=VARR(2,ITIME)
         VZC_EJEC=VARR(3,ITIME)
c Heliocentric distance and true anomaly at ejection time:
         RC_EJEC=DSQRT(XC_EJEC*XC_EJEC+YC_EJEC*YC_EJEC+ZC_EJEC*ZC_EJEC)
         THETAC_EJEC=TRUEANARR(ITIME)
C     IF IEJEC_MODE=3, WRITE OUT SUBSOLAR LATITUDE AS A
c     FUNCTION OF TIME, RC, AND THETAC
C          
          IF (IEJEC_MODE.EQ.3) THEN
            SUBSOLAT=DASIN(DSIN(NUC_INC)*DSIN(NUC_PHI+THETAC_EJEC))
            WRITE (29,24) TIME-PER_JD,RC_EJEC,THETAC_EJEC,SUBSOLAT/TORAD
          ENDIF        
          
C -----------------------          
C     FOR TRAIL COMPUTATION (SHORT-PERIOD COMETS ONLY), THE COMET ORBITAL PATH
C     IS ASSUMED TO BE THE SAME DURING ALL CONSIDERED PREVIOUS COMET ORBITS,
C     REPEATING THE SAME ORBITAL AND EJECTION PATTERNS FOR ALL ORBITS
         
         XDTIME=TIME-PER_JD
C THIS SECTION BECOMES ACTIVE FOR TRAIL COMPUTATION (SEVERAL COMET ORBITS):
C --- (FOR ELLIPTICAL COMETS ONLY) ----
         IF (EC.LT.1.D0) THEN
 505     IF (XDTIME.LT.-HALF_PER)  XDTIME=XDTIME+ORB_PER
         IF (XDTIME.GT. HALF_PER)  XDTIME=XDTIME-ORB_PER
         if (XDTIME.GT.-HALF_PER.AND.XDTIME.LT.HALF_PER) GOTO 510
         GOTO 505
         ENDIF
 510     CONTINUE 
         
C INTERPOLATE LINEARLY ALL AT ONCE, ONLY ONE CALL:
           CALL INTERP5 (NINPUTS,XDTIME,DMDTL,POWER,VFAC,RADMIN,RADMAX)
           DMDT=10.**DMDTL
          
          XMASS=DMDT*DELTAT*86400.D0 ! TOTAL EJECTED DUST MASS IN TIME INTERVAL 
           
         R1LOG=LOG10(RADMIN)
         R2LOG=LOG10(RADMAX)	
         STEPRLOG=(R2LOG-R1LOG)/DFLOAT(NSIZES)
         
C ----------------------------------------         
C -------------SIZE LOOP---------------
C ----------------------------------------         

        DO 3399 IRLOG=1,NSIZES
        RLOG=R1LOG+(IRLOG-1)*STEPRLOG
	R1=10.**RLOG
	R2=10.**(RLOG+STEPRLOG)
	RAD=DSQRT(R1*R2) ! RAD AT MID-INTERVAL,GEOMETRIC MEAN
C NUMBER OF PARTICLES EJECTED IN THE INTERVAL(R1,R2) WITHIN
C THE CORRESPONDING TIME INTERVAL:
	ENEINT=3.D0*XMASS*POWERINT(R1,R2,POWER)/	      
     &  (FOURPI*PDEN*POWERINT(RADMIN,RADMAX,POWER+3.D0))
C POSITION OF PARTICLE RELATIVE TO NUCLEUS AT EMISSION TIME: 
        XIT=0
        YIT=0
        ZIT=0
C INITIAL DUST SPEED (A FUNCTION OF THE PARTICLE RADIUS
C AND THE HELIOCENTRIC DISTANCE IN THIS EXAMPLE):
        BETA=1.191D-3/(2.D0*PDEN*RAD)
        VEJ=V0*VFAC*BETA**GAMMA*RC_EJEC**KAPPA


C IF IEJEC_MODE=3, PERFORM A PRE-SAMPLING FOR ONE SINGLE SIZE
C TO ESTIMATE THE PORTION OF THE AREA IN DARKNESS, E.G., COSZ<0
C THIS TEST IS NOT DONE IF ISUN.NE.1, AS IN THAT CASE EJECTION DOES
C NOT DEPEND ON SOLAR ILLUMINATION CONDITIONS
        

C ----- PRE-SAMPLING STARTS -----------           
        IF (IEJEC_MODE.EQ.3.AND.IRLOG.EQ.1.AND.ISUN.EQ.1) THEN

           NSHADOW=0
           DO KK=1,NEVENT
           
        CALL RANDOM_NUMBER (RANDOM1)
C LONGITUDE SAMPLING:
        AREA_LONGITUDE=LON_MIN+(LON_MAX-LON_MIN)*RANDOM1        
C LATITUDE SAMPLING:
        CL1=DCOS(HALFPI+LAT_MAX)
        CL2=DCOS(HALFPI+LAT_MIN)
        CALL RANDOM_NUMBER (RANDOM2)
        AREA_LATITUDE=DACOS((CL2-CL1)*DBLE(RANDOM2)+CL1)-HALFPI
c        write (*,*) AREA_LATITUDE/TORAD,AREA_LONGITUDE/TORAD
        CALL ANISOT_DIR2(NUC_PHI,NUC_INC,THETAC_EJEC,
     &       TIME,PER_JD,PERIOD,AREA_LATITUDE,AREA_LONGITUDE,
     &       THETA0,UR,UTHETA,UZ)

        COSZ=-UR

        
        IF (COSZ.LT.0.0D0) NSHADOW=NSHADOW+1

        

      ENDDO

      
      IF (NSHADOW.EQ.NEVENT) THEN ! All the area is in shadow
      WRITE (*,990) TIME-PER_JD
         DMDT=0.D0  ! No mass is ejected in this time step
         GOTO 7775  ! Go to next time step
      ELSE
C CONTINUE WITH A REDUCED MASS AND PARTICLE NUMBER
         RED_FACTOR=1.D0-DFLOAT(NSHADOW)/DFLOAT(NEVENT)
         DMDT=DMDT*RED_FACTOR
      ENDIF
      ENDIF ! Pre-sampling ends
      
C ---- PRE-SAMPLING ENDS ------

C REDUCE PARTICLE NUMBER (ONLY IF IEJEC_MODE=3 AND NSHADOW.NE.0, OTHERWISE
C THE REDUCTION FACTOR EQUALS UNITY)

          ENEINT=ENEINT*RED_FACTOR
      
C ----------------------------------------        
C --------DIRECTIONAL SAMPLING LOOP--------------
C ----------------------------------------        
        
        DO 1122 K=1,NEVENT

           
           
        IF (IEJEC_MODE.EQ.1) THEN
C ISOTROPIC:
        CALL RANDOM_NUMBER(RANDOM1)
        CALL RANDOM_NUMBER(RANDOM2)
        LAT=DACOS(2.D0*DBLE(RANDOM1)-1.D0)-HALFPI
        LON=TWOPI*DBLE(RANDOM2)
C DUST VELOCITY COMPONENTS:
        VXIT=DCOS(LAT)*DCOS(LON)        
        VYIT=DCOS(LAT)*DSIN(LON)        
        VZIT=DSIN(LAT)

        VXIT=VXIT*CTEVEL*VEJ
        VYIT=VYIT*CTEVEL*VEJ
        VZIT=VZIT*CTEVEL*VEJ
        GOTO 7799
        
        ELSE IF (IEJEC_MODE.EQ.2) THEN
C SUNWARD HEMISPHERICAL EMISSION (V $\propto$ COSZ**EXPOCOS):
 72     CALL RANDOM_NUMBER(RANDOM1)
        CALL RANDOM_NUMBER(RANDOM2)
        LAT=DACOS(2.D0*DBLE(RANDOM1)-1.D0)-HALFPI
        LON=TWOPI*DBLE(RANDOM2)
C DUST VELOCITY COMPONENTS:
        VXIT=DCOS(LAT)*DCOS(LON)        
        VYIT=DCOS(LAT)*DSIN(LON)        
        VZIT=DSIN(LAT)
        EJMOD=SQRT(XC_EJEC**2+YC_EJEC**2+ZC_EJEC**2)
        COSZ=-(VXIT*XC_EJEC+VYIT*YC_EJEC+VZIT*ZC_EJEC)/EJMOD

        IF (COSZ.LT.0.D0) GOTO 72 ! Only emission from illuminated areas
                                  ! so it goes back to 72 until cosz>0 
        
        VXIT=VXIT*CTEVEL*VEJ*(COSZ)**EXPOCOS
        VYIT=VYIT*CTEVEL*VEJ*(COSZ)**EXPOCOS
        VZIT=VZIT*CTEVEL*VEJ*(COSZ)**EXPOCOS
        GOTO 7799

        ELSE  ! IEJEC_MODE=3

 690    CALL RANDOM_NUMBER (RANDOM1)
C LONGITUDE SAMPLING:
        AREA_LONGITUDE=LON_MIN+(LON_MAX-LON_MIN)*RANDOM1        
        CL1=DCOS(HALFPI+LAT_MAX)
        CL2=DCOS(HALFPI+LAT_MIN)
C LATITUDE SAMPLING:
        CALL RANDOM_NUMBER (RANDOM2)
        AREA_LATITUDE=DACOS((CL2-CL1)*DBLE(RANDOM2)+CL1)-HALFPI

        CALL ANISOT_DIR2(NUC_PHI,NUC_INC,THETAC_EJEC,
     &       TIME,PER_JD,PERIOD,AREA_LATITUDE,AREA_LONGITUDE,
     &       THETA0,UR,UTHETA,UZ)
        COSZ=-UR                
        IF (COSZ.LT.0.0D0.AND.ISUN.EQ.1) goto 690 ! Some fraction of the
        ! area must be illuminated according to the above pre-sampling test
        ! If ISUN.NE.1 program simply continues (irrelevant illumination)
        VXOP=UR*DCOS(THETAC_EJEC)-UTHETA*DSIN(THETAC_EJEC)
        VYOP=UR*DSIN(THETAC_EJEC)+UTHETA*DCOS(THETAC_EJEC)
        VZOP=UZ
        
C TRANSFORM PLANE-OF-ORBIT TO HELIOCENTRIC PARTICLE VELOCITY 
        CALL HPO_TO_HE (VXOP,VYOP,VZOP,VXIT,VYIT,VZIT)

        VXIT=VXIT*CTEVEL*VEJ*(COSZ)**EXPOCOS
        VYIT=VYIT*CTEVEL*VEJ*(COSZ)**EXPOCOS
        VZIT=VZIT*CTEVEL*VEJ*(COSZ)**EXPOCOS


        GOTO 7799

        ENDIF  ! Ejection mode   
        
 7799   CALL HELIORBIT (PER_JD,TC,TAU,QP,PDEN,RAD,
     &  XC_EJEC,YC_EJEC,ZC_EJEC,VXC_EJEC,VYC_EJEC,VZC_EJEC,
     &  RCOBS,THETACOBS,
     &  XIT,YIT,ZIT,VXIT,VYIT,VZIT,NPAR,MPAR,LPAR)

        
c These NPAR,MPAR coordinates are rotated by angle PA counterclockwise to
c be transformed to coordinates in the RA,DEC axes:
          
          NROT=NPAR*DCOS(PA)-MPAR*DSIN(PA)
	  MROT=NPAR*DSIN(PA)+MPAR*DCOS(PA)
C Keep old notation for convenience:
	  NPAR=NROT
	  MPAR=MROT
          
          IF (IGRAPHO.EQ.1) THEN
             CALL RANDOM_NUMBER (RANDOM1)
             IF (RANDOM1.GT.PCP) THEN
             CALL PGPOINT (1,SNGL(NROT),SNGL(MROT),-1)
             ENDIF
          ENDIF
          
          IF (IPRN.EQ.1) THEN
          WRITE (16,19) NPAR,MPAR,LPAR,RAD,ITIME,TIMES(ITIME)
          CALL FLUSH(16)
          ENDIF
       
C  STORE IN DUST TAIL ARRAY BRIGHTNESS
          
          IF ( NPAR.GE.NMIN.AND.NPAR.LE.NMAX .AND.
     &         MPAR.GE.MMIN.AND.MPAR.LE.MMAX ) THEN
          II=1+IDINT((NPAR-NMIN)*AUGX) 
          JJ=1+IDINT((MPAR-MMIN)*AUGY)
C PARTICLE INTENSITY RATIO i/i_0:
          FLUX_OUT=CTE_PART*RAD**2
          FLUX(II,JJ)=FLUX(II,JJ)+
     &    SNGL(FLUX_OUT*ENEINT/DFLOAT(NEVENT))
C OPTICAL DEPTH SIMPLY USES FRAUNHOFFER APPROX, QEXT=2 
          OPT_DEPTH(II,JJ)=OPT_DEPTH(II,JJ)+
     &         SNGL(2.D0*ENEINT*PI*(RAD/(GRDSIZ*1.D3))**2)

C THE NUCLEUS BRIGHTNESS ATTENUATION IS ONLY AFFECTED BY PARTICLES
C HAVING LPAR > 0
       IF (LPAR.GT.0.D0.AND.II.EQ.INUC.AND.JJ.EQ.JNUC) THEN
          OPT_DEPTH_NUC=OPT_DEPTH_NUC+
     &     2.D0*ENEINT*PI*(RAD/(GRDSIZ*1.D3))**2
       ENDIF
          
       ENDIF

 1122  ENDDO   ! Directional sampling
 3399  ENDDO   ! Size sampling
 7775  TOTAL_MASS=TOTAL_MASS+DMDT*DELTAT*86400.D0 ! TOTAL DUST MASS EJECTED
       WRITE (77,995) TIME-PER_JD,RC_EJEC,DMDT
 7776  ENDDO    ! Time sampling

       WRITE (*,18) TOTAL_MASS

       
C ADD NUCLEUS CONTRIBUTION TO IMAGE, ATTENUATED BY OPTICAL THICKNESS:
          
           FLUX(INUC,JNUC)=FLUX(INUC,JNUC)+SNGL(BRNUCLEUS*
     &                     EXP(-OPT_DEPTH_NUC))
C ATTENUATES STAR FLUX BY OPTICAL DEPTH OF DUST TAIL
C -- ONLY WORKS FOR HIGHLY SAMPLED IMAGES, i.e. LARGE NTIMES x NSIZES x NEVENT)
          DO II=1,NX
              DO JJ=1,NY
               FLUX_STAR(II,JJ)=FLUX_STAR(II,JJ)*EXP(-OPT_DEPTH(II,JJ))
            ENDDO
         ENDDO
C TOTAL FLUX ON IMAGE:
C IF ANY STAR IS LOCATED AT NUCLEUS PIXEL, and THE DISTANCE FROM CENTER
C IS SMALLER THAN THE PHYSICAL NUCLEUS RADIUS THE STAR LIGHT WOULD BE
C BLOCKED, BUT THIS VERY UNLIKELY CASE IS NOT ACCOUNTED FOR  
           DO II=1,NX
              DO JJ=1,NY
               FLUX(II,JJ)=FLUX(II,JJ)+FLUX_STAR(II,JJ)
            ENDDO
         ENDDO
                  
C --- WRITE OUT RESULTS AS FITS FILES           
C -----------------------------------
         OUTIMAGE='tail_sdu.fits'
         IF (ICONV.EQ.1) THEN              
         CALL CONVOLY(SNGL(SFWHM),FLUX,CONVOLU,NX,NY)
         DO II=1,NX
         DO JJ=1,NY
         FLUX(II,JJ)=CONVOLU(II,JJ)
         ENDDO
         ENDDO
         ENDIF
            
C OUTPUT SYNTHETIC TAIL IMAGE IN SOLAR DISK INTENSITY UNITS (FITS FORMAT):      
         CALL WRITEIMAGE (OUTIMAGE,RECID,FLUX,NX,NY)
C COMPUTES AFRHO AND MAGNITUDE FOR SELECTED APERTURE AROUND THE PHOTOCENTER
         IF (IAP.EQ.1) RHO_AP=(RHO_AP/SCALE)*GRDSIZ  ! Aperture in km
         CALL AFRHOMAG (NX,NY,INUC,JNUC,MAGSUN,SCALE,DELTA,RCOBS,
     &        GRDSIZ,PHAS_ANG,RHO_AP,AFRHO,AFRHO_0,MAG)
         WRITE (*,23) RHO_AP,AFRHO,MAG
         WRITE (99,24) END_JD-PER_JD,AFRHO,AFRHO_0,MAG

C OUTPUT SYNTHETIC TAIL IMAGE IN MAG/ARCSEC**2 (FITS FORMAT):         
         OUTIMAGE='tail_mag.fits'         

         DO II=1,NX
         DO JJ=1,NY
         IF (FLUX(II,JJ).GT.0.0) THEN      
         FLUX(II,JJ)=SNGL(CTE_MAG)+MAGSUN-2.5*LOG10(FLUX(II,JJ))
         ELSE
         FLUX(II,JJ)=100.  ! mag/arcsec**2 (i.e., ~zero flux)
         ENDIF   
         ENDDO
         ENDDO

         CALL WRITEIMAGE (OUTIMAGE,RECID,FLUX,NX,NY)      

C OPTICAL DEPTH IMAGE (FITS FORMAT):      
       OUTIMAGE='OPT_DEPTH.fits'               
       CALL WRITEIMAGE (OUTIMAGE,RECID,OPT_DEPTH,NX,NY)      
        
       IF (IGRAPHO.EQ.1) CALL PGADVANCE
       
C ---END OF PROCEDURE ---------------        
         CALL CPU_TIME (TIME_END)
         WRITE (*,21) (TIME_END-TIME_START)/60.
C -------------------------------------------------         
 8       FORMAT (1X,
     &   ' Start activity date (Day,Month,Year): ',F7.3,1X,I2,1X,I5,
     & /,'  Observing date      (Day,Month,Year): ',F7.3,1X,I2,1X,I5) 
 9       FORMAT (1X,' JD at activation :',F14.5,'  ,i.e.,',F10.2,
     & ' days to perihelion',
     &        /'  JD at observation:',F14.5,'  ,i.e.,',F10.2,
     & ' days to perihelion')
 10      FORMAT (1X,' Orbital period: ',F10.3,' years ')
 11      FORMAT (1X,
     &   ' Helioc. coor. of the Earth at obs. date:',3(F12.7))
 12      FORMAT (1X,A20,1X,I3,I3,1X,F5.2)
 13      FORMAT (1X,' Phase angle = ',F7.3,' Pos. angle = ',F7.3,
     &    ' Plane angle= ',F7.3, ' (deg)') 
 14      FORMAT (1X,' Comet elements: EC= ',F9.6,
     &    ' QR= ',F9.6,' au    TP= ',F14.5)
 15      FORMAT (1X,' Comet elements: NODE= ',F7.3,
     &    ' W= ',F7.3,' INC= ',F7.3,' (deg)')

 16      FORMAT(1X,' Image scale=',F10.3,' km/px <==> ',
     &  F8.3,' arcsec/px ',' Field of view= ',F8.3,' x',F8.3,' arcmin')
 17      FORMAT (1X,' Target heliocentric distance=',F11.3,' au ',
     &  ' Distance to observer=',F11.3,' au')
 18      FORMAT (1X,' Total dust mass ejected= ',1PE10.3,' kg')
 19      FORMAT (1X,4(E15.5,1X),I6,1X,F18.4)
 20      FORMAT(A14)
 21      FORMAT (1X,' Elapsed CPU time: ',F12.5,' minutes')
 22      FORMAT(6X,F20.16,6X,F20.16,A37,F13.10,7X,F13.10,7X,F13.10)         
 23      FORMAT(1X,' Aperture(km)=', F10.2,
     & ' Afrho(m)=',F8.3,2X,' Mag=',F8.3)
 24      FORMAT (1X,E20.12,1X,6(E15.8,1X))
 880     FORMAT (1X,2(I5,1X),2(E12.4,1X))
 990     FORMAT (1X,' The entire prescribed area is in darkness at:',
     &         F10.3,' days to perihelion')
 995     FORMAT (1X,4(E15.5,1X))
C -----------------------------------------------------         
         END
                          

      SUBROUTINE HE_TO_HPO(XIN,YIN,ZIN,XOUT,YOUT,ZOUT)

C CONVERT HELIOCENTRIC ECLIPTIC COORDINATES (XIN,YIN,ZIN) TO 
C HELIOCENTRIC PLANE-OF-ORBIT (XOUT,YOU,ZOUT)

      IMPLICIT NONE
      REAL*8 XIN,YIN,ZIN,XOUT,YOUT,ZOUT
      REAL*8 XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
      COMMON /HELIOMATRIX/XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ

      XOUT=XX*XIN + XY*YIN + XZ*ZIN
      YOUT=YX*XIN + YY*YIN + YZ*ZIN
      ZOUT=ZX*XIN + ZY*YIN + ZZ*ZIN

      RETURN
      END

       SUBROUTINE SETHELIOMATRIX(INC,OM,WPER,XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ)
C SET ROTATION MATRIX ELEMENTS FROM HELIOCENTRIC ORBITAL PLANE
C TO HELIOCENTRIC ECLIPTIC AXES
       IMPLICIT NONE
       REAL*8 INC,OM,WPER,XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ
       REAL*8 CW,CO,CI
       REAL*8 SW,SO,SI

       CW=DCOS(WPER)
       CO=DCOS(OM)
       CI=DCOS(INC)        
       SW=DSIN(WPER)
       SO=DSIN(OM)
       SI=DSIN(INC)
       
       XX=CW*CO-CI*SO*SW
       YX=-SW*CO-CI*SO*CW
       ZX=SI*SO

       XY=CW*SO+CI*CO*SW
       YY=-SW*SO+CI*CO*CW
       ZY=-SI*CO

       XZ=SW*SI
       YZ=CW*SI
       ZZ=CI
       
       RETURN
       END

        SUBROUTINE WRITEIMAGE (OUTIMAGE,RECID,IMAGEFIT,NX,NY)
        IMPLICIT NONE
C  Create a FITS primary array containing a 2-D image
	include 'image.inc'
        INTEGER status,unit,blocksize,bitpix,naxis,naxes(2)
        INTEGER i,j,group,fpixel,nelements,nx,ny,ntotmc
        CHARACTER RECID*14,OUTIMAGE*80
        REAL*8 START_JD,END_JD,GRDSIZ
	REAL*4 IMAGEFIT(NXMAX,NYMAX),ARRAY(NX,NY)
        LOGICAL SIMPLE,EXTEND

        COMMON/FITSIMAGE/START_JD,END_JD,GRDSIZ,NTOTMC        
        
        status=0
C  Delete the file if it already exists, so we can then recreate it.
C  The deletefile subroutine is listed at the end of this file.
      call deletefile(OUTIMAGE,status)
C     Get an unused Logical Unit Number to use to create the FITS file
      call ftgiou(unit,status)
C     create the new empty FITS file
      blocksize=1
      call ftinit(unit,OUTIMAGE,blocksize,status)
C     initialize parameters about the FITS image (300 x 200 16-bit integers)
      simple=.true.
      bitpix=-32  
      naxis=2
      naxes(1)=nx
      naxes(2)=ny
      extend=.true.

      do j=1,naxes(2)
          do i=1,naxes(1)
              array(i,j)=imagefit(i,j)
          enddo
      enddo
C  Write the required header keywords
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
      do j=1,naxes(2)
          do i=1,naxes(1)
              array(i,j)=imagefit(i,j)
          enddo
      enddo


C  Write the array to the FITS file
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)
      call ftppre(unit,group,fpixel,nelements,array,status)
      
C  Write some more useful information to the header:
      
        call ftpkys(unit,'OBJECT',RECID,'Object id ',status)
 	call ftpkyf(unit,'START_TIME (JD)',SNGL(START_JD),5,
     &       'Activation date',status)
	call ftpkyf(unit,'OBS_TIME (JD)',SNGL(END_JD),5,
     &       'Observation date',status)
	call ftpkyf(unit,'PIXSIZ',SNGL(GRDSIZ),3,
     &       'Pixel size [km]',status)
	call ftpkyj(unit,'MC-Event',NTOTMC,
     &'Total MonteCarlo events',status)

C Close the file and free the unit number
      call ftclos(unit, status)
      call ftfiou(unit, status)
      return
      end


C *************************************************************************
      subroutine deletefile(filename,status)

C  A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*(*) filename

C  Simply return if status is greater than zero
      if (status .gt. 0)return

C  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)

C  Try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
C         file was opened;  so now delete it 
          call ftdelt(unit,status)
      else if (status .eq. 103)then
C         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
      else
C         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

C  Free the unit number for later reuse
      call ftfiou(unit, status)
      end

      SUBROUTINE IMAGE_GRID (NX,NY,SCALE,DELTA,GRDSIZ,ARCSEC_RAD,
     & NMIN,NMAX,MMIN,MMAX,NGRID,MGRID,INUC,JNUC)
      IMPLICIT NONE
      include 'image.inc'
      REAL*8 PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,ARCSEC_RAD,MU,CTE_MAG
      REAL*8 SCALE,DELTA,GRDSIZ,GR2,NGRID(NXMAX),MGRID(NYMAX)
      REAL*8 NMIN,NMAX,MMIN,MMAX
      INTEGER I,NX,NY,INUC,JNUC
      COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG

      GRDSIZ=SCALE*DELTA*AUKM/ARCSEC_RAD ! km/px
      GR2=GRDSIZ/2.D0

      INUC=NX/2
      JNUC=NY/2
      
      NMIN=(INUC-1)*GRDSIZ+GR2
      NMAX=(NX-INUC)*GRDSIZ-GR2
      NMIN=-NMIN

      MMIN=(JNUC-1)*GRDSIZ+GR2
      MMAX=(NY-JNUC)*GRDSIZ-GR2
      MMIN=-MMIN
      
        NGRID(1)=NMIN+GR2
        MGRID(1)=MMIN+GR2
        DO I=2,NX
           NGRID(I)=NGRID(I-1)+GRDSIZ
        ENDDO
        DO I=2,NY
         MGRID(I)=MGRID(I-1)+GRDSIZ
        ENDDO
        RETURN
        END
      
C --------------------------------------------------------


      REAL*8 FUNCTION POWERINT (XA,XB,ALPHA)
      IMPLICIT NONE
      REAL*8 XA,XB,ALPHA,UALP
C INTEGRAL OF A SINGLE POWER-LAW DISTRIBUTION OF INDEX ALPHA
C BETWEEN RADII XA AND XB
	IF (ALPHA.NE.-1.D0) THEN
	UALP=1.D0+ALPHA
	POWERINT=(XB**UALP-XA**UALP)/UALP
	ELSE
	POWERINT=DLOG(XB)-DLOG(XA)
	ENDIF
	RETURN
	END

 

      SUBROUTINE GET_COORVEL (RECID,END_JD,NTIMES,TIMES,
     &     XARR,VARR,XYZ_EARTH,QR,EC,INC,OM,WPER,TRUEANARR,PER_JD,
     &     RA,DEC,DELTA,PSANG,PLANG,PHAS_ANG)
      
C Query the JPL-Ephemeris generator to obtain observation point and target
C coordinates

      IMPLICIT NONE 
      INCLUDE 'timemax.inc'
      REAL*8 PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
      REAL*8 XARR(3,NTIMEMAX),VARR(3,NTIMEMAX),XYZ_EARTH(3)
      REAL*8 TIMES(NTIMEMAX),EC,INC,WPER,OM,PER_JD,QR
      REAL*8 TRUEANARR(NTIMEMAX),X,Y,Z,VX,VY,VZ
      REAL*8 XJD,END_JD
      REAL*8 RA,DEC,DELTA,DELDOT,PSANG,PSAMV,PLANG,PHAS_ANG
      REAL*8 GM,OM_RAD,INC_RAD,WPER_RAD
      REAL*8 TRUEANOMALY
      REAL*8 XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ      
      REAL*8 ORB_PER,SEMIAXIS,PER_D,PER_M,PER_Y
      CHARACTER*200 STRING
      CHARACTER*80 STARC
      CHARACTER*14 RECID
      CHARACTER*1 CONT
      INTEGER I,NTIMES

      COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
      COMMON /HELIOMATRIX/XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ

      CALL SET_OBSERVATION_POINT (END_JD)      
C QUERY THE JPL EPHEMERIS TO GET THE OBS. POINT POSITION AT OBS.TIME(END_JD)

      CALL EXECUTE_COMMAND_LINE('curl -s -F format=text -F 
     & input=@obs_point.txt 
     & https://ssd.jpl.nasa.gov/api/horizons_file.api > obs.txt')

      OPEN (1,FILE='obs.txt',STATUS='OLD')
      DO I=1,100000
         READ (1,10) STARC

         IF (STARC(1:5).EQ.'$$SOE') GOTO 800
      ENDDO
 800  READ (1,2) XJD
      READ (1,200) STRING
      CALL EXTRACTN(STRING,X,Y,Z)
      XYZ_EARTH(1)=X
      XYZ_EARTH(2)=Y         
      XYZ_EARTH(3)=Z
      CLOSE(1)
C GET  RA,DEC,PSANG,PHAS_ANG, AND PLANG OF COMET AT OBSERVATION DATE      
      
      CALL GET_COMET_RADEC(RECID,END_JD)
      CALL EXECUTE_COMMAND_LINE('curl -s -F format=text -F 
     & input=@comet_ra_dec_elements.txt 
     & https://ssd.jpl.nasa.gov/api/horizons_file.api > target.txt')
      OPEN (1,FILE='target.txt',STATUS='OLD')

      DO I=1,10000
         READ (1,10) STARC 
         IF (STARC(1:16).EQ.'Target body name') WRITE (*,*) ' '//STARC
         IF (STARC(1:7).EQ.'  EPOCH') THEN
            READ (1,200) STRING
            CALL EXTRACTN(STRING,EC,QR,PER_JD)
c  Correct PER_JD if it does not corresponds to the orbit
c  where observation occurs

            IF (EC.LT.1.D0) THEN
            SEMIAXIS=QR/(1.D0-EC) 
            ORB_PER=TWOPI*DSQRT(SEMIAXIS**3/MU) 
            IF (DABS(END_JD-PER_JD).GT.ORB_PER/2.D0) THEN
            WRITE (*,*)
     &  ' The orbit selected does not corresponds to the current one'
            WRITE (*,100) PER_JD
            CALL CALDATE (PER_JD,PER_D,PER_M,PER_Y)               
            WRITE (*,102) PER_D,INT(PER_M),INT(PER_Y) 
            PER_JD=PER_JD+ORB_PER
            WRITE (*,101) PER_JD
            CALL CALDATE (PER_JD,PER_D,PER_M,PER_Y)               
            WRITE (*,103) PER_D,INT(PER_M),INT(PER_Y) 
C   per_jd=2459361.6141d0 ! May 21, 2021, 7P/Pons-Winnecke nearest perihelion
            WRITE (*,*) ' Enter to continue .. '
            READ (*,1) CONT
 1          FORMAT(A1)
C            PAUSE
            ENDIF
            ENDIF

            READ (1,200) STRING
            CALL EXTRACTN(STRING,OM,WPER,INC)

         ENDIF   
         IF (STARC(1:5).EQ.'$$SOE') GOTO 750
      ENDDO
      STOP ' Error reading target.txt '


      
 750       READ (1,'(A)') STRING
           READ (STRING(22:200),*) RA,DEC,DELTA,DELDOT,
     &     PSANG,PSAMV,PLANG,PHAS_ANG
           CLOSE(1)


           
           
C DEFINES TRANSFORMATION MATRIX
      GM=MU !2.959122082855911D-4 !GM in au**3/day**2
      OM_RAD=OM*TORAD
      WPER_RAD=WPER*TORAD
      INC_RAD=INC*TORAD           
      CALL SETHELIOMATRIX (INC_RAD,OM_RAD,WPER_RAD,
     &                   XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ)        


C COMPUTE X,Y,Z,VX,VY,VZ, AND TRUE ANOMALY AS A FUNCTION OF TIME (time_ejec)   

      DO I=1,NTIMES
      CALL ELEMENTS_TO_XV (GM,TIMES(I),QR,EC,PER_JD,
     & X,Y,Z,VX,VY,VZ,TRUEANOMALY)

      

      
      TRUEANARR(I)=TRUEANOMALY


       XARR(1,I)=X
       XARR(2,I)=Y
       XARR(3,I)=Z
       VARR(1,I)=VX
       VARR(2,I)=VY
       VARR(3,I)=VZ
       ENDDO
       
 2    FORMAT(F17.9)
 100  FORMAT (1X,' Perihelion date from requested orbit, JD=',f20.10) 
 101  FORMAT (1X,' Updated perihelion date ( Check!), JD=',f20.10) 
 102  FORMAT (1X,' Perihelion date (Day,Month,Year): ',F7.3,1X,I2,1X,I5) 
 103  FORMAT (1X,' Updated perihelion date (Day,Month,Year): ',
     & F7.3,1X,I2,1X,I5) 
            
 10   FORMAT(A80)
 200  FORMAT(A200)
      RETURN
      END

C ---------------      
      SUBROUTINE SET_OBSERVATION_POINT (END_JD)
      
      IMPLICIT NONE 
      REAL*8 END_JD
      CHARACTER*40 CLINES(9),CE
C---      CHARACTER*28 CTIME1
      CHARACTER*47 CTIME1
      
      INTEGER I

C SET EARTH FILE TO GET HELIOCENTRIC POSITIONS AND VELOCITIES
C AT THE REQUESTED TIME SPAN

      OPEN (1,FILE='obs_point.txt',STATUS='UNKNOWN')
      
      CLINES(1)="!$$SOF"
      CLINES(2)="COMMAND='399'"  ! Set to Earth 
      CLINES(3)="OBJ_DATA='NO'"   
      CLINES(4)="MAKE_EPHEM='YES'"
      CLINES(5)="REF_PLANE='ECLIPTIC'"
      CLINES(6)="TABLE_TYPE='VECTORS'"
      CLINES(7)="CENTER='500@10'"  ! Center=Sun
      CLINES(8)="OUT_UNITS='AU-D'"
      CLINES(9)="VEC_TABLE='2'"
      WRITE (CE,100)END_JD

      DO I=1,9
         WRITE (1,14) CLINES(I)
      ENDDO
      
      CTIME1='TLIST='//"'"//CE
      WRITE (1,15) CTIME1//"'"
      CLOSE(1)

 14   FORMAT(g0,A40)
 15   FORMAT(g0,A28)
      
 100  FORMAT(F14.6)

      RETURN
      END
      

      SUBROUTINE GET_COMET_RADEC (RECID,END_JD)
      IMPLICIT NONE 
      REAL*8 END_JD
      CHARACTER*40 CLINES(10),CE
C ---      CHARACTER*28 CTIME1
      CHARACTER*47 CTIME1
      
      CHARACTER*14 RECID
      INTEGER I

C SET EARTH FILE TO GET HELIOCENTRIC POSITIONS AND VELOCITIES
C AT THE REQUESTED TIME SPAN

      OPEN (1,FILE='comet_ra_dec_elements.txt',STATUS='UNKNOWN')
      CLINES(1)="!$$SOF"
      CLINES(2)="COMMAND="//"'"//RECID//";'"
      CLINES(3)="OBJ_DATA='NO'"   
      CLINES(4)="EPHEM_TYPE='OBS'"
      CLINES(5)="REF_PLANE='ECLIPTIC'"
      CLINES(6)="ANG_FORMAT='DEG'"
      CLINES(7)="CENTER='399'"  ! geocentric - Change for another viewpoint
      CLINES(8)="OUT_UNITS='AU-D'"
      CLINES(9)="CAL_FORMAT='JD'"
      CLINES(10)="QUANTITIES='1,20,27,28,43'" ! RA,DEC (Astrometric)
                                              ! PsAng,PlAng,Phase angle
      WRITE (CE,100)END_JD
      CTIME1='TLIST='//"'"//CE
      DO I=1,10
         WRITE (1,14) CLINES(I)
      ENDDO
      WRITE (1,15) CTIME1//"'"
      CLOSE(1)

 14   FORMAT(g0,A40)
 15   FORMAT(g0,A28)
 100  FORMAT(F14.6)
      RETURN
      END

      
      SUBROUTINE SET_COMET_POSVEL (RECID,START_JD,END_JD,NTIMES)
      IMPLICIT NONE
      REAL*8 START_JD,END_JD

      CHARACTER*40 CLINES(9),CS,CE
      CHARACTER*28 CTIME1,CTIME2,CTIME3
      CHARACTER*6 SPAN
      CHARACTER*14 RECID
      INTEGER I,NTIMES

      
C SET COMET FILE TO GET HELIOCENTRIC POSITIONS AND VELOCITIES
C AT THE REQUESTED TIME SPAN

      OPEN (1,FILE='comet.txt',STATUS='UNKNOWN')
      
      CLINES(1)="!$$SOF"
      CLINES(2)="COMMAND="//"'"//RECID//";'"
      CLINES(3)="OBJ_DATA='NO'"   
      CLINES(4)="MAKE_EPHEM='YES'"
      CLINES(5)="REF_PLANE='ECLIPTIC'"
      CLINES(6)="TABLE_TYPE='VECTORS'"
      CLINES(7)="CENTER='500@10'" ! Center=Sun
      CLINES(8)="OUT_UNITS='AU-D'"
      CLINES(9)="VEC_TABLE='2'"

      WRITE (CS,100)START_JD
      WRITE (CE,100)END_JD
      WRITE (SPAN,101)NTIMES

      
      CTIME1='START_TIME='//"'JD "//CS
      CTIME2='STOP_TIME='//"'JD "//CE
      CTIME3='STEP_SIZE='//"'"//SPAN//"'"
      CTIME1=CTIME1//"'"

      DO I=1,9
         WRITE (1,14) CLINES(I)
      ENDDO
      WRITE (1,15) CTIME1//"'"
      WRITE (1,15) CTIME2//"'"
      WRITE (1,15) CTIME3
      CLOSE(1)

      

 14   FORMAT(g0,A40)
 15   FORMAT(g0,A28)
      
 100  FORMAT(F14.6)
 101  FORMAT(I6)
      RETURN
      END
      
      

      SUBROUTINE EXTRACTN(STRING,X,Y,Z)
      
C EXTRACTS 3 REAL*8 NUMBERS FROM A STRING CHAIN
C USEFUL FOR READ IN FROM JPL EPHEMERIS GENERATOR
C POSITION/VELOCITY OUTPUT FILE (VECTOR FORMAT) OR
C ORBITAL ELEMENT

      IMPLICIT NONE
      REAL*8 X,Y,Z
      CHARACTER*200 STRING,STRINGN,CNUM1,CNUM2,CNUM3
      INTEGER I,J,K,L,N,N1,N2,N3

      N=1
      DO I=1,200

         IF (  STRING(I:I).EQ.'0'.OR.STRING(I:I).EQ.'1'.OR.
     &         STRING(I:I).EQ.'2'.OR.STRING(I:I).EQ.'3'.OR.
     &         STRING(I:I).EQ.'4'.OR.STRING(I:I).EQ.'5'.OR.
     &         STRING(I:I).EQ.'6'.OR.STRING(I:I).EQ.'7'.OR.
     &         STRING(I:I).EQ.'8'.OR.STRING(I:I).EQ.'9'.OR.
     &         STRING(I:I).EQ.'.'.OR.(STRING(I:I).EQ.'E'.
     &         AND.STRING(I+1:I+1).NE.'C').OR.
     &         STRING(I:I).EQ.'+'.OR.STRING(I:I).EQ.'-' )THEN    
            STRINGN(N:N)=STRING(I:I)
            N=N+1
         ELSE
             STRINGN(N:N)=' '
             N=N+1
        ENDIF
      ENDDO


      N=1
      DO I=1,200
         IF (STRINGN(I:I).NE.' ') THEN
            CNUM1(N:N)=STRINGN(I:I)
            N=N+1
         ENDIF
            IF (N.GT.1.AND.STRINGN(I:I).EQ.' ') GOTO 1000
       ENDDO
         
 1000  N1=N-1

       N=1
         DO J=I,200

         IF (STRINGN(J:J).NE.' ') THEN
            CNUM2(1:1)=STRINGN(J:J)
            GOTO 138
         ENDIF
      ENDDO
 138  CONTINUE
      N=2
      DO K=J+1,200
         IF (STRINGN(K:K).NE.' ') THEN
         CNUM2(N:N)=STRINGN(K:K)   
         N=N+1
      ELSE
         GOTO 139
         ENDIF
      ENDDO
 139  N2=N-1
      N=1
      DO L=K+1,200
         IF (STRINGN(L:L).NE.' ') THEN
            CNUM3(N:N)=STRINGN(L:L)
            N=N+1
         ENDIF
      ENDDO
      N3=N-1

         READ (CNUM1(1:N1),3) X
         READ (CNUM2(1:N2),3) Y
         READ (CNUM3(1:N3),3) Z
      
 3       FORMAT (E25.14)

         RETURN
         END




      subroutine convoly (FWHM,IMAGE,CONVOLU,NX,NY)
c Returns convolved image (CONVOLU) of array IMAGE by a Gaussian function
c having width=FWHM (pixels). Separable convolution is applied.
c Borders of width wbar (see below) will appear in convolved image --      
      IMPLICIT NONE
      include 'image.inc'
      integer NX,NY
      integer i,w,wbar,ix,iy,its,itest
      real*4 fwhm,sigma,s2,val,TS
      REAL*8 PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
      real*4 image(nxmax,nymax)
      real*4 aux(nxmax,nymax)
      real*4 convolu(nxmax,nymax),garr(200)
      COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG


      sigma=FWHM/2.35482
      s2=sigma**2
      TS=10.*sigma !5.*SIGMA  ! KERNEL WIDTH ~ 5 sigma


      
      ITS=INT(TS)
      ITEST=MOD(ITS,2)
      W=ITS
      IF (ITEST.EQ.0) W=ITS+1  ! Making the kernel with odd dimension
      IF (W.GT.200) STOP ' REDIM ARRAY garr '
      wbar=(w-1)/2

c Convolution kernel (Gaussian):
      do i=1,w
         garr(i)=SNGL( (1.D0/(DSQRT(TWOPI)*DBLE(SIGMA)))*
     &           DBLE(exp(-(i-wbar-1)**2/(2.*s2))) )
      enddo
c Convolve image along horizontal direction:       

      do iy=1,ny
         do ix=wbar+1,nx-wbar
            val=0
            do i=1,w
               val=val+garr(i)*image(ix+wbar-i+1,iy)
            enddo
            aux(ix,iy)=val            
         enddo
      enddo

c Convolve array aux along vertical direction and build output: 
       do iy=wbar+1,ny-wbar
         do ix=1,nx
            val=0
            do i=1,w
               val=val+garr(i)*aux(ix,iy+wbar-i+1)
            enddo
           convolu(ix,iy)=val
         enddo
      enddo

      return
      end
c-------------------------------

       SUBROUTINE ELEMENTS_TO_XV (GM,T,QR,EC,PER_JD,
     & XH,YH,ZH,VXH,VYH,VZH,THETA)
       IMPLICIT NONE
       REAL*8 PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
       REAL*8 XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
       REAL*8 GM,T,QR,EC,PER_JD
       REAL*8 XH,YH,ZH,VXH,VYH,VZH,R,THETA
       REAL*8 AX,ENE,EME,EE,H
       REAL*8 X,Y,Z,VC,VX,VY,VZ
       REAL*8 EKEPL2,HKEPLER
       REAL*8 CC,WPAR,ZPAR,ONETHIRD
       COMMON /HELIOMATRIX/XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
       COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
       
       IF (EC.EQ.1.D0) GOTO 10
       AX=QR/(1.D0-EC)
 10    CONTINUE
       
       IF (EC.LT.1.D0) THEN ! ELLIPTIC TRAJECTORY
          ENE=DSQRT(GM/AX**3)
            EME=ENE*(T-PER_JD)
            EE=EKEPL2(EME,EC)
            THETA=2.D0*DATAN(DSQRT((EC+1.D0)/(1.D0-EC))*DTAN(EE/2.D0))
            R=AX*(1.D0-EC*EC)/(1.D0+EC*DCOS(THETA))
       ELSE IF (EC.GT.1.D0) THEN  ! HYPERBOLIC TRAJECTORY
           ENE=DSQRT(-GM/AX**3)
           EME=ENE*(T-PER_JD)
           EE=HKEPLER(EME,EC)
           R = -AX*(EC*DCOSH(EE)-1.D0)
	   THETA = 2.*DATAN(DSQRT((EC+1.D0)/(EC-1.D0))*DTANH(EE/2.D0))
        ELSE                    !  PARABOLIC TRAJECTORY 
          ONETHIRD=1.D0/3.D0
          CC=3.D0*DSQRT(GM/(2.D0*QR**3))*(T-PER_JD)
          WPAR=4.D0*CC+DSQRT(64.D0+16.D0*CC**2)
          ZPAR=0.5*WPAR**ONETHIRD-2.D0*WPAR**(-ONETHIRD)
          THETA=2.*DATAN(ZPAR)
          R=QR*(1.D0+ZPAR**2)
       ENDIF

      X=R*DCOS(THETA)
      Y=R*DSIN(THETA)
      Z=0.0D0

C TRANSFORM FROM PLANE OF ORBIT TO HELIOCENTRIC COORDINATES
C (ELEMENTS OF TRANSFORMATION MATRIX ARE ALREADY IN COMMON HELIOMATRIX)
      CALL HPO_TO_HE(X,Y,Z,XH,YH,ZH)
C VELOCITIES
      IF (EC.LT.1.D0) THEN    !ELLIPTIC COMET
         VC=DSQRT(GM*AX)/R
         VX=-VC*DSIN(EE)
         VY=VC*DSQRT(1.D0-EC*EC)*DCOS(EE)
         VZ=0.0D0
      ELSE IF (EC.GT.1.D0) THEN !HYPERBOLIC COMET
         VC=-DSQRT(-GM*AX)/R
         H=2.D0*DATANH(DSQRT((EC-1.D0)/(EC+1.D0))*DTAN(THETA/2.D0))
         VX=VC*DSINH(H)
         VY=-VC*DSQRT(EC*EC-1.D0)*DCOSH(H)
         VZ=0.0D0
      ELSE                       !PARABOLIC COMET
          VY = DSQRT(2.D0*GM/QR)/(1.D0+ZPAR**2)
          VX = -ZPAR * VY
          VZ=0.D0
       ENDIF
       
C TRANSFORM FROM PLANE OF ORBIT TO HELIOCENTRIC VELOCITY VECTOR

      CALL HPO_TO_HE(VX,VY,VZ,VXH,VYH,VZH)

      RETURN
      END
C -----------
      SUBROUTINE CALDATE (JD,DD,MM,YY)
C ALGOR. BASED ON "Practical Astronomy with your calculator" by P.Duffett-Smith 
      IMPLICIT NONE
      REAL*8 JD,XJD,AR,F,DD,MM,YY
      INTEGER I,A,B,C,D,E,G
      XJD=JD+0.5D0
      I=INT(XJD)
      F=XJD-I
      IF (I.GT.2299160) THEN
         AR=((I-1867216.25D0)/36524.25D0)
         A=INT(AR)
         B=I+1+A-INT(FLOAT(A)/4.)
      ELSE
         B=I
      ENDIF
      C=B+1524
      D=INT((DFLOAT(C)-122.1D0)/365.25D0)
      E=INT(365.25D0*FLOAT(D))
      G=INT(DFLOAT(C-E)/30.6001D0)
      DD=C-E+F-INT(30.6001D0*DFLOAT(G))
      IF (G.LT.13.5) THEN
         MM=G-1
      ELSE
         MM=G-13
      ENDIF
      IF (MM.GT.2.5) THEN
         YY=D-4716
      ELSE
         YY=D-4715
      ENDIF
      
      RETURN
      END
      
      
c------------------------------------------------------------------
      
      SUBROUTINE STDCOOR(ALPHA0_IN,DELTA0_IN,ALPHA_IN,DELTA_IN,X,Y)
C CONVERT (RA,DEC) TO STANDARD COORDINATES (X,Y)
C (RA,DEC) IN DEGREE
      IMPLICIT NONE
      REAL*8 PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
      REAL*8 ALPHA0_IN,DELTA0_IN,ALPHA_IN,DELTA_IN
      REAL*8 ALPHA0,DELTA0,ALPHA,DELTA,X,Y,DENO
      COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
      


C CONVERT ALL INPUTS TO RADIANS      
      ALPHA0=ALPHA0_IN*TORAD
      DELTA0=DELTA0_IN*TORAD
      ALPHA=ALPHA_IN*TORAD
      DELTA=DELTA_IN*TORAD
      
      DENO=DCOS(DELTA0)*DCOS(DELTA)*DCOS(ALPHA-ALPHA0)+
     &     DSIN(DELTA)*DSIN(DELTA0)
      X=DCOS(DELTA)*DSIN(ALPHA-ALPHA0)/DENO
      Y=(DSIN(DELTA0)*DCOS(DELTA)*DCOS(ALPHA-ALPHA0)-
     &     DCOS(DELTA0)*DSIN(DELTA))/DENO
      RETURN
      END

      SUBROUTINE AFRHOMAG (NX,NY,INUC,JNUC,MAGSUN,SCALE,DELTA,RCOBS,
     & GRDSIZ,PHAS_ANG,RHO_AP,AFRHO,AFRHO_0,MAG)

c COMPUTES MAGNITUDE AND AFRHO OF MODELED TAIL WITHIN APERTURE RHO_AP (in km) 
C OUTPUT AFRHO is in metres
        IMPLICIT NONE
        INCLUDE 'image.inc'
        REAL*8 FCOMET,FSUN,RADION,GRDSIZ,PHAS_ANG,RHO_AP,FIMAG
        REAL*8 XMAG,SCALE,CTE,DELTA,RCOBS,AFRHO,AFRHO_0,MAG
        INTEGER I,J,NX,NY,INUC,JNUC
        REAL*8 PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
        REAL*4 MAGSUN,IMAGE(NXMAX,NYMAX)
        REAL*8 XP,CORREC,C(7)
C SIXTH ORDER POLYNOMIAL FIT TO SCHLEICHER PHASE FUNCTION
        DATA C/-7.4004978D-03,-1.6492566D-02,1.0950353D-04,
     & 8.3640600D-07,1.0157539D-09,-9.6882641D-11,4.4184372D-13/
        
        COMMON/OUTPUT_IMAGE/IMAGE
        COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
        
C m_star= CTE_MAG + MAGSUN -2.5*log10(i/i_0) (mag/arcsec**2)
        
        Fcomet=0.0
	DO I=1,NX
	   DO J=1,NY
	     RADION=GRDSIZ*DSQRT(DFLOAT(INUC-I)**2+DFLOAT(JNUC-J)**2)
	     IF (RADION.LE.RHO_AP) THEN
             FIMAG=IMAGE(I,J)
             IF (FIMAG.LE.0.0) FIMAG=1.E-25
           XMAG=CTE_MAG+MAGSUN-2.5*DLOG10(FIMAG)-2.5*DLOG10(SCALE*SCALE)
             Fcomet=Fcomet+10**(-XMAG/2.5)
             ENDIF
           ENDDO
        ENDDO
        MAG=-2.5D0*DLOG10(Fcomet)
        Fsun=10**(-MAGSUN/2.5D0)
        CTE=(2.D0*DELTA*AUKM*RCOBS/RHO_AP)**2
	AFRHO=CTE*(Fcomet/Fsun)*RHO_AP*1.D3
        XP=PHAS_ANG
        CORREC=10**(C(1)+C(2)*XP+C(3)*XP**2+
     &  C(4)*XP**3+C(5)*XP**4+C(6)*XP**5+C(7)*XP**6)
        AFRHO_0=AFRHO/CORREC ! Afrho at 0 deg phase angle
	 RETURN
	 END

      SUBROUTINE INTERP5 (NINPUTS,XDTIME,DMDTL,POWER,VFAC,RADMIN,RADMAX)
C SIMPLE LINEAR INTERPOLATION ON THE dmdt_vel_power_rmin_rmax.dat FILE
C TIME IS ASSUMED TO BE MONOTONICALLY INCREASING             
      IMPLICIT NONE      
      include 'ninputs.inc'

       REAL*8 DTIME(NINPUTSMAX),DMDTLOG(NINPUTSMAX),VELFAC(NINPUTSMAX)
       REAL*8 POWERA(NINPUTSMAX),RADIOMIN(NINPUTSMAX)
       REAL*8 RADIOMAX(NINPUTSMAX),M,B,DELTA
       REAL*8 XDTIME,DMDTL,POWER,VFAC,RADMIN,RADMAX
       INTEGER I,NINPUTS
       COMMON/INTERP/DTIME,DMDTLOG,VELFAC,POWERA,RADIOMIN,RADIOMAX

       IF (XDTIME.LT.DTIME(1).OR.XDTIME.GT.DTIME(NINPUTS)) THEN
       WRITE (*,*) ' Time to perihelion=',XDTIME
       STOP ' Requested time out of bounds - revise input file'
       ENDIF
       
       IF (XDTIME.EQ.DTIME(NINPUTS)) THEN
          DMDTL=DMDTLOG(NINPUTS)
          VFAC=VELFAC(NINPUTS)
          POWER=POWERA(NINPUTS)
          RADMIN=RADIOMIN(NINPUTS)
          RADMAX=RADIOMAX(NINPUTS)
         RETURN
       END IF      
        
      DO I=1,NINPUTS
         IF (DTIME(I).GT.XDTIME) THEN

            DELTA=DTIME(I)-DTIME(I-1)

            M=(DMDTLOG(I)-DMDTLOG(I-1))/DELTA
            B=DMDTLOG(I)-M*DTIME(I)
            DMDTL=M*XDTIME+B

            M=(VELFAC(I)-VELFAC(I-1))/DELTA
            B=VELFAC(I)-M*DTIME(I)
            VFAC=M*XDTIME+B

            M=(POWERA(I)-POWERA(I-1))/DELTA
            B=POWERA(I)-M*DTIME(I)
            POWER=M*XDTIME+B

            M=(RADIOMIN(I)-RADIOMIN(I-1))/DELTA
            B=RADIOMIN(I)-M*DTIME(I)
            RADMIN=M*XDTIME+B

             M=(RADIOMAX(I)-RADIOMAX(I-1))/DELTA
             B=RADIOMAX(I)-M*DTIME(I)
             RADMAX=M*XDTIME+B
           
            RETURN
         ENDIF
      ENDDO

      END

         
      SUBROUTINE ANISOT_DIR2(NUC_PHI,NUC_INC,THETAC_EJEC,
     &       TIME,PER_JD,PERIOD,AREA_LATITUDE,AREA_LONGITUDE,
     &       THETA0,UR,UTHETA,UZ)
      IMPLICIT NONE
      REAL*8 PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
      REAL*8 NUC_PHI,NUC_INC,THETAC_EJEC,TIME,PER_JD,PERIOD
      REAL*8 AREA_LATITUDE,AREA_LONGITUDE,THETA0,THETA0_T,THETA
      REAL*8 FINU,COSFINU,SINFINU,SINI,COSI,UR,UTHETA,UZ
      REAL*8 A11,A12,A13,A21,A22,A23,A31,A32,A33,V1,V2,V3
      
c Based on formulae by Sekanina (1981) and Sekanina & Larson (1984)      
c Gives particle velocity components: ur,utheta,uz, as a function of:
c Time of ejection: TIME (JD,DAYS)       
c NUC_INC:(rad)obliquity (angle between the rotation axis and comet orbit plane)
c NUC_PHI: (rad)argument of the subsolar meridian at perihelion
c thetac_ejec: (rad)true anomaly of the comet at time TIME
c area_longitude: (rad) active area longitude
c area_latitude:  (rad) active area latitude
c period: Rotation period of the comet, in days

        COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG

	FINU=NUC_PHI+THETAC_EJEC
	COSFINU=DCOS(FINU)
	SINFINU=DSIN(FINU)
	COSI=DCOS(NUC_INC)
	SINI=DSIN(NUC_INC)

	A11=COSFINU
	A12=COSI*SINFINU
	A13=SINI*SINFINU

	A21=-SINFINU
	A22=COSI*COSFINU
	A23=SINI*COSFINU

	A31=0.D0
	A32=SINI
	A33=-COSI

        THETA0=DATAN(DTAN(THETAC_EJEC+NUC_PHI)*COSI)
        THETA0_T=DATAN(DTAN(NUC_PHI)*COSI)

       THETA=(TWOPI/PERIOD)*(TIME-PER_JD)+THETA0_T+THETA0-AREA_LONGITUDE
        
 	V1=DCOS(AREA_LATITUDE)*DCOS(THETA+THETA0)
	V2=DCOS(AREA_LATITUDE)*DSIN(THETA+THETA0)
	V3=DSIN(AREA_LATITUDE)

c Output vector (ur,utheta,uz):
	UR =    -(A11*V1+A12*V2+A13*V3)
	UTHETA= -(A21*V1+A22*V2+A23*V3)
	UZ=     -(A31*V1+A32*V2+A33*V3)

        RETURN
        END
