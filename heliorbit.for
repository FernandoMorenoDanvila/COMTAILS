

      SUBROUTINE HELIORBIT (PER_JD,TC,TAU,QP,PDEN,RAD,
     &  XC_EJEC,YC_EJEC,ZC_EJEC,VXC_EJEC,VYC_EJEC,VZC_EJEC,
     &  RC_OBS,THETAC_OBS,
     &  X,Y,Z,VX,VY,VZ,NPAR,MPAR,LPAR)
      
C PER_JD: perihelion day      
C TC: OBSERVATION TIME RELATIVE TO PERIHELION (days)      
C TAU: EJECTION TIME RELATIVE TO TC (days)
C E,Q: Comet eccentricity and perihelion distance (au)
C QP: Scattering efficiency for radiation pressure (QP=1 usually)
C PDEN: particle density, MKS
C RAD: particle radius (meters)            
C XC_EJEC,YC_EJEC,ZC_EJEC Heliocentric position of comet
C VXC_EJEC,VYC_EJEC,VZC_EJEC Heliocentric velocity of comet      
C X,Y,Z,VX,VY,VZ Heliocentric position & velocity of particle rel. to comet
C Units: position, au; velocity, au/day      
C RC_EJEC,THETAC_EJEC: Comet helioc. dist (au) and true anom (rad) at ejec time
C XC_EJEC,YC_EJEC,ZC_EJEC: Comet helioc. coordinates, au
C OUTPUT:(N,M) SKY PLANE COORDINATES, and L, all in km (Finson-Probstein)

      IMPLICIT NONE
      REAL*8 QD,ED,IND,WD,ND
      REAL*8 XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
      REAL*8 XXP,XYP,XZP,YXP,YYP,YZP,ZXP,ZYP,ZZP
      REAL*8 DELTA,NMPAR,NMPAR1,NMPAR2,NMPAR3,NMPAR4,NMPAR5
      REAL*8 NMPAR6,NMPAR7,NMPAR8
      REAL*8 EME,EPS1
      REAL*8 PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
      REAL*8 UMU,GM,TIEM0,TC,TAU,PER_JD
      REAL*8 X,Y,Z,XC_EJEC,YC_EJEC,ZC_EJEC,XOUTP,YOUTP,ZOUTP
      REAL*8 VX,VY,VZ,VXC_EJEC,VYC_EJEC,VZC_EJEC,VXOUTP,VYOUTP,VZOUTP
      REAL*8 RC_OBS,THETAC_OBS,THETAD,T0D,AD,NPAR,MPAR,LPAR
      REAL*8 QP,PDEN,RAD,ENE

      COMMON/ORBITAL_ELEMENTS_PARTICLE/QD,ED,IND,WD,ND
      COMMON /HELIOMATRIX/XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ 
      COMMON /HELIOMATRIX_PARTICLE/XXP,XYP,XZP,YXP,YYP,YZP,ZXP,ZYP,ZZP
      COMMON/PARAM_NM/DELTA,NMPAR,NMPAR1,NMPAR2,NMPAR3,NMPAR4,
     &NMPAR5,NMPAR6,NMPAR7,NMPAR8
      COMMON/KEPLER/EME,EPS1
      COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG

C      DATA MU/2.959122082855911D-4/ ! GM 
c      write (*,*) PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
c      pause

      
      TIEM0=TC-TAU+PER_JD   ! Ejection time, t0 (JD)
C Parameters (1-mu) and GM for the given radius radius (RAD):

      UMU=1.191D-3*QP/(2.D0*PDEN*RAD)
      GM=(1.D0-UMU)*MU          ! REDUCED mu*GM


      
C HELIOCENTRIC ECLIPTIC POSITION AND VELOCITY (TOTAL) OF PARTICLE:
        
        XOUTP=X+XC_EJEC
        YOUTP=Y+YC_EJEC
        ZOUTP=Z+ZC_EJEC

        VXOUTP=VX+VXC_EJEC
        VYOUTP=VY+VYC_EJEC
        VZOUTP=VZ+VZC_EJEC

        
        IF (UMU.GT.1.D0) THEN     ! REPULSIVE HYPERBOLA

      CALL HE_TO_ORBITAL_ELEMENTS_R (-GM,TIEM0,XOUTP,
     &  YOUTP,ZOUTP,VXOUTP,VYOUTP,VZOUTP,
     &  QD,ED,IND,ND,WD,AD,THETAD,ENE,T0D)
      ELSE                  ! EITHER ELLIPTIC OR ATTRACTIVE HYPERBOLA
      CALL HE_TO_ORBITAL_ELEMENTS(UMU,GM,TIEM0,
     &  XOUTP,YOUTP,ZOUTP,VXOUTP,VYOUTP,VZOUTP,
     &  QD,ED,IND,ND,WD,AD,THETAD,ENE,T0D)
      ENDIF

      
C COMPUTE (N,M) COORDINATES (SKY PLANE COORDINATES): 
      CALL NM(GM,TIEM0,TC,TAU,
     &     RC_OBS,THETAC_OBS,THETAD,T0D,UMU,AD,NPAR,MPAR,LPAR)
C CONVERT TO KM:
      NPAR=NPAR*AUKM
      MPAR=MPAR*AUKM
      LPAR=LPAR*AUKM
      RETURN
      END
      
      SUBROUTINE HE_TO_ORBITAL_ELEMENTS (UMU,GM,TIEM0,X,Y,Z,
     &      VX,VY,VZ,Q,E,I,OM,W,A,ANOM,ENE,T0D)
C Follows Sterne, An Introduction to Celestial Mechanics pp 56-58 
      IMPLICIT NONE
      REAL*8 PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
      REAL*8 GM,X,Y,Z,R,VX,VY,VZ,V,V2,Q,E,I,W,OM,A,ZZ
      REAL*8 HX,HY,HZ,U,RDOT,H,H2,SINI,COSI,RV,TIEM0,T0D
      REAL*8 ANOM,ENE,TEMP1,TEMP2,UMU

      COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG

      
C ANGULAR MOMENTUM VECTOR:
      CALL VECTORIAL (X,Y,Z,VX,VY,VZ,HX,HY,HZ)
      H2=HX*HX+HY*HY+HZ*HZ
      H=DSQRT(H2)
      R=DSQRT(X*X+Y*Y+Z*Z)
      V2=VX*VX+VY*VY+VZ*VZ
      V=DSQRT(V2)
      RV=X*VX+Y*VY+Z*VZ
C SEMI-MAJOR AXIS FROM VIS VIVA EQUATION:
      A=1.D0/(2.D0/R-V2/GM)
C ECCENTRICITY:      
      E=DSQRT(1.D0-H2/(GM*A))
C  PERIHELION DISTANCE:
      IF (E.NE.1.D0) THEN
         Q=A*(1.D0-E)
      ELSE
      WRITE (*,*) ' PARABOLIC DUST (1) - Possible error '
         Q=H2/(2.D0*GM)  ! E=1, Parabolic dust
      ENDIF
     
c ----------------
C NODE:
      OM=DATAN2(HX,-HY)
      IF(OM.LT.0.D0) OM=OM+TWOPI      
C INCLINATION:                 
      COSI=HZ/H
      SINI=(HX/H)/DSIN(OM)
      I=DATAN2(SINI,COSI)
      


      
      U=DATAN2(Z/SINI,X*DCOS(OM)+Y*DSIN(OM))
      RDOT=(X*VX+Y*VY+Z*VZ)/R
C TRUE ANOMALY:      
      ANOM=DATAN2(H*RDOT/GM,H2/(R*GM)-1.D0)
C ARGUMENT OF PERIHELION:
      W=U-ANOM
      IF (W.LT.0.D0) W=W+TWOPI
      IF (W.GT.TWOPI) W=MOD(W,TWOPI)
      IF (E.LT.1.D0.AND.UMU.LT.1.D0) THEN ! ELLIPTIC TRAJECTORY
      ENE=DSQRT(GM/A**3)
      TEMP1=DSQRT((1.D0+E)/(1.D0-E))
      TEMP2=2.*DATAN((1.D0/TEMP1)*DTAN(ANOM/2.D0)) !eccentric anomaly
      T0D=TIEM0-(1.D0/ENE)*(TEMP2-E*DSIN(TEMP2))
      ELSE IF (E.GT.1.D0) THEN            ! HYPERBOLIC (ATTRACTIVE)
      TEMP1 = (E-1.D0)/(E+1.D0)
      TEMP2=DSQRT(TEMP1)*DTAN(ANOM/2.D0)
      TEMP1=2.*DATANH(TEMP2)
      ENE=A**3/GM
      T0D=TIEM0-DSQRT(-ENE)*(E*DSINH(TEMP1)-TEMP1)
      ENE=DSQRT(-GM/A**3)
      ELSE             ! E=1, PARABOLIC DUST (VERY UNLIKELY CASE)
         WRITE (*,*) ' PARABOLIC DUST (2) -- Posible NaN error: '
         WRITE (*,*) ' UMU,GM,VELOCITY: ',UMU,GM,vx,vy,vz
      ZZ=DTAN(ANOM/2.D0)
      ENE=DSQRT(GM/(2.D0*Q**3))
      T0D=TIEM0-(1.D0/ENE)*(ZZ+ZZ**3/3.D0)         
      ENDIF    
      
      RETURN
      END

      
      SUBROUTINE HE_TO_ORBITAL_ELEMENTS_R (GM,TIEM0,X,Y,Z,
     & VX,VY,VZ,Q,E,I,OM,W,A,ANOM,ENE,T0D)
C Follows Sterne, BUT FOR REPULSIVE HYPERBOLAE
      IMPLICIT NONE
      REAL*8 GM,X,Y,Z,R,VX,VY,VZ,V,V2,Q,E,I,W,OM,A
      REAL*8 HX,HY,HZ,U,H,H2,COSI,SINI,RV,ENE,F
      REAL*8 ANOM,RDOT,TIEM0,T0D
      REAL*8 PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
      COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG

C ANGULAR MOMENTUM VECTOR:
      CALL VECTORIAL (X,Y,Z,VX,VY,VZ,HX,HY,HZ)
      H2=HX*HX+HY*HY+HZ*HZ
      H=DSQRT(H2)
      R=DSQRT(X*X+Y*Y+Z*Z)
      V2=VX*VX+VY*VY+VZ*VZ
      V=DSQRT(V2)
      RV=X*VX+Y*VY+Z*VZ
C SEMI-MAJOR AXIS FROM VIS-VIVA EQUATION:
      A=1.D0/(V**2/GM+2.D0/R)
      A=-DABS(A)
C ECCENTRICITY:      
      E=DSQRT(1.D0-H2/(GM*A))  
C PERIHELION DISTANCE:
      Q=DABS(A)*(1.D0+E)
C NODE:
      OM=DATAN2(HX,-HY)
      IF(OM.LT.0.D0) OM=OM+TWOPI      
C INCLINATION:                 
      COSI=HZ/H
      SINI=(HX/H)/DSIN(OM)
      I=DATAN2(SINI,COSI)

      U=DATAN2(Z/SINI,X*DCOS(OM)+Y*DSIN(OM))      
      RDOT=(X*VX+Y*VY+Z*VZ)/R
      
C TRUE ANOMALY:      
      ANOM=DATAN2(H*RDOT/GM,H2/(R*GM)+1.D0)

C ARGUMENT OF PERIHELION:
      W=U-ANOM
      IF (W.LT.0.D0) W=W+TWOPI
      IF (W.GT.TWOPI) W=MOD(W,TWOPI)

c  ECCENTRIC ANOMALY:
      F=2.*DATANH(TAN(ANOM/2.D0)*DSQRT((E+1.D0)/(E-1.D0)))
      ENE=DSQRT(GM/DABS(A))

C  PERIHELION TIME:
      T0D=TIEM0-(E*DSINH(F)+F)/DSQRT(GM/DABS(A)**3)

      RETURN
      END
      
      SUBROUTINE VECTORIAL (A1,A2,A3,B1,B2,B3,C1,C2,C3)
      IMPLICIT NONE
      REAL*8 A1,A2,A3,B1,B2,B3,C1,C2,C3
C GIVES CROSS PRODUCT C=AxB
      C1=A2*B3-A3*B2
      C2=A3*B1-A1*B3
      C3=A1*B2-A2*B1
      RETURN
      END
         

         
      SUBROUTINE TRANSFORM_MATRIX_PARTICLE (XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ)
      IMPLICIT NONE
      REAL*8 XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
      REAL*8 QD,ED,IND,WD,ND
      REAL*8 CW,SW,CI,SI,CN,SN
C COMPUTES TRANSFORM_MATRIX 
        COMMON/ORBITAL_ELEMENTS_PARTICLE/QD,ED,IND,WD,ND
        
        CW=DCOS(WD)
        SW=DSIN(WD)
        CI=DCOS(IND)
        SI=DSIN(IND)
        CN=DCOS(ND)
        SN=DSIN(ND)
        
      XX=CW*CN - CI*SN*SW 
      YX=-SW*CN - CI*SN*CW
      ZX=SI*SN

      XY=CW*SN + CI*CN*SW
      YY=-SW*SN + CI*CN*CW
      ZY=-SI*CN

      XZ=SW*SI
      YZ=CW*SI
      ZZ=CI

      
      
      RETURN
      END
      
       SUBROUTINE NM(GM,TIEM0,TC_IN,TAU_IN,RC,THETAC,THETAD,
     &     T0DEX,UMU,AD,NPAR,MPAR,LPAR)
C L COORDINATE ADDED TO COMPUTE OPTICAL DEPTH ALONG NUCLEUS       
       IMPLICIT NONE
C GIVES (N,M) SKY PLANE COORDINATES in AU
       REAL*8 TC_IN,TAU_IN
       REAL*8 GM,TC,TAU,RC,THETAC,THETAD,UMU,AD,NPAR,MPAR,LPAR
       REAL*8 QD,ED,IND,WD,ND,EPS1,ENE,EME,EEP,RTBIS,RD
       real*8 XD,YD,ZD
       REAL*8 EKEPL2,HKEPLER
       REAL*8 XXP,XYP,XZP,YXP,YYP,YZP,ZXP,ZYP,ZZP
       REAL*8 XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ
       REAL*8 XDE,YDE,ZDE,XOUT,YOUT,ZOUT
c       REAL*8 ATANH
       REAL*8 DELTA,NMPAR,NMPAR1,NMPAR2,NMPAR3,NMPAR4,NMPAR5
       REAL*8 NMPAR6,NMPAR7,NMPAR8
       REAL*8 CHITA,ETA,GITA
       REAL*8 PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
       REAL*8 SINT,COST,TIEM0,T0DEX
       
       
       EXTERNAL RTBIS,EKEPL2,HKEPLER  !,ATANH
       
       COMMON/ORBITAL_ELEMENTS_PARTICLE/QD,ED,IND,WD,ND
       COMMON /HELIOMATRIX/XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ 
       COMMON /HELIOMATRIX_PARTICLE/XXP,XYP,XZP,YXP,YYP,YZP,ZXP,ZYP,ZZP
       COMMON/PARAM_NM/DELTA,NMPAR,NMPAR1,NMPAR2,NMPAR3,NMPAR4,NMPAR5,
     & NMPAR6,NMPAR7,NMPAR8
       COMMON/KEPLER/EME,EPS1
       COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
       

        TC=TC_IN
        TAU=TAU_IN


        EPS1=ED

        
       IF (ED.LT.1.D0.AND.UMU.LT.1.D0) THEN ! ELLIPTIC TRAJECTORY
            ENE=DSQRT(GM/AD**3)
            EME=ENE*(tiem0-t0dex+tau)
            EEP=EKEPL2(EME,ED)
            THETAD=2.D0*DATAN(DSQRT((ED+1.D0)/(1.D0-ED))*DTAN(EEP/2.D0))
            RD=AD*(1.D0-ED*ED)/(1.D0+ED*DCOS(THETAD))
            
        ELSE                            ! HYPERBOLIC TRAJECTORIES (2 TYPES):

           IF ((1.D0-UMU).GT.0.D0) THEN ! Attractive hiperbola
           ENE=DSQRT(-GM/AD**3)
           EME=ENE*(tiem0-t0dex+tau)
           EEP=HKEPLER(EME,ED)
           RD = -AD*(ED*DCOSH(EEP)-1.D0)
	   THETAD = 2.*DATAN(DSQRT((ED+1.D0)/(ED-1.D0))*DTANH(EEP/2.D0))

        ELSE                             ! Repulsive hyperbola, 1-mu > 0  

           ENE=DSQRT(GM/AD**3)
           EME=ENE*(tiem0-t0dex+tau)

           EPS1=ED
           EEP=HKEPLER(-EME,-ED)
           thetad=2.*datan( dsqrt( (ed-1.d0)/(ed+1.d0) )*tanh(EEP/2.d0))
           rd=abs(ad)*(ed*cosh(EEP)+1.)
        ENDIF

      ENDIF





       XD=RD*DCOS(THETAD)
       YD=RD*DSIN(THETAD)
       ZD=0.0
      
C CALCULATE TRANSFORMATION MATRIX FROM DUST ORBIT PLANE TO HELIOC. ECLIPTIC
        CALL TRANSFORM_MATRIX_PARTICLE(XXP,XYP,XZP,YXP,YYP,YZP,
     &     ZXP,ZYP,ZZP)
        
C COOR HELIOCENTRIC ECLIPTIC OF PARTICLE
        CALL HPO_TO_HE_PARTICLE(XD,YD,ZD,XDE,YDE,ZDE)

C TRANSFORM ECLIPTIC HELIOCENTRIC COORDINATES TO PLANE-OF-ORBIT
C COORDINATES, WHERE NOW, THE PLANE-OF-ORBIT IS THE COMET ORBIT PLANE 
        CALL HE_TO_HPO(XDE,YDE,ZDE,XOUT,YOUT,ZOUT)

C COMPUTE COORDINATES IN THE COMETOCENTRIC FRAME BY FINSON & PROBSTEIN (1968)

        SINT=DSIN(THETAC)
        COST=DCOS(THETAC)
        CHITA=XOUT*COST+YOUT*SINT-RC
        ETA=  XOUT*SINT-YOUT*COST
        GITA=ZOUT

c Convert to plane-of-sky (L,M,N) coordinates (Finson-Probstein):
	MPAR=NMPAR1*CHITA - NMPAR2*ETA - NMPAR3*GITA
	NPAR=NMPAR4*ETA -   NMPAR5*GITA
        LPAR=NMPAR6*CHITA + NMPAR7*ETA + NMPAR8*GITA 

        
        RETURN
        END

      SUBROUTINE HPO_TO_HE_PARTICLE(XIN,YIN,ZIN,XOUT,YOUT,ZOUT)

C CONVERT HELIOCENTRIC PLANE-OF-ORBIT COORDINATES (XIN,YIN,ZIN) TO 
C HELIOCENTRIC ECLIPTIC COORDINATES (XOUT,YOU,ZOUT)
C MEDIANTE LA MATRIZ DE TRANSFORMACION

      IMPLICIT NONE
      REAL*8 XIN,YIN,ZIN,XOUT,YOUT,ZOUT
      REAL*8 XXP,XYP,XZP,YXP,YYP,YZP,ZXP,ZYP,ZZP
      COMMON /HELIOMATRIX_PARTICLE/XXP,XYP,XZP,YXP,YYP,YZP,ZXP,ZYP,ZZP

      XOUT=XXP*XIN + YXP*YIN + ZXP*ZIN
      YOUT=XYP*XIN + YYP*YIN + ZYP*ZIN
      ZOUT=XZP*XIN + YZP*YIN + ZZP*ZIN

      RETURN
      END


      REAL*8 FUNCTION EKEPL2(EM,E)
C SOLVES KEPLER EQUATION EM=EKEPL2-E*DSIN(EKEPL2)
C WITH E IN THE RANGE [0,1]
C CODE BY A.W. ODELL AND R.H. GOODING, Cel. Mech., 38, 307, (1986).
      IMPLICIT NONE
      REAL*8 EM,E,PINEG,SW,AHALF,ASIXTH,ATHIRD,A,B,EMR,EE
      REAL*8 W,E1,FDD,FDDD,F,FD,DEE,EMKEPL
      REAL*8 PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG
      INTEGER ITER
      LOGICAL L
      COMMON/CONSTANTS/PI,HALFPI,TWOPI,AUKM,CTEVEL,TORAD,MU,CTE_MAG

      PINEG=-PI
      SW=0.1D0
      AHALF=0.5D0
      ASIXTH=AHALF/3.D0
      ATHIRD=ASIXTH*2.D0
      A=(PI-1.D0)**2/(PI+2.D0/3.D0)
      B=2.D0*(PI-ASIXTH)**2/(PI+2.D0/3.D0)
      EMR=DMOD(EM,TWOPI)
      IF(EMR.LT.PINEG) EMR=EMR+TWOPI
      IF(EMR.GT.PI) EMR=EMR-TWOPI
      EE=EMR
      
C      IF (EE) 1,4,2

      IF (EE.LT.0.D0) GOTO 1
      IF (EE.EQ.0.D0) GOTO 4
      IF (EE.GT.0.D0) GOTO 2
C ------------      
 1    EE=-EE
 2    IF(EE.LT.ASIXTH) THEN
         EE=(6.D0*EE)**ATHIRD
      ELSE
         W=PI-EE
         EE=PI-A*W/(B-W)
      ENDIF
      IF (EMR.LT.0.D0) EE=-EE
      EE=EMR+(EE-EMR)*E
      E1=1.D0-E
      L=(E1+EE*EE/6.D0).GE.SW
      DO 3 ITER=1,2
         FDD=E*DSIN(EE)
         FDDD=E*DCOS(EE)
         IF (L) THEN
            F=(EE-FDD)-EMR
            FD=1.D0-FDDD
            ELSE
            F=EMKEPL(E,EE)-EMR
            FD=E1+2.D0*E*DSIN(AHALF*EE)**2
          ENDIF
          DEE=F*FD/(AHALF*F*FD-FD*FD)
          W=FD+AHALF*DEE*(FDD+ATHIRD*DEE*FDDD)
          FD=FD+DEE*(FDD+AHALF*DEE*FDDD)
          EE=EE-(F-DEE*(FD-W))/FD
 3        CONTINUE
 4        EKEPL2=EE+(EM-EMR)
          RETURN
          END

      REAL*8 FUNCTION EMKEPL(E,EE)
C AUXILIARY FUNCTION FOR KEPLER SOLVER EKEPL2
      IMPLICIT NONE
      REAL*8 X,E,EE,EE2,TERM,D,X0
      X=(1.D0-E)*DSIN(EE)
      EE2=-EE*EE
      TERM=EE
      D=0.D0
 1    D=D+2.D0
      TERM=TERM*EE2/(D*(D+1.D0))
      X0=X
      X=X-TERM
      IF (X.NE.X0) GOTO 1
      EMKEPL=X
      RETURN
      END


      REAL*8 FUNCTION HKEPLER (EME,E)
c Solves hyperbolic Kepler equation M=e*sinh(H)-H, i.e., gives H given M and e
C CODE BY A.W. ODELL AND R.H. GOODING, Cel. Mech., 44, 267, (1988).
      IMPLICIT REAL*8 (A-H,O-Z)
      EL=EME/E
      G1=1.D0-1.D0/E
	HKEPLER=DASINH(SHKEPL(EL,G1))
      RETURN
      END


      REAL*8 FUNCTION SHKEPL(EL,G1)
C SOLVES EQUATION: EL=SHKEPL+(G1-1)*ASINH(SHKEPL)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (SW=0.5D0,AHALF=0.5D0,
     & ASIXTH=AHALF/3.D0,ATHIRD=ASIXTH*2.D0)

      S=EL
      IF (EL.EQ.0.D0) GOTO 2
      G=1.D0-G1
      CL=DSQRT(1.D0+EL**2)
      AL=DASINH(EL)
      

        w = g**2*al/cl**3
        s = 1.D0 - g/cl
        s = el + g*al/dcubrt(s**3+w*el*(1.5D0 - g/0.75D0)) !**(1.d0/3.d0))

C Two iterations (at most) of halley-then-newton process

        do 1 iter=1,2 
            s0 = s*s
            s1 = s0 + 1.D0
            s2 = dsqrt(s1)
            s3 = s1*s2
            fdd = g*s/s3
            fddd = g*(1.D0 - 2.D0*s0)/(s1*s3)
            if (asixth*s0 + g1 .GE. sw) then
                f = (s - g*dasinh(s)) - el
                fd = 1.D0 - g/s2
            else
                f = shmkep(g1, s) - el
                fd = (s0/(s2 + 1.D0) + g1)/s2
            end if
            ds = f*fd/(ahalf*f*fdd - fd*fd)
            stemp = s + ds
            if (stemp.EQ.s) GOTO 2
            f = f + ds*(fd + ahalf*ds*(fdd + athird*ds*fddd))
            fd = fd + ds*(fdd + ahalf*ds*fddd)
            s = stemp - f/fd
 1          continue
 2          shkepl = s
            RETURN
            END
      
C THIS FUNCTION TO COMPUTE ACCURATELY THE CUBE ROOT      
      REAL*8 FUNCTION DCUBRT(X)
      IMPLICIT NONE
      REAL*8 ATHIRD,X,C,Y
C      IMPLICIT REAL*8 (A-Z)
C      PARAMETER (SW=0.5D0,AHALF=0.5D0,
C     & ASIXTH=AHALF/3.D0,ATHIRD=ASIXTH*2.D0)

      ATHIRD=1.D0/3.D0
      
      IF (X.EQ.0.0D0) THEN
         C=0.D0
      ELSE
         Y=DABS(X)
         C=Y**ATHIRD
         C=C-ATHIRD*(C-Y/C**2)
         C=SIGN(C,X)
      ENDIF
      DCUBRT=C
      RETURN
      END


      
        REAL*8 function shmkep (g1, s)

        implicit none

    
        REAL*8 g1
        REAL*8 s,one,two

        REAL*8 g,t,tsq,x,term,twoi1,x0
        one=1.d0
        two=2.d0
        g = one - g1
        t = s/(one + dsqrt(one + s*s))
        tsq = t*t
        x = s*(g1 + g*tsq)
        term = two*g*t
        twoi1 = one
 1      twoi1 = twoi1 + two
        term = term*tsq
        x0 = x
        x = x - term/twoi1
        if (x.ne.x0) goto 1
        shmkep = x
        return
        end


      SUBROUTINE HPO_TO_HE(XIN,YIN,ZIN,XOUT,YOUT,ZOUT)

C CONVERT HELIOCENTRIC PLANE-OF-ORBIT COORDINATES (XIN,YIN,ZIN) TO 
C HELIOCENTRIC ECLIPTIC COORDINATES (XOUT,YOU,ZOUT)

      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /HELIOMATRIX/XXD,XYD,XZD,YXD,YYD,YZD,ZXD,ZYD,ZZD

      XOUT=XXD*XIN + YXD*YIN + ZXD*ZIN
      YOUT=XYD*XIN + YYD*YIN + ZYD*ZIN
      ZOUT=XZD*XIN + YZD*YIN + ZZD*ZIN

      RETURN
      END
