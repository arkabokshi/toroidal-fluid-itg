PROGRAM itg
USE inputdata,&
ONLY:m0,length,dt,NumSteps,RunTime,xx,dp,NumModes,InitialMode,FinalMode,GenIso_Transition,		&
    eta,qprofile,sigma,tau,c,shear,PtDensity,x0,FlowShear,ShearingRate,dGammaE_dt,restart,		&
    MyRank,MySize,ierror,NumModes_by_2,ModesPerProc,nqp,rs,ktheta,NumSteps,MinRad,			&
    ci,DelPrint,FlowOnOff,Init_gammaE,Final_gammaE,idelta_m, etag, etac, m0, n0, qedge, NumTheta,	&
    Quad_gammaE, TaylorFlow, noiseStart, calcThetaMaxima, gammatol, navg, runpath
IMPLICIT NONE
INCLUDE 'mpif.h'
EXTERNAL Evolve,OdeSolver

DOUBLE COMPLEX, DIMENSION( 3*length,NumModes )  ::   u1 = 0.0_dp, u0 = 0.0_dp
DOUBLE COMPLEX, DIMENSION( NumSteps,NumModes )  ::   Omega_t  = 0.0_dp
DOUBLE COMPLEX, DIMENSION( NumSteps )           ::   PolPot   = 0.0_dp
DOUBLE PRECISION,DIMENSION(NumSteps )		::   GlobalGamma = 0.0, GlobalOmega = 0.0, ThetaMaxima = 0.0
DOUBLE COMPLEX, DIMENSION( length )             ::   ax,bx,cx
REAL(KIND=dp) , DIMENSION( NumSteps )           ::   gammaE_t = -1.0_dp
DOUBLE COMPLEX					::   Omega_0

INTEGER						::   TimeStep,mode,q_m,NumFiles,ii,jj,mm,MaxTheta, &
						     begtime, endtime
REAL(KIND=dp)					::   mu,delm,delx,t0,t1,Amp,tt,wm,gm,RealMax,ImagMax
 CHARACTER(LEN=20)				::   t

DOUBLE PRECISION				::   ABS_PHIm_t0, ABS_PHIm_t1, ABS_Gm_t1, theta
DOUBLE PRECISION, DIMENSION(NumTheta)		::   ABS_Potential
DOUBLE PRECISION, PARAMETER			::   pi = ATAN(1.0_dp) * 4.0_dp
DOUBLE COMPLEX, DIMENSION(length)		::   ComplexPotential
REAL(KIND=dp), DIMENSION( length)   ::   ureal, uimag
REAL(KIND=dp)					::   old_gamma, new_gamma, delta_gamma


! ---------------- !
! Initialising MPI !
! ---------------- !

  CALL MPI_INIT( ierror )
  CALL MPI_COMM_SIZE( MPI_COMM_WORLD,MySize,ierror )
  CALL MPI_COMM_RANK( MPI_COMM_WORLD,MyRank,ierror )
  

  ModesPerProc = NumModes / MySize

  IF (MOD(NumModes,MySize).NE.0) THEN
      PRINT*, 'EXITING: Number of processors should exactly divide number of modes'
      CALL EXIT()
  END IF


!--------------------!
! Setting parameters !
!--------------------!

  eta = etag * ( 1.0_dp - etac * (xx/nqp)**2 )

!--------------------------!
! Create runpath directory !
!--------------------------!

  CALL SYSTEM('mkdir ' // runpath)

! ------------------- !
! Initialising fields !
! ------------------- !

  Omega_0 = DCMPLX(-0.5,0.5 )			! Arbitrary / or would have a relation between them help converge faster
  mu  = 0.1					! No noticable effect
  Amp = 1.0

  DO mode = 1,NumModes
    delm = DBLE( InitialMode-m0-1+mode )
    u0(:length,mode) = ((1)**mode) * Amp * DCMPLX( EXP(-mu*((xx-delm)**2 )),EXP(-mu*((xx-delm)**2)) )
  END DO

  IF (noiseStart) THEN
      CALL RANDOM_NUMBER(ureal)
      CALL RANDOM_NUMBER(uimag)
      u0(:length,1) = Amp * DCMPLX( ureal,uimag )
      DO mode = 2,NumModes
          u0(:length,mode) = u0(:length,mode-1) * ci * mode/NumModes
      END DO
  END IF

  u0 ( 1*length+1 : 2*length, : ) = u0 ( :length , : ) * Omega_0**1
  u0 ( 2*length+1 : 3*length, : ) = u0 ( :length , : ) * Omega_0**2

  
  IF (restart) THEN
    
      OPEN (5,FILE=TRIM(runpath)//'/FinalFields.dat',FORM='unformatted',ACCESS='stream')
      DO mode=1,NumModes
	READ (5) u0(:,mode)
      END DO
      CLOSE (5)

  END IF



! --------------------------- !
! Performing matrix inversion !
! for use in alphainverse     !
! --------------------------- !

  ax =  ( c*(shear**2) )
  bx =  ( 0.0_dp )
  cx = -( c + 1.0_dp + idelta_m + TaylorFlow )

  delx = xx(2)-xx(1)

  CALL odesolver( length,delx,ax,bx,cx )


! ---------------------------------- !
! Writing initial conditions to file !
! ---------------------------------- !

  IF (MyRank.EQ.0) THEN

	OPEN (5,FILE=TRIM(runpath)//'/parameters.txt')
	WRITE (5,*) NumSteps,dt,length,NumModes,InitialMode,FinalMode,tau
	CLOSE (5)
	OPEN (5,FILE=TRIM(runpath)//'/xx.txt')
	WRITE (5,*) xx
	CLOSE (5)
	OPEN (5,FILE=TRIM(runpath)//'/profiles.dat',FORM='unformatted',ACCESS='stream')
	WRITE (5) eta,qprofile,sigma,Final_gammaE*xx + Quad_gammaE*(xx**2)
	CLOSE (5)
	OPEN (5,FILE=TRIM(runpath)//'/100000.dat',FORM='unformatted',ACCESS='stream')
	WRITE (5) REAL( u0(:length,:) )
	CLOSE (5)
	OPEN (5,FILE=TRIM(runpath)//'/200000.dat',FORM='unformatted',ACCESS='stream')
	WRITE (5) AIMAG( u0(:length,:) )
	CLOSE (5)

  END IF



! ----------------------- !
! Evolving system in time !
! ----------------------- !

  old_gamma = -1.0
  NumFiles = 1
  t0 = MPI_WTIME()
  
  DO  TimeStep = 1,NumSteps


	! -------------------------------- !
	! Dynamically adjusting flow-shear !
	! -------------------------------- !

	!  IF (TimeStep*dt.LT.GenIso_Transition) THEN
	!    ShearingRate = Init_gammaE
	!  ELSE IF (TimeStep*dt.GE.GenIso_Transition) THEN
	!    ShearingRate = ShearingRate + dt*dGammaE_dt
	!    IF (ShearingRate.GE.Final_gammaE) ShearingRate = Final_gammaE
	!  END IF

	!  ShearingRate = Init_gammaE
	!  gammaE_t(TimeStep) = ShearingRate

	!IF ( TimeStep*dt.LT.GenIso_Transition ) THEN
	!	FlowShear = Init_gammaE * xx
	!ELSE
	!	FlowShear = (Final_gammaE*xx) + (Quad_gammaE * (xx**2))
	!END IF  


	! ------------- !
	! Evolve fields !
	! ------------- !

	  CALL Evolve(u0,u1)


	! -------- !
	! Analysis !
	! -------- !
	
	  ABS_PHIm_t0  =   0.0_dp
	  ABS_PHIm_t1  =   0.0_dp
	  ABS_Gm_t1    =   0.0_dp

	  DO mode = 1,NumModes

	  ! (1)
	  ! Omega_m vs. Time
	    q_m = x0 + PtDensity*(InitialMode-m0-1+mode)
	    Omega_t ( TimeStep,mode ) = ( u1(q_m+length,mode) / u1(q_m,mode) )


	  ! (2)
	  ! Global [ Omega,Gamma ]
	    ABS_PHIm_t0 = ABS_PHIm_t0  +  SUM( (ABS(u0( 0*length+1: 1*length ,mode )))**2 )
	    ABS_PHIm_t1 = ABS_PHIm_t1  +  SUM( (ABS(u1( 0*length+1: 1*length ,mode )))**2 )
	    ABS_Gm_t1   = ABS_Gm_t1    +  SUM( (ABS(u1( 1*length+1: 2*length ,mode )))**2 )


	  END DO

	IF (calcThetaMaxima) THEN
	  ! (3)
	  ! Max [ Theta ]	  
	  DO jj = 1,NumTheta
		Theta = 2.0d0 * pi * DBLE(jj-1)/DBLE(NumTheta-1)
		ComplexPotential = DCMPLX(0.0d0,0.0d0)
			DO mm = InitialMode,FinalMode
				mode = mm-InitialMode+1
				ComplexPotential = ComplexPotential + u1(:length,mode )*exp(-ci*mode*Theta)
			END DO
		ABS_Potential(jj) = maxval( ABS(ComplexPotential) )
	  END DO
	  ThetaMaxima(TimeStep) = 2d0 * DBLE( maxloc(ABS_Potential,DIM=1) ) / DBLE(NumTheta-1)	! In units of pi radians
	ELSE
	  ThetaMaxima(TimeStep) = 0
	END IF


	  GlobalOmega(TimeStep) = SQRT( ABS_Gm_t1 / ABS_PHIm_t1 )
	  GlobalGamma(TimeStep) = ( DLOG(ABS_PHIm_t1)-DLOG(ABS_PHIm_t0) ) / dt		! NOTE: Not square-rooting, so factor 2 in analysisV2.py needed

	! -------- !
	! Printing !
	! -------- !

	  IF (MyRank.EQ.0) THEN

	  
	      IF ( MOD(TimeStep,DelPrint).EQ.0 ) THEN


		    NumFiles = NumFiles + 1
		    WRITE(t,'(I5.5)') TimeStep

		    OPEN (5,FILE=TRIM(runpath)//'/1'//TRIM(t)//'.dat',FORM='unformatted',ACCESS='stream')
		    WRITE (5) REAL(u1(:length,:))
		    CLOSE (5)
		    OPEN (5,FILE=TRIM(runpath)//'/2'//TRIM(t)//'.dat',FORM='unformatted',ACCESS='stream')
		    WRITE (5) AIMAG(u1(:length,:))
		    CLOSE (5)
		    OPEN (5,FILE=TRIM(runpath)//'/frequency.dat',FORM='unformatted',ACCESS='stream')
		    WRITE (5) REAL(Omega_t)
		    CLOSE (5)
		    OPEN (5,FILE=TRIM(runpath)//'/growthrate.dat',FORM='unformatted',ACCESS='stream')
		    WRITE (5) AIMAG(Omega_t)
		    CLOSE (5)
		    OPEN (5,FILE=TRIM(runpath)//'/gammaE_t.dat',FORM='unformatted',ACCESS='stream')
		    WRITE (5) gammaE_t
		    CLOSE (5)
		    OPEN (5,FILE=TRIM(runpath)//'/globalgamma.dat',FORM='unformatted',ACCESS='stream')
		    WRITE (5) GlobalGamma
		    CLOSE (5)
		    OPEN (5,FILE=TRIM(runpath)//'/globalomega.dat',FORM='unformatted',ACCESS='stream')
		    WRITE (5) GlobalOmega
		    CLOSE (5)
		    OPEN (5,FILE=TRIM(runpath)//'/thetamaxima.dat',FORM='unformatted',ACCESS='stream')
		    WRITE (5) ThetaMaxima
		    CLOSE (5)

	      END IF

	  END IF


	!------------------!
	! NORMALISE FIELDS !
	!------------------!

	  RealMax = MAXVAL(REAL(u1))
	  ImagMax = MAXVAL(AIMAG(u1))
	  u0 = u1 / MAX(RealMax,ImagMax)


	  CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

	IF (MOD(TimeStep,navg).EQ.0) THEN

		! take mean over window with width navg
		begtime = TimeStep+1 - navg
		endtime = TimeStep
		new_gamma = SUM( GlobalGamma(begtime:endtime) ) / navg
		delta_gamma = ABS(new_gamma-old_gamma) / ABS(new_gamma)
		old_gamma = new_gamma

	        PRINT*, TimeStep*dt,'/',NumSteps*dt
		PRINT*, delta_gamma, begtime, endtime

		IF (delta_gamma .LE. gammatol) EXIT

	END IF


  END DO



! ------------------------- !
! Printing final conditions !
! ------------------------- !

  IF (MyRank.EQ.0) THEN
  
      OPEN (5,FILE=TRIM(runpath)//'/videoparam.txt')
      WRITE (5,*) NumFiles,DelPrint
      CLOSE (5)
      OPEN (5,FILE=TRIM(runpath)//'/realfieldend.dat',FORM='unformatted',ACCESS='stream')
      WRITE (5) REAL(u1(:length,:))
      CLOSE (5)
      OPEN (5,FILE=TRIM(runpath)//'/imagfieldend.dat',FORM='unformatted',ACCESS='stream')
      WRITE (5) AIMAG(u1(:length,:))
      CLOSE (5)
      OPEN (5,FILE=TRIM(runpath)//'/frequency.dat',FORM='unformatted',ACCESS='stream')
      WRITE (5) REAL(Omega_t)
      CLOSE (5)
      OPEN (5,FILE=TRIM(runpath)//'/growthrate.dat',FORM='unformatted',ACCESS='stream')
      WRITE (5) AIMAG(Omega_t)
      CLOSE (5)
      OPEN (5,FILE=TRIM(runpath)//'/gammaE_t.dat',FORM='unformatted',ACCESS='stream')
      WRITE (5) gammaE_t
      CLOSE (5)
      OPEN (5,FILE=TRIM(runpath)//'/globalgamma.dat',FORM='unformatted',ACCESS='stream')
      WRITE (5) GlobalGamma
      CLOSE (5)
      OPEN (5,FILE=TRIM(runpath)//'/globalomega.dat',FORM='unformatted',ACCESS='stream')
      WRITE (5) GlobalOmega
      CLOSE (5)
      OPEN (5,FILE=TRIM(runpath)//'/thetamaxima.dat',FORM='unformatted',ACCESS='stream')
      WRITE (5) ThetaMaxima
      CLOSE (5)

      OPEN (5,FILE=TRIM(runpath)//'/FinalFields.dat',FORM='unformatted',ACCESS='stream')
      DO mode=1,NumModes
	WRITE (5)  u1(:,mode)
      END DO


  END IF



! ------------ !
! Finalize MPI !
! ------------ !

  t1 = MPI_WTIME()
  IF (MyRank.EQ.0) PRINT*, 'Run-time(s):',t1-t0

  CALL MPI_FINALIZE( ierror )



END PROGRAM itg
