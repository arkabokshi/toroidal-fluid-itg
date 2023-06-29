PROGRAM FastestGrowingMode
IMPLICIT NONE

INTEGER,PARAMETER::dp=SELECTED_REAL_KIND(15,300)
DOUBLE COMPLEX::aa,bb,cc,omega_sub

REAL(KIND=dp),PARAMETER::		MinRad = 0.5_dp,					&
					MajRad = 5.0_dp,					& 
					m0 = 70.0_dp,						&
					n0 = 50.0_dp,						&
					shear = 25.0_dp,					&
					delta = -1.0_dp,					&	! keep NEGATIVE (not too high or low)
					krho = SQRT(0.04_dp),					&	! LOWER krho
					b = krho**2,						&
					q0 = (m0)/n0							! LOWER q

INTEGER,PARAMETER::			MaxOrder = 100

INTEGER::				i,j,k,e,eta_size,tau_size,MaxGrowth_Index,eps_size
REAL(KIND=dp)::				eta_min,eta_max,tau_min,tau_max,eps_min,eps_max,d_tau,d_eta,d_eps,		&
					sigma,eta,eps,tau,mode,c
REAL(KIND=dp),DIMENSION(2)::		eta_range,tau_range,eps_range

REAL(KIND=dp),DIMENSION(:),ALLOCATABLE::eta_space,tau_space,eps_space,AnalyGamma,AnalyOmega,Gamma_eps,Omega_eps,most_unstable_mode
INTEGER,DIMENSION(:,:),ALLOCATABLE::ModeSpace
REAL(KIND=dp),DIMENSION(:,:),ALLOCATABLE::GammaSpace,OmegaSpace


! NOTE: eta*tau >> 1
  eta_range = (/ 2.0_dp, 2.0_dp /)		! I am looking at the first values...
  tau_range = (/ 1.0_dp, 1.0_dp /)		! ...
  eps_range = (/ 0.02_dp, 0.15_dp /)

! Array size
  eta_size = 2
  tau_size = 2
  eps_size = 1000


! Allocatable arrays
! 1-D 
  ALLOCATE(AnalyGamma(0:MaxOrder),AnalyOmega(0:MaxOrder))
  ALLOCATE(eta_space(1:eta_size),tau_space(1:tau_size))
  ALLOCATE(eps_space(1:eps_size),Gamma_eps(1:eps_size),Omega_eps(1:eps_size),most_unstable_mode(1:eps_size))
! 2-D
  ALLOCATE(ModeSpace(1:eta_size,1:tau_size))
  ALLOCATE(GammaSpace(1:eta_size,1:tau_size)) 
  ALLOCATE(OmegaSpace(1:eta_size,1:tau_size))


! ---------------------------- !
! Creating eta, tau, eps space !
! ---------------------------- !
  eta_min = eta_range(1); eta_max = eta_range(2)
  tau_min = tau_range(1); tau_max = tau_range(2)
  eps_min = eps_range(1); eps_max = eps_range(2)

  d_eps = (eps_max - eps_min) / (eps_size - 1)
  d_eta = (eta_max - eta_min) / (eta_size - 1)
  d_tau = (tau_max - tau_min) / (tau_size - 1)

  DO e = 1,eps_size
    eps_space(e) = eps_min + d_eps*(e-1)
  END DO
  DO i = 1,eta_size
    eta_space(i) = eta_min + d_eta*(i-1)
  END DO
  DO j = 1,tau_size
    tau_space(j) = tau_min + d_tau*(j-1)
  END DO


! --------- !
! Main loop !
! --------- !

  DO e = 1,eps_size

    DO i = 1,eta_size							! Looping over -- eta
      DO j = 1,tau_size							! Looping over -- tau

	DO k = 0,MaxOrder						! Looping over -- Hermite-mode number

	    eta = eta_space(i)
	    tau = tau_space(j)
	    eps = eps_space(e)

	    c = b * tau
	    mode = DBLE( 2*k + 1 )
	    sigma = eps / ( SQRT(c)*q0 )

	    aa = DCMPLX( 1d0 + c , delta )
	    bb = DCMPLX( c*eta - 1d0 , mode*sigma*shear*SQRT(c) )
	    cc = DCMPLX( 0d0 , mode*eta*sigma*shear*SQRT(c) )

	    omega_sub = ( -bb - SQRT(bb**2 - 4d0*aa*cc) ) / ( 2d0*aa )	

	    AnalyGamma(k) = AIMAG(omega_sub)				! Growth-rate for each mode
	    AnalyOmega(k) = REAL(omega_sub)				! Frequency of each mode

	END DO

	MaxGrowth_Index = MAXLOC(AnalyGamma,DIM=1)-1			! Maximum growth-rate for (eta,tau)

	ModeSpace(i,j)  = MaxGrowth_Index		
	GammaSpace(i,j) = AnalyGamma( MaxGrowth_Index )		! corresponding growth rate...
	OmegaSpace(i,j) = AnalyOmega( MaxGrowth_Index )		! ...and mode frequency

      END DO
    END DO

    Gamma_eps(e) = GammaSpace(1,1)
    Omega_eps(e) = OmegaSpace(1,1)
    most_unstable_mode(e) = ModeSpace(1,1)

  END DO

  !PRINT*,'Most unstable mode:',most_unstable_mode(1)
  !PRINT*,'Gamma:',Gamma_eps(1)
  !PRINT*,'Omega:',Omega_eps(1)


  OPEN (5,FILE='gammaeps.txt',ACTION='WRITE')
  WRITE (5,*) Gamma_eps
  CLOSE (5)
  OPEN (5,FILE='omegaeps.txt',ACTION='WRITE')
  WRITE (5,*) Omega_eps
  CLOSE (5)
  OPEN (5,FILE='epsspace.txt',ACTION='WRITE')
  WRITE (5,*) eps_space
  CLOSE (5)
  OPEN (5,FILE='mostunstablemode.txt',ACTION='WRITE')
  WRITE (5,*) most_unstable_mode
  CLOSE (5)

  OPEN (5,FILE='etaspace.txt',ACTION='WRITE')
  WRITE (5,*) eta_space
  CLOSE (5)
  OPEN (5,FILE='tauspace.txt',ACTION='WRITE')
  WRITE (5,*) tau_space
  CLOSE (5)
  OPEN (5,FILE='modespace.txt',ACTION='WRITE')
    DO i=1,eta_size
    WRITE (5,*) (ModeSpace(i,j),j=1,tau_size)
    END DO
  CLOSE (5)
  OPEN (5,FILE='gammaspace.txt',ACTION='WRITE')
    DO i=1,eta_size
    WRITE (5,*) (GammaSpace(i,j),j=1,tau_size)
    END DO
  CLOSE (5)
  OPEN (5,FILE='omegaspace.txt',ACTION='WRITE')
    DO i=1,eta_size
    WRITE (5,*) (OmegaSpace(i,j),j=1,tau_size)
    END DO
  CLOSE (5)


  DEALLOCATE( AnalyGamma,AnalyOmega,eta_space,tau_space )
  DEALLOCATE( ModeSpace,GammaSpace,OmegaSpace )

END PROGRAM FastestGrowingMode
