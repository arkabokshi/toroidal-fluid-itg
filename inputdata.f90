MODULE inputdata
IMPLICIT NONE

SAVE
INTEGER::					u
INTEGER,PARAMETER::				initial=-800,final=800,							&	! Total grid points
					   
						length = (final-initial+1),						&
						x0 = (length+1)/2,							&
						dp = SELECTED_REAL_KIND(15,300),					&

						m0 = 70,								&
						n0 = 50,								&
						NumModes_by_2 = 40,							&	! Num. modes on either side of central m0

						InitialMode = m0-NumModes_by_2,						&
						FinalMode   = m0+NumModes_by_2,						&
						NumModes = FinalMode-InitialMode + 1


REAL(KIND=dp),PARAMETER::			MinRad = 0.5_dp,							&	
						MajRad = 5.0_dp,							& 	! Not used explicitly anywhere in the code

						krhoi  =  0.2,								&	! Parameters explained in the paper
						shear  =  25.0_dp,							&	! Bokshi PPCF 2016
						CURV   =  1.0_dp,							&
						etag   =  2.0_dp,							&
						qedge  =  3.45_dp,							&

						etac  =  1062.5_dp,							&
						tau   =  1.0_dp,							&
						epsilonn  = 0.08_dp,							&
						delta_m   = 0.0_dp,							&	! I don't think this is used anywhere...

						b = krhoi**2,								&
						c = b * tau,								&
				
						rs = MinRad * ( DBLE(m0)/(qedge*DBLE(n0)) )**(1/shear),			&
						ktheta = DBLE(m0) / rs,							&
						nqp = ktheta * shear,							&

						dt = 0.1_dp,								&	! Time-step for RK4
						RunTime = 500.1_dp,							&
						GenIso_Transition = 0.0_dp,						&	! Useful for transition studies, eg. Fig 7 in paper
						Init_gammaE = 0.0e-3_dp,						&
						Final_gammaE =0.0e-3_dp,						&
						Quad_gammaE  = 0.0e-4_dp							! Effectively k1/n term in Fig. 6

INTEGER,PARAMETER::				NumSteps = INT(RunTime/dt),						&
						PtDensity = 10,								&	! Determines domain in radial coordinate xx
						DelPrint  = 100,							&	! How frequently to output to file
						NumTheta = 360									! Points in poloidal direction
					
REAL(KIND=dp)::					FlowOnOff = 1.0_dp,							&
						ShearingRate,								&
						dGammaE_dt = 1.0e-1_dp								! If dynamically adjusting flow-shear, uncomment
																! line 134 onwards in itg.f90 and make appropriate 																	! minor changes to the file
REAL(KIND=dp),DIMENSION(length),PARAMETER::	xx = (/ (u,u=initial,final,1) /) / DBLE(PtDensity)


REAL(KIND=dp),DIMENSION(length),PARAMETER::	qprofile = DBLE(m0)/DBLE(n0),						&
 						!qprofile = qedge * ( (xx/nqp + rs)/MinRad )**shear,			&	! Can give a q-profile
						!eta = etag * ( 1.0_dp - etac * (xx/nqp)**2 ),				&	! can give a drive-profile
						sigma = epsilonn / ( SQRT(b*tau)*qprofile ),				&
						TaylorFlow = 0.0_dp

REAL(KIND=dp),DIMENSION(length)::		FlowShear = (Init_gammaE*xx) + (Quad_gammaE*(xx**2)),			&
						eta										! Defined in itg.f90 - as are other parameters
																! not defined here
DOUBLE COMPLEX,PARAMETER::			ci = DCMPLX( 0.0_dp,1.0_dp ),						&
						idelta_m = DCMPLX( 0.0_dp,delta_m )

INTEGER::					MyRank,MySize,ierror,ModesPerProc
DOUBLE COMPLEX,DIMENSION(length,length)::	inv_tridiag_matrix

COMPLEX*16,DIMENSION(length) ::			diagonal
COMPLEX*16,DIMENSION(length-1)::		upp_diag,low_diag,super_diag
COMPLEX*16,DIMENSION(length,NumModes)::		FULLtoinvert


LOGICAL,PARAMETER::				restart  = .FALSE.,							&	! If restarting simulation, define location of 
						with_MPI = .TRUE.,                                                      &! FinalFields.txt from previous run
                                                noiseStart = .TRUE.
						calcThetaMaxima = .false.

END MODULE inputdata


!------!
! NOTE !
!------!

! (1) Running 'make' should compile and run the code. Check 'makefile' for the softwares needed. 
!     f90 scripts need MPI, LAPACK and BLAS
! (2) All python scripts are for post-processing. Analysis scripts can be run before simulation has finished.








