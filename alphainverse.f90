SUBROUTINE alphainverse(PHIminus,PHI,PHIplus,Gminus,G,Gplus,Hminus,H,Hplus,delPolMode,F)
USE inputdata,&
ONLY:length,eta,sigma,xx,dp,epsilonn,shear,c,CURV,Init_gammaE, 	&
     inv_tridiag_matrix,FlowShear,idelta_m,low_diag,diagonal,	&
     upp_diag,super_diag, TaylorFlow
IMPLICIT NONE
INCLUDE 'mpif.h'
REAL(KIND=dp),INTENT(IN)::delPolMode
DOUBLE COMPLEX,DIMENSION(length),INTENT(IN)::PHIminus,PHI,PHIplus,Gminus,G,Gplus,Hminus,H,Hplus
DOUBLE COMPLEX,DIMENSION(length),INTENT(OUT)::F
DOUBLE COMPLEX,DIMENSION(length)::betaterm,gammaterm,deltaterm,couplingterm,DField,toinvert
REAL(KIND=dp),PARAMETER::Delxx=ABS(xx(2)-xx(1))
INTEGER::info, ipiv

COMPLEX*16,DIMENSION(length)::cen_diag
COMPLEX*16,DIMENSION(length-1)::lft_diag,rgt_diag

! DELTA-TERM
  deltaterm =	O_delta(PHI) 				+ 	&
		O_gamma(PHI) * FlowShear	 	+ 	&
		O_beta(PHI)  * (FlowShear**2) 		+ 	&
		O_alpha(PHI) * (FlowShear**3)


! GAMMA-TERM
  gammaterm = 	O_gamma(G)				+	&
		O_beta(G)  * 2.0_dp * FlowShear		+	&
		O_alpha(G) * 3.0_dp * (FlowShear**2) 


! BETA-TERM
  betaterm =	O_beta(H)				+	&
		O_alpha(H) * 3.0_dp * FlowShear




! COUPLING TERM
  couplingterm =	O_couple(Hminus,Hplus)					+	&
			O_couple(Gminus,Gplus) * ( eta + 2.0_dp*FlowShear )	+	&
			O_couple(PHIminus,PHIplus) * ( FlowShear**2 + eta*FlowShear)


! TO-INVERT
  toinvert = -(gammaterm + betaterm + deltaterm) + CURV*epsilonn*couplingterm



! INVERTED FIELD Fm
! the inverted tridiagonal matrix is calculated once at the beginning
  toinvert(1)       =   DCMPLX( 0.0_dp,0.0_dp )
  toinvert(length)  =   DCMPLX( 0.0_dp,0.0_dp )
  !F = MATMUL(inv_tridiag_matrix,toinvert)
  
  lft_diag = low_diag
  rgt_diag = upp_diag
  cen_diag = diagonal
  CALL ZGTSV( length,1,lft_diag,cen_diag,rgt_diag,toinvert,length,INFO  )
  F = toinvert

  CONTAINS


! GAMMA Operator
  FUNCTION O_gamma(field)
  DOUBLE COMPLEX,DIMENSION(length) :: O_gamma
  DOUBLE COMPLEX,DIMENSION(length),INTENT(IN)::field
    O_gamma = (sigma**2) * ((xx-delPolMode)**2) * field
  END FUNCTION O_gamma


! DELTA Operator
  FUNCTION O_delta(field)
  DOUBLE COMPLEX,DIMENSION(length) :: O_delta
  DOUBLE COMPLEX,DIMENSION(length),INTENT(IN)::field
    O_delta = (sigma**2) * ((xx-delPolMode)**2) * eta* field
  END FUNCTION O_delta


! ALPHA Operator
  FUNCTION O_alpha(field)
  DOUBLE COMPLEX,DIMENSION(length) :: O_alpha
  DOUBLE COMPLEX,DIMENSION(length),INTENT(IN)::field
    CALL Deriv( 2, field,DField,Delxx,length )
    O_alpha = c*(shear**2)*DField - ( c+1.0_dp+idelta_m+TaylorFlow )*field
  END FUNCTION O_alpha


! BETA Operator
  FUNCTION O_beta(field)
  DOUBLE COMPLEX,DIMENSION(length) :: O_beta
  DOUBLE COMPLEX,DIMENSION(length),INTENT(IN)::field
    CALL Deriv( 2, field,DField,Delxx,length )
    O_beta = eta*c*(shear**2)*DField + ( 1.0_dp - eta*(c+TaylorFlow) )*field
  END FUNCTION O_beta


! COUPLING Operator
  FUNCTION O_couple(field_minus,field_plus)
  DOUBLE COMPLEX,DIMENSION(length) :: O_couple
  DOUBLE COMPLEX,DIMENSION(length),INTENT(IN)::field_minus,field_plus
    CALL Deriv( 1, field_plus-field_minus,DField,Delxx,length )
    O_couple = ( field_plus+field_minus ) + shear*DField
  END FUNCTION O_couple




END SUBROUTINE alphainverse
