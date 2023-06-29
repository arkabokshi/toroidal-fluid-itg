SUBROUTINE plasmaprofiles()
USE inputdata,&
ONLY: eta0,ShearingRate,epsilonn,c,qedge,shear,MinRad,dp,&
      xx,eta,qprofile,sigma,FlowShear,&
      InitialMode,FinalMode,NumModes,m0,n0,length
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER::ii,mm
REAL(KIND=dp),DIMENSION(NumModes)::qq_temp,xx_temp,workarray
REAL(KIND=dp),DIMENSION(length)::d_qprofile
REAL(KIND=dp)::Delta_x,xx_start,nqp,qp

! CHECK ONCE FOR ERRORS

      Delta_x = 1.0_dp
      xx_start = InitialMode-m0

      DO ii = 0,NumModes-1
	  mm = InitialMode + ii
	  qq_temp(ii+1) = DBLE(mm)/DBLE(n0)
	  xx_temp(ii+1) = xx_start + ii*Delta_x
      END DO

      ! qprofile = DBLE(m0)/DBLE(n0)
      ! qprofile = qedge * ( (xx/nqp + rs)/MinRad )**shear				! rs commented out in inputdata.f90
      CALL spline1d( qprofile,xx,length,qq_temp,xx_temp,NumModes,workarray )		! shear must be = 1.0


!       eta = eta0 * ( 1.0_dp - 75.0_dp*(xx/nqp)**2 )					! NEED to change to variables in inputdata.f90
!       sigma = epsilonn / ( SQRT(c)*qprofile )
!       FlowShear = ShearingRate * xx
	  


END SUBROUTINE plasmaprofiles