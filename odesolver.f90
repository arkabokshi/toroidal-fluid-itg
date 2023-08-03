module odesolver_mod
  use iso_fortran_env, only: real64

  complex(kind = real64), parameter :: zero = 0.0, one = 1.0
  contains
    SUBROUTINE odesolver(nn,dh,ax,bx,cx, low_diag,diagonal,upp_diag)
      use iso_fortran_env, only: real64
    IMPLICIT NONE

! INPUTS
    INTEGER,INTENT(IN)::nn
    REAL(KIND=real64),INTENT(IN)::dh
    DOUBLE COMPLEX,DIMENSION(nn),INTENT(IN) :: ax,bx,cx
    DOUBLE COMPLEX, dimension(nn), intent(out) :: diagonal
    DOUBLE COMPLEX, dimension(nn-1), intent(out) :: low_diag, upp_diag
    INTEGER::i


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constructing A,B,C,D of the form          !
! A(i)f(i-1) + B(i)f(i) + C(i)f(i+1) = D(i) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    upp_diag(1) = zero
    diagonal(1) = one
    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i) &
    !$OMP SHARED(nn, low_diag, diagonal, upp_diag, ax, dh, bx, cx) &
    !$OMP SCHEDULE(static)
    DO i=2,nn-1
        low_diag(i-1) = ( ax(i)/(dh**2) - bx(i)/(2d0*dh) )
        diagonal(i) = ( -2d0*ax(i)/(dh**2) + cx(i) )
        upp_diag(i) = ( ax(i)/(dh**2) + bx(i)/(2d0*dh) )
    END DO
    !$OMP END PARALLEL DO
    diagonal(nn) = one
    low_diag(nn) = zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Boundary conditions                              !
! NOTE: 1. this will be specific to the problem    !
!	2. ALWAYS USE 2nd order b.c. for dy/dx...  !
!	   code extremely sensitive - refer notes  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! ! 1. DIRICHLET - setting boundary to fixed value
!   BB(1) = 1d0; CC(1) = 0d0; DD(1) = 0d0
! ! 2. NEUMANN - setting slope at the boundary to a fixed value
!   AA(nn) = -1d0/dh; BB(nn) = 1d0/dh; DD(nn) = 0d0

END SUBROUTINE odesolver
end module odesolver_mod
