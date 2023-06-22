SUBROUTINE odesolver(nn,dh,ax,bx,cx)
USE inputdata,ONLY:inv_tridiag_matrix,low_diag,diagonal,upp_diag
IMPLICIT NONE
INCLUDE 'mpif.h'
INTEGER,PARAMETER::dp=SELECTED_REAL_KIND(P=15,R=300)

! INPUTS
  INTEGER,INTENT(IN)::nn
  REAL(KIND=dp),INTENT(IN)::dh
  !REAL(KIND=dp),DIMENSION(nn),INTENT(IN) :: ax,bx,cx
  DOUBLE COMPLEX,DIMENSION(nn),INTENT(IN) :: ax,bx,cx

! LOCAL VARIABLES
  !REAL(KIND=dp),DIMENSION(nn,nn) :: tridiag_matrix
  !REAL(KIND=dp),DIMENSION(nn)::AA,BB,CC
  DOUBLE COMPLEX,DIMENSION(nn,nn) :: tridiag_matrix
  DOUBLE COMPLEX,DIMENSION(nn)	  :: AA, BB, CC
  INTEGER::i,m,n


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constructing A,B,C,D of the form          !
! A(i)f(i-1) + B(i)f(i) + C(i)f(i+1) = D(i) !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO i=1,nn
    AA(i) = ( ax(i)/(dh**2) - bx(i)/(2d0*dh) )
    BB(i) = ( -2d0*ax(i)/(dh**2) + cx(i) )
    CC(i) = ( ax(i)/(dh**2) + bx(i)/(2d0*dh) )
!     DD(i) = ( dx(i) ) 
  END DO



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

  BB(1)  = DCMPLX(1d0,1d0);  CC(1)  = DCMPLX(0D0,0d0)   !; DD(1) = 0D0
  BB(nn) = DCMPLX(1d0,1d0);  AA(nn) = DCMPLX(0D0,0d0)   !; DD(nn) = 0D0



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Constructing the tridiagonal matrix ! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    tridiag_matrix = DCMPLX(0d0,0d0)
  ! DIAGONAL TERM
    DO i = 1,nn
      tridiag_matrix(i,i) = BB(i)
	diagonal(i)  = BB(i) 
    END DO
  ! LOW-DIAG TERM
    DO i = 2,nn
      tridiag_matrix(i,i-1) = AA(i)
	low_diag(i-1) = AA(i) 
    END DO
  ! UPP_DIAG TERM
    DO i = 1,nn-1
      tridiag_matrix(i,i+1) = CC(i)
	upp_diag(i) = CC(i)
    END DO


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Inverting the matrix using external library !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !CALL matrixinvert(nn,tridiag_matrix,inv_tridiag_matrix)


!!!!!!!!!!!
! MODULES !
!!!!!!!!!!!

  CONTAINS

  SUBROUTINE matrixinvert(MM,MATRIX,INV_MATRIX)
  IMPLICIT NONE
  EXTERNAL DGETRF,DGETRI

  INTEGER,INTENT(IN)::MM
  !REAL(KIND=dp),DIMENSION(MM,MM),INTENT(IN)::MATRIX
  !REAL(KIND=dp),DIMENSION(MM,MM),INTENT(OUT)::INV_MATRIX
  !REAL(KIND=dp),DIMENSION(MM)::WORK
  DOUBLE COMPLEX,DIMENSION(MM,MM),INTENT(IN)::MATRIX
  DOUBLE COMPLEX,DIMENSION(MM,MM),INTENT(OUT)::INV_MATRIX
  DOUBLE COMPLEX,DIMENSION(MM)::WORK
  
  INTEGER,DIMENSION(MM)::IPIV
  INTEGER::INFO1,INFO2

      INV_MATRIX = MATRIX
      !CALL DGETRF(MM,MM,INV_MATRIX,MM,IPIV,INFO1)
      !CALL DGETRI(MM,INV_MATRIX,MM,IPIV,WORK,MM,INFO2)
      CALL ZGETRF(MM,MM,INV_MATRIX,MM,IPIV,INFO1)
      CALL ZGETRI(MM,INV_MATRIX,MM,IPIV,WORK,MM,INFO2)

  END SUBROUTINE matrixinvert



END SUBROUTINE odesolver


!-------------------------------!
!  DOCUMENTATION		!
!-------------------------------!

!  DGETRI computes the inverse of a matrix using the LU factorization
!  computed by DGETRF.
!
!  This method inverts U and then computes inv(A) by solving the system
!  inv(A)*L = inv(U) for inv(A).
!
!  Arguments
!  =========
!
!  N       (input) INTEGER
!          The order of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the factors L and U from the factorization
!          A = P*L*U as computed by DGETRF.
!          On exit, if INFO = 0, the inverse of the original matrix A.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
!  IPIV    (input) INTEGER array, dimension (N)
!          The pivot indices from DGETRF; for 1<=i<=N, row i of the
!          matrix was interchanged with row IPIV(i).
!
!  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
!
!  LWORK   (input) INTEGER
!          The dimension of the array WORK.  LWORK >= max(1,N).
!          For optimal performance LWORK >= N*NB, where NB is
!          the optimal blocksize returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
!                singular and its inverse could not be computed.
!
!  =====================================================================

!  DGETRF computes an LU factorization of a general M-by-N matrix A
!  using partial pivoting with row interchanges.
!
!  The factorization has the form
!     A = P * L * U
!  where P is a permutation matrix, L is lower triangular with unit
!  diagonal elements (lower trapezoidal if m > n), and U is upper
!  triangular (upper trapezoidal if m < n).
!
!  This is the right-looking Level 3 BLAS version of the algorithm.
!
!  Arguments
!  =========
!
!  M       (input) INTEGER
!          The number of rows of the matrix A.  M >= 0.
!
!  N       (input) INTEGER
!          The number of columns of the matrix A.  N >= 0.
!
!  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix to be factored.
!          On exit, the factors L and U from the factorization
!          A = P*L*U; the unit diagonal elements of L are not stored.
!
!  LDA     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!  IPIV    (output) INTEGER array, dimension (min(M,N))
!          The pivot indices; for 1 <= i <= min(M,N), row i of the
!          matrix was interchanged with row IPIV(i).
!
!  INFO    (output) INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
!                has been completed, but the factor U is exactly
!                singular, and division by zero will occur if it is used
!                to solve a system of equations.
!
!  =====================================================================
