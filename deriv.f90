SUBROUTINE Deriv( Order,InArray,OutArray,Delh,ArrayLength )
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    INTEGER,PARAMETER::dp=SELECTED_REAL_KIND(15,300)
    INTEGER,INTENT(IN)::Order,ArrayLength
    REAL(KIND=dp),INTENT(IN)::Delh
    DOUBLE COMPLEX,DIMENSION(ArrayLength),INTENT(IN)::InArray
    DOUBLE COMPLEX,DIMENSION(ArrayLength),INTENT(OUT)::OutArray
    INTEGER::ii


    IF (Order.EQ.2) THEN
!!!!!!!!!!!!!!!!!!!!!
! SECOND DERIVATIVE !
!!!!!!!!!!!!!!!!!!!!!

        ii = 1
        OutArray(ii) = ( InArray(ii+2) - 2.0_dp*InArray(ii+1) + InArray(ii) ) / ( Delh**2 )
        ii = ArrayLength
        OutArray(ii) = ( InArray(ii) - 2.0_dp*InArray(ii-1) + InArray(ii-2) ) / ( Delh**2 )

        DO ii=2,ArrayLength-1
            OutArray(ii) = ( InArray(ii+1) - 2.0_dp*InArray(ii) + InArray(ii-1) ) / ( Delh**2 )
        END DO


    ELSE IF (Order.EQ.1) THEN
!!!!!!!!!!!!!!!!!!!!
! FIRST DERIVATIVE !
!!!!!!!!!!!!!!!!!!!!

        ii = 1
        OutArray(ii) = ( InArray(ii+1)-InArray(ii) ) / ( Delh )
        ii = ArrayLength
        OutArray(ii) = ( InArray(ii)-InArray(ii-1) ) / ( Delh )

        DO ii=2,ArrayLength-1
            OutArray(ii) = ( InArray(ii+1) - InArray(ii-1) ) / ( 2.0_dp*Delh )
        END DO


    ELSE

        PRINT*,'Can only compute first or second derivatives'
        CALL EXIT()

    END IF



END SUBROUTINE Deriv
