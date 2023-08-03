module deriv_mod
contains

  pure SUBROUTINE Deriv( Order,InArray,OutArray,ArrayLength )
    use inputdata, only: xx
    IMPLICIT NONE
    INTEGER,PARAMETER::dp=SELECTED_REAL_KIND(15,300)
    INTEGER,INTENT(IN)::Order,ArrayLength
    REAL(KIND=dp),PARAMETER::Delh=ABS(xx(2)-xx(1))
    DOUBLE COMPLEX,DIMENSION(ArrayLength),INTENT(IN)::InArray
    DOUBLE COMPLEX,DIMENSION(ArrayLength),INTENT(OUT)::OutArray
    INTEGER::ii


    IF (Order.EQ.2) THEN
!!!!!!!!!!!!!!!!!!!!!
! SECOND DERIVATIVE !
!!!!!!!!!!!!!!!!!!!!!

        ii = 1
        OutArray(ii) = ( InArray(ii+2) - 2.0_dp*InArray(ii+1) + InArray(ii) ) / ( Delh**2 )
        DO ii=2,ArrayLength-1
            OutArray(ii) = ( InArray(ii+1) - 2.0_dp*InArray(ii) + InArray(ii-1) ) / ( Delh**2 )
        END DO

        ii = ArrayLength
        OutArray(ii) = ( InArray(ii) - 2.0_dp*InArray(ii-1) + InArray(ii-2) ) / ( Delh**2 )


    ELSE IF (Order.EQ.1) THEN
!!!!!!!!!!!!!!!!!!!!
! FIRST DERIVATIVE !
!!!!!!!!!!!!!!!!!!!!

        ii = 1
        OutArray(ii) = ( InArray(ii+1)-InArray(ii) ) / ( Delh )
        DO ii=2,ArrayLength-1
            OutArray(ii) = ( InArray(ii+1) - InArray(ii-1) ) / ( 2.0_dp*Delh )
        END DO

        ii = ArrayLength
        OutArray(ii) = ( InArray(ii)-InArray(ii-1) ) / ( Delh )


    ELSE

!        PRINT*,'Can only compute first or second derivatives'
 !       CALL EXIT()
       error stop 'Can only compute first or second derivatives'
    END IF



END SUBROUTINE Deriv
end module deriv_mod
