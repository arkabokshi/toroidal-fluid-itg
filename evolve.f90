module evolve_mod
  contains
SUBROUTINE Evolve(u0,u1)
    USE inputdata,&
        ONLY:length,dp,dt,NumModes,InitialMode,m0,xx,FlowShear,ci,FlowOnOff,	&
        MyRank,MySize,ierror,ModesPerProc
    IMPLICIT NONE

    DOUBLE COMPLEX,DIMENSION(3*length,NumModes), INTENT(IN):: u0
    DOUBLE COMPLEX,DIMENSION(3*length,NumModes),INTENT(OUT):: u1
    DOUBLE COMPLEX,DIMENSION(3*length,NumModes):: k1,k2,k3,k4
    integer :: i

    k1 = df ( u0 )
    k2 = df ( u0 + dt/2.0_dp*k1 )
    k3 = df ( u0 + dt/2.0_dp*k2 )
    k4 = df ( u0 + dt*k3 )

    !$OMP PARALLEL DO DEFAULT(none) &
    !$OMP PRIVATE(i) &
    !$OMP SHARED(u1, u0, k1, k2, k3, k4) &
    !$OMP SCHEDULE(static)
    do i = 1, NumModes
       u1(:, i) = u0(:, i) + dt/6.0_dp*( k1(:, i) + &
            2.0_dp*(k2(:,i) + k3(:,i)) + k4(:,i) )
    end do
    !$OMP END PARALLEL DO
  end SUBROUTINE Evolve


    ! ------------------------------------ !
    ! Determining time-derivatives for RK4 !
    ! ------------------------------------ !

    FUNCTION df(statevector)
    use alphainverse_mod, only: alphainverse
    USE inputdata,&
        ONLY:length,dp,dt,NumModes,InitialMode,m0,xx,FlowShear,ci,FlowOnOff,	&
        MyRank,MySize,ierror,ModesPerProc
        IMPLICIT NONE
        DOUBLE COMPLEX,DIMENSION(3*length,NumModes)::df
        DOUBLE COMPLEX,DIMENSION(3*length,NumModes),INTENT(IN)::statevector
        DOUBLE COMPLEX,DIMENSION(length)::PHI,G,H,F,PHIminus,PHIplus,Gplus,Gminus,Hplus,Hminus
        REAL(KIND=dp)::delPolMode
        INTEGER::ModeNumber
        integer :: phi_start, phi_end, g_start, g_end, h_start, h_end
        phi_start = 1 ; phi_end = phi_start + length - 1
        g_start = phi_end + 1 ; g_end = g_start + length - 1
        h_start = g_end + 1 ; h_end = h_start + length - 1

        ! ------------------------------ !
        ! Looping over the various modes !
        ! ------------------------------ !

        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(ModeNumber, delPolMode, phi, g, h, phiplus, gplus, &
        !$OMP hplus, phiminus, gminus, hminus, F) &
        !$OMP SHARED(df, statevector, phi_start, phi_end, g_start, g_end, &
        !$OMP h_start, h_end) &
        !$OMP SCHEDULE(static)
        DO ModeNumber = 1, NumModes


            delPolMode = DBLE( InitialMode-m0-1+ModeNumber )

            PHI  =  statevector( phi_start:phi_end , ModeNumber )
            G    =  statevector( g_start:g_end , ModeNumber )
            H    =  statevector( h_start:h_end , ModeNumber )


            ! ------------------ !
            ! Setting neighbours !
            ! ------------------ !

            IF (ModeNumber.EQ.1) THEN

                PHIplus  =   statevector( phi_start:phi_end , ModeNumber+1 )
                Gplus    =   statevector( g_start:g_end , ModeNumber+1 )
                Hplus    =   statevector( h_start:h_end , ModeNumber+1 )
                PHIminus =   0.0_dp
                Gminus   =   0.0_dp
                Hminus   =   0.0_dp

            ELSE IF (ModeNumber.EQ.NumModes) THEN

                PHIplus  =   0.0_dp
                Gplus    =   0.0_dp
                Hplus    =   0.0_dp
                PHIminus =   statevector( phi_start:phi_end , NumModes-1 )
                Gminus   =   statevector( g_start:g_end , NumModes-1 )
                Hminus   =   statevector( h_start:h_end , NumModes-1 )

            ELSE

                PHIplus  =   statevector( phi_start:phi_end , ModeNumber+1 )
                PHIminus =   statevector( phi_start:phi_end , ModeNumber-1 )
                Gplus    =   statevector( g_start:g_end , ModeNumber+1 )
                Gminus   =   statevector( g_start:g_end , ModeNumber-1 )
                Hplus    =   statevector( h_start:h_end , ModeNumber+1 )
                Hminus   =   statevector( h_start:h_end , ModeNumber-1 )

            END IF

            ! --------------------------- !
            ! Constructing F by inversion !
            ! and statevector to evolve   !
            ! --------------------------- !

            CALL alphainverse(PHIminus,PHI,PHIplus,Gminus,G,Gplus,Hminus,H,Hplus,delPolMode,F)

            df( phi_start:phi_end, ModeNumber ) = -ci * G
            df( g_start:g_end, ModeNumber) = -ci * H
            df( h_start:h_end, ModeNumber) = -ci * F

        END DO
        !$OMP END PARALLEL DO
    END FUNCTION df


end module evolve_mod
