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


    k1 = df ( u0 )
    k2 = df ( u0 + dt/2.0_dp*k1 )
    k3 = df ( u0 + dt/2.0_dp*k2 )
    k4 = df ( u0 + dt*k3 )

    u1 = u0 + dt/6.0_dp*( k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4 )

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
        INTEGER::ModeNumber,ll


        ! ------------------------------ !
        ! Looping over the various modes !
        ! ------------------------------ !

        !$OMP PARALLEL DO DEFAULT(none) &
        !$OMP PRIVATE(ModeNumber, delPolMode, phi, g, h, phiplus, gplus, &
        !$OMP hplus, phiminus, gminus, hminus, F) &
        !$OMP SHARED(ModesPerProc, MyRank, df, statevector) &
        !$OMP SCHEDULE(static)
        DO ModeNumber = 1+ModesPerProc*MyRank,(MyRank+1)*ModesPerProc


            delPolMode = DBLE( InitialMode-m0-1+ModeNumber )

            PHI  =  statevector( 0*length+1:1*length , ModeNumber )
            G    =  statevector( 1*length+1:2*length , ModeNumber )
            H    =  statevector( 2*length+1:3*length , ModeNumber )


            ! ------------------ !
            ! Setting neighbours !
            ! ------------------ !

            IF (ModeNumber.EQ.1) THEN

                PHIplus  =   statevector( 0*length+1:1*length , ModeNumber+1 )
                Gplus    =   statevector( 1*length+1:2*length , ModeNumber+1 )
                Hplus    =   statevector( 2*length+1:3*length , ModeNumber+1 )
                PHIminus =   0.0_dp
                Gminus   =   0.0_dp
                Hminus   =   0.0_dp

            ELSE IF (ModeNumber.EQ.NumModes) THEN

                PHIplus  =   0.0_dp
                Gplus    =   0.0_dp
                Hplus    =   0.0_dp
                PHIminus =   statevector( 0*length+1:1*length , NumModes-1 )
                Gminus   =   statevector( 1*length+1:2*length , NumModes-1 )
                Hminus   =   statevector( 2*length+1:3*length , NumModes-1 )

            ELSE

                PHIplus  =   statevector( 0*length+1:1*length , ModeNumber+1 )
                PHIminus =   statevector( 0*length+1:1*length , ModeNumber-1 )
                Gplus    =   statevector( 1*length+1:2*length , ModeNumber+1 )
                Gminus   =   statevector( 1*length+1:2*length , ModeNumber-1 )
                Hplus    =   statevector( 2*length+1:3*length , ModeNumber+1 )
                Hminus   =   statevector( 2*length+1:3*length , ModeNumber-1 )

            END IF


            ! --------------------------- !
            ! Constructing F by inversion !
            ! and statevector to evolve   !
            ! --------------------------- !

            CALL alphainverse(PHIminus,PHI,PHIplus,Gminus,G,Gplus,Hminus,H,Hplus,delPolMode,F)

!            ll = ll+1
            !            TempData( :, ll ) = -ci * [ G, H, F ]
            df( :, ModeNumber ) = -ci * [ G, H, F ]


        END DO
        !$OMP END PARALLEL DO
    END FUNCTION df


end module evolve_mod
