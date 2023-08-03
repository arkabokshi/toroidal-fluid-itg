module alphainverse_mod
  contains
SUBROUTINE alphainverse(PHIminus,PHI,PHIplus,Gminus,G,Gplus,Hminus,H,Hplus,delPolMode,F)
    USE inputdata,&
        ONLY:length,eta,sigma,xx,dp,epsilonn,shear,c,CURV,Init_gammaE,  &
        FlowShear,idelta_m,low_diag,diagonal,	&
        upp_diag, TaylorFlow, has_flow
    IMPLICIT NONE
    REAL(KIND=dp),INTENT(IN)::delPolMode
    DOUBLE COMPLEX,DIMENSION(length),INTENT(IN)::PHIminus,PHI,PHIplus,Gminus,G,Gplus,Hminus,H,Hplus
    DOUBLE COMPLEX,DIMENSION(length),INTENT(OUT)::F
    DOUBLE COMPLEX,DIMENSION(length)::betaterm,gammaterm,deltaterm,couplingterm,toinvert

    INTEGER::info

    DOUBLE COMPLEX,DIMENSION(length)::cen_diag
    DOUBLE COMPLEX,DIMENSION(length-1)::lft_diag,rgt_diag

! DELTA-TERM
    deltaterm =	O_delta(PHI, delPolMode)
    if (has_flow) then
       deltaterm = deltaterm +       &
        O_gamma(PHI, delPolMode) * FlowShear            +       &
        O_beta(PHI)  * (FlowShear**2)           +       &
        O_alpha(PHI) * (FlowShear**3)
    end if


! GAMMA-TERM
    gammaterm =         O_gamma(G, delPolMode)
    if (has_flow) then
       gammaterm = gammaterm + &
            O_beta(G)  * 2.0_dp * FlowShear		+	&
            O_alpha(G) * 3.0_dp * (FlowShear**2)
    end if


! BETA-TERM
    betaterm =	O_beta(H)
    if(has_flow) then
       betaterm = betaterm + &
        O_alpha(H) * 3.0_dp * FlowShear
    end if

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

  END SUBROUTINE alphainverse

! GAMMA Operator
    pure FUNCTION O_gamma(field, delPolMode)
    USE inputdata,&
        ONLY:length,eta,sigma,xx,dp,epsilonn,shear,c,CURV,Init_gammaE,  &
        FlowShear,idelta_m,low_diag,diagonal,	&
        upp_diag, TaylorFlow
    IMPLICIT none
        DOUBLE COMPLEX,DIMENSION(length) :: O_gamma
        DOUBLE COMPLEX,DIMENSION(length),INTENT(IN)::field
        REAL(KIND=dp),INTENT(IN)::delPolMode

        O_gamma = (sigma**2) * ((xx-delPolMode)**2) * field
    END FUNCTION O_gamma


! DELTA Operator
    pure FUNCTION O_delta(field, delPolMode)
    USE inputdata,&
        ONLY:length,eta,sigma,xx,dp,epsilonn,shear,c,CURV,Init_gammaE,  &
        FlowShear,idelta_m,low_diag,diagonal,	&
        upp_diag, TaylorFlow
    IMPLICIT none
        DOUBLE COMPLEX,DIMENSION(length) :: O_delta
        DOUBLE COMPLEX,DIMENSION(length),INTENT(IN)::field
        REAL(KIND=dp),INTENT(IN)::delPolMode

        O_delta = (sigma**2) * ((xx-delPolMode)**2) * eta* field
    END FUNCTION O_delta


! ALPHA Operator
    FUNCTION O_alpha(field)
    use deriv_mod, only: deriv
    USE inputdata,&
        ONLY:length,eta,sigma,xx,dp,epsilonn,shear,c,CURV,Init_gammaE,  &
        FlowShear,idelta_m,low_diag,diagonal,	&
        upp_diag, TaylorFlow
    IMPLICIT none
        DOUBLE COMPLEX,DIMENSION(length) :: O_alpha
        DOUBLE COMPLEX,DIMENSION(length),INTENT(IN)::field
        DOUBLE COMPLEX,DIMENSION(length)::DField
        CALL Deriv( 2, field,DField,length )
        O_alpha = c*(shear**2)*DField - ( c+1.0_dp+idelta_m+TaylorFlow )*field
    END FUNCTION O_alpha


! BETA Operator
    pure FUNCTION O_beta(field)
    use deriv_mod, only: deriv
    USE inputdata,&
        ONLY:length,eta,sigma,xx,dp,epsilonn,shear,c,CURV,Init_gammaE,  &
        FlowShear,idelta_m,low_diag,diagonal,	&
        upp_diag, TaylorFlow
    IMPLICIT none
        DOUBLE COMPLEX,DIMENSION(length) :: O_beta
        DOUBLE COMPLEX,DIMENSION(length),INTENT(IN)::field
        DOUBLE COMPLEX,DIMENSION(length)::DField
        CALL Deriv( 2, field,DField,length )
        O_beta = eta*c*(shear**2)*DField + ( 1.0_dp - eta*(c+TaylorFlow) )*field
    END FUNCTION O_beta


! COUPLING Operator
    pure FUNCTION O_couple(field_minus,field_plus)
    use deriv_mod, only: deriv
    USE inputdata,&
        ONLY:length,eta,sigma,xx,dp,epsilonn,shear,c,CURV,Init_gammaE,  &
        FlowShear,idelta_m,low_diag,diagonal,	&
        upp_diag, TaylorFlow
    IMPLICIT none
        DOUBLE COMPLEX,DIMENSION(length) :: O_couple
        DOUBLE COMPLEX,DIMENSION(length),INTENT(IN)::field_minus,field_plus
        DOUBLE COMPLEX,DIMENSION(length)::tmp, DField
        tmp = field_plus - field_minus
        CALL Deriv( 1, tmp,DField,length )
        O_couple = ( field_plus+field_minus ) + shear*DField
    END FUNCTION O_couple

end module alphainverse_mod
