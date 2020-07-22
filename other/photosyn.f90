

!**********************************************************************

SUBROUTINE PSILTUZET(PAR,TLEAF,CS,RH,VPD,PATM, &
        JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN, &
        THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW, &
        MODELGS,EMAXLEAF,KTOT,WEIGHTEDSWP,PSILIN,PSIL,GSMIN, &
	    G0,D0L,GAMMA,G1,GK,TUZFUN,SF,PSIV,HMSHAPE, &
        FSOIL,GS,ALEAF,ALEAFHM,RD)

   IMPLICIT NONE
   
   DOUBLE PRECISION PATM,FSOIL,ETEST
   DOUBLE PRECISION PAR,TLEAF,CS,RH,PSIL,KTOT
   DOUBLE PRECISION JMAX25,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC
   DOUBLE PRECISION THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW,TVJUP,TVJDN
   DOUBLE PRECISION G0,D0L,GAMMA,G1,GK,GMESO
   DOUBLE PRECISION GS,ALEAF,RD,MINLEAFWP,RD0ACC
   DOUBLE PRECISION GAMMASTAR,KM,JMAX,VCMAX,J,VJ
   DOUBLE PRECISION A,B,C,AC,AJ,GSDIVA,CIC,CIJ
   DOUBLE PRECISION VPD,VPDG,ALEAFHM,HMSHAPE
   DOUBLE PRECISION EMAXLEAF,WEIGHTEDSWP,GSV,GSMIN
   DOUBLE PRECISION GSVGSC,RESPLIGHTFRAC,SF,PSIV,PSILIN
   DOUBLE PRECISION PSILCUR,MINLWP,LWPDELTA,LASTDIFF,CI
   INTEGER MODELGS,IQERROR,IECO,I,TUZFUN
   LOGICAL SOLVED

   PSILCUR = WEIGHTEDSWP
   MINLWP = -20
   LWPDELTA = 0.01
   SOLVED = .FALSE.
   LASTDIFF = 1E09
   CI = -999

    ! ITERATE IN STEPS OF 0.01, STARTING AT PSIS DOWNWARD
    ! DO WHILE (PSILCUR .GT. MINLWP .AND. SOLVED .EQV. .FALSE.)
        
        ! CALL PHOTOSYN(PAR,TLEAF,CS,RH,VPD,PATM, &
        ! JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN, &
        ! THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW, &
        ! MODELGS,EMAXLEAF,KTOT,WEIGHTEDSWP,PSILCUR,PSIL,GSMIN, &
	    ! G0,D0L,GAMMA,G1,GK,TUZFUN,SF,PSIV,HMSHAPE, &
        ! FSOIL,GS,ALEAF,ALEAFHM,RD,CI,)
        
        ! If current
        ! IF (ABS(PSILCUR - PSIL).GT.LASTDIFF) THEN
        !     SOLVED = .TRUE.
        !     PSILCUR = PSILCUR + LWPDELTA  ! Value for previous iteration
        ! ELSE
        !     PSILCUR = PSILCUR - LWPDELTA  ! Reduce PSIL, try again.
        !     LASTDIFF = PSILCUR - PSIL
        ! ENDIF
    
    ! ENDDO
    
    ! Now find a finer resolution solution.
    !LWPDELTA = 0.001
    !SOLVED = .FALSE.
    !PSILCUR = PSILCUR + LWPDELTA ! extra step above
    !LASTDIFF = 1E09
    !DO WHILE (PSILCUR .GT. MINLWP .AND. SOLVED .EQV. .FALSE.)
        
     !   CALL PHOTOSYN(PAR,TLEAF,CS,RH,VPD,PATM, &
     !   JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN, &
     !   THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW, &
     !   MODELGS,EMAXLEAF,KTOT,WEIGHTEDSWP,PSILCUR,PSIL,GSMIN, &
	 !   G0,D0L,GAMMA,G1,GK,SF,PSIV,HMSHAPE, &
     !   FSOIL,GS,ALEAF,ALEAFHM,RD)
        
        ! If current
      !  IF (ABS(PSILCUR - PSIL).GT.LASTDIFF) THEN
      !      SOLVED = .TRUE.
      !      PSILCUR = PSILCUR + LWPDELTA  ! Value for previous iteration
     !   ELSE
     !       PSILCUR = PSILCUR - LWPDELTA  ! Reduce PSIL, try again.
     !       LASTDIFF = PSILCUR - PSIL
     !   ENDIF
    
    !ENDDO
    
    ! This subroutine should just be called for LWP calculation,
    ! because all other parameters are calculated for the next step LWP,
    ! and are therefore a bit different from the best solution.
    PSILIN = PSILCUR
    PSIL = PSILCUR
    
    
END SUBROUTINE PSILTUZET



!**********************************************************************

DOUBLE PRECISION FUNCTION FPSIL(PSIL,SF,PSIV,FUNTYPE)

IMPLICIT NONE
DOUBLE PRECISION PSIL,SF,PSIV
INTEGER FUNTYPE

IF(FUNTYPE.EQ.1)THEN
	FPSIL = (1 + EXP(SF*PSIV) )/ (1 + EXP(SF*(PSIV-PSIL)))
ENDIF

IF(FUNTYPE.EQ.2)THEN
	FPSIL = 1 - PSIL/SF
ENDIF

END FUNCTION FPSIL



!**********************************************************************
SUBROUTINE PHOTOSYN(PAR,TLEAF,CS,RH,VPD,PATM, &
        JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN, &
        THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW, &
        MODELGS,EMAXLEAF,KTOT,WEIGHTEDSWP,PSILIN,PSIL,GSMIN, &
	    G0,D0L,GAMMA,G1,GK,TUZFUN,SF,PSIV,HMSHAPE, &
        FSOIL,GS,ALEAF,ALEAFHM,RD,CI,GSDIVA,AJ,AC,VJ,VCMAX,JMAX,KM,GAMMASTAR)
        
! This subroutine calculates photosynthesis according to the ECOCRAFT
! agreed formulation of the Farquhar & von Caemmerer (1982) equations.
! Stomatal conductance may be calculated according to the Jarvis,
! Ball-Berry or BB-Leuning models.
! NB ALEAF is NET leaf photosynthesis.
! FSOIL is now output (RAD 2008).
!**********************************************************************

    IMPLICIT NONE

    DOUBLE PRECISION PATM,FSOIL,ETEST
    DOUBLE PRECISION PAR,TLEAF,CS,RH,PSIL,KTOT
    DOUBLE PRECISION JMAX25,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC
    DOUBLE PRECISION THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW,TVJUP,TVJDN
    DOUBLE PRECISION G0,D0L,GAMMA,G1,GK
    DOUBLE PRECISION GS,ALEAF,RD,MINLEAFWP,RD0ACC
    DOUBLE PRECISION GAMMASTAR,KM,JMAX,VCMAX,J,VJ
    DOUBLE PRECISION A,B,C,AC,AJ,GSDIVA,CIC,CIJ
    DOUBLE PRECISION VPD,VPDG,ALEAFHM,HMSHAPE,GMESO
    DOUBLE PRECISION EMAXLEAF,WEIGHTEDSWP,GSV,GSMIN
    DOUBLE PRECISION GSVGSC,RESPLIGHTFRAC,SF,PSIV,PSILIN,CI
    DOUBLE PRECISION, EXTERNAL :: GAMMAFN
    DOUBLE PRECISION, EXTERNAL :: KMFN
    DOUBLE PRECISION, EXTERNAL :: JMAXTFN
    DOUBLE PRECISION, EXTERNAL :: VCMAXTFN
    DOUBLE PRECISION, EXTERNAL :: RESP
    DOUBLE PRECISION, EXTERNAL :: QUADM
    DOUBLE PRECISION, EXTERNAL :: QUADP
    DOUBLE PRECISION, EXTERNAL :: FPSIL
        
    INTEGER*4 MODELGS,IQERROR,IECO,TUZFUN

	IF(G0.LT.1E-09)G0 = 1E-09
	
	
    ! Calculate photosynthetic parameters from leaf temperature.
    GAMMASTAR = GAMMAFN(TLEAF,IECO)                   ! CO2 compensati
    KM = KMFN(TLEAF,IECO)                             ! Michaelis-Ment
    JMAX = JMAXTFN(JMAX25,TLEAF,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN)      ! Po
    VCMAX = VCMAXTFN(VCMAX25,TLEAF,EAVC,EDVC,DELSC,TVJUP,TVJDN)   ! Ma
	
    ! Make sure light suppression of dark respiration only occurs when it is light.
    ! See Atkin et al. 1998 (or 2001?)
    IF(PAR.LT.100.0)THEN
	    RESPLIGHTFRAC = 1.0
    ELSE
	    RESPLIGHTFRAC = DAYRESP
    ENDIF 

    RD = RESP(RD0,RD0ACC,TLEAF,Q10F,RTEMP,RESPLIGHTFRAC,TBELOW)
    J = QUADM(THETA,-(AJQ*PAR+JMAX),AJQ*PAR*JMAX,IQERROR) ! Actual e-
    ! VJ = J/DBLE(4.0)                                      ! RuBP-regen rat

    ! CONSTANTS
    GSVGSC = DBLE(1.57)
	
    ! Deal with extreme cases
    IF ((JMAX.LE.0.0).OR.(VCMAX.LE.0.0)) THEN
        ALEAF = -RD
        GS = G0
        RETURN
    END IF

        IF (MODELGS.EQ.2) THEN
            ! Ball-Berry model
            GSDIVA = G1 * RH / (CS - GAMMA) * FSOIL
        ELSE IF (MODELGS.EQ.3) THEN
            ! Ball-Berry-Leuning model
            GSDIVA = G1 / (CS - GAMMA) / (1 + VPD/D0L) * FSOIL
        ELSE IF (MODELGS.EQ.4) THEN
            IF(VPD.LT.50)THEN
                VPDG = 50.0/1000.0
            ELSE
                VPDG =VPD/1000.0
            ENDIF
            GSDIVA = G1 / (CS - GAMMA) / SQRT(VPDG)
        ELSE IF (MODELGS.EQ.5) THEN
            IF(VPD.LT.50)THEN
                VPDG = 50.0/1000.0
            ELSE
                VPDG =VPD/1000.0
            ENDIF
            GSDIVA = G1 / (CS - GAMMA) / VPDG**GK    
        ELSE IF (MODELGS.EQ.6) THEN
            IF(VPD.LT.50)THEN
                VPDG = 50.0/1000.0
            ELSE
                VPDG =VPD/1000.0
            ENDIF
            GSDIVA = (G1 / (CS - GAMMA)) * FPSIL(PSILIN,SF,PSIV,TUZFUN)    
        END IF

        ! Following calculations are used for both BB & BBL models.
        ! Solution when Rubisco activity is limiting
        A = G0 + GSDIVA * (VCMAX - RD)
        B = (1. - CS*GSDIVA) * (VCMAX - RD) + G0 * (KM - CS)- GSDIVA * (VCMAX*GAMMASTAR + KM*RD)
        C = -(1. - CS*GSDIVA) * (VCMAX*GAMMASTAR + KM*RD) - G0*KM*CS

        CIC = QUADP(A,B,C,IQERROR)

        IF ((IQERROR.EQ.1).OR.(CIC.LE.0.0).OR.(CIC.GT.CS)) THEN
            AC = 0.0
        ELSE
            AC = VCMAX * (CIC - GAMMASTAR) / (CIC + KM)
        END IF

        ! Solution when electron transport rate is limiting
        A = G0 + GSDIVA * (VJ - RD)
        B = (1. - CS*GSDIVA) * (VJ - RD) + G0 * (2.*GAMMASTAR - CS) &
            - GSDIVA * (VJ*GAMMASTAR + 2.*GAMMASTAR*RD)
        C = -(1. - CS*GSDIVA) * GAMMASTAR * (VJ + 2.*RD) &
            - G0*2.*GAMMASTAR*CS
        CIJ = QUADP(A,B,C,IQERROR)

        AJ = VJ * (CIJ - GAMMASTAR) / (CIJ + 2.*GAMMASTAR)
        IF (AJ-RD.LT.1E-6) THEN        ! Below light compensation point
            CIJ = CS
            AJ = VJ * (CIJ - GAMMASTAR) / (CIJ + 2.*GAMMASTAR)
        END IF

        ALEAF = DMIN1(AC,AJ) - RD  
        
        ! Hyperbolic minimum; use for optimization to avoid discontinuity.
        ! HMSHAPE determines how smooth the transition is.
        ALEAFHM = (AC+AJ - sqrt((AC+AJ)**2-4*HMSHAPE*AC*AJ))/(2*HMSHAPE) - RD
        
        ! Solution for Ball-Berry model
        IF(HMSHAPE.LT.1) THEN
        GS = G0 + GSDIVA*ALEAFHM
        ELSE
        GS = G0 + GSDIVA*ALEAF
        ENDIF
        
        IF (GS.LT.G0) GS = G0


        ! If E > Emax, set E to Emax, and corresponding gs and A.
        ! Leaf transpiration in mmol m-2s-1.
        ETEST = 1000 * (VPD/PATM) * GS * GSVGSC

        ! Leaf water potential
        PSIL = WEIGHTEDSWP - ETEST/KTOT

		! Return CI.
        IF(GS.GT.0.AND.ALEAF.GT.0)THEN
            CI = CS - ALEAF/GS
        ELSE
            CI = CS
        ENDIF	
			
			
        IF(MODELGS.LT.6)THEN
            
            IF(ETEST > EMAXLEAF)THEN

                ! Just for output:
                FSOIL = EMAXLEAF / ETEST

                ! Gs in mol m-2 s-1
                GSV = 1E-03 * EMAXLEAF / (VPD/PATM)
               
				if(gsv.lt.g0)then
					gsv = g0
				endif
				
				GS = GSV / GSVGSC
				
				
                ! Minimum leaf water potential reached
                PSIL = WEIGHTEDSWP - EMAXLEAF/KTOT

                !IF(GS.LT.1E-06) THEN
                !    GS = 0.0
                !    ALEAF = 0.0
                !ELSE

				
                    ! Now that GS is known, solve for CI and A as in the Jar
                    ! Photosynthesis when Rubisco is limiting
					A = 1./GS
                    B = (RD - VCMAX)/GS - CS - KM
                    C = VCMAX * (CS - GAMMASTAR) - RD * (CS + KM)

                    AC = QUADM(A,B,C,IQERROR)

                    !IF (IQERROR.EQ.1) THEN
                    !    GS = GSMIN
                    !    AC = - RD
                    !END IF

                    ! Photosynthesis when electron transport is limiting
                    B = (RD - VJ)/GS - CS - 2*GAMMASTAR
                    C = VJ * (CS - GAMMASTAR) - RD * (CS + 2*GAMMASTAR)
                    AJ = QUADM(A,B,C,IQERROR)

                    !IF (IQERROR.EQ.1) THEN
                    !    GS = GSMIN
                    !    AJ = - RD
                    !END IF

                    ALEAF = DMIN1(AC,AJ)  ! - RD       ! Emax model solution.

					! For droughted situation, calculate CI with gross photosynthesis,
					! to avoid strange behavior at very low GS:
					! Return CI.
					IF(GS.GT.0.AND.ALEAF.GT.0)THEN
						CI = CS - ALEAF/GS
					ELSE
						CI = CS
					ENDIF	
										
                !END IF ! if gs < 1e-05
            END IF ! if (E>EMAX)
        END IF ! if(modelgs.lt.6)
        

        
        
    RETURN
END SUBROUTINE PHOTOSYN

!**********************************************************************
SUBROUTINE ACI(CI, PAR,TLEAF,PATM, &
        JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN, &
        THETA,AJQ,GMESO,RD0,Q10F,RTEMP,DAYRESP,TBELOW,HMSHAPE, &
        AC,AJ,ALEAF,ALEAFHM,RD,GAMMASTAR,VJ,KM,VCMAX)

! This subroutine calculates photosynthesis according to the ECOCRAFT
! agreed formulation of the Farquhar & von Caemmerer (1982) equations.
! Stomatal conductance may be calculated according to the Jarvis,
! Ball-Berry or BB-Leuning models.
! NB ALEAF is NET leaf photosynthesis.
! FSOIL is now output (RAD 2008).
!**********************************************************************

    IMPLICIT NONE

    DOUBLE PRECISION PATM,FSOIL,ETEST,CI
    DOUBLE PRECISION PAR,TLEAF,CS,RH,PSIL,KTOT
    DOUBLE PRECISION JMAX25,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC
    DOUBLE PRECISION THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW,TVJUP,TVJDN
    DOUBLE PRECISION G0,D0L,GAMMA,G1,GK,GMESO
    DOUBLE PRECISION GS,ALEAF,RD,MINLEAFWP,RD0ACC
    DOUBLE PRECISION GAMMASTAR,KM,JMAX,VCMAX,J,VJ
    DOUBLE PRECISION A,B,C,AC,AJ,GSDIVA,CIC,CIJ
    DOUBLE PRECISION VPD,VPDG,HMSHAPE,ALEAFHM,AP,BP
    DOUBLE PRECISION EMAXLEAF,WEIGHTEDSWP,GSV,GSMIN
    DOUBLE PRECISION GSVGSC,RESPLIGHTFRAC,SF,PSIV,PSILIN
    DOUBLE PRECISION, EXTERNAL :: GAMMAFN
    DOUBLE PRECISION, EXTERNAL :: KMFN
    DOUBLE PRECISION, EXTERNAL :: JMAXTFN
    DOUBLE PRECISION, EXTERNAL :: VCMAXTFN
    DOUBLE PRECISION, EXTERNAL :: RESP
    DOUBLE PRECISION, EXTERNAL :: QUADM
    DOUBLE PRECISION, EXTERNAL :: QUADP
    DOUBLE PRECISION, EXTERNAL :: FPSIL
        
    INTEGER MODELGS,IQERROR,IECO

    ! Calculate photosynthetic parameters from leaf temperature.
    GAMMASTAR = GAMMAFN(TLEAF,IECO)                   ! CO2 compensati
    KM = KMFN(TLEAF,IECO)                             ! Michaelis-Ment
    JMAX = JMAXTFN(JMAX25,TLEAF,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN)      ! Po
    VCMAX = VCMAXTFN(VCMAX25,TLEAF,EAVC,EDVC,DELSC,TVJUP,TVJDN)   ! Ma
	
    ! Make sure light suppression of dark respiration only occurs when it is light.
    ! See Atkin et al. 1998.
    IF(PAR.LT.100.0)THEN
	     RESPLIGHTFRAC = 1.0
    ELSE
	     RESPLIGHTFRAC = DAYRESP
    ENDIF 

    RD = RESP(RD0,RD0ACC,TLEAF,Q10F,RTEMP,RESPLIGHTFRAC,TBELOW)
    
    IF(GMESO.LT.0)THEN
        J = QUADM(THETA,-(AJQ*PAR+JMAX),AJQ*PAR*JMAX,IQERROR) ! Actual e-
        VJ = J/DBLE(4.0)                                      ! RuBP-regen rat

	    ! CO2-limited rate of photosynthesis.
	    AC = VCMAX * (CI - GAMMASTAR) / (CI + KM)
    	
	    ! Electron-transport limited rate.
	    AJ = VJ * (CI - GAMMASTAR) / (CI + 2.*GAMMASTAR)

    ELSE
        ! Use Niinemets et al. 2009 formulation of A-Cc curve. We now assume that
        ! VCMAX and JMAX are the rates given CC, not CI! This puts constraints on how they
        ! should be estimated from data.
        
!        ! Eq. 6.
!        AC = QUADM(1/GMESO, (VCMAX-RD)/GMESO-CI-KM, VCMAX*(CI-GAMMASTAR)-RD*(CI+KM) )
!        ! Eq. 7
!        AJ = QUADM(4/GMESO, -(J-4*RD)/GMESO - 4*CI - 8*GAMMASTAR, J*(CI-GAMMASTAR)-4*RD*(CI+2*GAMMASTAR))
!        
!        AC = AC + RD
!        AJ = AJ + RD
        
        ! Got nonsense with the above; try using Ethier and Livingston (2004) (Equation 10).
        AP = -1/GMESO
        BP = (VCMAX - RD)/GMESO + CI + KM
        AC = QUADP(AP,BP, -(VCMAX-RD)*(CI-GAMMASTAR) , IQERROR)
        AJ = QUADP(AP,BP, RD*(CI+KM)-VCMAX*(CI-GAMMASTAR), IQERROR )

        AC = AC + RD
        AJ = AJ + RD
    ENDIF
        	
	    ALEAF = DMIN1(AC,AJ) - RD  
    	
        !  Hyperbolic minimum; use for optimization to avoid discontinuity.
        ! HMSHAPE determines how smooth the transition is (0.999 is a good value).
        ALEAFHM = (AC+AJ - sqrt((AC+AJ)**2-4*HMSHAPE*AC*AJ))/(2*HMSHAPE) - RD
    
    
    RETURN
END SUBROUTINE ACI


!**********************************************************************
      DOUBLE PRECISION FUNCTION TK(TCELSIUS)
! Converts Celsius temperature to Kelvin.
!**********************************************************************

      IMPLICIT NONE
      DOUBLE PRECISION TCELSIUS
      
      TK = TCELSIUS + 273.15 !ABSZERO
      RETURN
      END !TCelsius


!**********************************************************************
DOUBLE PRECISION FUNCTION GAMMAFN(TLEAF, IECO)
! This subroutine calculates Gamma(star), or the CO2 compensation point
! in the absence of non-photorespiratory respiration.
! This is the ECOCRAFT-agreed formulation of this function.
!**********************************************************************

    IMPLICIT NONE
    
    DOUBLE PRECISION, EXTERNAL :: ARRH
    DOUBLE PRECISION :: TLEAF
    INTEGER :: IECO
    
    IF (IECO.EQ.1) THEN
        ! Ecocraft fomulation; based on Brooks & Farquhar and von Caemmerer et al.
        ! If TLEAF < -1.0 then calculate Gamma for T = -1 (quadratic not applicable)
        IF (TLEAF.LT.-1.0) THEN
            GAMMAFN = 36.9 + 1.88*(-26.0) + 0.036*(-26.0)*(-26.0)
        ELSE
            GAMMAFN = 36.9 + 1.88*(TLEAF-25) + 0.036*(TLEAF-25)*(TLEAF-25)
        END IF
    ELSE      ! Bernacchi et al 2001 PCE 24: 253-260
            GAMMAFN = ARRH(DBLE(42.75),DBLE(37830.0),TLEAF,DBLE(25.0))
    ENDIF
    RETURN
END FUNCTION GAMMAFN

!**********************************************************************
DOUBLE PRECISION FUNCTION KMFN(TLEAF,IECO)
! This subroutine calculates Km, or the effective Michaelis-Menten
! coefficient of Rubisco activity.
! This is the ECOCRAFT-agreed formulation of this function.
!**********************************************************************

    !USE maestcom
    IMPLICIT NONE
    
    INTEGER :: IECO
    DOUBLE PRECISION OI,KC25,KO25,KCEA,KOEA,KC,KO, TLEAF
    DOUBLE PRECISION, EXTERNAL :: ARRH
    
    OI = 205000         ! Oxygen partial pressure (umol mol-1)
    IF (IECO.EQ.1) THEN
    ! Physiological constants - values agreed by Ecocraft - Badger & Collatz values
        KC25 = 404          ! MM coefft of Rubisco for CO2 (umol mol-1)
        KO25 = 248000       ! MM coefft of Rubisco for O2 (umol mol-1)
        KCEA = 59400        ! Temp. response of Kc (J mol-1)
        KOEA = 36000        ! Temp. response of Ko (J mol-1)
    ELSE !  Bernacchi et al 2001 PCE 24: 253-260
        KC25 = 404.9        ! MM coefft of Rubisco for CO2 (umol mol-1)
        KO25 = 278400       ! MM coefft of Rubisco for O2 (umol mol-1)
        KCEA = 79430        ! Temp. response of Kc (J mol-1)
        KOEA = 36380        ! Temp. response of Ko (J mol-1)
    END IF

    ! This function is well-behaved for TLEAF < 0.0
    KC = ARRH(KC25,KCEA,TLEAF,DBLE(25.0))
    KO = ARRH(KO25,KOEA,TLEAF,DBLE(25.0))
    KMFN = KC * (1. + OI/KO)

    RETURN
END FUNCTION KMFN


!**********************************************************************
DOUBLE PRECISION FUNCTION JMAXTFN(JMAX25,TLEAF,EAVJ,EDVJ,DELSJ,TVJUP,TVJDN)
! This subroutine calculates the potential electron transport rate
! (Jmax) at the leaf temperature.
! This is the ECOCRAFT-agreed formulation of this function.
!**********************************************************************

    !USE maestcom
    IMPLICIT NONE
   
    DOUBLE PRECISION :: JMAX25,TLEAF,EAVJ,EDVJ,DELSJ, TVJUP, TVJDN, RCONST
    DOUBLE PRECISION, EXTERNAL :: TK
	
    RCONST = 8.314
   
    ! This function is well-behaved for TLEAF < 0.0
    JMAXTFN = JMAX25 * EXP((TLEAF-25)*EAVJ/(RCONST*TK(TLEAF)*TK(DBLE(25.))))  &
                    * (1.+EXP((DELSJ*TK(DBLE(25.))-EDVJ)/(RCONST*TK(DBLE(25.)))))   &
                    / (1.+EXP((DELSJ*TK(TLEAF)-EDVJ)/(RCONST*TK(TLEAF))))

    ! Function allowing Vcmax to be forced linearly to zero at low T -
    ! introduced for Duke data
    IF (TLEAF.LT.TVJDN) THEN
        JMAXTFN = 0.0
    ELSE IF (TLEAF.LT.TVJUP) THEN
        JMAXTFN = (TLEAF - TVJDN)/(TVJUP - TVJDN)*JMAXTFN
    END IF

    RETURN
END FUNCTION JMAXTFN

!**********************************************************************
DOUBLE PRECISION FUNCTION VCMAXTFN(VCMAX25,TLEAF,EAVC,EDVC,DELSC,TVJUP,TVJDN)
! This subroutine calculates the maximum Rubisco activity
! (Vcmax) at the leaf temperature.
! This is the ECOCRAFT-agreed formulation of this function.
!**********************************************************************

    !USE maestcom
    IMPLICIT NONE
    DOUBLE PRECISION VCMAX25,TLEAF,EAVC,EDVC,DELSC,TVJUP,TVJDN,RCONST
    DOUBLE PRECISION, EXTERNAL :: TK
    RCONST = 8.314
    
    ! There is still disagreement as to whether this function has an
    ! optimum or not. Both forms are available here. If EDVC <= 0 (default 0)
    ! then the no-optimum form is used.
    ! Both functions are well-behaved for TLEAF < 0.0

    IF (EDVC.LE.0.0) THEN
        VCMAXTFN = VCMAX25 * EXP(EAVC*(TLEAF - 25)/ (TK(DBLE(25.))*RCONST*TK(TLEAF)))
    ELSE
        VCMAXTFN = VCMAX25 * EXP((TLEAF-25.)*EAVC/(RCONST*TK(TLEAF)*TK(DBLE(25.)))) &
                    * (1.+EXP((DELSC*TK(DBLE(25.))-EDVC)/(RCONST*TK(DBLE(25.))))) &
                    / (1.+EXP((DELSC*TK(TLEAF)-EDVC)/(RCONST*TK(TLEAF))))
    END IF

    ! Function allowing Vcmax to be forced linearly to zero at low T -
    ! introduced for Duke data
    IF (TLEAF.LT.TVJDN) THEN
        VCMAXTFN = 0.0
    ELSE IF (TLEAF.LT.TVJUP) THEN
        VCMAXTFN = (TLEAF - TVJDN)/(TVJUP - TVJDN)*VCMAXTFN
    END IF

    RETURN
END FUNCTION VCMAXTFN


!**********************************************************************
DOUBLE PRECISION FUNCTION RESP(RD0,RD0ACC,TLEAF,Q10F,RTEMP,DAYRESP,TBELOW)
! This function calculates respiration from temperature
! using a Q10 (exponential) formulation.
!**********************************************************************

    IMPLICIT NONE      
    DOUBLE PRECISION RD0,TLEAF,Q10F,RTEMP,DAYRESP,RD0ACC,TBELOW

    IF (TLEAF.GE.TBELOW) THEN
        RESP =  RD0 * EXP(Q10F * (TLEAF-RTEMP)/10) * DAYRESP
    ELSE
        RESP = 0.0
    END IF

    RETURN
END FUNCTION RESP


!**********************************************************************
DOUBLE PRECISION FUNCTION ARRH(KT,EA,T,TREF)
! The Arrhenius function.
! KT is the value at Tref deg C; Ea the activation energy (J mol-1) and T the temp (deg C).
!**********************************************************************

    !USE maestcom
    IMPLICIT NONE
    DOUBLE PRECISION KT, EA, T, TREF, RCONST, ABSZERO
    
    ABSZERO = -273.15
    RCONST = 8.314
	
    ARRH = KT*EXP(EA*(T-TREF)/(RCONST*(T-ABSZERO)*(TREF-ABSZERO)))
    RETURN
END FUNCTION ARRH


!**********************************************************************
DOUBLE PRECISION FUNCTION QUADM(A, B, C, IQERROR)
! Solves the quadratic equation - finds smaller root.
!**********************************************************************

    !USE maestcom
    IMPLICIT NONE
    DOUBLE PRECISION A, B, C
    INTEGER :: IQERROR
    
    IF ((B*B - 4.*A*C).LT.0.0) THEN
        !CALL SUBERROR('WARNING:IMAGINARY ROOTS IN QUADRATIC',IWARN,0)
        IQERROR = 1
        QUADM = 0.0
    ELSE
        IF (A.EQ.0.0) THEN
            IF (B.EQ.0.0) THEN
                QUADM = 0.0
                !IF (C.NE.0.0) CALL SUBERROR('ERROR: CANT SOLVE QUADRATIC',IWARN,0)
            ELSE
                QUADM = -C/B
            END IF
        ELSE
            QUADM = (- B - SQRT(B*B - 4*A*C)) / (2.*A)
        END IF
    END IF

    RETURN
END FUNCTION QUADM


!**********************************************************************
DOUBLE PRECISION FUNCTION QUADP(A,B,C,IQERROR)
! Solves the quadratic equation - finds larger root.
!**********************************************************************

    !USE maestcom
    IMPLICIT NONE
    DOUBLE PRECISION A,B,C
    INTEGER IQERROR
   

    IF ((B*B - 4.*A*C).LT.0.0) THEN
        !CALL SUBERROR('WARNING:IMAGINARY ROOTS IN QUADRATIC',IWARN,0)
        IQERROR = 1
        QUADP = 0.0
    ELSE
        IF (A.EQ.0.0) THEN
            IF (B.EQ.0.0) THEN
                QUADP = 0.0
                !IF (C.NE.0.0) CALL SUBERROR('ERROR: CANT SOLVE QUADRATIC',IWARN,0)
            ELSE
                QUADP = -C/B
            END IF
        ELSE
            QUADP = (- B + SQRT(B*B - 4*A*C)) / (2.*A)
        END IF
    END IF

    RETURN
END FUNCTION QUADP





SUBROUTINE PSILFIND(PAR,TLEAF,CS,RH,VPD,PATM, &
        JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN, &
        THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW, &
        MODELGS,EMAXLEAF,KTOT,WEIGHTEDSWP,PSIL,GSMIN, &
	    G0,D0L,GAMMA,G1,GK,TUZFUN,SF,PSIV,HMSHAPE, &
        FSOIL,GS,ALEAF,ALEAFHM,RD)

        IMPLICIT NONE

        INTEGER IECO,MODELGS
        DOUBLE PRECISION PAR,TLEAF,CS,RH,VPD,PATM,JMAX25,EAVJ,EDVJ
        DOUBLE PRECISION DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN
        DOUBLE PRECISION THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW
        DOUBLE PRECISION EMAXLEAF,KTOT,WEIGHTEDSWP,PSILCUR,PSIL,GSMIN
        DOUBLE PRECISION G0,D0L,GAMMA,G1,GK,SF,PSIV,HMSHAPE
        DOUBLE PRECISION FSOIL,GS,ALEAF,ALEAFHM,RD,T1,T2,XACC
        DOUBLE PRECISION EXTRAPARS(50),PSILIN
        INTEGER EXTRAINT(10)
		INTEGER TUZFUN
        DOUBLE PRECISION, EXTERNAL :: ZBRENT
        DOUBLE PRECISION, EXTERNAL :: PSILOBJFUN
        
        ! INIT
        PSILIN = -0.1
        PSIL = -0.1
        
        
        EXTRAINT(1) = IECO
        EXTRAINT(2) = MODELGS
		EXTRAINT(3) = TUZFUN
		
        EXTRAPARS(1) = PAR
        EXTRAPARS(2) = TLEAF  
        EXTRAPARS(3) = CS
        EXTRAPARS(4) = RH 
        EXTRAPARS(5) = VPD 
        EXTRAPARS(6) = PATM  
        EXTRAPARS(7) = JMAX25 
        EXTRAPARS(8) = EAVJ 
        EXTRAPARS(9) = EDVJ 
        EXTRAPARS(10) = DELSJ 
        EXTRAPARS(11) = VCMAX25 
        EXTRAPARS(12) = EAVC 
        EXTRAPARS(13) = EDVC 
        EXTRAPARS(14) = DELSC 
        EXTRAPARS(15) = TVJUP 
        EXTRAPARS(16) = TVJDN 
        EXTRAPARS(17) = THETA 
        EXTRAPARS(18) = AJQ 
        EXTRAPARS(19) = RD0 
        EXTRAPARS(20) = Q10F
        EXTRAPARS(21) = RTEMP
        EXTRAPARS(22) = DAYRESP 
        EXTRAPARS(23) = TBELOW 
        EXTRAPARS(24) = EMAXLEAF 
        EXTRAPARS(25) = KTOT 
        EXTRAPARS(26) = WEIGHTEDSWP 
        !EXTRAPARS(27) = PSIL 
        EXTRAPARS(28) = PSILIN 
        EXTRAPARS(27) = -999
        !EXTRAPARS(28) = -999
        EXTRAPARS(29) = GSMIN 
        EXTRAPARS(30) = G0 
        EXTRAPARS(31) = D0L
        EXTRAPARS(32) = GAMMA
        EXTRAPARS(33) = G1 
        EXTRAPARS(34) = GK 
        EXTRAPARS(35) = SF 
        EXTRAPARS(36) = PSIV 
        EXTRAPARS(37) = HMSHAPE 
        EXTRAPARS(38) = FSOIL 
        EXTRAPARS(39) = GS
        EXTRAPARS(40) = ALEAF 
        EXTRAPARS(41) = ALEAFHM 
        EXTRAPARS(42) = RD 
        
        ! Set bounds for root-finding
        T1 = -100.0
        T2 = 0.0

        ! Error tolerance.
        XACC = 1E-03
        
        PSIL = ZBRENT(PSILOBJFUN,T1,T2,XACC,EXTRAPARS,EXTRAINT)


END




!**********************************************************************
DOUBLE PRECISION FUNCTION PSILOBJFUN(PSILIN, EXTRAPARS, EXTRAINT)

! New PSIL finder:
! Make wrapper function that takes PSIL as first argument, EXTRAPARS an array
! second argument, and returns the squared difference in PSILIN and PSIL
! Then : PSILOPT = PSIFLFINDER(PSIL, EXTRAPARS)
! Load and unload parameters from EXTRAPARS...
!**********************************************************************
        IMPLICIT NONE

        INTEGER IECO,MODELGS
        DOUBLE PRECISION PAR,TLEAF,CS,RH,VPD,PATM,JMAX25,EAVJ,EDVJ
        DOUBLE PRECISION DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN
        DOUBLE PRECISION THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW
        DOUBLE PRECISION EMAXLEAF,KTOT,WEIGHTEDSWP,PSIL,GSMIN
        DOUBLE PRECISION G0,D0L,GAMMA,G1,GK,SF,PSIV,HMSHAPE
        DOUBLE PRECISION FSOIL,GS,ALEAF,ALEAFHM,RD
        DOUBLE PRECISION EXTRAPARS(50),PSILIN,CI
        INTEGER EXTRAINT(10)
		INTEGER TUZFUN

        IECO = EXTRAINT(1)
        MODELGS = EXTRAINT(2)
		TUZFUN = EXTRAINT(3)
		
        PAR = EXTRAPARS(1) 
        TLEAF = EXTRAPARS(2) 
        CS = EXTRAPARS(3)
        RH = EXTRAPARS(4)
        VPD = EXTRAPARS(5)
        PATM = EXTRAPARS(6)
        JMAX25 = EXTRAPARS(7)
        EAVJ = EXTRAPARS(8) 
        EDVJ = EXTRAPARS(9)
        DELSJ = EXTRAPARS(10)
        VCMAX25 = EXTRAPARS(11)
        EAVC = EXTRAPARS(12)
        EDVC = EXTRAPARS(13)
        DELSC = EXTRAPARS(14)
        TVJUP = EXTRAPARS(15)
        TVJDN = EXTRAPARS(16)
        THETA = EXTRAPARS(17)
        AJQ = EXTRAPARS(18)
        RD0 = EXTRAPARS(19)
        Q10F = EXTRAPARS(20)
        RTEMP = EXTRAPARS(21)
        DAYRESP = EXTRAPARS(22)
        TBELOW = EXTRAPARS(23)
        EMAXLEAF = EXTRAPARS(24)
        KTOT = EXTRAPARS(25)
        WEIGHTEDSWP = EXTRAPARS(26)
        !PSIL = EXTRAPARS(27)
        !PSILIN = EXTRAPARS(28)
        GSMIN = EXTRAPARS(29)
        G0 = EXTRAPARS(30)
        D0L = EXTRAPARS(31)
        GAMMA = EXTRAPARS(32)
        G1 = EXTRAPARS(33)
        GK = EXTRAPARS(34)
        SF = EXTRAPARS(35)
        PSIV = EXTRAPARS(36)
        HMSHAPE = EXTRAPARS(37)
        FSOIL = EXTRAPARS(38)
        GS = EXTRAPARS(39)
        ALEAF = EXTRAPARS(40)
        ALEAFHM = EXTRAPARS(41)
        RD = EXTRAPARS(42)
		
        CI = -999 ! output, not used here
        
        ! CALL PHOTOSYN(PAR,TLEAF,CS,RH,VPD,PATM, &
        ! JMAX25,IECO,EAVJ,EDVJ,DELSJ,VCMAX25,EAVC,EDVC,DELSC,TVJUP,TVJDN, &
        ! THETA,AJQ,RD0,Q10F,RTEMP,DAYRESP,TBELOW, &
        ! MODELGS,EMAXLEAF,KTOT,WEIGHTEDSWP,PSILIN,PSIL,GSMIN, &
	    ! G0,D0L,GAMMA,G1,GK,TUZFUN,SF,PSIV,HMSHAPE, &
        ! FSOIL,GS,ALEAF,ALEAFHM,RD,CI)
        
        PSILOBJFUN = PSILIN - PSIL

END

!**********************************************************************
DOUBLE PRECISION FUNCTION ZBRENT(FUNC,X1,X2,TOL,EXTRAPARS,EXTRAINT)
! Function taken from SPA.
! Bisection routine; finds the root of FUNC in the interval (X1,X2)
!**********************************************************************
    IMPLICIT NONE
    INTEGER ITMAX,ITER
    DOUBLE PRECISION TOL,X1,X2,EPS
    DOUBLE PRECISION EXTRAPARS(50)  ! Parameters passed to FUNC.
    INTEGER EXTRAINT(10)
    PARAMETER (ITMAX=30,EPS=3.E-8)
    DOUBLE PRECISION A,B,C,D,E,FA,FB,FC,P,Q,R,S,TOL1,XM
    DOUBLE PRECISION , EXTERNAL :: FUNC
    
    A = X1
    B = X2
    FA = FUNC(A, EXTRAPARS, EXTRAINT)
    FB = FUNC(B, EXTRAPARS, EXTRAINT)
    IF((FA.GT.0..AND.FB.GT.0.).OR.(FA.LT.0..AND.FB.LT.0.))THEN
        FA=FUNC(A, EXTRAPARS, EXTRAINT)
        FB=FUNC(B, EXTRAPARS, EXTRAINT)
        ! WRITE(*,*)' 	    FA		  FB		  X1		    X2'
        ! WRITE(*,*)FA,FB,X1,X2
        ! WRITE(*,*)'ROOT MUST BE BRACKETED FOR ZBRENT'
        FA=FUNC(A, EXTRAPARS, EXTRAINT)
        FB=FUNC(B, EXTRAPARS, EXTRAINT)
    ENDIF
    C=B
    FC=FB
    DO ITER=1,ITMAX
        IF((FB.GT.0..AND.FC.GT.0.).OR.(FB.LT.0..AND.FC.LT.0.))THEN
            C=A
            FC=FA
            D=B-A
            E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
            A=B
            B=C
            C=A
            FA=FB
            FB=FC
            FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
            ZBRENT=B      
            RETURN        
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
            S=FB/FA           
            IF(A.EQ.C) THEN   
                P=2.*XM*S
                Q=1.-S
            ELSE
                Q=FA/FC
                R=FB/FC
                P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
                Q=(Q-1.)*(R-1.)*(S-1.)
            ENDIF
            IF(P.GT.0.) Q=-Q
            P=ABS(P)
            IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
                E=D
                D=P/Q
            ELSE
                D=XM
                E=D
            ENDIF
        ELSE
            D=XM
            E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
            B=B+D
        ELSE
            B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B, EXTRAPARS, EXTRAINT)
    END DO
    !WRITE(*,*) 'ZBRENT EXCEEDING MAXIMUM ITERATIONS'
    ZBRENT=B
    RETURN
END FUNCTION ZBRENT
         
