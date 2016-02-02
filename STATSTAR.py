#
#   General Description:
#   ====================
#       StatStar computes a ZAMS model using a number of simplifying assumptions about the physics.  The code
#       is not designed to produce precise, research-quality models of ZAMS stars; rather, it is meant to
#       illustrate many of the fundamental physical ideas discussed in
#
#           "An Introduction to Modern Astrophysics"
#           Bradley W. Carroll and Dale A. Ostlie
#           Second Edition, Addison Wesley,   Copyright 2007.
#
#       StatStar performs an inward integration of the stellar structure equations, given values for the
#       star's mass, luminosity, effective temperature, and composition.
#
#       The simplifying assumptions made here include:
#           (a) A constant composition throughout (characteristic of a ZAMS star).
#           (b) The use of the ideal gas law throughout.
#           (c) The gas is assumed to be completely ionized throughout.
#           (d) Radiation pressure is incorporated.
#           (e) Convection, when present, is taken to be purely adiabatic.
#           (f) Opacity is computed by using approximation formulae:
#               1.  Bound-free processes via Kramer's formula.
#               2.  Free-free processes via Kramer's formula.
#               3.  Electron scattering via Thomson formula.
#               4.  H- ion via fitting function.
#           (g) Surface is assumed to have P = 0, T = 0, rho = 0.
#           (h) Outermost (optically thin) zone is assumed to be radiative.
#           (i) No attempt is made to satisfy the Eddington approximation by
#               adjusting the outside boundary condition.
#
#---------------------------------------------------------------------

    IMPLICIT NONE
    REAL(dp)                                ::  Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, mu
    REAL(dp)                                ::  r, P, T, M_r, L_r
    REAL(dp)                                ::  P_0, T_0, rho_0, epsilon_0, kappa_0
    REAL(dp)                                ::  dr
    REAL(dp),   DIMENSION(n)                ::  PMLT0, PMLT
    REAL(dp),                   PARAMETER   ::  dr_over_r = 1.0E-03         #Initial fractional step size

    REAL(dp),                   PARAMETER   ::  M_fraction_limit = 0.01     #Mass fraction stop condition
    REAL(dp),                   PARAMETER   ::  L_fraction_limit = 0.10     #Luminosity stop condition
    REAL(dp),                   PARAMETER   ::  r_fraction_limit = 0.02     #radius stop condition
    INTEGER,                    PARAMETER   ::  maximum_zones = 10000       #Maximum number of zones allowed

    INTEGER,                    PARAMETER   ::  n_surface = 1               #Number of surface boundary zones
    INTEGER                                 ::  i                           #Zone counter
    INTEGER                                 ::  ios                         #I/O status flag
    CHARACTER(1)                            ::  all_new = "Y"               #Select new model parameter

    LOGICAL                                 ::  ok_surface, ok_core, ok_Runge
    LOGICAL,                    PARAMETER   ::  adjust_step_size = .TRUE.   #Allow variable step size

    CHARACTER(26)                           ::  format_table = '(I5, 9ES11.3, 2X, A, F5.1)'

#---------------------------------------------------------------------

    New_model: DO
        i = 0       #Initialize zone number

        CALL Initial_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new)
        IF (all_new == "E" .OR. all_new == "e") EXIT New_model

        OPEN(UNIT = 10, FILE = "ZAMSmodel.txt", STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ios)
        IF (ios /= 0) THEN
            WRITE (*,*) "Unable to open output file 'ZAMSmodel.txt' --- terminating calculation"
            EXIT New_model
        END IF

#       Write input data to the output file
        WRITE (10,'(T46,"A ZAMS Stellar Model")')
        WRITE (10,'(T46,"--------------------",/)')
        WRITE (10,'(T46,"M    = ", F11.5, " solar")') Msolar
        WRITE (10,'(T46,"L    = ", F11.5, " solar")') Lsolar
        WRITE (10,'(T46,"R    = ", F11.5, " solar")') Rsolar
        WRITE (10,'(T46,"Teff = ", F11.5, " K")')     Teff
        WRITE (10,'(T46,"X    = ", F11.5)')           X
        WRITE (10,'(T46,"Y    = ", F11.5)')           Y
        WRITE (10,'(T46,"Z    = ", F11.5)')           Z
        WRITE (10,'(//)')

#       Set up previous zone values
        Pm   = 0
        Tm   = 0
        Xm   = X
        Zm   = Z
        rm   = Rs
        taum = 0
        rhom = 0
        kappam = 0
        epsilonm = 0
        dlnPdlnT = 99.9
        rc_flag = "r"

        WRITE (10,'(" zone      r         tau     1-M_r/Ms      L_r         T          P         rho        &
                    &kap        eps    dlnPdlnT")')

        WRITE (10,format_table) i, rm, 0.0, 0.0, Ls, Tm, Pm, rhom, kappam, epsilonm, rc_flag, dlnPdlnT

#       Compute surface zones and step size
        Mm = Ms
        Lm = Ls
        rm = Rs
        dr = -dr_over_r*Rs
        step_size_condition = 0
        Surface_boundary: DO
            i = i + 1

#           Update last zone values
            IF (i > 1) THEN
                Mm = M_r
                Lm = L_r
                rm = r
                Pm = P
                Tm = T
                Xm = X
                Zm = Z
                taum = tau
                rhom = rho
                kappam = kappa
                epsilonm = epsilon
            END IF

            CALL Surface(i, Mm, Lm, rm, X, Z, dr, r, P, T, M_r, L_r, rho, kappa, epsilon, ok_surface)
            IF (.NOT. ok_surface) EXIT Surface_boundary

            tau = taum + Optical_Depth_Change(kappa, kappam, rho, rhom, r, rm)
            WRITE (10,format_table) i, r, tau, 1 - M_r/Ms, L_r, T, P, rho, kappa, epsilon, rc_flag, dlnPdlnT

            IF (i == n_surface) EXIT Surface_boundary
        END DO Surface_boundary

        Satisfactory_Surface: IF (ok_surface) THEN
#           Load array of first derivatives to start the general inward integration
            Y        = Helium(X, Z)
            mu       = Mean_Molecular_Weight(X, Y, Z)
            gamma    = Specific_Heat_Ratio()
            dlnPdlnT = PTgradient(Pm, P, Tm, T)

            dfdr0(1) = dPdr(M_r, rho, r)
            dfdr0(2) = dMdr(r, rho)
            dfdr0(3) = dLdr(r, rho, epsilon)
            dfdr0(4) = dTdr(kappa, rho, T, L_r, r, mu, M_r, gamma, dlnPdlnT)

#           Main inward integration loop
            Main_Loop: DO
                i = i + 1

#               Update last zone values
                Mm = M_r
                Lm = L_r
                Pm = P
                Tm = T
                Xm = X
                Zm = Z
                rm = r
                taum = tau
                rhom = rho
                kappam = kappa
                epsilonm = epsilon

                PMLT0 = (/ Pm, Mm, Lm, Tm /)
                CALL RK_4(n, rm, dr, PMLT0, PMLT, dfdr0, Structure_Eqns, ok_Runge)
                IF (.NOT. ok_Runge) EXIT Main_Loop

#               Results from the current step
                P   = PMLT(1)
                M_r = PMLT(2)
                L_r = PMLT(3)
                T   = PMLT(4)

                tau = taum + Optical_Depth_Change(kappa, kappam, rho, rhom, rm + dr, rm)

                WRITE (10,format_table) i, r, tau, 1-M_r/Ms, L_r, T, P, rho, kappa, epsilon, rc_flag, dlnPdlnT
                IF ((M_r/Ms < M_fraction_limit .AND. L_r/Ls < L_fraction_limit .AND. r/Rs < r_fraction_limit) &
                        .OR. T < 0 .OR. P < 0) EXIT Main_Loop
                IF (i > maximum_zones) THEN
                    CALL Too_Many_Zones(i, Msolar, Ms, M_r, Lsolar, Ls, L_r, r, Rs, Rsolar, Teff, X, Y, Z, P_0, T_0, &
                        rho_0, kappa_0, epsilon_0, rc_flag)
                    ok_Runge = .FALSE.
                    EXIT Main_Loop
                END IF

#               Is it time to change step size?
                IF (adjust_step_size) THEN
                    SELECT CASE (step_size_condition)
                        CASE(0)
                            IF (M_r < 0.99*Ms) THEN
                                dr = -Rs/100
                                step_size_condition = 1
                            END IF
                        CASE(1)
                            IF (ABS(dr) > 5*r) THEN
                                dr = dr/10
                                step_size_condition = 2
                            END IF
                    END SELECT
                END IF
                r = r + dr
            END DO Main_Loop

            Core_Extrapolation: IF (ok_Runge) THEN
#               Determine core conditions
                i = i + 1
                CALL Core(M_r, L_r, P, T, X, Z, r, P_0, T_0, rho_0, kappa_0, epsilon_0, rc_flag, dlnPdlnT, ok_core)
                IF (.NOT. ok_core) THEN
                    WRITE (*,'(/,"WARNING:  There was a problem with the core extrapolation routine",/)')
                END IF

                tau = tau + Optical_Depth_Change(kappa_0, kappa, rho_0, rho, 0.0_dp, r)
                WRITE (10,format_table) i, 0, tau, 1-M_r/Ms, L_r, T_0, P_0, rho_0, kappa_0, epsilon_0, rc_flag, dlnPdlnT

#               Write initial and final conditions to the screen
                CALL Final_Results(i, Msolar, Ms, M_r, Lsolar, Ls, L_r, r, Rs, Rsolar, Teff, X, Y, Z, &
                        P, T, rho, kappa, epsilon, P_0, T_0, rho_0, kappa_0, epsilon_0, rc_flag)
            END IF Core_Extrapolation
        END IF Satisfactory_Surface

#       Does the user want to compute a new model?
        all_new = "Y"
        CALL Change_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new)
        CLOSE (UNIT = 10, IOSTAT = ios)
        IF (ios /= 0) THEN
            WRITE (*,*) "Unable to close the output file - the new model calculation is being aborted"
            EXIT New_model
        END IF
        IF (all_new == "E") Exit New_model
        WRITE (*, '(//)')
    END DO New_model

#   General Description:
#   ====================
#
#       This module contains the basic equations of stellar structure.  The module also
#       contains a driver function that selects among the required equations for the
#       Runge Kutta routines.
#---------------------------------------------------------------------

    USE Constants, ONLY         :   dp, G, four_pi, a => a_rad, four_ac_o3, c, m_H, k => k_B
    USE Zone_Quantities, ONLY   :   n
    PRIVATE                     ::  dp, G, four_pi, a, four_ac_o3, c, m_H, k

    CONTAINS

#   Driver for stellar structure equations
    REAL(dp) FUNCTION Structure_Eqns(i, r, S, ok)   RESULT(dfdr)

        USE Composition
        USE Physics
        USE Zone_Quantities,  X => Xm, Z => Zm

        IMPLICIT NONE
        INTEGER,                    INTENT(IN)  ::  i
        REAL(dp),                   INTENT(IN)  ::  r       #independent variable
        REAL(dp),   DIMENSION(n),   INTENT(IN)  ::  S       #dependent variables
        LOGICAL,                    INTENT(OUT) ::  ok      #Returns .TRUE. if the derivative calculation was successful

        REAL(dp)                                ::  P, M_r, L_r, T
        REAL(dp)                                ::  Y, XCNO, mu

        ok = .TRUE.

        P   = S(1)
        M_r = S(2)
        L_r = S(3)
        T   = S(4)

        Y   = Helium(X, Z)
        mu  = Mean_Molecular_Weight(X, Y, Z)
        rho = Density(T, P, mu)
        IF (rho < 0) THEN
            WRITE (*, '("Density calculation error in FUNCTION Structure_Eqns")')
            ok = .FALSE.
        END IF

        SELECT CASE (i)
            CASE (1)
                if (ok):
                    dfdr = dPdr(M_r, rho, r)
                else:
                    dfdr = 0
                dfdr0(1) = dfdr                     #Save result for next zone start

            CASE (2)
                if (ok):
                    dfdr = dMdr(r, rho)
                else:
                    dfdr = 0
                dfdr0(2) = dfdr                     #Save result for next zone start

            CASE (3)
                if (ok):
                    epsilon  = Nuclear(T, rho, X, Z)
                    dfdr     = dLdr(r, rho, epsilon)
                else
                    dfdr     = 0
                dfdr0(3) = dfdr                     #Save result for next zone start

            CASE (4)
                if (ok):
                    kappa    = Opacity(T, rho, X, Z)
                    gamma    = Specific_Heat_Ratio()
                    dlnPdlnT = PTgradient(Pm, P, Tm, T)
                    dfdr     = dTdr(kappa, rho, T, L_r, r, mu, M_r, gamma, dlnPdlnT)
                else:
                    dfdr     = 0
                dfdr0(4) = dfdr                     #Save result for next zone start
        END SELECT
    END FUNCTION Structure_Eqns
#---------------------------------------------------------------------

#   Hydrostatic Equilibrium
    REAL(dp) FUNCTION dPdr(M_r, rho, r)
        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  M_r, rho, r

        dPdr = -G*M_r*rho/r**2              #Eq. (10.6)
    END FUNCTION dPdr
#---------------------------------------------------------------------

#   Mass Conservation
    REAL(dp) FUNCTION dMdr(r, rho)
        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  r, rho

        dMdr = four_pi*r**2*rho             #Eq. (10.7)
    END FUNCTION dMdr
#---------------------------------------------------------------------

#   Luminosity Gradient
    REAL(dp) FUNCTION dLdr(r, rho, epsilon)
        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  r, rho, epsilon

        dLdr = four_pi*r**2*rho*epsilon     #Eq. (10.36)
    END FUNCTION dLdr
#---------------------------------------------------------------------

#   Temperature Gradient
    REAL(dp) FUNCTION dTdr(kappa, rho, T, L_r, r, mu, M_r, gamma, dlnPdlnT)

        USE Zone_Quantities, ONLY   : rc_flag

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  kappa, rho, T, L_r, r, mu, M_r, gamma, dlnPdlnT
        REAL(dp)                ::  gamma_ratio

        gamma_ratio = gamma/(gamma - 1)
        IF (dlnPdlnT > gamma_ratio) THEN                                #radiation criterion,   Eq. (10.95)
            dTdr = -(kappa*rho/T**3)*(L_r/(four_pi*r**2))/four_ac_o3    #radiation,             Eq. (10.68)
            rc_flag = "r"
        ELSE
            dTdr = -(1/gamma_ratio)*(mu*m_H/k)*(G*M_r/r**2)             #adiabatic convection,  Eq. (10.89)
            rc_flag = "c"
        END IF
    END FUNCTION dTdr

#   General Description:
#   ====================
#       This module takes care of all the user input requests
#---------------------------------------------------------------------

    USE Constants, ONLY :   dp, M_Sun, L_Sun, R_Sun, four_pi, sigma
    PRIVATE             ::  dp, M_Sun, L_Sun, R_Sun, four_pi, sigma

    CONTAINS
    SUBROUTINE Initial_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new)

#       General Description:
#       ====================
#           Gather the initial model data

        USE Composition

        IMPLICIT NONE
        REAL(dp),       INTENT(OUT)     ::  Msolar, Lsolar, Rsolar
        REAL(dp),       INTENT(OUT)     ::  Ms, Ls, Rs, Teff, X, Y, Z
        CHARACTER,      INTENT(INOUT)   ::  all_new
        CHARACTER(1)                    ::  y_n
        INTEGER                         ::  ios

#       Write introductory information for user
        WRITE (*,'(   T15, "StatStar is designed to build a ZAMS star")')

        WRITE (*,'(//,T15, "Details of the code are described in:")')
        WRITE (*,'(   T15, "   An Introduction to Modern Astrophysics")')
        WRITE (*,'(   T15, "   Bradley W. Carroll and Dale A. Ostlie")')
        WRITE (*,'(   T15, "      Second Edition, Addison Wesley")')
        WRITE (*,'(   T15, "      copyright 2007",//)')

        WRITE (*,'(   T15, "The user will be asked to enter the following quantities:")')
        WRITE (*,'(   T15, "   Mass of the star       (in solar units)")')
        WRITE (*,'(   T15, "   Luminosity of the star (in solar units)")')
        WRITE (*,'(   T15, "   Effective Temperature  (in K)")')
        WRITE (*,'(   T15, "   Hydrogen mass fraction (X)")')
        WRITE (*,'(   T15, "   Metals mass fraction   (Z)",//)')

        All_parameters: DO
            IF (all_new == "Y" .OR. all_new == "y" .OR. all_new == "A" .OR. all_new == "a") THEN
                M_loop: DO
                    WRITE (*,'("Enter the mass (in solar units)       :  Msolar = ")', ADVANCE = 'NO')
                    READ  (*,*, iostat = ios) Msolar
                    IF (ios == 0 .AND. Msolar > 0) EXIT M_loop
                    WRITE (*,'("Invalid value entered - please try again",/)')
                END DO M_loop

                T_loop: DO
                    WRITE (*,'("Enter the effective temperature (in K):  Teff   = ")', ADVANCE = 'NO')
                    READ  (*,*, iostat = ios) Teff
                    IF (ios == 0 .AND. Teff > 0) EXIT T_loop
                    WRITE (*,'("Invalid value entered - please try again",/)')
                END DO T_loop

                L_loop: DO
                    WRITE (*,'("Enter the luminosity (in solar units) :  Lsolar = ")', ADVANCE = 'NO')
                    READ  (*,*, iostat = ios) Lsolar
                    IF (ios == 0 .AND. Lsolar > 0) EXIT L_loop
                    WRITE (*,'("Invalid value entered - please try again",/)')
                END DO L_loop

                XYZ: DO
                    X_loop: DO
                        WRITE (*,'("Enter the mass fraction of hydrogen   :  X      = ")', ADVANCE = 'NO')
                        READ  (*,*, iostat = ios) X
                        IF (ios == 0) THEN
                            IF (X >= 0 .AND. X <= 1) THEN
                                EXIT X_loop
                            ELSE
                                WRITE (*,'(/,"0 <= X <= 1 is required")')
                            END IF
                        ELSE
                            WRITE (*,'(/)')
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO X_loop

                    Z_loop: DO
                        WRITE (*,'("Enter the mass fraction of metals     :  Z      = ")', ADVANCE = 'NO')
                        READ  (*,*, iostat = ios) Z
                        IF (ios == 0) THEN
                            IF (Z >= 0 .AND. Z <= 1) THEN
                                EXIT Z_loop
                            ELSE
                                WRITE (*,'(/,"0 <= Z <= 1 is required")')
                            END IF
                        ELSE
                            WRITE (*,'(/)')
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO Z_loop

                    Y = Helium(X, Z)
                    IF (Y < 0 .OR. Y > 1) THEN
                        WRITE (*,'("Note that 0 <= X + Z <= 1 is required", /, "  Please reenter composition values",/)')
                    ELSE
                        EXIT XYZ
                    END IF
                END DO XYZ
            END IF

#           Compute SI values
            Ms = Msolar*M_Sun
            Ls = Lsolar*L_Sun
            Rs = SQRT(Ls/(four_pi*sigma*Teff**4))   #Eq. (3.17)
            Rsolar = Rs/R_Sun                       #Solar radius from SI value

#           Allow the user the opportunity to change values as needed
            Fix_parameters: DO
                WRITE (*,'(/,"Your model parameters are:", &
                    &      T42, "M      = ", F11.5, " solar")') Msolar
                WRITE (*,'(T42, "Teff   = ", F11.5, " K    ")') Teff
                WRITE (*,'(T42, "L      = ", F11.5, " solar")') Lsolar
                WRITE (*,'(T42, "R      = ", F11.5, " solar")') Rsolar
                WRITE (*,'(T42, "X      = ", F11.5)') X
                WRITE (*,'(T42, "Y      = ", F11.5)') Y
                WRITE (*,'(T42, "Z      = ", F11.5)') Z

                WRITE (*,'(/,"Are these values ok (y/n)? ")', ADVANCE = 'NO')
                Yes_No: DO
                    READ  (*,*) y_n
                    IF (y_n == "Y" .OR. y_n == "y" .OR. y_n == "N" .OR. y_n == "n") EXIT Yes_No
                    WRITE (*,'("Please answer Yes (y) or No (n):  ")', ADVANCE = 'NO')
                END DO Yes_No
                IF (y_n == "Y" .OR. y_n == "y") EXIT All_parameters

                all_new = "N"
                CALL Change_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new)
                IF (all_new == "E" .OR. all_new == "e") EXIT All_parameters
                IF (all_new == "A" .OR. all_new == "a") EXIT Fix_parameters
            END DO Fix_parameters
        END DO All_parameters
    END SUBROUTINE Initial_Model
#---------------------------------------------------------------------

    SUBROUTINE Change_Model(Msolar, Lsolar, Rsolar, Ms, Ls, Rs, Teff, X, Y, Z, all_new)

#       General Description:
#       ====================
#           Get updated model input data

        USE Composition

        IMPLICIT NONE
        REAL(dp),       INTENT(INOUT)   ::  Msolar, Lsolar, Rsolar
        REAL(dp),       INTENT(INOUT)   ::  Ms, Ls, Rs, Teff, X, Y, Z
        CHARACTER(1),   INTENT(INOUT)   ::  all_new
        CHARACTER(1)                    ::  y_n
        INTEGER                         ::  ios

        IF (all_new == "Y" .OR. all_new == "y") THEN
            WRITE (*,'(/,"Would you like to run another model?", /, &
                    &    "Your previous results will be overwritten in the output file. (y/n): ")', ADVANCE = 'NO')
            Yes_No: DO
                READ  (*,*) y_n
                IF (y_n == "Y" .OR. y_n == "y" .OR. y_n == "N" .OR. y_n == "n") EXIT Yes_No
                WRITE (*,'("Please answer Yes (y) or No (n):  ")', ADVANCE = 'NO')
            END DO Yes_No

            IF (y_n == "Y" .OR. y_n == "y") THEN
                WRITE (*, '(/,"Which variable would you like to change?")')
                WRITE (*, '("     M = Mass",                       T40, "Current value = ", F11.5, " solar")') Msolar
                WRITE (*, '("     T = effective Temperature",      T40, "Current value = ", F11.5, " K    ")') Teff
                WRITE (*, '("     L = Luminosity",                 T40, "Current value = ", F11.5, " solar")') Lsolar
                WRITE (*, '("     X = hydrogen mass mraction (X)", T40, "Current value = ", F11.5, "      ")') X
                WRITE (*, '("     Z = metals mass fraction (Z)",   T40, "Current value = ", F11.5, "      ")') Z
                WRITE (*, '("     A = select an All new set of model parameters")')
                WRITE (*, '("     E = Exit the calculation")')
                WRITE (*, '(/, "Select a letter: ")', ADVANCE = 'NO')
                All_New1: DO
                    READ  (*,*) all_new
                    IF (all_new == "M" .OR. all_new == "m" .OR. &
                        all_new == "T" .OR. all_new == "t" .OR. &
                        all_new == "L" .OR. all_new == "l" .OR. &
                        all_new == "X" .OR. all_new == "x" .OR. &
                        all_new == "Z" .OR. all_new == "z" .OR. &
                        all_new == "A" .OR. all_new == "a" .OR. &
                        all_new == "E" .OR. all_new == "e") EXIT All_New1
                    WRITE (*,'("Please respond with one of the options listed:  ")', ADVANCE = 'NO')
                END DO All_New1
            END IF
        ELSE
            y_n = "Y"
            WRITE (*, '(/,"Which variable would you like to change?")')
            WRITE (*, '("     M = Mass")')
            WRITE (*, '("     T = effective Temperature")')
            WRITE (*, '("     L = Luminosity")')
            WRITE (*, '("     X = hydrogen mass fraction (X)")')
            WRITE (*, '("     Z = metals mass fraction (Z)")')
            WRITE (*, '("     A = select an All new set of model parameters")')
            WRITE (*, '("     E = Exit the calculation")')
            WRITE (*, '(/, "Select a letter: ")', ADVANCE = 'NO')
            All_New2: DO
                READ  (*,*) all_new
                IF (all_new == "M" .OR. all_new == "m" .OR. &
                    all_new == "T" .OR. all_new == "t" .OR. &
                    all_new == "L" .OR. all_new == "l" .OR. &
                    all_new == "X" .OR. all_new == "x" .OR. &
                    all_new == "Z" .OR. all_new == "z" .OR. &
                    all_new == "A" .OR. all_new == "a" .OR. &
                    all_new == "E" .OR. all_new == "e") EXIT All_New2
                WRITE (*,'("Please respond with one of the options listed:  ")', ADVANCE = 'NO')
            END DO All_New2
            IF (all_new == "E" .OR. all_new == "e") y_n = "N"
        END IF

        IF (y_n == "Y" .OR. y_n == "y") THEN
            SELECT CASE (all_new)
                CASE ("M", "m")
                    M_loop: DO
                        WRITE (*,'("Enter the mass (in solar units)       :  Msolar = ")', ADVANCE = 'NO')
                        READ  (*,*, iostat = ios) Msolar
                        IF (ios == 0 .AND. Msolar > 0) THEN
                            Ms = Msolar*M_Sun
                            EXIT M_loop
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO M_loop
                    all_new = "n"
                CASE ("T", "t")
                    T_loop: DO
                        WRITE (*,'("Enter the effective temperature (in K):  Teff   = ")', ADVANCE = 'NO')
                        READ  (*,*, iostat = ios) Teff
                        IF (ios == 0 .AND. Teff > 0) THEN
                            Rs = SQRT(Ls/(four_pi*sigma*Teff**4))   #Eq. (3.17)
                            Rsolar = Rs/R_Sun                       #Solar radius from SI value
                            EXIT T_loop
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO T_loop
                    all_new = "n"
                CASE ("L", "l")
                    L_loop: DO
                        WRITE (*,'("Enter the luminosity (in solar units) :  Lsolar = ")', ADVANCE = 'NO')
                        READ  (*,*, iostat = ios) Lsolar
                        IF (ios == 0 .AND. Lsolar > 0) THEN
                            Ls = Lsolar*L_Sun
                            Rs = SQRT(Ls/(four_pi*sigma*Teff**4))   #Eq. (3.17)
                            Rsolar = Rs/R_Sun                       #Solar radius from SI value
                            EXIT L_loop
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO L_loop
                    all_new = "n"
                CASE ("X", "x")
                    X_loop: DO
                        WRITE (*,'("Enter the mass fraction of hydrogen   :  X      = ")', ADVANCE = 'NO')
                        READ  (*,*, iostat = ios) X
                        IF (ios == 0) THEN
                            IF (X >= 0 .AND. X <= 1) THEN
                                EXIT X_loop
                            ELSE
                                WRITE (*,'(/,"0 <= X <= 1 is required")')
                            END IF
                        ELSE
                            WRITE (*,'(/)')
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO X_loop
                    X_again: DO
                        Y = Helium(X, Z)
                        IF (Y >= 0 .AND. Y <= 1) EXIT X_again
                        WRITE (*,'("Note that 0 <= X + Z < = 1 is required")')
                        X_repeat: DO
                            WRITE (*,'("  Please try again:                      X      = ")', ADVANCE = 'NO')
                            READ (*,*, iostat = ios) X
                            IF (ios == 0) THEN
                                IF (X >= 0 .AND. X <= 1) THEN
                                    EXIT X_repeat
                                ELSE
                                    WRITE (*,'(/,"0 <= X <= 1 is required")')
                                END IF
                            ELSE
                                WRITE (*,'(/)')
                            END IF
                            WRITE (*,'("Invalid value entered",/)')
                        END DO X_repeat
                    END DO X_again
                    all_new = "n"
                CASE ("Z", "z")
                    Z_loop: DO
                        WRITE (*,'("Enter the mass fraction of metals     :  Z      = ")', ADVANCE = 'NO')
                        READ  (*,*, iostat = ios) Z
                        IF (ios == 0) THEN
                            IF (Z >= 0 .AND. Z <= 1) THEN
                                EXIT Z_loop
                            ELSE
                                WRITE (*,'(/,"0 <= Z <= 1 is required")')
                            END IF
                        ELSE
                            WRITE (*,'(/)')
                        END IF
                        WRITE (*,'("Invalid value entered - please try again",/)')
                    END DO Z_loop
                    Z_again: DO
                        Y = Helium(X, Z)
                        IF (Y >= 0 .AND. Y <= 1) EXIT Z_again
                        WRITE (*,'("Note that 0 <= X + Z < = 1 is required")')
                        Z_repeat: DO
                            WRITE (*,'("  Please try again:                      Z      = ")', ADVANCE = 'NO')
                            READ (*,*, iostat = ios) Z
                            IF (ios == 0) THEN
                                IF (Z >= 0 .AND. Z <= 1) THEN
                                    EXIT Z_repeat
                                ELSE
                                    WRITE (*,'(/,"0 <= Z <= 1 is required")')
                                END IF
                            ELSE
                                WRITE (*,'(/)')
                            END IF
                            WRITE (*,'("Invalid value entered",/)')
                        END DO Z_repeat
                    END DO Z_again
                    all_new = "n"
                CASE ("E", "e")
                    y_n = "N"
                CASE ("A", "a")
                    y_n = "Y"
                CASE DEFAULT
                    all_new = "y"
            END SELECT
        ELSE
            all_new = "E"           #Exit calculations
        END IF
    END SUBROUTINE Change_Model
#---------------------------------------------------------------------

    SUBROUTINE Too_Many_Zones(i, Msolar, Ms, M_r, Lsolar, Ls, L_r, r, Rs, Rsolar, Teff, X, Y, Z, P_0, T_0, &
                  & rho_0, kappa_0, epsilon_0, rc_flag)

#       General Description:
#       ====================
#           Tell user that the maximum number of zones has been exceeded

        IMPLICIT NONE
        INTEGER,        INTENT(IN)  ::  i
        REAL(dp),       INTENT(IN)  ::  Msolar, Lsolar, Rsolar
        REAL(dp),       INTENT(IN)  ::  Rs, Ms, Ls, Teff, X, Y, Z
        REAL(dp),       INTENT(IN)  ::  r, M_r, L_r
        REAL(dp),       INTENT(IN)  ::  P_0, T_0, rho_0, kappa_0, epsilon_0
        CHARACTER(1),   INTENT(IN)  ::  rc_flag

        WRITE (*,'(/,"The maximum number of zones has been exceeded for this model - Sorry#")')
        WRITE (*,'("The conditions at the time the model was terminated were: ")')
        WRITE (*,'(T5,"Surface Conditions:",           T35, "Last Zone Calculated:")')
        WRITE (*,'(T5,"-------------------",           T35, "---------------------")')
        WRITE (*,'(T5,"M    = ", F12.6, " solar",      T35, "M_r/Ms  = ", F12.6)')                 Msolar, M_r/Ms
        WRITE (*,'(T5,"Teff = ", F12.6, " K",          T35, "L_r/LS  = ", F12.6)')                 Teff,   L_r/Ls
        WRITE (*,'(T5,"L    = ", F12.6, " solar",      T35, "r/Rs    = ", F12.6)')                 Lsolar, r/Rs
        WRITE (*,'(T5,"R    = ", F12.6, " solar",      T35, "P       = ", ES12.5, " N/m^2")')      Rsolar, P_0
        WRITE (*,'(T5,"X    = ", F12.6,                T35, "T       = ", ES12.5, " K")')          X,      T_0
        WRITE (*,'(T5,"Y    = ", F12.6,                T35, "rho     = ", ES12.5, " kg/m^3")')     Y,      rho_0
        WRITE (*,'(T5,"Z    = ", F12.6,                T35, "kappa   = ", ES12.5, " m^2/kg")')     Z,      kappa_0
        WRITE (*,'(                                    T35, "epsilon = ", ES12.5, " W/kg")')               epsilon_0
        IF (rc_flag == "r") THEN
            WRITE (*,'(T35, "The core is RADIATIVE")')
        ELSE
            WRITE (*,'(T35, "The core is CONVECTIVE")')
        END IF
    END SUBROUTINE Too_Many_Zones
#---------------------------------------------------------------------

    SUBROUTINE Final_Results(i, Msolar, Ms, M_r, Lsolar, Ls, L_r, r, Rs, Rsolar, Teff, X, Y, Z, &
                 & P, T, rho, kappa, epsilon, P_0, T_0, rho_0, kappa_0, epsilon_0, rc_flag)

#       General Description:
#       ====================
#           Tell the user the conditions at the surface and core of the completed model

        IMPLICIT NONE
        INTEGER,        INTENT(IN)  ::  i
        REAL(dp),       INTENT(IN)  ::  Msolar, Lsolar, Rsolar
        REAL(dp),       INTENT(IN)  ::  Ms, Ls, Rs, Teff, X, Y, Z
        REAL(dp),       INTENT(IN)  ::  r, M_r, L_r
        REAL(dp),       INTENT(IN)  ::  p, T, rho, kappa, epsilon
        REAL(dp),       INTENT(IN)  ::  P_0, T_0, rho_0, kappa_0, epsilon_0
        CHARACTER(1),   INTENT(IN)  ::  rc_flag

        WRITE (*,'(//,"*********************THE INTEGRATION HAS BEEN COMPLETED*********************")')
        WRITE (*,'(T5,"Surface Conditions:",           T35, "Core Conditions:")')
        WRITE (*,'(T5,"-------------------",           T35, "----------------")')
        WRITE (*,'(T5,"M    = ", F12.6, " solar",      T35, "M_r/Ms  = ", F12.6)')                 Msolar, M_r/Ms
        WRITE (*,'(T5,"Teff = ", F12.6, " K",          T35, "L_r/LS  = ", F12.6)')                 Teff,   L_r/Ls
        WRITE (*,'(T5,"L    = ", F12.6, " solar",      T35, "r/Rs    = ", F12.6)')                 Lsolar, r/Rs
        WRITE (*,'(T5,"R    = ", F12.6, " solar",      T35, "P       = ", ES12.5, " N/m^2")')      Rsolar, P_0
        WRITE (*,'(T5,"X    = ", F12.6,                T35, "T       = ", ES12.5, " K")')          X,      T_0
        WRITE (*,'(T5,"Y    = ", F12.6,                T35, "rho     = ", ES12.5, " kg/m^3")')     Y,      rho_0
        WRITE (*,'(T5,"Z    = ", F12.6,                T35, "kappa   = ", ES12.5, " m^2/kg")')     Z,      kappa_0
        WRITE (*,'(                                    T35, "epsilon = ", ES12.5, " W/kg")')               epsilon_0
        IF (rc_flag == "r") THEN
            WRITE (*,'(T35, "The core is RADIATIVE")')
        ELSE
            WRITE (*,'(T35, "The core is CONVECTIVE")')
        END IF

        WRITE (*,'(/,"For your information, the conditions in the last zone above the core are:")')
        WRITE (*,'(T35, "P       = ", ES12.5, " N/m^2")') P
        WRITE (*,'(T35, "T       = ", ES12.5, " K")') T
        WRITE (*,'(T35, "rho     = ", ES12.5, " kg/m^3")') rho
        WRITE (*,'(T35, "kappa   = ", ES12.5, " m^2/kg")') kappa
        WRITE (*,'(T35, "epsilon = ", ES12.5, " W/kg")') epsilon

        WRITE (*,'(/,"The number of mass shells in this model: ", I5)') i
        WRITE (*,'("The details of the model are available in ZAMSmodel.txt")')
    END SUBROUTINE Final_Results
#
#   General Description:
#   ====================
#       This module holds required data from the previous and current zones
#
#---------------------------------------------------------------------

    USE Constants, ONLY :   dp

    IMPLICIT NONE
    PRIVATE             ::  dp

#   Previous zone data
    REAL(dp)                                ::  Mm, Lm, rm
    REAL(dp)                                ::  Pm, Tm
    REAL(dp)                                ::  Xm, Zm
    REAL(dp)                                ::  rhom, kappam, taum, epsilonm

#   Current zone data
    REAL(dp)                                ::  rho, kappa, tau, epsilon, gamma, dlnPdlnT
    CHARACTER(1)                            ::  rc_flag

#   Number of stellar structure equations
    INTEGER,                    PARAMETER   ::  n = 4

#   Current step size flag
    INTEGER                                 ::  step_size_condition

#   The first derivatives from the stellar structure equations to be used by Runge Kutta routines
    REAL(dp),   DIMENSION(n)                ::  dfdr0
#   General Description:
#   ====================
#       This module contains the physics routines required for the construction of stellar models.
#
#---------------------------------------------------------------------

    USE Constants, ONLY :   dp, a_o3 => a_rad_o3, k => k_B, m_H
    PRIVATE             ::  dp, a_o3, k, m_H

    CONTAINS
    REAL(dp) FUNCTION PTgradient(Pm, P, Tm, T) RESULT(dlnPdlnT)

#       General Description:
#       ====================
#           Compute the pressure gradient with respect to temperature to determine whether convection
#           is required. Limit value of dlnPdlnT for output purposes.

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  Pm, P, Tm, T

        dlnPdlnT = ((Tm + T)/(Pm + P))*((Pm - P)/(Tm - T))
        IF (dlnPdlnT > 99.9) dlnPdlnT = 99.9
    END FUNCTION PTgradient
#---------------------------------------------------------------------

    REAL(dp) FUNCTION Specific_Heat_Ratio() RESULT(Gamma)

#       General Description:
#       ====================
#           Compute the ratio C_P/C_V

        IMPLICIT NONE
        REAL(dp),   PARAMETER   ::  monatomic = 5/3.0_dp

        gamma = monatomic                               #Assume a purely monatomic gas, Eq. (10.80)
    END FUNCTION Specific_Heat_Ratio
#---------------------------------------------------------------------

    REAL(dp) FUNCTION Density(T, P, mu) RESULT(rho)

#       General Description:
#       ====================
#           Density computes the density of the gas, assuming the ideal gas law and radiation pressure
#           A negative value for the density indicates that an error was detected in the routine

        USE Zone_Quantities, ONLY : step_size_condition

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  T, P, mu
        REAL(dp)                ::  P_gas

        P_gas = P - a_o3*T**4                           #Eq. (10.20)
        IF (P_gas <= 0 .AND. T > 0) THEN                #Do something desperate
            SELECT CASE (step_size_condition)
                CASE (0)
                    P_gas = P
                CASE (1)
                    P_gas = 0.001*P
                CASE (2)
                    P_gas = 0.0001*P
            END SELECT
        END IF

        IF (T > 0 .AND. P_gas > 0) THEN
            rho = P_gas*mu*m_H/(k*T)                    #Eq. (10.11)
        ELSE
            rho = -1
        END IF
        IF (rho < 0) THEN
            WRITE (*,'("A negative density was computed#", /, &
                     & "Sorry but I am not programmed to handle this new physics :-)", /, &
                     & "Terminating calculation with: ", /, &
                     & "         T     = ", ES13.6, /, &
                     & "         P     = ", ES13.6, /, &
                     & "         P_gas = ", ES13.6)') T, P, P_gas
        END IF
    END FUNCTION Density
#---------------------------------------------------------------------

    REAL(dp) FUNCTION Opacity(T, rho, X, Z) RESULT(kappa)

#       General Description:
#       ====================
#           Opacity computes an approximation of the Rosseland Mean Opacity, based on approximation formulae

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  T, rho, X, Z
        REAL(dp)                ::  kappa_bf, kappa_ff, kappa_es, kappa_Hminus
        REAL(dp)                ::  tog_bf
        REAL(dp),   PARAMETER   ::  g_ff = 1                    #the free-free Gaunt factor is on the order of unity
        REAL(dp),   PARAMETER   ::  A_bf = 4.34E21, A_ff = 3.68E18, A_es = 0.02, A_Hm = 7.9E-34

        tog_bf = 0.708*(rho*(1 + X))**0.2                       #Taken from Novotny (1973), p. 469

        kappa_bf = (A_bf/tog_bf)*Z*(1 + X)*rho/T**3.5           #Eq. (9.22)
        kappa_ff = A_ff*g_ff*(1 - Z)*(1 + X)*rho/T**3.5         #Eq. (9.23)
        kappa_es = A_es*(1 + X)                                 #Eq. (9.27)

        IF ((T > 3000 .AND. T < 6000) .AND. (rho > 1E-10 .AND. rho < 1E-5) .AND. (Z > 0.001 .AND. Z < 0.03)) THEN
            kappa_Hminus = A_Hm*(Z/0.02)*SQRT(rho)*T**9         #Eq. (9.28)
        ELSE
            kappa_Hminus = 0
        END IF

        kappa = kappa_bf + kappa_ff + kappa_es + kappa_Hminus
    END FUNCTION Opacity
#---------------------------------------------------------------------

    REAL(dp) FUNCTION Optical_Depth_Change(kappa, kappam, rho, rhom, r, rm) RESULT(dtau)

#       General Description:
#       ====================
#           Compute the change in optical depth across the zone

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  kappa, kappam, rho, rhom, r, rm

        dtau = -(kappa*rho + kappam*rhom)*(r - rm)/2            #Eq. (9.15)
    END FUNCTION Optical_Depth_Change
#---------------------------------------------------------------------

    REAL(dp) FUNCTION Nuclear(T, rho, X, Z) RESULT(epsilon)

#       General Description:
#       ====================
#           Nuclear computes the nuclear energy generation rates for the proton-proton chains, the CNO cycle,
#           and helium burning.

        USE Composition

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  T, rho, X, Z
        REAL(dp)                ::  psipp, Cpp, XCNO, CCNO, Y
        REAL(dp),   PARAMETER   ::  fpp = 1, f3a = 1                                 #screening factors
        REAL(dp)                ::  eps_pp, eps_CNO, eps_He
        REAL(dp)                ::  T6, T8
        REAL(dp),   PARAMETER   ::  onethird = 1/3.0_dp, twothirds = 2*onethird
        REAL(dp),   PARAMETER   ::  fourthirds = 4*onethird, fivethirds = 5*onethird
        REAL(dp),   PARAMETER   ::  A_pp = 0.241, A_CNO = 8.67E20, A_He = 50.9      #reaction rate coefficients

        T6 = T*1.0E-06
        T8 = T*1.0E-08

#       PP chains (see Hansen and Kawaler, Eq. 6.65, 6.73, and 6.74)
        psipp = 1 + 1.412E8*(1/X - 1)*EXP(-49.98*T6**(-onethird))
        Cpp = 1 + 0.0123*T6**onethird + 0.0109*T6**twothirds + 0.000938*T6
        eps_pp = A_pp*rho*X*X*fpp*psipp*Cpp*T6**(-twothirds)*EXP(-33.80*T6**(-onethird))    #Eq. (10.46)

#       CNO cycle (Kippenhahn and Weigert, Eq. 18.65)
        XCNO = CNO(Z)
        CCNO = 1 + 0.0027*T6**onethird - 0.00778*T6**twothirds - 0.000149*T6
        eps_CNO = A_CNO*rho*X*XCNO*CCNO*T6**(-twothirds)*EXP(-152.28*T6**(-onethird))       #Eq. (10.58)

#       Helium burning (Kippenhahn and Weigert, Eq. 18.67)
        Y = Helium(X, Z)
        eps_He = A_He*rho**2*Y**3/T8**3*f3a*EXP(-44.027/T8)                                 #Eq. (10.62)

#       Combined energy generation rate
        epsilon = eps_pp + eps_CNO + eps_He
    END FUNCTION Nuclear
#
#   General Description:
#   ====================
#
#       Module ODE_Integrator contains the numerical integration routine used to integrate a set
#       of n first order linear  differential equations that depend on one independent variable, x.
#
#       This module accepts n ODEs with vectors of initial conditions.
#       A user-defined EXTERNAL function that returns the derivative of one of n specified functions is
#       required to have the form:
#
#           REAL(8) FUNCTION Deriv(i, x, y, ok)
#
#       where
#           (a)  all function derivative definitions are included in the single routine,
#           (b)  i designates which function value is to be returned
#           (c)  x is a scalar representing the independent variable
#           (d)  y is a vector of dimension n representing the set of dependent variables at x
#           (e)  ok is a logical flag.  If ok == .TRUE. the derivative was computed successfully.
#
#---------------------------------------------------------------------

    USE Constants,   ONLY       :   dp
    PRIVATE                     ::  dp

    CONTAINS
    SUBROUTINE RK_4(n, x0, h, y0, y4, f0, f, ok)

#       General Description:
#       ====================
#           This routine uses the fourth order Runge Kutta scheme to move the solution forward
#           by one step size, h.
#
#           For details of the method, see
#
#               Press, Teukolsky, Vetterling, and Flannery, "Numerical Recipes in Fortran 77:  The Art
#                   of Scientific Computing", second edition, Cambridge University Press, Cambridge, 1996.
#
#---------------------------------------------------------------------

        IMPLICIT NONE
        INTEGER,                    INTENT(IN)  ::  n       #Number of ODEs
        REAL(dp),                   INTENT(IN)  ::  x0, h   #The independent variable and step size
        REAL(dp),   DIMENSION(:),   INTENT(IN)  ::  y0      #The array of initial y values at the start of the interval
        REAL(dp),   DIMENSION(:),   INTENT(OUT) ::  y4      #The array of results to 4th order at x + h
        REAL(dp),   DIMENSION(:),   INTENT(IN)  ::  f0      #The array of first derivatives at the start of the interval
        REAL(dp),                   EXTERNAL    ::  f       #The name of the user-defined function that returns n derivatives
        LOGICAL,                    INTENT(OUT) ::  ok      #Reports if function evaluation was successful

        INTEGER                                 ::  i

#       Temporary work arrays containing values at intermediate steps
        REAL(dp),   DIMENSION(n)                ::  k1, k2, k3, k4

#---------------------------------------------------------------------

        ok = .TRUE.

#       Calcualtion intermediate derivatives using the user-defined external function
        k1 = h*f0

        DO i = 1, n
            k2(i) = h*f(i, x0 + h/2, y0 + k1/2, ok)
            IF (.NOT. ok) RETURN
        END DO

        DO i = 1, n
            k3(i) = h*f(i, x0 + h/2, y0 + k2/2, ok)
            IF (.NOT. ok) RETURN
        END DO

        DO i = 1, n
            k4(i) = h*f(i, x0 + h,   y0 + k3,   ok)
            IF (.NOT. ok) RETURN
        END DO

#       Compute the variables for the next shell using the 4th order Runge-Kutta formula
        y4 = y0 + k1/6 + k2/3 + k3/3 + k4/6
    END SUBROUTINE RK_4
#
#   General Description:
#   ====================
#
#       This module contains the most up-to-date physical and
#       astronomical constants in SI units.  This module also identifies
#       the correct kind parameters for the current machine.
#
#       "An Introduction to Modern Astrophysics", Appendix I
#       Bradley W. Carroll and Dale A. Ostlie
#       Addison Wesley, 2007
#
#       Weber State University
#       Ogden, UT
#       modastro@weber.edu
#-------------------------------------------------------------------

    IMPLICIT NONE

#The system's precision and range
    INTEGER,    PARAMETER   ::  sp          = SELECTED_REAL_KIND(p = 6, r = 37)
    INTEGER,    PARAMETER   ::  dp          = SELECTED_REAL_KIND(p = 15, r = 307)
    INTEGER,    PARAMETER   ::  qp          = SELECTED_REAL_KIND(p = 33, r = 4931)

    INTEGER,    PARAMETER   ::  i1          = SELECTED_INT_KIND(2)
    INTEGER,    PARAMETER   ::  i2          = SELECTED_INT_KIND(4)
    INTEGER,    PARAMETER   ::  i4          = SELECTED_INT_KIND(8)
    INTEGER,    PARAMETER   ::  i8          = SELECTED_INT_KIND(15)

#The smallest non-zero number and the number of significant figures
    REAL(sp),   PARAMETER   ::  tiny_sp     = TINY(1.0_sp)
    REAL(dp),   PARAMETER   ::  tiny_dp     = TINY(1.0_dp)
    REAL(qp),   PARAMETER   ::  tiny_qp     = TINY(1.0_qp)
    INTEGER,    PARAMETER   ::  sig_fig_sp  = PRECISION(1.0_sp)
    INTEGER,    PARAMETER   ::  sig_fig_dp  = PRECISION(1.0_dp)
    INTEGER,    PARAMETER   ::  sig_fig_qp  = PRECISION(1.0_qp)
    REAL(sp),   PARAMETER   ::  eps_sp      = 10.0_sp**(-sig_fig_sp)
    REAL(dp),   PARAMETER   ::  eps_dp      = 10.0_dp**(-sig_fig_dp)
    REAL(qp),   PARAMETER   ::  eps_qp      = 10.0_qp**(-sig_fig_qp)

#The largest number for given precision
    REAL(sp),   PARAMETER   ::  biggest_sp  = HUGE(1.0_sp)
    REAL(dp),   PARAMETER   ::  biggest_dp  = HUGE(1.0_dp)
    REAL(qp),   PARAMETER   ::  biggest_qp  = HUGE(1.0_qp)
    INTEGER(i1),PARAMETER   ::  biggest_i1  = HUGE(1_i1)
    INTEGER(i2),PARAMETER   ::  biggest_i2  = HUGE(1_i2)
    INTEGER(i4),PARAMETER   ::  biggest_i4  = HUGE(1_i4)
    INTEGER(i8),PARAMETER   ::  biggest_i8  = HUGE(1_i8)

#Values related to pi and e
    INTEGER,    PARAMETER,  PRIVATE ::  rpi = SELECTED_REAL_KIND(p = 33, r = 2)
    REAL(rpi),  PARAMETER   ::  pi          = 3.14159265358979323846264338327950_rpi
    REAL(rpi),  PARAMETER   ::  two_pi      = 2*pi
    REAL(rpi),  PARAMETER   ::  four_pi     = 4*pi
    REAL(rpi),  PARAMETER   ::  four_pi_o3  = four_pi/3
    REAL(rpi),  PARAMETER   ::  pi_over_2   = pi/2

    REAL(rpi),  PARAMETER   ::  natural_e   = 2.71828182845904523536028747135266_rpi

#Conversions for radians to degrees and degrees to radians
    REAL(rpi),  PARAMETER   ::  degrees_to_radians = pi/180
    REAL(rpi),  PARAMETER   ::  radians_to_degrees = 180/pi

#Physical constants
    INTEGER,    PARAMETER,  PRIVATE :: rG   = SELECTED_REAL_KIND(p = 4, r = 11)
    REAL(rG),   PARAMETER   ::  G           = 6.673e-11_rG
    REAL(qp),   PARAMETER   ::  c           = 2.99792458e08_qp
    REAL(qp),   PARAMETER   ::  mu_0        = four_pi*1e-07_qp
    REAL(qp),   PARAMETER   ::  epsilon_0   = 1/(mu_0*c**2)

    INTEGER,    PARAMETER, PRIVATE  :: reC  = SELECTED_REAL_KIND(p = 10, r = 19)
    REAL(reC),  PARAMETER   ::  e_C         = 1.602176462e-19_reC
    REAL(reC),  PARAMETER   ::  eV          = e_C
    REAL(reC),  PARAMETER   ::  keV         = eV*1.0e3_reC
    REAL(reC),  PARAMETER   ::  MeV         = eV*1.0e6_reC
    REAL(reC),  PARAMETER   ::  GeV         = eV*1.0e9_reC

    INTEGER,    PARAMETER,  PRIVATE :: rh   = SELECTED_REAL_KIND(p = 9, r = 34)
    REAL(rh),   PARAMETER   ::  h           = 6.62606876e-34_rh
    REAL(rh),   PARAMETER   ::  hbar        = h/two_pi

    INTEGER,    PARAMETER,  PRIVATE :: rkB  = SELECTED_REAL_KIND(p = 8, r = 23)
    REAL(rkB),  PARAMETER   ::  k_B         = 1.3806503e-23_rkB

    INTEGER,    PARAMETER,  PRIVATE :: rsig = SELECTED_REAL_KIND(p = 8, r = 25)
    REAL(rsig), PARAMETER   ::  sigma       = 2*pi**5*k_B**4/(15*c**2*h**3)
    REAL(rsig), PARAMETER   ::  a_rad       = 4*sigma/c
    REAL(rsig), PARAMETER   ::  a_rad_o3    = a_rad/3
    REAL(rsig), PARAMETER   ::  four_ac_o3  = 4*a_rad_o3*c

    INTEGER,    PARAMETER, PRIVATE  :: rme  = SELECTED_REAL_KIND(p = 9, r = 31)
    REAL(rme),  PARAMETER   ::  m_e         = 9.10938188e-31_rme

    INTEGER,    PARAMETER,  PRIVATE :: rmp  = SELECTED_REAL_KIND(p = 9, r = 27)
    REAL(rmp),  PARAMETER   ::  m_p         = 1.67262158e-27_rmp

    INTEGER,    PARAMETER, PRIVATE  :: rmn  = SELECTED_REAL_KIND(p = 9, r = 27)
    REAL(rmn),  PARAMETER   ::  m_n         = 1.67492716e-27_rmn

    INTEGER,    PARAMETER, PRIVATE  :: rmH  = SELECTED_REAL_KIND(p = 10, r = 27)
    REAL(rmH),  PARAMETER   ::  m_H         = 1.673532499e-27_rmH

    INTEGER,    PARAMETER, PRIVATE  :: ru   = SELECTED_REAL_KIND(p = 9, r = 27)
    REAL(ru),   PARAMETER   ::  u           = 1.66053873e-27_ru

    INTEGER,    PARAMETER, PRIVATE  :: rNA  = SELECTED_REAL_KIND(p = 9, r = 23)
    REAL(rNA),  PARAMETER   ::  N_A         = 6.02214199e23_rNA

    INTEGER,    PARAMETER, PRIVATE  :: rR   = SELECTED_REAL_KIND(p = 7, r = 1)
    REAL(rR),   PARAMETER   ::  R_gas       = 8.314472_rR

    INTEGER,    PARAMETER, PRIVATE  :: ra0  = SELECTED_REAL_KIND(p = 10, r = 11)
    REAL(ra0),  PARAMETER   ::  a_0         = four_pi*epsilon_0*hbar**2/(m_e*e_C**2)

    INTEGER,    PARAMETER, PRIVATE  :: rRH  = SELECTED_REAL_KIND(p = 14, r = 7)
    REAL(rRH),  PARAMETER   ::  R_infty     = m_e*e_C**4/(64*pi**3*epsilon_0**2*hbar**3*c)
    REAL(rRH),  PARAMETER   ::  R_H         = m_p/(m_e + m_p)*R_infty

#Time constants
    INTEGER(i2),PARAMETER   ::  hr          = 3600
    INTEGER(i4),PARAMETER   ::  day         = 24*hr
    REAL(qp),   PARAMETER   ::  J_yr        = 365.25_qp*day

    INTEGER,    PARAMETER, PRIVATE  :: ryr  = SELECTED_REAL_KIND(p = 9, r = 7)
    REAL(ryr),  PARAMETER   ::  yr          = 3.15581450e7_ryr

    INTEGER,    PARAMETER, PRIVATE  :: rTyr = SELECTED_REAL_KIND(p = 10, r = 7)
    REAL(rTyr), PARAMETER   ::  T_yr        = 3.155692519e7_rTyr

    INTEGER,    PARAMETER, PRIVATE  :: rGyr = SELECTED_REAL_KIND(p = 8, r = 7)
    REAL(rGyr), PARAMETER   ::  G_yr        = 3.1556952e7_rGyr

#Astronomical length constants
    INTEGER,    PARAMETER, PRIVATE  :: rAU  = SELECTED_REAL_KIND(p = 11, r = 11)
    REAL(rAU),  PARAMETER   ::  AU          = 1.4959787066e11_dp

    INTEGER,    PARAMETER, PRIVATE  :: rpc  = SELECTED_REAL_KIND(p = 11, r = 17)
    REAL(rpc),  PARAMETER   ::  pc          = 206264.806_rpc*AU

    REAL(qp),   PARAMETER   ::  ly          = c*J_yr

#Solar constants
    INTEGER,    PARAMETER, PRIVATE  :: rMs  = SELECTED_REAL_KIND(p = 5, r = 30)
    REAL(rMs),  PARAMETER   ::  M_Sun       = 1.9891e30_rMs

    INTEGER,    PARAMETER, PRIVATE  :: rSs  = SELECTED_REAL_KIND(p = 4, r = 3)
    REAL(rSs),  PARAMETER   ::  S_Sun       = 1.365e3_rSs

    INTEGER,    PARAMETER, PRIVATE  :: rLs  = SELECTED_REAL_KIND(p = 4, r = 26)
    REAL(rLs),  PARAMETER   ::  L_Sun       = four_pi*AU**2*S_Sun

    INTEGER,    PARAMETER, PRIVATE  :: rRs  = SELECTED_REAL_KIND(p = 6, r = 6)
    REAL(rMs),  PARAMETER   ::  R_Sun       = 6.95508e8_rRs

    INTEGER,    PARAMETER, PRIVATE  :: rTs  = SELECTED_REAL_KIND(p = 4, r = 26)
    REAL(rTs),  PARAMETER   ::  Te_Sun      = (L_Sun/(four_pi*R_Sun**2*sigma))**0.25_qp

#Solar magnitudes
    REAL(sp),   PARAMETER   ::  Mbol_Sun    =   4.74
    REAL(sp),   PARAMETER   ::  MU_Sun      =   5.67
    REAL(sp),   PARAMETER   ::  MB_Sun      =   5.47
    REAL(sp),   PARAMETER   ::  MV_Sun      =   4.82
    REAL(sp),   PARAMETER   ::  Mbol_Sun_ap = -26.83
    REAL(sp),   PARAMETER   ::  MU_Sun_ap   = -25.91
    REAL(sp),   PARAMETER   ::  MB_Sun_ap   = -26.10
    REAL(sp),   PARAMETER   ::  MV_Sun_ap   = -26.75
    REAL(sp),   PARAMETER   ::  BC_Sun      =  -0.08

#Earth constants
    INTEGER,    PARAMETER, PRIVATE  :: rMea = SELECTED_REAL_KIND(p = 5, r = 24)
    REAL(rMea), PARAMETER   ::  M_Earth     = 5.9736e24_rMea

    INTEGER,    PARAMETER, PRIVATE  :: rRea = SELECTED_REAL_KIND(p = 7, r = 6)
    REAL(rRea), PARAMETER   ::  R_Earth     = 6.378136e6_rRea

#Unit Conversions
    REAL(sp),   PARAMETER   ::  cm          = 1e-2
    REAL(sp),   PARAMETER   ::  gram        = 1e-3
    REAL(sp),   PARAMETER   ::  erg         = 1e-7
    REAL(sp),   PARAMETER   ::  dyne        = 1e-5
    REAL(dp),   PARAMETER   ::  esu         = 3.335640952e-10
    REAL(dp),   PARAMETER   ::  statvolt    = 2.997924580e2
    REAL(sp),   PARAMETER   ::  gauss       = 1e-4
    REAL(sp),   PARAMETER   ::  angstrom    = 1e-10
    REAL(sp),   PARAMETER   ::  jansky      = 1e-26
#
#   General Description:
#   ====================
#       This module contains the information about the composition of the gas
#
#---------------------------------------------------------------------

    USE Constants, ONLY :   dp
    PRIVATE             ::  dp

    CONTAINS
    REAL(dp) FUNCTION Mean_Molecular_Weight(X, Y, Z)   RESULT(mu)

#       General Description:
#       ====================
#           Calculate the mean molecular weight of the gas

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  X, Y, Z

        mu = 1/(2*X + 3*Y/4 + Z/2)           #Assume complete ionization, Eq. (10.16)
    END FUNCTION Mean_Molecular_Weight
#---------------------------------------------------------------------

    REAL(dp) FUNCTION Helium(X, Z) RESULT(Y)

#       General Description:
#       ====================
#           Calculate the amount of Helium-4 in the mixture

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  X, Z

        Y = 1 - X - Z
    END FUNCTION Helium
#---------------------------------------------------------------------

    REAL(dp) FUNCTION CNO(Z) RESULT(XCNO)

#       General Description:
#       ====================
#           Calculate the mass fraction of C, N, and O in the mixture

        IMPLICIT NONE
        REAL(dp),   INTENT(IN)  ::  Z

        XCNO = Z/2
    END FUNCTION CNO
#
#   General Description:
#   ====================
#       This module contains:
#           (a) the starting conditions for the surface inward integration,
#               assuming that the surface pressure, temperature, and density are all zero.
#           (b) extrapolation formuale to estimate central conditions from the last computed zone
#
#---------------------------------------------------------------------

    USE Constants, ONLY :   dp, pi, two_pi, four_pi_o3, G, a => a_rad, a_o3 => a_rad_o3, c, k => k_B, m_H
    PRIVATE             ::  dp, pi, two_pi, four_pi_o3, G, a, a_o3, c, k, m_H

    CONTAINS
    SUBROUTINE Surface(i, Ms, Ls, rm, X, Z, dr, r, P, T, M_r, L_r, rho, kappa, epsilon, good_surface)

#       General Description:
#       ====================
#           Estimate the temperature and pressure of the outermost zone from the zero boundary condition.
#           Electron scattering and H- ion contributions to the opacity are neglected for simplification.

        USE Composition
        USE Physics
        USE Stellar_Structure_Equations
        USE Zone_Quantities, ONLY   :   Mm, Lm, Pm, Tm, dlnPdlnT, gamma, rc_flag

        IMPLICIT NONE
        INTEGER,    INTENT(IN)      ::  i                                   #zone number
        REAL(dp),   INTENT(IN)      ::  Ms, Ls, rm, X, Z
        REAL(dp),   INTENT(INOUT)   ::  dr
        REAL(dp),   INTENT(OUT)     ::  r, P, T, M_r, L_r, rho, kappa, epsilon
        LOGICAL,    INTENT(OUT)     ::  good_surface

        REAL(dp)                    ::  Y, mu, Aop, tog_bf, XCNO
        REAL(dp),   PARAMETER       ::  g_ff = 1                                    #the free-free Gaunt factor is on the order of unity
        REAL(dp),   PARAMETER       ::  A_bf = 4.34E21, A_ff = 3.68E18              #Bound-free and free-free coefficients
        REAL(dp)                    ::  gamma_ratio, kPadiabatic

        REAL(dp),   PARAMETER       ::  maximum   = 1.0E-8                          #Maximum change in Ms and Ls over surface zone

        INTEGER                     ::  j
        INTEGER,    PARAMETER       ::  j_max = 50

        r = rm + dr

        Y  = Helium(X, Z)
        mu = Mean_Molecular_Weight(X, Y, Z)
        gamma = Specific_Heat_Ratio()
        gamma_ratio = gamma/(gamma - 1)

        j = 0
        outerzone: DO
#           Compute the temperature and pressure for the radiative boundary condition
            rc_flag = "r"
            T = G*Ms*(mu*m_H/(4.25*k))*(1/r - 1/rm)                                 #Eq. (L.2); radiative assumption

            IF (i < 2) THEN
                tog_bf = 0.01                                                       #Assume small value for surface
            ELSE
                tog_bf = 2.82*(rho*(1 + X))**0.2                                    #Taken from Novotny (1973), p. 469
            END IF
            Aop = (A_bf/tog_bf)*Z*(1+X) + A_ff*g_ff*(1-Z)*(1+X)                     #From Eq. (9.22) and (9.23)
            P = SQRT((1/4.25)*(16*pi/3)*(G*Ms/Ls)*(a*c*k/(Aop*mu*m_H)))*T**4.25     #Eq. (L.1)

#           If the zone is convective, recompute the adiabatic temperature and pressure
            dlnPdlnT = PTgradient(Pm, P, Tm, T)
            IF (dlnPdlnT < gamma_ratio .AND. i > 2) THEN
                rc_flag = "c"
                kPadiabatic = Pm/Tm**gamma_ratio
                T = G*Ms*(mu*m_H/(k*gamma_ratio))*(1/r - 1/rm)                      #Eq. (L.3)
                P = kPadiabatic*T**gamma_ratio                                      #Eq. (10.83)
            END IF

#           Compute remaining surface quantities
            rho = Density(T, P, mu)
            IF (rho < 0) THEN
                good_surface = .FALSE.
                EXIT outerzone
            END IF
            kappa   = Opacity(T, rho, X, Z)
            XCNO    = CNO(Z)
            epsilon = Nuclear(T, rho, X, Z)

#           Test to be sure that variations in M_r and L_r are not too large
            M_r = Mm + dMdr(r, rho)*dr
            L_r = Lm + dLdr(r, rho, epsilon)*dr
            IF (ABS((Ms - M_r)/Ms) < maximum .AND. ABS((Ls - L_r)/Ls) < maximum) THEN
                good_surface = .TRUE.
                EXIT outerzone
            END IF

#           If changes in M_r and L_r were too large, repeat with one-half the step size
            j = j + 1
            IF (j > j_max) THEN
                WRITE (*,*) "Unable to converge in SUBROUTINE Surface --- Exiting"
                good_surface = .FALSE.
                EXIT outerzone
            END IF
            dr = dr/2
            r = rm + dr
        END DO outerzone

        IF (.NOT. good_surface) THEN
            WRITE (*,'("The last values obtained by SUBROUTINE Surface were: ", /, &
                    "     M_r = ", ES13.6, "   dM_r/Ms = ", ES13.6, /, &
                    "     L_r = ", ES13.6, "   dL_r/Ls = ", ES13.6)')  M_r, (Ms - M_r)/Ms, L_r, (Ls - L_r)/Ls
        END IF
    END SUBROUTINE Surface
#---------------------------------------------------------------------

    SUBROUTINE Core(M_r, L_r, P, T, X, Z, r, P_0, T_0, rho_0, kappa_0, epsilon_0, rc_flag, dlnPdlnT, good_T)

#       General Description:
#       ====================
#           This routine extrapolates from the inner-most zone to obtain estimates of core conditions in the star

        USE Composition
        USE Physics

        IMPLICIT NONE
        REAL(dp),       INTENT(IN)  ::  M_r, L_r, P, T, X, Z, r
        REAL(dp),       INTENT(OUT) ::  P_0, T_0, rho_0, kappa_0, epsilon_0
        REAL(dp),       INTENT(OUT) ::  dlnPdlnT
        LOGICAL,        INTENT(OUT) ::  good_T
        CHARACTER(1),   INTENT(OUT) ::  rc_flag

        REAL(dp)                    ::  Y, mu, dT
        REAL(dp)                    ::  gamma
        REAL(dp),       PARAMETER   ::  converged = 1.0E-8
        INTEGER                     ::  i
        INTEGER,        PARAMETER   ::  i_max = 50

        rho_0     = M_r/(four_pi_o3*r**3)             #Average density of the central ball
        P_0       = P + (two_pi/3)*G*rho_0**2*r**2    #Central pressure, Eq. (L.4)
        epsilon_0 = L_r/M_r                           #Average energy generation rate of the central ball

#       Find core temperature by Newton-Raphson method (including radiation pressure)
        Y   = Helium(X, Z)
        mu  = Mean_Molecular_Weight(X, Y, Z)

        IF (rho_0 > 0) THEN
            i = 0
            T_0 = T
            good_T = .TRUE.
            Find_T_0: DO
                i = i + 1
                dT = -f(T_0)/dfdT(T_0)
                IF (ABS(dT/T_0) < converged) EXIT Find_T_0
                T_0 = T_0 + dT
                IF (i > i_max) THEN
                    WRITE (*,*) "Unable to converge on core temperature in SUBROUTINE Core --- Exiting"
                    good_T = .FALSE.
                    EXIT FIND_T_0
                END IF
            END DO Find_T_0
        ELSE
            T_0 = -T
            good_T = .FALSE.
        END IF

        IF (good_T) THEN
            kappa_0  = Opacity(T_0, rho_0, X, Z)
            dlnPdlnT = PTgradient(P, P_0, T, T_0)
            gamma    = Specific_Heat_Ratio()
            IF (dlnPdlnT < (gamma/(gamma - 1))) THEN
                rc_flag = "c"
            ELSE
                rc_flag = "r"
            END IF
        ELSE
            kappa_0  = -99.9
            dlnPdlnT = -99.9
            rc_flag  = "*"
        END IF

        CONTAINS
            REAL(dp) FUNCTION f(T)
                IMPLICIT NONE
                REAL(dp),   INTENT(IN)  ::  T

                f = rho_0*k*T/(mu*m_H) + a_o3*T**4 - P_0    #f = Ideal Gas Law + Radiation Pressure - core P = 0
            END FUNCTION f

            REAL(dp) FUNCTION dfdT(T)
                IMPLICIT NONE
                REAL(dp),   INTENT(IN)  ::  T

                dfdT = rho_0*k/(mu*m_H) + 4*a_o3*T**3       #df/dT
            END FUNCTION dfdT
    END SUBROUTINE Core
