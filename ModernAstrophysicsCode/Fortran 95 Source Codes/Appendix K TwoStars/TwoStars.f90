PROGRAM TwoStars
!
!   General Description:
!   ====================
!       TwoStars computes the radial velocities and light curves of a 
!       binary star system.  The code is not designed to produce 
!       high-precision, state-of-the-art solutions to the binary star 
!       problem.  Rather, this code is designed to be illustrative of 
!       the fundamental ideas discussed in
!
!           "An Introduction to Modern Astrophysics", Appendix K
!           Bradley W. Carroll and Dale A. Ostlie
!           Second Edition, Addison Wesley,   Copyright 2007.
!
!           Weber State University
!           Ogden, UT
!           modastro@weber.edu
!
!       The model used here assumes or includes the following:
!           1.  The stars are spherically symmetric, implying
!               a.	tidal distortions and centrifugal effects are 
!                   neglected.
!               b.	star spots or reflective heating are not 
!                   included.
!           2.  A simple formula for limb darkening is included.
!           3.  Elliptical orbits are allowed.
!           4.  An arbitrary inclination angle is allowed.
!           5.  An arbitrary periastron angle is allowed.
!           6.  The user may enter the velocity vector of the center 
!               of mass through space.
!
!       The plane of the sky is the (y',z') plane, with x' being along 
!       the line of sight.  The plane of the binary star orbit is the 
!       (x,y) plane, which is tilted at an inclination angle of i 
!       relative to the (y',z') plane.  Note that y' = y. The center 
!       of mass is initally assumed to be at the common origin of the 
!       two coordinate systems.  However, the center of mass can also 
!       be given an arbitrary velocity vector through space.
!
!       If phi =  0 degrees, the major axis is aligned with the x axis,
!       with periastron in +x direction.
!       If phi = 90 degrees, the major axis is aligned with y axis, 
!       with periastron in +y direction.
!
!             *         | z'     *
!              *        |       *
!               * z     |      * 
!                *      |     * <--- plane of orbit
!                 *     |    * 
!                  *    | i *
!                   *   |  *
!                    *  | *
!       x'            * |*      (y, y' out of paper)
!       <---------------+---------------
!        to            *|*              
!        observer     * | *
!                    *  |  *
!                   *   |   *
!                x *    |    *
!                 *     |     *
!                *      |      *
!               *       | <--- plane of sky
!
!       The transformation between coordinate systems is:
!           x' = z*COS(i) + x*SIN(i)
!           y' = y
!           z' = z*SIN(i) - x*COS(i)
!
!---------------------------------------------------------------

    USE Constants, ONLY     :   i1, dp, M_Sun, R_Sun, L_Sun, Mbol_Sun, G, sigma, &
                                two_pi, four_pi, AU, day, degrees_to_radians

    IMPLICIT NONE
    REAL(dp)                ::  m1, m2, T1, T2, R1, R2, L1, L2, Mbol
    REAL(dp)                ::  P, e, i, phi
    REAL(dp)                ::  a, a1, a2
    REAL(dp)                ::  mu, M, L, Lt, dS, r, v, vr, x, y 
    REAL(dp)                ::  x1, y1, x2, y2, v1r, v2r, y1p, z1p, y2p, z2p
    REAL(dp)                ::  vcmxp, vcmyp, vcmzp
    REAL(dp)                ::  L_ang, dAdt
    REAL(dp)                ::  t, dt, theta, dtheta
    REAL(dp)                ::  rS1, rS2, drS1, drS2, S1, S2, S

    REAL(dp)                ::  Mbol_max = -99999
    REAL(dp)                ::  t_max

    CHARACTER(80)           ::  file_name
    INTEGER(i1)             ::  ierr
    LOGICAL                 ::  open_file
    
    INTEGER,    PARAMETER   ::  N  = 1001       !time steps per orbit
    INTEGER,    PARAMETER   ::  Nr = 100        !surface rings
    INTEGER                 ::  j
    
    CHARACTER               ::  xpause

    INTERFACE
        SUBROUTINE Eclipse (x1, y1, R1, T1, x2, y2, R2, T2, i, dS, y1p, z1p, y2p, z2p)
            USE Constants, ONLY     :   dp
            REAL(dp), INTENT(IN)    ::  x1, y1, R1, T1, x2, y2, R2, T2
            REAL(dp), INTENT(IN)    ::  i
            REAL(dp), INTENT(OUT)   ::  dS, y1p, z1p, y2p, z2p
        END SUBROUTINE Eclipse

        REAL(dp) FUNCTION F(r_prime, R, T)
            USE Constants, ONLY     :   dp
            REAL(dp),   INTENT(IN)  ::  r_prime, R, T
        END FUNCTION F
    END INTERFACE

!---------------------------------------------------------------

    new_file: DO
        WRITE (*, '(A)', ADVANCE = 'NO')    &
            "Specify the name of your output file: "
        READ  (*,*) file_name
        OPEN  (10, file = file_name, FORM = 'FORMATTED', &
            STATUS = 'REPLACE', ACTION = 'WRITE', IOSTAT = ierr)
        IF (ierr == 0) EXIT new_file
        WRITE (*, '(/,A)')  "An error has occurred in opening&
            & the output file, please try again"
        WRITE (*, '(A)')    "  or hit CNTRL-C to quit"
    END DO new_file

        
    WRITE (*, '(//,A)')              "Enter the data for Star #1"
    WRITE (*, '(A)', ADVANCE = 'NO') "   Mass (solar masses):             "
    READ  (*,*) m1
    WRITE (*, '(A)', ADVANCE = 'NO') "   Radius (solar radii):            "
    READ  (*,*) R1
    WRITE (*, '(A)', ADVANCE = 'NO') "   Effective Temperature (K):       "
    READ  (*,*) T1

    WRITE (*,'(/,A)')                "Enter the data for Star #2"
    WRITE (*, '(A)', ADVANCE = 'NO') "   Mass (solar masses):             "
    READ  (*,*) m2
    WRITE (*, '(A)', ADVANCE = 'NO') "   Radius (solar radii):            "
    READ  (*,*) R2
    WRITE (*, '(A)', ADVANCE = 'NO') "   Effective Temperature (K):       "
    READ  (*,*) T2

    WRITE (*,'(/,A)')                "Enter the desired orbital parameters"
    WRITE (*,'(A)', ADVANCE = 'NO')  "   Orbital Period (days):           "
    READ  (*,*) P
    WRITE (*,'(A)', ADVANCE = 'NO')  "   Orbital Eccentricity:            "
    READ  (*,*) e
    WRITE (*,'(A)', ADVANCE = 'NO')  "   Orbital Inclination (deg):       "
    READ  (*,*) i
    WRITE (*,'(A)', ADVANCE = 'NO')  "   Orientation of Periastron (deg): "
    READ  (*,*) phi

    WRITE (*,'(/,A)')          "Enter the x', y', and z' components of the&
        & center of mass velocity vector:"
    WRITE (*,'(A)')     "Notes:  (1)  The plane of the sky is (y',z')"
    WRITE (*,'(A,/)')   "        (2)  If v_x' < 0, then the center of mass&
        & is blueshifted"
    
    WRITE (*,'(A)', ADVANCE = 'NO')     &
        "   v_x' (km/s)                      "
    READ  (*,*) vcmxp
    WRITE (*,'(A)', ADVANCE = 'NO')     &
        "   v_y' (km/s)                      "
    READ  (*,*) vcmyp
    WRITE (*,'(A)', ADVANCE = 'NO')     &
        "   v_z' (km/s)                      "
    READ  (*,*) vcmzp

!   Convert values to conventional SI units and radians
    m1  = m1*M_Sun
    m2  = m2*M_Sun
    M   = m1 + m2                   ! Total mass of system
    R1  = R1*R_Sun
    R2  = R2*R_Sun
    P   = P*day
    i   = i*degrees_to_radians
    phi = phi*degrees_to_radians

    vcmxp = vcmxp*1000
    vcmyp = vcmyp*1000
    vcmzp = vcmzp*1000

!   Compute the semimajor axes of the orbits
    mu = m1*m2/M                            ! Reduced mass, Eq. (2.22)
    a = (P**2*G*M/two_pi**2)**(1.0_dp/3)    ! Kepler's 3rd Law, Eq. (2.37)
    a1 = (mu/m1)*a                          ! Semimajor axis, Eq. (2.23)
    a2 = (mu/m2)*a                          ! Semimajor axis, Eq. (2.24)

    WRITE (*,'(//,A,F13.6,A)') "The semimajor axis of the reduced mass is ",a/AU, " AU"
    WRITE (*,'(A,F13.6,A)')    "   a1 = ", a1/AU, " AU"
    WRITE (*,'(A,F13.6,A,//)') "   a2 = ", a2/AU, " AU"
    
    WRITE (*,'(//,A)', ADVANCE = 'NO') "Enter any character to perform calculation:  "
    READ (*,*) xpause


    IF (a < R1 + R2) THEN
        WRITE (*,'(//,A)')      "Your two stars are in contact!"
        WRITE (*,'(A, F13.6, A)')  " R1 + R2 = ", (R1 + R2)/AU, " AU"
        WRITE (*,'(//A)')       "The spherically symmetric approximation is clearly invalid."
        WRITE (*,'(A)')         "****** Terminating Calculation ******"
    ELSE

!       Compute the luminosity of each star and the total luminosity
        L1   = four_pi*R1**2*sigma*T1**4    !Stefan-Boltzmann, Eq. (3.17)
        L2   = four_pi*R2**2*sigma*T2**4    !Star #2
        L    = L1 + L2                      !Total uneclipsed luminosity

!       Determine the energy produced by the uneclipsed, projected disks
        drS1 = R1/Nr
        drS2 = R2/Nr
        rS1  = 0
        rS2  = 0
        S1   = 0
        S2   = 0
        DO j = 1, Nr            ! Numerical integration loop (Eq. J.7)
            rS1 = rS1 + drS1
            rS2 = rS2 + drS2
            S1  = S1 + two_pi*F(rS1 - drS1/2, R1, T1)*(rS1 - drS1/2)*drS1
            S2  = S2 + two_pi*F(rS2 - drS2/2, R2, T2)*(rS2 - drS2/2)*drS2
        END DO
        S = S1 + S2

!       Initial orbit loop conditions
        t       = 0
        dt      = P/N                           ! time step
        theta   = 0
        L_ang   = mu*SQRT(G*M*a*(1 - e**2))     ! L, Eq. (2.30)
        dAdt    = L_ang/(2*mu)                  ! 2nd Law, Eq. (2.32)

!       Write output headers
        WRITE (*,'(4X,  5(5X, A))') "t/P   ", "v1r (km/s)", "v2r (km/s)", "    Mbol  ", " dS (W) "
        WRITE(10,'(4X, 11(5X, A))') "t/P   ", "v1r (km/s)", "v2r (km/s)", "    Mbol  ", " dS (W) ", &
            " y1p (AU) ", " z1p (AU) ", " y2p (AU) ", " z2p (AU) ", "ycmp (AU) ", "zcmp (AU) "

!       Reduced mass orbit loop
        DO
            r   =  a*(1 - e**2)/(1 + e*COS(theta))  ! position, Eq. (2.3)
            v   =  SQRT(G*M*(2/r - 1/a))            ! velocity, Eq. (2.36)
            vr  = -v*SIN(i)*SIN(theta + phi)        ! radial velocity
            v1r =  (mu/m1)*vr                       ! Eq. (2.23)
            v2r = -(mu/m2)*vr                       ! Eq. (2.24)

!           Determine (x,y) positions of centers of stars
            x   = r*COS(theta + phi)                    ! reduced mass
            y   = r*SIN(theta + phi)                    ! reduced mass
            x1  =  (mu/m1)*x
            y1  =  (mu/m1)*y
            x2  = -(mu/m2)*x
            y2  = -(mu/m2)*y

            CALL Eclipse (x1, y1, R1, T1, x2, y2, R2, T2, i, dS, y1p, z1p, y2p, z2p)
            Lt = L*(1 - dS/S)
            Mbol = Mbol_Sun - 5*LOG10(Lt/L_Sun)/2   ! Mbol, Eq. (3.8)
            IF (Mbol > Mbol_max) THEN
                Mbol_max = Mbol
                t_max = t
            END IF

!           Print results with radial velocities in km/s
            WRITE (*,'(4F15.6,  ES15.6)')   t/P, (v1r + vcmxp)/1000, (v2r + vcmxp)/1000, Mbol, Lt*dS/S
            WRITE(10,'(4F15.6, 7ES15.6)')   t/P, (v1r + vcmxp)/1000, (v2r + vcmxp)/1000, Mbol, Lt*dS/S, &
                (y1p + vcmyp*t)/AU, (z1p + vcmzp*t)/AU, (y2p + vcmyp*t)/AU, (z2p + vcmzp*t)/AU, vcmyp*t/AU, vcmzp*t/AU
                
            IF (t > P) EXIT

            dtheta = (2*dAdt/r**2)*dt                   ! Eq. (2.31)
            theta  = theta + dtheta
            t      = t + dt
        END DO

        WRITE (*,'(//,A,F13.6)')   "The deepest minimum in the light curve occurred at t/P = ", t_max/P
        WRITE (*,'(A,F13.6)')      "with a value of Mbol = ", Mbol_max
    END IF
    
    WRITE (*,'(//,A)', ADVANCE = 'NO') "Your calculation has finished; enter a character and hit <enter> to exit: "
    READ (*,*) xpause
    
END PROGRAM TwoStars

SUBROUTINE Eclipse (x1, y1, R1, T1, x2, y2, R2, T2, i, dS, y1p, z1p, y2p, z2p)
!
!   General Description:
!   ====================
!       This routine computes the change in observed luminosity due to
!       an eclipse.
!
!---------------------------------------------------------------

    USE Constants, ONLY     :   dp, pi

    IMPLICIT NONE
    REAL(dp),   INTENT(IN)  ::  x1, y1, R1, T1, x2, y2, R2, T2
    REAL(dp),   INTENT(IN)  ::  i
    REAL(dp),   INTENT(OUT) ::  dS, y1p, z1p, y2p, z2p

    REAL(dp)                ::  Rb, Tb, Rf, Tf
    REAL(dp)                ::  xfp, yfp, zfp, xbp, ybp, zbp
    REAL(dp)                ::  x1p, x2p
    REAL(dp)                ::  r_prime, dr_prime, theta0_prime, theta_prime
    REAL(dp)                ::  yp_dA, zp_dA, d, r_stop

!   Number of steps for r', theta' integrations
    INTEGER,    PARAMETER   ::  Nr = 100, Ntheta = 500
    REAL(dp),   PARAMETER   ::  dtheta_prime = pi/Ntheta

    INTERFACE
        SUBROUTINE Transformation (x, y, xp, yp, zp, i)
            USE Constants, ONLY     :   dp
            REAL(dp),   INTENT(IN)  ::  x, y
            REAL(dp),   INTENT(IN)  ::  i
            REAL(dp),   INTENT(OUT) ::  xp, yp, zp
        END SUBROUTINE Transformation

        REAL(dp) FUNCTION F(r_prime, R, T)
            USE Constants, ONLY     :   dp
            REAL(dp),   INTENT(IN)  ::  r_prime, R, T
        END FUNCTION F
    END INTERFACE
!---------------------------------------------------------------

    dS = 0

!   Perform coordinate transformations
    CALL Transformation (x1, y1, x1p, y1p, z1p, i)
    CALL Transformation (x2, y2, x2p, y2p, z2p, i)

!   Determine which star is in front (f) and which star is in back (b)
    IF (x1p > 0) THEN
        xfp = x1p; yfp = y1p; zfp = z1p; Rf = R1; Tf = T1
        xbp = x2p; ybp = y2p; zbp = z2p; Rb = R2; Tb = T2
    ELSE
        xfp = x2p; yfp = y2p; zfp = z2p; Rf = R2; Tf = T2
        xbp = x1p; ybp = y1p; zbp = z1p; Rb = R1; Tb = T1
    END IF
    
!   Are the two stars close enough for an eclipse?
    d = SQRT((yfp - ybp)**2 + (zfp - zbp)**2)           ! Eq. (J.4)
    IF (d <= Rf + Rb) THEN

!       Find the angle between y' and the projected line between the 
!       centers of the stars in the (y',z') plane.  The polar coordinate
!       integration will be centered on this line to take advantage 
!       of spherical symmetry.
        theta0_prime = ATAN2((zfp - zbp), (yfp - ybp))  ! Eq. (J.5)

!       Determine the starting radius for the integration
        IF (d < Rb - Rf) THEN
            r_prime = d + Rf    !Foreground star disk entirely 
            r_stop  = d - Rf    !inside background star disk
            IF (r_stop < 0) r_stop = 0
        ELSE
            r_prime = Rb
            r_stop  = 0
        END IF
        dr_prime = r_prime/Nr

!       The surface integration loop
        r_loop: DO

!           Determine the limits of the angular integration 
!           for the current r_prime
            theta_prime = theta0_prime
            theta_loop: DO
                yp_dA = r_prime*COS(theta_prime + dtheta_prime) + ybp
                zp_dA = r_prime*SIN(theta_prime + dtheta_prime) + zbp
                IF (SQRT((yp_dA - yfp)**2 + (zp_dA - zfp)**2) > Rf) EXIT theta_loop         ! Eq. (J.6)

                theta_prime = theta_prime + dtheta_prime
                IF ((theta_prime - theta0_prime) > pi) EXIT theta_loop
            END DO theta_loop

!           Add the luminosity change for differential area  (Eq. J.7)
            dS = dS + 2*F(r_prime - dr_prime/2, Rb, Tb)*(r_prime - dr_prime/2)&
                *dr_prime*(theta_prime - theta0_prime)

!           Check to see that there is no remaining overlap or if center 
!           of disk has been reached
            r_prime = r_prime - dr_prime
            IF (r_prime < r_stop) EXIT r_loop
        END DO r_loop
    END IF
END SUBROUTINE Eclipse

SUBROUTINE Transformation (x, y, xp, yp, zp, i)
!
!   General Description:
!   ====================
!       This routine performs the coordinate transformation between
!       the orbital plane coordinates (x,y) and the plane of the sky
!       coordinates (xp, yp, zp), based on the angle of inclination, i.
!
!---------------------------------------------------------------

    USE Constants, ONLY     :   dp

    IMPLICIT NONE
    REAL(dp),   INTENT(IN)  ::  x, y
    REAL(dp),   INTENT(IN)  ::  i
    REAL(dp),   INTENT(OUT) ::  xp, yp, zp
!---------------------------------------------------------------

    xp =  x*SIN(i)          ! Eq. (J.1)
    yp =  y                 ! Eq. (J.2)
    zp = -x*COS(i)          ! Eq. (J.3)
END SUBROUTINE Transformation

REAL(dp) FUNCTION F(r_prime, R, T) RESULT(flux)
!
!   General Description:
!   ====================
!       This routine calculates the flux for a specific value
!       of the distance from the center of the stellar disk, r_prime.
!
!---------------------------------------------------------------

    USE Constants, ONLY     :   dp, sigma

    IMPLICIT NONE
    REAL(dp),   INTENT(IN)  ::  r_prime, R, T

    INTERFACE
        REAL(dp) FUNCTION Limb_Darkening (r_prime, R)
            USE Constants, ONLY     :   dp
            REAL(dp),   INTENT(IN)  ::  r_prime, R
        END FUNCTION Limb_Darkening
    END INTERFACE
!---------------------------------------------------------------

    flux = sigma*T**4*Limb_Darkening (r_prime, R)   ! Eq. (3.18)
END FUNCTION F

REAL(dp) FUNCTION Limb_Darkening (r_prime, R)
!
!   General Description:
!   ====================
!       This routine calculates limb darkening for a specific value 
!       of the distance from the center of the stellar disk, r_prime.
!
!       Data are due to Van Hamme, W., Astronomical Journal, 
!           Vol. 106, 1096, 1993.
!
!---------------------------------------------------------------

    USE Constants, ONLY     :   dp, sigma, pi_over_2

    IMPLICIT NONE
    REAL(dp),   INTENT(IN)  ::  r_prime, R

!   Van Hamme model #109 (solar-like)
    REAL(dp),   PARAMETER   ::  x = 0.648, y = 0.207
    REAL(dp)                ::  mu
!---------------------------------------------------------------

    mu = COS(r_prime/R*pi_over_2)
    Limb_Darkening = 1 - x*(1 - mu) - y*mu*LOG(mu)  !Van Hamme Eq. (3)
END FUNCTION Limb_Darkening
