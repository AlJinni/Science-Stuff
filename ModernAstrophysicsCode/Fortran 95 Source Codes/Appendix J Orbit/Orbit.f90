PROGRAM Orbit
!
!   General Description:
!   ====================
!       Orbit computes the orbit of a small mass about a much larger mass, 
!       or it can be considered as computing the motion of the reduced mass 
!       about the center of mass.
!
!       "An Introduction to Modern Astrophysics", Appendix J
!       Bradley W. Carroll and Dale A. Ostlie
!       Second Edition, Addison Wesley, 2007
!
!       Weber State University
!       Ogden, UT
!       modastro@weber.edu
!-------------------------------------------------------------------

    USE Constants,  ONLY    :   i1, dp, G, AU, M_Sun, pi, two_pi, yr, &
                                radians_to_degrees, eps_dp

    IMPLICIT NONE
    REAL(dp)                ::  a, aAU, e, t, dt, LoM, P
    REAL(dp)                ::  Mstrsun, Mstar, theta, dtheta, r
    INTEGER                 ::  n, k, kmax
    INTEGER(i1)             ::  ios         !I/O error flag
    REAL(dp)                ::  delta       !error range at end of period
    CHARACTER               ::  xpause

!   Open the output file
    OPEN (UNIT = 10, FILE = "Orbit.txt", STATUS = 'REPLACE', ACTION = 'WRITE', &
        IOSTAT = ios)
    IF (ios /= 0) THEN
        WRITE (*,'(" Unable to open Orbit.txt.  --- Terminating calculation")')
        STOP
    END IF

!   Write introductory information for user
    WRITE (*,'(   T5,  "Orbit computes the orbit of a small mass about a much &
        &larger mass.")')

    WRITE (*,'(//,T15, "Details of the code are described in:")')
    WRITE (*,'(   T15, "An Introduction to Modern Astrophysics")')
    WRITE (*,'(   T15, " Bradley W. Carroll and Dale A. Ostlie")')
    WRITE (*,'(   T15, "         Addison Wesley")')
    WRITE (*,'(   T15, "         copyright 2007",//)')

!   Get input from user
    WRITE (*,'(" Enter the mass of the parent star (in solar masses):  ")', &
        ADVANCE = 'NO')
    READ  (*,*)   Mstrsun

    WRITE (*,'(" Enter the semimajor axis of the orbit (in AU):        ")', &
        ADVANCE = 'NO')
    READ  (*,*)   aAU

    WRITE (*,'(" Enter the orbital eccentricity:                       ")', &
        ADVANCE = 'NO')  
    READ  (*,*)   e

!   Convert entered values to conventional SI units
    Mstar = Mstrsun*M_Sun
    a = aAU*AU

!   Calculate the orbital period in seconds using Kepler's Third Law (Eq. 2.37)
    P = SQRT(4*pi**2*a**3/(G*Mstar))

!   Convert the orbital period to years and print the result
    WRITE (*,'(/," The period of this orbit is ", F10.3, " yr")') P/yr

!   Enter the number of time steps and the time interval to be printed
    WRITE (*,'(//," Please enter the number of time steps to be calculated and the", /, &
                  " frequency with which you want time steps printed.", /, &
                /," Note that taking too large a time step during the calculation", /, &
                  " will produce inaccurate results.",/)') 
    WRITE (*,'(   " Enter the number of time steps desired for the calculation: ")', &
        ADVANCE = 'NO')
    READ  (*,*) n
    n = n + 1   !increment to include t=0 (initial) point
    WRITE (*,'(/, " How often do you want time steps to be printed?", /, &
                  "            1 = every time step", /, &
                  "            2 = every second time step", /, &
                  "            3 = every third time step", /, &
                  "                        etc.", /)')
    WRITE (*,'(   " Frequency: ")', ADVANCE = 'NO')
    READ  (*,*) kmax

!   Print header information for output file
    WRITE (10,'(T27,"Elliptical Orbit Data",//,&
                T27,"Mstar = ", F10.3, " Msun",/,&
                T27,"a     = ", F10.3, " AU",/,&
                T27,"P     = ", F10.3, " yr",/, &
                T27,"e     = ", F10.3,///,     &
                T9, "t (yr)", 12X, "r (AU)", 12X, "theta (deg)", &
                    12X, "x (AU)", 12X, "y (AU)",/, &
                T9, "------", 12X, "------", 12X, "-----------", &
                    12X, "------", 12X, "------")') &
                Mstrsun, aAU, P/yr, e

!   Initialize print counter, angle, elapsed time, and time step.
    k = 1                   !printer counter
    theta = 0               !angle from direction to perihelion (radians)
    t = 0                   !elapsed time (s)
    dt = P/(n-1)            !time step (s)
    delta = 2*eps_dp        !allowable error at end of period

!   Start main time step loop
    DO
!       Calculate the distance from the principal focus using Eq. (2.3); Kepler's First Law.
        r = a*(1 - e**2)/(1 + e*COS(theta))

!       If time to print, convert to cartesian coordinates.  Be sure to print last point also.
        IF (k == 1) THEN
            WRITE (10, '(4X, F10.3, 8X, F10.3, 13X, F10.3, 2(8X, F10.3))') &
                t/yr, r/AU, theta*radians_to_degrees, r*COS(theta)/AU, r*SIN(theta)/AU
        END IF

!       Prepare for the next time step:  Update the elapsed time.
        t = t + dt

!       Calculate the angular momentum per unit mass, L/m (Eq. 2.30).
        LoM = SQRT(G*Mstar*a*(1 - e**2))

!       Compute the next value for theta using the fixed time step by combining
!           Eq. (2.31) with Eq. (2.32), which is Kepler's Second Law.
        dtheta = LoM/r**2*dt
        theta = theta + dtheta

!       Reset the print counter if necessary
        k = k + 1
        IF (k > kmax .OR. (theta - two_pi)/two_pi > delta .OR. (t - P)/P > delta) k = 1

!       Exit the loop if the calculation has gone more than 2pi radians, or more than one orbital period.
        IF ((theta - two_pi) > dtheta/2 .OR. (t - P) > dt/2) EXIT
    END DO

    WRITE (*,'(/,"The calculation is finished and the data are in Orbit.txt")')
    
    WRITE (*,'(//,"Enter any character and press <enter> to exit:  ")', ADVANCE = 'NO')
    READ  (*,*) xpause
END PROGRAM Orbit