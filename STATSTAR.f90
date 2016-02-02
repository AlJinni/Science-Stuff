Program STATSTAR
!
! This program is listed in Appendix H of both "An Introduction
! to Modern Astrophysics," Bradley W. Carroll and Dale A. Ostlie,
! Addison-Wesley Publishing Company, copyright 1996, and "An
! Introduction to Modern Stellar Astrophysics," Dale A. Ostlie
! and Bradley W. Carroll, Addison-Wesley Publishing Company,
! copyright 1996.
!
!
! This program will calculate a static stellar model using the
! equations developed in the text.  The user is expected to supply the
! star's mass, luminosity, effective temperature, and composition
! (X and Z).  If the choices for these quantities are not consistent
! with the central boundary conditions, an error message will be
! generated and the user will then need to supply a different set of
! initial values.
!
real*8 r(999), P(999), M_r(999), L_r(999), T(999), deltar, Te
real*8 rho(999), kappa(999), epslon(999), dlPdlT(999)
real*8 X, Y, Z, XCNO, mu
real*8 Ms, Ls, Rs, T0, P0
real*8 Pcore, Tcore, rhocor, epscor, rhomax
real*8 Rsolar, Msolar, Lsolar, Qm, Rcrat, Mcrat, Lcrat
real*8 deltam, dlPlim
real*8 Rsun, Msun, Lsun
real*8 f_im1(4), f_i(4), dfdr(4)
real*8 dMdr, dPdr, dLdr, dTdr
real*8 sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad,  tog_bf, g_ff
character clim, rcf
common /cnstnt/ sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad, g_ff
!
! deltar = radius integration step
! idrflg = set size flag
!        = 0 (initial surface step size of Rs/1000.)
!        = 1 (standard step size of Rs/100.)
!        = 2 (core step size of Rs/5000.)
c
! Nstart = number of steps for which starting equations are to be used
!          (the outermost zone is assumed to be radiative)
! Nstop = maximum number of allowed zones in the star
! Igoof = final model condition flag
!       = -1 (number of zones exceeded; also the initial value)
!       =  0 (good model)
!       =  1 (core density was extreme)
!       =  2 (core luminosity was extreme)
!       =  3 (extrapolated core temperature is too low)
!       =  4 (mass became negative before center was reached)
!       =  5 (luminosity became negative before center was reached)
! X, Y, Z = mass fractions of hydrogen, helium, and metals
! T0, P0 = surface temperature and pressure (T0 = P0 = 0 is assumed)
! Ms, Ls, Rs = mass, luminosity, and radius of the star (cgs units)
!
data Nstart, Nstop, Igoof, ierr / 10, 999, -1, 0 /
data P0, T0, dlPlim / 0.0d0, 0.0d0, 99.9d0 /
data Rsun, Msun, Lsun / 6.9599d+10, 1.989d+33, 3.826d+33 /
data sigma, c, a, G, k_B, m_H, pi, gamma, tog_bf, g_ff / 5.67051d-5, 2.99792458d+10, 7.56591d-15, 6.67259d-8, 1.380658d-16, 1.673534d-24, 3.141592654d0, 1.6666667d0, 0.01d0, 1.0d0 /
!
! Assign values to constants (cgs units)
!    Rsun = radius of the Sun
!    Msun = mass of the Sun
!    Lsun = luminosity of the Sun
!    sigma = Stefan-Boltzmann constant
!    != speed of light in vacuum
!    a = 4*sigma/!(radiation pressure constant)
!    G = universal gravitational constant
!    K_B = Boltzmann constant
!    m_H = mass of hydrogen atom
!    pi = 3.141592654
!    gamma = 5/3 (adiabati!gamma for a monatomic gas)
!    gamrat = gamma/(gamma-1)
!    kPad = P/T**(gamma/(gamma-1)) (adiabatic gas law constant)
!    tog_bf = bound-free opacity constant
!            (ratio of guillotine to gaunt factors)
!    g_ff = free-free opacity gaunt factor (assumed to be unity)
!
open(unit=20,file='starmodl.dat',form='formatted',
1 status='unknown')
!
! Enter desired stellar parameters
!
write(*,*) ' Enter the mass of the star (in solar units):'
read(*,*) Msolar
write(*,*) ' Enter the luminosity of the star (in solar units):'
read(*,*) Lsolar
write(*,*) ' Enter the effective temperature of the star (in K):'
read(*,*) Te
10 continue
    write(*,*) ' Enter the mass fraction of hydrogen (X):'
    read(*,*) X
    write(*,*) ' Enter the mass fraction of metals (Z):'
    read(*,*) Z
    Y = 1.d0 - X - Z
    if (Y.lt.0.0d0) then
        write(*,100)
        go to 10
    end if
!
! Select the mass fraction CNO to be 50% of Z.
!
XCNO = Z/2.0d0
!
! Calculate the mass, luminosity, and radius of the star.
! The radius is calculated from Eq. (3.17).
!
Ms = Msolar*Msun
Ls = Lsolar*Lsun
Rs = sqrt(Ls/(4.d0*pi*sigma))/Te**2
Rsolar = Rs/Rsun
!
! Begin with a very small step size since surface conditions vary
! rapidly.
!
deltar = -Rs/1000.0d0
idrflg = 0
!
! Calculate mean molecular weight mu assuming complete ionization
! (see Eq. 10.21).
!
mu = 1.0d0/(2.0d0*X + 0.75d0*Y + 0.5d0*Z)
!
! Calculate the delimiter between adiabatic convection and radiation
! (see Eq. 10.87).
!
gamrat = gamma/(gamma - 1.0d0)
!
! Initialize values of r, P, M_r, L_r, T, rho, kappa, and epslon at
! the surface.  The outermost zone is assumed to be zone 1.  The zone
! number increases toward the center.
!
r(1)   = Rs
M_r(1) = Ms
L_r(1) = Ls
T(1)   = T0
P(1)   = P0
if (P0.le.0.0d0 .or. T0.le.0.0d0) then
    rho(1)    = 0.0d0
    kappa(1)  = 0.0d0
    epslon(1) = 0.0d0
else
    call EOS(X, Z, XCNO, mu, P(1), T(1), rho(1), kappa(1), epslon(1), tog_bf, 1 ,ierr)
    if (ierr.ne.0) stop
end if
!
! Apply approximate surface solutions to begin the integration,
! assuming radiation transport in the outermost zone (the do 20 loop).
! irc = 0 for radiation, irc = 1 for convection.
! Assume arbitrary initial values for kPad, and dlPdlT.
! dlPdlT = dlnP/dlnT (see Eq. 10.87)
!
kPad = 0.3d0
irc = 0
dlPdlT(1) = 4.25d0
do 20 i = 1, Nstart
    ip1 = i + 1
    call STARTMDL(deltar, X, Z, mu, Rs, r(i), M_r(i), L_r(i), r(ip1), P(ip1), M_r(ip1), L_r(ip1), T(ip1), tog_bf, irc)
    call EOS(X, Z, XCNO, mu, P(ip1), T(ip1), rho(ip1), kappa(ip1), epslon(ip1), tog_bf, ip1, ierr)
    if (ierr.ne.0) then
        write(*,400) r(i)/Rs, rho(i), M_r(i)/Ms, kappa(i), T(i), epslon(i), P(i), L_r(i)/Ls
        stop
    end if
!
! Determine whether convection will be operating in the next zone by
! calculating dlnP/dlnT numerically between zones i and i+1 (ip1).
! Update the adiabatic gas constant if necessary.
!
    if (i.gt.1) then
        dlPdlT(ip1) = log(P(ip1)/P(i))/log(T(ip1)/T(i))
    else
        dlPdlT(ip1) = dlPdlT(i)
    end if
    if (dlPdlT(ip1).lt.gamrat) then
        irc = 1
    else
        irc = 0
        kPad = P(ip1)/T(ip1)**gamrat
    end if
!
! Test to see whether the surface assumption of constant mass is still
! valid.
!
    deltaM = deltar*dMdr(r(ip1), rho(ip1))
    M_r(ip1) = M_r(i) + deltaM
    if (abs(deltaM).gt.0.001d0*Ms) then
        write(*,200)
        if (ip1.gt.2) ip1 = ip1 - 1
        go to 30
    end if
20 continue
!
! This is the main integration loop.  The assumptions of constant
! interior mass and luminosity are no longer applied.
!
30 Nsrtp1 = ip1 + 1
do 40 i = Nsrtp1, Nstop
    im1 = i - 1
!
! Initialize the Runge-Kutta routine by specifying zone i-1 quantities
! and their derivatives.  Note that the pressure, mass, luminosity,
! and temperature are stored in the memory locations f_im1(1),
! f_im1(2), f_im1(3), and f_im1(4), respectively.  The derivatives of
! those quantities with respect to radius are stored in dfdr(1),
! dfdr(2), dfdr(3), and dfdr(4).  Finally, the resulting values for
! P, M_r, L_r, and T are returned from the Runge-Kutta routine in
! f_i(1), f_i(2), f_i(3), and f_i(4), respectively.
!
! The stellar structure equations dPdr (Eq. 10.7), dMdr (Eq. 10.8),
! dLdr (Eq. 10.45), and dTdr (Eq. 10.61 or Eq. 10.81) are calculated
! in function calls, defined later in the code.
!
    f_im1(1) = P(im1)
    f_im1(2) = M_r(im1)
    f_im1(3) = L_r(im1)
    f_im1(4) = T(im1)
    dfdr(1)  = dPdr(r(im1), M_r(im1), rho(im1))
    dfdr(2)  = dMdr(r(im1), rho(im1))
    dfdr(3)  = dLdr(r(im1), rho(im1), epslon(im1))
    dfdr(4)  = dTdr(r(im1), M_r(im1), L_r(im1), T(im1), rho(im1), kappa(im1), mu, irc)
    call RUNGE(f_im1, dfdr, f_i, r(im1), deltar, irc, X, Z, XCNO, mu, i, ierr)
    if (ierr.ne.0) then
        write(*,300)
        write(*,400) r(im1)/Rs, rho(im1), M_r(im1)/Ms, kappa(im1), T(im1), epslon(im1), P(im1), L_r(im1)/Ls
        stop
    end if
!
! Update stellar parameters for the next zone, including adding
! dr to the old radius (note that dr <  0 since the integration is
! inward).
!
    r(i)   = r(im1) + deltar
    P(i)   = f_i(1)
    M_r(i) = f_i(2)
    L_r(i) = f_i(3)
    T(i)   = f_i(4)
!
! Calculate the density, opacity, and energy generation rate for
! this zone.
!
    call EOS(X, Z, XCNO, mu, P(i), T(i), rho(i), kappa(i), epslon(i), tog_bf, i, ierr)
    if (ierr.ne.0) then
        write(*,400) r(im1)/Rs, rho(im1), M_r(im1)/Ms, kappa(im1), T(im1), epslon(im1), P(im1), L_r(im1)/Ls
        stop
    end if
!
! Determine whether convection will be operating in the next zone by
! calculating dlnP/dlnT and comparing it to gamma/(gamma-1)
! (see Eq. 10.87).  Set the convection flag appropriately.
!
    dlPdlT(i) = log(P(i)/P(im1))/log(T(i)/T(im1))
    if (dlPdlT(i).lt.gamrat) then
        irc = 1
    else
        irc = 0
    end if
!
! Check to see whether the center has been reached.  If so, set Igoof and
! estimate the central conditions rhocor, epscor, Pcore, and Tcore.
! The central density is estimated to be the average density of the
! remaining central ball, the central pressure is determined by
! applying Eq. (H.4), and the central value for the energy
! generation rate is calculated to be the remaining interior
! luminosity divided by the mass of the central ball.  Finally, the
! central temperature is computed by applying the ideal gas law
! (where radiation pressure is neglected).
!
    if (r(i).le.abs(deltar) .and. (L_r(i).ge.0.1d0*Ls .or. M_r(i).ge.0.01d0*Ms)) then
        Igoof = 6
    else if (L_r(i).le.0.0d0) then
        Igoof = 5
        rhocor = M_r(i)/(4.0d0/3.0d0*pi*r(i)**3)
        if (M_r(i).ne.0.0d0) then
            epscor = L_r(i)/M_r(i)
        else
            epscor = 0.0d0
        end if
        Pcore = P(i) + 2.0d0/3.0d0*pi*G*rhocor**2*r(i)**2
        Tcore = Pcore*mu*m_H/(rhocor*k_B)
    else if (M_r(i).le.0.0d0) then
        Igoof  = 4
        Rhocor = 0.0d0
        epscor = 0.0d0
        Pcore  = 0.0d0
        Tcore  = 0.0d0
    else if (r(i).lt.0.02d0*Rs .and. M_r(i).lt.0.01d0*Ms .and. L_r(i).lt.0.1d0*Ls) then
        rhocor = m_r(i)/(4./3.*pi*r(i)**3)
        rhomax = 10.0d0*(rho(i)/rho(im1))*rho(i)
        epscor = L_r(i)/M_r(i)
        Pcore  = P(i) + 2.0d0/3.0d0*pi*G*rhocor**2*r(i)**2
        Tcore  = Pcore*mu*m_H/(rhocor*k_B)
        if (rhocor.lt.rho(i) .or. rhocor.gt.rhomax) then
            Igoof = 1
        else if (epscor.lt.epslon(i)) then
            Igoof = 2
        else if (Tcore.lt.T(i)) then
            Igoof = 3
        else
            Igoof = 0
        end if
    end if
    if (Igoof.ne.-1) then
        istop = i
        go to 50
    end if
!
! Is it time to change the step size?
!
    if (idrflg.eq.0 .and. M_r(i).lt.0.99d0*Ms) then
        deltar = -Rs/100.0d0
        idrflg = 1
    end if
    if (idrflg.eq.1 .and.deltar.ge.0.5*r(i)) then
        deltar = -Rs/5000.0d0
        idrflg = 2
    end if
    istop = i
40 continue
!
! Generate warning messages for the central conditions.
!
rhocor = M_r(istop)/(4.0d0/3.0d0*pi*r(istop)**3)
epscor = L_r(istop)/M_r(istop)
Pcore  = P(istop) + 2.0d0/3.0d0*pi*G*rhocor**2*r(istop)**2
Tcore  = Pcore*mu*m_H/(rhocor*k_B)
50 continue
if (Igoof.ne.0) then
    if (Igoof.eq.-1) then
        write(*,5000)
        write(*,5100)
    else if (Igoof.eq.1) then
        write(*,6000)
        write(*,5200) rho(istop)
        if (rhocor.gt.1.0d10) write(*,5300)
    else if (Igoof.eq.2) then
        write(*,6000)
        write(*,5400) epslon(istop)
    else if (Igoof.eq.3) then
        write(*,6000)
        write(*,5500) T(istop)
    else if (Igoof.eq.4) then
        write(*,5000)
        write(*,5600)
    else if (Igoof.eq.5) then
        write(*,5000)
        write(*,5700)
    else if (Igoof.eq.6) then
        write(*,5000)
        write(*,5800)
    end if
else
    write(*,7000)
end if
!
! Print the central conditions.  If necessary, set limits for the
! central radius, mass, and luminosity if necessary, to avoid format
! field overflows.
!
Rcrat = r(istop)/Rs
if (Rcrat.lt.-9.999d0) Rcrat = -9.999d0
Mcrat = M_r(istop)/Ms
if (Mcrat.lt.-9.999d0) Mcrat = -9.999d0
Lcrat = L_r(istop)/Ls
if (Lcrat.lt.-9.999d0) Lcrat = -9.999d0
write( *,2000) Msolar, Mcrat, Rsolar, Rcrat, Lsolar, Lcrat, Te,
1 rhocor, X, Tcore, Y, Pcore, Z, epscor, dlPdlT(istop)
write(20,1000)
write(20,2000) Msolar, Mcrat, Rsolar, Rcrat, Lsolar, Lcrat, Te,
1 rhocor, X, Tcore, Y, Pcore, Z, epscor, dlPdlT(istop)
write(20,2500) Ms
!
! Print data from the center of the star outward, labeling convective
! or radiative zones by c or r, respectively.  If abs(dlnP/dlnT)
! exceeds 99.9, set a print warning flag (*) and set the output limit
! to +99.9 or -99.9 as appropriate to avoid format field overflows.
!
write(20,3000)
do 60 ic = 1, istop
    i = istop - ic + 1
    Qm = 1.0d0 - M_r(i)/Ms
    if (dlPdlT(i).lt.gamrat) then
        rcf = 'c'
    else
        rcf = 'r'
    end if
    if (abs(dlPdlT(i)).gt.dlPlim) then
        dlPdlT(i) = sign(dlPlim,dlPdlT(i))
        clim = '*'
    else
        clim = ' '
    end if
    write(20,4000) r(i), Qm, L_r(i), T(i), P(i), rho(i), kappa(i), epslon(i), clim, rcf, dlPdlT(i)
60 continue
write(*,9000)
!
! Format statements
!
100 format(' ',/,' You must have X + Z <= 1',/,' please reenter composition',/)
200 format(' ',/,' The variation in mass has become larger than 0.001*Ms',/,' leaving the approximation loop before Nstart was reached',/)
300 format(' ',/,' The problem occurred in the Runge-Kutta routine',/)
400 format(' Values from the previous zone are:',/,10x,'r/Rs   = ',1pe12.5,'  ',12x,'rho     = ',1pe12.5,' g/cm**3',/,10x,'M_r/Ms = ',1pe12.5,'  ',12x,'kappa   = ',1pe12.5,' cm**2/g',/,10x,'T      = ',1pe12.5,' K',12x,'epsilon = ',1pe12.5,' ergs/g/s',/,10x,'P      = ',1pe12.5,' dynes/cm**2',/,10x,'L_r/Ls = ',1pe12.5)
1000 format(' ',15x,'A Homogeneous Main-Sequence Model',/)
2000 format(' ',/,' The surface conditions are:',10x,'The central conditions are:',//,' Mtot = ',0pf13.6,' Msun',12x,'Mc/Mtot     = ',1pe12.5,/,' Rtot = ',0pf13.6,' Rsun',12x,'Rc/Rtot     = ',1pe12.5,/,' Ltot = ',0pf13.6,' Lsun',12x,'Lc/Ltot     = ',1pe12.5,/,' Teff = ',0pf13.6,' K   ',12x,'Density     = ',1pe12.5,' g/cm**3',/,' X    = ',0pf13.6,'     ',12x,'Temperature = ',1pe12.5,' K',/,' Y    = ',0pf13.6,'     ',12x,'Pressure    = ',1pe12.5,' dynes/cm**2',/,' Z    = ',0pf13.6,'     ',12x,'epsilon     = ',1pe12.5,' ergs/s/g',/,'        ',   13x ,'     ',12x,'dlnP/dlnT   = ',1pe12.5,//)
2500 format(' ',/,' Notes: ',/,' (1) Mass is listed as Qm = 1.0 - M_r/Mtot, where Mtot = ',1pe13.6,' g',/,' (2) Convective zones are indicated by c, radiative zones by r',/,' (3) dlnP/dlnT may be limited to +99.9 or -99.9; if so it is',' labeled by *',//)
3000 format(' ',5x,'r',7x,'Qm',7x,'L_r',7x,'T',8x,'P',7x,'rho',6x,'kap',6x,'eps',3x,'dlPdlT')
4000 format(' ',1p8e9.2,2a1,0pf5.1)
5000 format(' ',/,15x,'Sorry to be the bearer of bad news, but...',/,15x,'       Your model has some problems',/)
5100 format(' ',   8x,'The number of allowed shells has been exceeded',/)
5200 format(' ', 14x,'The core density seems a bit off,'/,5x,' density should increase smoothly toward the center.',/,5x,' The density of the last zone calculated was rho = ',1pe10.3,' gm/cm**3',/)
5300 format(' ',   1x,'It looks like you will need a degenerate',' neutron gas and general relativity',/,' to solve this core.  Who do you think I am, Einstein?',/)
5400 format(' ',  14x,'The core epsilon seems a bit off,',/,9x,' epsilon should vary smoothly near the center.',/,9x,' The value calculated for the last zone was eps =',1pe10.3,' ergs/g/s',/)
5500 format(' ',8x,' Your extrapolated central temperature is too low',/,8x,' a little more fine tuning ought to do it.',/,8x,' The value calculated for the last zone was T = ',1pe10.3,' K',/)
5600 format(' ',  10x,'You created a star with a hole in the center!',/)
5700 format(' ',  10x,'This star has a negative central luminosity!',/)
5800 format(' ',  5x,'You hit the center before the mass and/or ','luminosity were depleted!',/)
6000 format(' ',///,15x,'It looks like you are getting close,',/,12x,'however, there are still a few minor errors',/)
7000 format(' ',///,15x,'CONGRATULATIONS, I THINK YOU FOUND IT!',/,9x,'However, be sure to look at your model carefully.',//)
9000 format(' ',10x,'***** The integration has been completed *****',/,10x,'   The model has been stored in starmodl.dat',/)
stop
end
!
! Subroutine STARTMDL computes values of M_r, L_r, P, and T, near the
! surface of the star using expansions of the stellar structure
! equations (M_r and L_r are assumed to be constant).
!
Subroutine STARTMDL(deltar, X, Z, mu, Rs, r_i, M_ri, L_ri, r, P_ip1, M_rip1, L_rip1, T_ip1, tog_bf, irc)
real*8 deltar, X, Z, mu, Rs, r_i, M_ri, L_ri, r, P_ip1, M_rip1, L_rip1, T_ip1, tog_bf
real*8 A_bf, A_ff, Afac
real*8 sigma, c, a, G , k_B, m_H, pi, gamma, gamrat, kPad, g_ff
common /cnstnt/ sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad, g_ff
!
r = r_i + deltar
M_rip1 = M_ri
L_rip1 = L_ri
!
! This is the radiative approximation (neglect radiation pressure
! and electron scattering opacity); see Eqs. (H.1), (H.2), (9.19),
! and (9.20).
!
if (irc.eq.0) then
    T_ip1 = G*M_rip1*mu*m_H/(4.25d0*k_B)*(1.0d0/r - 1.0d0/Rs)
    A_bf = 4.34d25*Z*(1.0d0 + X)/tog_bf
    A_ff = 3.68d22*g_ff*(1.0d0 - Z)*(1.0d0 + X)
    Afac = A_bf + A_ff
    P_ip1 = sqrt((1.0d0/4.25d0)*(16.0d0/3.0d0*pi*a*c)*(G*M_rip1/L_rip1)*(k_B/(afac*mu*m_H)))*T_ip1**4.25d0
!
! This is the convective approximation; see Eqs. (H.3) and (10.75).
!
else
    T_ip1 = G*M_rip1*mu*m_H/k_B*(1.0d0/r - 1.0d0/Rs)/gamrat
    P_ip1 = kPad*T_ip1**gamrat
end if
return
end
!
! Subroutine EOS calculates the values of density, opacity, the
! guillotine-to-gaunt factor ratio, and the energy generation rate at
! the radius r.
!
Subroutine EOS(X, Z, XCNO, mu, P, T, rho, kappa, epslon, tog_bf, izone, ierr)
real*8 X, Z, mu, P, T, rho, kappa, epslon, tog_bf, Prad, Pgas
real*8 k_bf, k_ff, k_e
real*8 T6, fx, fpp, epspp, epsCNO, XCNO, oneo3, twoo3, psipp, Cpp
real*8 sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad, g_ff
common /cnstnt/ sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad, g_ff
data oneo3, twoo3 / 0.333333333d0, 0.666666667d0 /
!
! Solve for density from the ideal gas law (remove radiation pressure);
! see Eq. (10.26).
!
if (T.le.0.0d0 .or. P.le.0.0d0) then
  ierr = 1
  write(*,100) izone, T, P
  return
end if
Prad = a*T**4/3.0d0
Pgas = P - Prad
rho=(mu*m_H/k_B)*(Pgas/T)
if (rho.lt.0.0d0) then
  ierr = 1
  write(*,200) izone, T, P, Prad, Pgas, rho
  return
end if
!
! Calculate opacity, including the guillotine-to-gaunt factor ratio;
! see Novotny (1973), p. 469. k_bf, k_ff, and k_e are the bound-free,
! free-free, and electron scattering opacities, given by Eqs. (9.19),
! (9.20), and (9.21), respectively.
!
tog_bf = 2.82d0*(rho*(1.0d0 + X))**0.2d0
k_bf = 4.34d25/tog_bf*Z*(1.0d0 + X)*rho/T**3.5d0
k_ff = 3.68d22*g_ff*(1.0d0 - Z)*(1.0d0 + X)*rho/T**3.5d0
k_e = 0.2d0*(1.0d0 + X)
kappa = k_bf + k_ff + k_e
!
! Compute energy generation by the pp chain and the CNO cycle.  These
! are calculated using Eqs. (10.49) and (10.53), which come from
! Fowler, Caughlan, and Zimmerman (1975). The screening factor for the
! pp chain is calculated as fpp; see Clayton (1968), p. 359ff.
!
T6 = T*1.0d-06
fx = 0.133d0*X*sqrt((3.0d0 + X)*rho)/T6**1.5d0
fpp = 1.0d0 + fx*X
psipp = 1.0d0 + 1.412d8*(1.0d0/X - 1.0d0)*exp(-49.98*T6**(-oneo3))
Cpp = 1.0d0 + 0.0123d0*T6**oneo3 + 0.0109d0*T6**twoo3 + 0.000938d0*T6
epspp = 2.38d6*rho*X*X*fpp*psipp*Cpp*T6**(-twoo3)*exp(-33.80d0*T6**(-oneo3))
CCNO = 1.0d0 + 0.0027d0*T6**oneo3 - 0.00778d0*T6**twoo3 - 0.000149d0*T6
epsCNO = 8.67d27*rho*X*XCNO*CCNO*T6**(-twoo3)*exp(-152.28d0*T6**(-oneo3))
epslon = epspp + epsCNO
!
! Formats
!
100 format(' ',/,' Something is a little wrong here.',/,' You are asking me to deal with either a negative temperature'/,' or a negative pressure.  I am sorry but that is not in my contract!',/,' You will have to try again with different',' initial conditions.',/,' In case it helps, I detected the problem in zone ',i3,' with the following',/,' conditions:',/,10x,'T = ',1pe10.3,' K',/,10x,'P = ',1pe10.3,' dynes/cm**2')
200 format(' ',/,' I am sorry, but a negative density was detected.',/,' my equation-of-state routine is a bit baffled by this new',/,' physical system you have created.  The radiation pressure',/,' is probably too great, implying that the star is unstable.'/,' Please try something a little less radical next time.',/,' In case it helps, I detected the problem in zone ',i3,' with the following',/,' conditions:',/,10x,'T       = ',1pe10.3,' K',/,10x,'P_total = ',1pe10.3,' dynes/cm**2',/,10x,'P_rad   = ',1pe10.3,' dynes/cm**2',/,10x,'P_gas   = ',1pe10.3,' dynes/cm**2',/,10x,'rho     = ',1pe10.3,' g/cm**3')
return
end
!
! The following four function subprograms calculate the gradients of
! pressure, mass, luminosity, and temperature at r.
!
real*8 Function dPdr(r, M_r, rho)
real*8 r, M_r, rho
real*8 sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad, g_ff
common /cnstnt/ sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad, g_ff
! Eq. (10.7)
dPdr = -G*rho*M_r/r**2
return
end
!
real*8 Function dMdr(r, rho)
real*8 r, rho
real*8 sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad, g_ff
common /cnstnt/ sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad, g_ff
! Eq. (10.8)
dMdr = 4.0d0*pi*rho*r**2
return
end
!
real*8 Function dLdr(r, rho, epslon)
real*8 r, rho, epslon
real*8 sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad, g_ff
common /cnstnt/ sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad, g_ff
! Eq. (10.45)
dLdr = 4.0d0*pi*rho*epslon*r**2
return
end
!
real*8 Function dTdr(r, M_r, L_r, T, rho, kappa, mu, irc)
real*8 r, M_r, L_r, T, rho, kappa, mu
real*8 sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad, g_ff
common /cnstnt/ sigma, c, a, G, k_B, m_H, pi, gamma, gamrat, kPad, g_ff
! This is the radiative temperature gradient (Eq. 10.61).
if (irc.eq.0) then
  dTdr = -(3.0d0/(16.0d0*pi*a*c))*kappa*rho/T**3*L_r/r**2
! This is the adiabatic convective temperature gradient (Eq. 10.81).
else
  dTdr = -1.0d0/gamrat*G*M_r/r**2*mu*m_H/k_B
end if
return
end
!
! This is a fourth-order Runge-Kutta integration routine.
!
Subroutine RUNGE(f_im1, dfdr, f_i, r_im1, deltar, irc, X, Z, XCNO, mu, izone, ierr)
real*8 f_im1(4), dfdr(4), f_i(4), df1(4), df2(4), df3(4), f_temp(4)
real*8 r_im1, r_i, deltar, dr12, dr16, r12, X, Z, XCNO, mu
!
dr12 = deltar/2.0d0
dr16 = deltar/6.0d0
r12  = r_im1 + dr12
r_i  = r_im1 + deltar
!
! Calculate intermediate derivatives from the fundamental stellar
! structure equations found in Subroutine FUNDEQ.
!
do 10 i = 1, 4
    f_temp(i) = f_im1(i) + dr12*dfdr(i)
10 continue
call FUNDEQ(r12, f_temp, df1, irc, X, Z, XCNO, mu, izone, ierr)
if (ierr.ne.0) return
!
do 20 i = 1, 4
    f_temp(i) = f_im1(i) + dr12*df1(i)
20 continue
call FUNDEQ(r12, f_temp, df2, irc, X, Z, XCNO, mu, izone, ierr)
if (ierr.ne.0) return
!
do 30 i = 1, 4
    f_temp(i) = f_im1(i) + deltar*df2(i)
30 continue
call FUNDEQ(r_i, f_temp, df3, irc, X, Z, XCNO, mu, izone, ierr)
if (ierr.ne.0) return
!
! Calculate the variables at the next shell (i + 1).
!
do 40 i = 1, 4
    f_i(i) = f_im1(i) + dr16*(dfdr(i) + 2.0d0*df1(i) + 2.0d0*df2(i) + df3(i))
40 continue
return
end
!
! This subroutine returns the required derivatives for RUNGE, the
! Runge-Kutta integration routine.
!
Subroutine FUNDEQ(r, f, dfdr, irc, X, Z, XCNO, mu, izone, ierr)
real*8 f(4), dfdr(4), X, Z, XCNO, r, M_r, P, T, L_r, rho, kappa, epslon, tog_bf, mu
real*8 dPdr, dMdr, dLdr, dTdr
!
P   = f(1)
M_r = f(2)
L_r = f(3)
T   = f(4)
call EOS(X, Z, XCNO, mu, P, T, rho, kappa, epslon, tog_bf, izone, ierr)
dfdr(1) = dPdr(r, M_r, rho)
dfdr(2) = dMdr(r, rho)
dfdr(3) = dLdr(r, rho, epslon)
dfdr(4) = dTdr(r, M_r, L_r, T, rho, kappa, mu, irc)
return
end
