/*PROGRAM TwoStars
    General Description:
    ====================
        TwoStars computes the radial velocities and light curves of a 
        binary star system.  The code is not designed to produce 
        high-precision, state-of-the-art solutions to the binary star 
        problem.  Rather, this code is designed to be illustrative of 
        the fundamental ideas discussed in
 
            "An Introduction to Modern Astrophysics", Appendix K
            Bradley W. Carroll and Dale A. Ostlie
            Second Edition, Addison Wesley,   Copyright 2007.
 
            Weber State University
            Ogden, UT
            modastro@weber.edu
 
        The model used here assumes or includes the following:
            1.  The stars are spherically symmetric, implying
                a.	tidal distortions and centrifugal effects are 
                    neglected.
                b.	star spots or reflective heating are not 
                    included.
            2.  A simple formula for limb darkening is included.
            3.  Elliptical orbits are allowed.
            4.  An arbitrary inclination angle is allowed.
            5.  An arbitrary periastron angle is allowed.
            6.  The user may enter the velocity vector of the center 
                of mass through space.
 
        The plane of the sky is the (y',z') plane, with x' being along 
        the line of sight.  The plane of the binary star orbit is the 
        (x,y) plane, which is tilted at an inclination angle of i 
        relative to the (y',z') plane.  Note that y' = y. The center 
        of mass is initally assumed to be at the common origin of the 
        two coordinate systems.  However, the center of mass can also 
        be given an arbitrary velocity vector through space.
 
        If phi =  0 degrees, the major axis is aligned with the x axis,
        with periastron in +x direction.
        If phi = 90 degrees, the major axis is aligned with y axis, 
        with periastron in +y direction.
 
              *         | z'     *
               *        |       *
                * z     |      * 
                 *      |     * <--- plane of orbit
                  *     |    * 
                   *    | i *
                    *   |  *
                     *  | *
        x'            * |*      (y, y' out of paper)
        <---------------+---------------
         to            *|*              
         observer     * | *
                     *  |  *
                    *   |   *
                 x *    |    *
                  *     |     *
                 *      |      *
                *       | <--- plane of sky
 
        The transformation between coordinate systems is:
            x' = z*COS(i) + x*SIN(i)
            y' = y
            z' = z*SIN(i) - x*COS(i)
 
----------------------------------------------------------------*/

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include "Constants.h"

using namespace std;

int main(){

    using ModAstroConstants::M_Sun;
    using ModAstroConstants::R_Sun;
    using ModAstroConstants::L_Sun;
    using ModAstroConstants::Mbol_Sun;
    using ModAstroConstants::G;
    using ModAstroConstants::sigma;
    using ModAstroConstants::two_pi;
    using ModAstroConstants::four_pi;
    using ModAstroConstants::AU;
    using ModAstroConstants::day;
    using ModAstroConstants::degrees_to_radians;
    
    void Eclipse (double x1, double y1, double R1, double T1, double x2, double y2, double R2, double T2, double i, 
              double& dS, double& y1p, double& z1p, double& y2p, double& z2p);
    double F (double r_prime, double R, double T);

//open output file
    char xpause;
    string file_name;

    cout << "Specify the name of your output file: ";
    getline(cin, file_name);

    ofstream outdata;
    outdata.open(file_name.c_str(), ios::out);
    if (outdata.fail()){
        cout << " Unable to open Orbit.txt --- Terminating calculation" 
            << "\n\n"
            << "Enter any character to quit: ";
        cin  >> xpause;   
        exit(1);
    }


//read input data
    double m1, R1, T1;
    cout << "\n\n"
         << "Enter the data for Star #1"
         << "\n   Mass (solar masses):             ";
    cin  >> m1;
    cout <<   "   Radius (solar radii):            ";
    cin  >> R1;
    cout <<   "   Effective Temperature (K):       ";
    cin  >> T1;

    double m2, R2, T2;
    cout << "\n\n"
         << "Enter the data for Star #2"
         << "\n   Mass (solar masses):             ";
    cin  >> m2;
    cout <<   "   Radius (solar radii):            ";
    cin  >> R2;
    cout <<   "   Effective Temperature (K):       ";
    cin  >> T2;

    double P, e, i, phi;
    cout << "\n\n"
         << "Enter the desired orbital parameters"
         << "\n   Orbital Period (days):           ";
    cin  >> P;
    cout <<   "   Orbital Eccentricity:            ";
    cin  >> e;
    cout <<   "   Orbital Inclination (deg):       ";
    cin  >> i;
    cout <<   "   Orientation of Periastron (deg): ";
    cin  >> phi;

    double vcmxp, vcmyp, vcmzp;
    cout << "\n\n"
         << "Enter the x', y', and z' components of the center of mass velocity vector:\n"
         << "Notes:  (1)  The plane of the sky is (y',z')\n"
         << "        (2)  If v_x' < 0, then the center of mass is blueshifted" << endl;
    cout << "\n\n"
         <<   "   v_x' (km/s)                      ";
    cin  >> vcmxp;
    cout <<   "   v_y' (km/s)                      ";
    cin  >> vcmyp;
    cout <<   "   v_z' (km/s)                      ";
    cin  >> vcmzp;

//Convert values to conventional SI units and radians
    double M;
    m1  = m1*M_Sun;
    m2  = m2*M_Sun;
    M   = m1 + m2;                              //Total mass of system
    R1  = R1*R_Sun;
    R2  = R2*R_Sun;
    P   = P*day;
    i   = i*degrees_to_radians;
    phi = phi*degrees_to_radians;

    vcmxp = vcmxp*1000;
    vcmyp = vcmyp*1000;
    vcmzp = vcmzp*1000;

//Compute the semimajor axes of the orbits
    double mu, a, a1, a2;
    mu = m1*m2/M;                                                   //Reduced mass, Eq. (2.22)
    a  = pow(static_cast<double>(P*P*G*M/(two_pi*two_pi)), 1/3.);   //Kepler's 3rd Law, Eq. (2.37)
    a1 = (mu/m1)*a;                                                 //Semimajor axis, Eq. (2.23)
    a2 = (mu/m2)*a;                                                 //Semimajor axis, Eq. (2.24)

    cout << fixed << showpoint << setprecision(6) << "\n\n"
         << "The semimajor axis of the reduced mass is " << setw(13) << a/AU << " AU\n"
         << "   a1 = " << setw(13) << a1/AU << " AU\n"
         << "   a2 = " << setw(13) << a2/AU << " AU" << endl;

    cout << "\n\n"
         << "Enter any character to perform calculation: ";
    cin  >> xpause;   

    if (a < R1 + R2) {
        cout << "\n\n"
             << "Your two stars are in contact!\n"
             << " R1 + R2 = " << setw(13) << (R1 + R2)/AU << " AU\n"
             << "The spherically symmetric approximation is clearly invalid.\n\n"
             << "****** Terminating Calculation ******" << endl;}
    else {

    //Compute the luminosity of each star and the total luminosity
        double L1, L2, L;
        L1   = four_pi*R1*R1*sigma*pow(T1,4);   //Stefan-Boltzmann, Eq. (3.17)
        L2   = four_pi*R2*R2*sigma*pow(T2,4);   //Star #2
        L    = L1 + L2;                         //Total uneclipsed luminosity

    //Determine the energy produced by the uneclipsed, projected disks
        const int Nr = 100;
        double drS1, drS2, rS1, rS2, S1, S2;
        drS1 = R1/Nr;
        drS2 = R2/Nr;
        rS1  = 0;
        rS2  = 0;
        S1   = 0;
        S2   = 0;

    //Numerical integration loop (Eq. J.7)
        for (int j = 1; j <= Nr; j++) {
            rS1 += drS1;
            rS2 += drS2;
            S1  += two_pi*F(rS1 - drS1/2, R1, T1)*(rS1 - drS1/2)*drS1;
            S2  += two_pi*F(rS2 - drS2/2, R2, T2)*(rS2 - drS2/2)*drS2;
        }
        double S = S1 + S2;

    //Initial orbit loop conditions
        const int N = 1001;
        double t, dt, theta, L_ang, dAdt;
        t       = 0;
        dt      = P/N;                          //time step
        theta   = 0;
        L_ang   = mu*sqrt(G*M*a*(1 - e*e));     //L, Eq. (2.30)
        dAdt    = L_ang/(2*mu);                 //2nd Law, Eq. (2.32)

    //Write output headers
        cout << "\n\n"
             << "         t/P        v1r (km/s)     v2r (km/s)         Mbol        dS (W) " << endl;
        
        outdata << "         t/P        v1r (km/s)     v2r (km/s)         Mbol        dS (W) "
                << "      y1p (AU)       z1p (AU)       y2p (AU)       z2p (AU)      ycmp (AU)      zcmp (AU) "
                << endl;

        //Reduced mass orbit loop
        double r, v, vr, v1r, v2r, x, y, x1, y1, x2, y2;
        double dS, y1p, z1p, y2p, z2p, Lt, Mbol, t_max, dtheta;
        double Mbol_max = -999999;
        while (t < P + dt/2) {
            r   =  a*(1 - e*e)/(1 + e*cos(theta));  //position, Eq. (2.3)
            v   =  sqrt(G*M*(2/r - 1/a));           //velocity, Eq. (2.36)
            vr  = -v*sin(i)*sin(theta + phi);       //radial velocity
            v1r =  (mu/m1)*vr;                      //Eq. (2.23)
            v2r = -(mu/m2)*vr;                      //Eq. (2.24)

        //Determine (x,y) positions of centers of stars
            x   = r*cos(theta + phi);               //reduced mass
            y   = r*sin(theta + phi);               //reduced mass
            x1  =  (mu/m1)*x;
            y1  =  (mu/m1)*y;
            x2  = -(mu/m2)*x;
            y2  = -(mu/m2)*y;

            Eclipse (x1, y1, R1, T1, x2, y2, R2, T2, i, dS, y1p, z1p, y2p, z2p);
            Lt = L*(1 - dS/S);
            Mbol = Mbol_Sun - 5*log10(Lt/L_Sun)/2;  //Mbol, Eq. (3.8)
            if (Mbol > Mbol_max) {
                Mbol_max = Mbol;
                t_max = t;}

        //Print results to the screen with radial velocities in km/s
            cout << fixed << showpoint << right << setprecision(6);
            cout << setw(15) << t/P << setw(15) << (v1r + vcmxp)/1000 << setw(15) 
                << (v2r + vcmxp)/1000 << setw(15) << Mbol;
            cout << scientific;
            cout << setw(15) << Lt*dS/S << endl;

        //Print results to the output file with radial velocities in km/s
            outdata << fixed << showpoint << right << setprecision(6);
            outdata << setw(15) << t/P << setw(15) << (v1r + vcmxp)/1000 << setw(15) 
                << (v2r + vcmxp)/1000 << setw(15) << Mbol;
            outdata << scientific;
            outdata << setw(15) << Lt*dS/S << setw(15) << (y1p + vcmyp*t)/AU 
                << setw(15) << (z1p + vcmzp*t)/AU << setw(15) << (y2p + vcmyp*t)/AU 
                << setw(15) << (z2p + vcmzp*t)/AU << setw(15) << vcmyp*t/AU 
                << setw(15) << vcmzp*t/AU << endl;
                
            dtheta = (2*dAdt/(r*r))*dt;             //Eq. (2.31)
            theta  += dtheta;
            t      += dt;
        }

        cout << fixed << setprecision(6) << "\n\n"
             << "The deepest minimum in the light curve occurred at t/P = " << setw(13) << t_max/P << endl;
        cout << "with a value of Mbol = " << setw(13) << Mbol_max << endl;
    };

    cout << "\n\n"
         << "Your calculation has finished; enter a character and hit <enter> to exit: ";
    cin  >> xpause;
}


void Eclipse (double x1, double y1, double R1, double T1, double x2, double y2, double R2, double T2, double i, 
              double& dS, double& y1p, double& z1p, double& y2p, double& z2p)
{
//   General Description:
//   ====================
//      This routine computes the change in observed luminosity due to an eclipse.
//---------------------------------------------------------------

    using ModAstroConstants::pi;

    void Transformation (double x, double y, double& xp, double& yp, double& zp, double i);
    double F (double r_prime, double R, double T);

//  Number of steps for r', theta' integrations
    const int       Nr = 100, Ntheta = 500;
    const double    dtheta_prime = pi/Ntheta;

    dS = 0;

    double x1p, x2p;

//  Perform coordinate transformations
    Transformation (x1, y1, x1p, y1p, z1p, i);
    Transformation (x2, y2, x2p, y2p, z2p, i);

//  Determine which star is in front (f) and which star is in back (b)
    double xfp, yfp, zfp, Rf, Tf;
    double xbp, ybp, zbp, Rb, Tb;
    if (x1p > 0) {
        xfp = x1p; yfp = y1p; zfp = z1p; Rf = R1; Tf = T1;
        xbp = x2p; ybp = y2p; zbp = z2p; Rb = R2; Tb = T2;}
    else {
        xfp = x2p; yfp = y2p; zfp = z2p; Rf = R2; Tf = T2;
        xbp = x1p; ybp = y1p; zbp = z1p; Rb = R1; Tb = T1;}
    
//  Are the two stars close enough for an eclipse?
    double d = sqrt(pow((yfp - ybp), 2.) + pow((zfp - zbp), 2.));               //Eq. (J.4)
    if (d <= Rf + Rb) {
        //Find the angle between y' and the projected line between the 
        //centers of the stars in the (y',z') plane.  The polar coordinate
        //integration will be centered on this line to take advantage 
        //of spherical symmetry.

        double theta0_prime = atan2((zfp - zbp), (yfp - ybp));      //Eq. (J.5)

        //Determine the starting radius for the integration
        double r_prime, r_stop;
        if (d < Rb - Rf) {
            r_prime = d + Rf;       //Foreground star disk entirely 
            r_stop  = d - Rf;       //inside background star disk
            if (r_stop < 0) r_stop = 0;}
        else {
            r_prime = Rb;
            r_stop  = 0;}
        double dr_prime = r_prime/Nr;
        
        //The surface integration loop
        do {
            
            //Determine the limits of the angular integration for the current r_prime
            double theta_prime = theta0_prime;
            bool finished = false;
            while (!finished) {
                double yp_dA = r_prime*cos(theta_prime + dtheta_prime) + ybp;
                double zp_dA = r_prime*sin(theta_prime + dtheta_prime) + zbp;
                if (sqrt(pow((yp_dA - yfp), 2.) + pow((zp_dA - zfp), 2.)) > Rf) finished = true;  //Eq. (J.6)

                theta_prime += dtheta_prime;
                if ((theta_prime - theta0_prime) > pi) finished = true;
            }

            //Add the luminosity change for differential area  (Eq. J.7)
            dS += 2*F(r_prime - dr_prime/2, Rb, Tb)*(r_prime - dr_prime/2)*dr_prime*(theta_prime - theta0_prime);

            //Check to see that there is no remaining overlap or if center of disk has been reached
            r_prime -= dr_prime;
        } while (r_prime >= r_stop);
    }
}


void Transformation (double x, double y, double& xp, double& yp, double& zp, double i)
{
//  General Description:
//  ====================
//      This routine performs the coordinate transformation between
//      the orbital plane coordinates (x,y) and the plane of the sky
//      coordinates (xp, yp, zp), based on the angle of inclination, i.
//--------------------------------------------------------------

    xp =  x*sin(i);         //Eq. (J.1)
    yp =  y;                //Eq. (J.2)
    zp = -x*cos(i);         //Eq. (J.3)
}

double F(double r_prime, double R, double T)
{
//  General Description:
//  ====================
//      This routine calculates the flux for a specific value
//      of the distance from the center of the stellar disk, r_prime.
//--------------------------------------------------------------

    using ModAstroConstants::sigma;

    double Limb_Darkening (double r_prime, double R);

    return sigma*pow(T, 4.)*Limb_Darkening (r_prime, R);              //Eq. (3.18)
}

double Limb_Darkening (double r_prime, double R)
{
//  General Description:
//  ====================
//      This routine calculates limb darkening for a specific value 
//      of the distance from the center of the stellar disk, r_prime.
//
//      Data are due to Van Hamme, W., Astronomical Journal, 
//          Vol. 106, 1096, 1993.
//--------------------------------------------------------------

    using ModAstroConstants::pi_over_2;
    
//  Van Hamme model #109 (solar-like)
    const double x = 0.648, y = 0.207;

    double mu = cos(r_prime/R*pi_over_2);
    return 1 - x*(1 - mu) - y*mu*log(mu);                       //Van Hamme Eq. (3)
}