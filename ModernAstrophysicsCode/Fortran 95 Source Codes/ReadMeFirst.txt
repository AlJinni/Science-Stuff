This folder contains all of the source files for the Fortran 95 codes presented in the text

	An Introduction to Modern Astrophysics
	Bradley W. Carroll and Dale A. Ostlie
	Addison Wesley, copyright 2007

and in the associated texts

	An Introduction to Modern Stellar Astrophysics
	Bradley W. Carroll and Dale A. Ostlie
	Addison Wesley, copyright 2007

	An Introduction to Modern Galactic Astronomy and Cosmology
	Bradley W. Carroll and Dale A. Ostlie
	Addison Wesley, copyright 2007

The subfolders contained here are:

(a)  Constants (Appendix I):  
	This folder contains an extensive list of astronomical, physical, and machine constants in the file Constants.F90.  The constants are presented to the highest accuracy known or possible for the machine being used.

(b)  Orbit (Appendix J):
	Orbit.F90 is a simple code that computes planetary orbits, or orbits of the reduced mass of a binary system.  Orbit.F90 requires Constants.F90 in order to be compiled.

(c)  TwoStars (Appendix K):
	TwoStars.F90 computes light curves and radial velocity curves for eclipsing, detached binary systems under the assumption that the stars are spherically symmetric.  TwoStars.F90 does incorporate a simple expression for limb darkening.  TwoStars.F90 requires Constants.F90 in order to be compiled.

(d)  StatStar.F90 (Appendix L):
	StatStar computes ZAMS models of stars using the physics presented in the text.  Simplifying assumptions made in this code include using Kramers opacity laws, adiabatic convection, and unsophisticated surface boundary conditions.  As such StatStar is not a research-grade stellar structure code, but is meant to introduce the student to the basic concepts of stellar model building.  StatStar is composed of a number of Fortran 95 modules, each containing one or more subroutines or functions.  StatStar also requires Constants.F90 in order to be compiled.