import matplotlib.pyplot as plt
import numpy
import scipy.constants

v=numpy.linspace(1,3000,1000)
T=float(raw_input('Temperature (in K): '))

h=scipy.constants.h
c=scipy.constants.c
k=scipy.constants.k

f = c/(v*10**-9)

data = ((2*h*f**3)/(c**2))*(1/(numpy.e**((h*f)/(k*T))-1))

plt.plot(v, data)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Power Density (W/s)')
plt.show()
