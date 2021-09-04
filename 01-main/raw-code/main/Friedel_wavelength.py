import numpy as np
def Fermi_vector(mu, t, omega):
    return np.sqrt((mu+omega)/t+4)
#
def Friedel_wavelength(k_F):
    return np.pi/k_F
#
def Main():
    k_F=Fermi_vector(mu, t, omega)
    x=Friedel_wavelength(k_F)
    print('Fermi vector: '+str(k_F))
    print('Friedel wavelength: '+str(x))
mu=-3.88
t=1.0
omega=0.0
Main()