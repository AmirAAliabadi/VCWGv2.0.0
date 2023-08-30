import numpy
import math

"""
Calculate the turbulent diffusion coefficient
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: March 2019
Originally developed by Scott Krayenhoff
"""

def TurbCoeff(nz,Ck,tke,dlk):

    """
    ------
    INPUT:
    nz: Number of grid points in the urban area
    Ck: Model constant
    tke: Turbulent kinetic energy profile [m^2 s^-2]
    dlk: Turbulent length scale [m]
    -------
    OUTPUT:
    Km: Turbulent diffusion coefficient [m^2 s^-1]
    """

    # Define turbulent diffusion coefficient [m^2 s^-1]
    Km = numpy.zeros(nz+1)

    # Km should be zero at street level
    Km[0] = 0.0

    # Calculate turbulent diffusion coefficient [m^2 s^-1] (eq. 4.8, Krayenhoff, PhD thesis)
    for i in range(1,nz-1):
        # Discretize TKE and length scale (vertical resolution (dz) is kept constant)
        tke_m = (tke[i-1]+tke[i])/2
        dlk_m = (dlk[i-1]+dlk[i])/2

        Km[i] = Ck*dlk_m * (math.sqrt(tke_m))

    Km[nz-1] = Km[nz-2]
    Km[nz] = Km[nz-2]

    return Km