import numpy
import math

"""
This Class is used to calculate buoyancy term for the TKE equation. (eq. 5.2, term IX, Krayenhoff, PhD thesis)
Developed by Scott Krayenhoff
University of Guelph, Guelph, Canada
Last update: March 2017
"""

def BuoProd(cdmin,nz,dz,th,Km,th0,prandtl):
    """
    ------
    INPUT:
    cdmin: minimum drag
    nz: Number of grid points in vertical column
    dz: Grid resolution [m]
    th: Potential temperature [K]
    Km: Turbulent diffusion coefficient [m^2 s^-1]
    th0: Reference potential temperature [K]
    prandtl: Turbulent Prandtl number
    -------
    OUTPUT:
    bu: Buoyant production [m^2 s^-3]
    """

    # Define buoyancy term [m^2 s^-3]
    bu = numpy.zeros(nz)
    # Set the buoyancy value at street level to zero
    bu[0] = 0
    # Calculate buoyancy using (eq. 5.2, term IX, Krayenhoff, PhD thesis)
    # buoyant production [m^2 s^-3] = (g/th0)*(Km/prandtl)*(dth/dz)
    for i in range(1,nz-1):
        # Discretize potential temperature gradient
        dthdz1 = (th[i]-th[i-1])/dz
        dthdz2 = (th[i+1]-th[i])/dz
        cdm = max(0.5*(Km[i]+Km[i+1])/prandtl,cdmin)

        dthmdz = 0.5*(dthdz1+dthdz2)

        bu[i] = -9.81*cdm*dthmdz/th0[i]

    bu[nz-1] = 0

    return bu