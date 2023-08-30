import numpy
import math

"""
Shear production calculation
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: March 2019
Originally developed by Alberto Martilli, Scott Krayenhoff, and Negin Nazarian
"""

def ShearProd(cdmin,nz,dz,vx,vy,km):
    '''
    :param:
    cdmin: minimum drag
    nz: Number of grid points in vertical column
    dz: Grid resolution [m]
    vx: x component of horizontal wind speed [m s^-1]
    vy: y component of horizontal wind speed [m s^-1]
    km: Turbulent diffusion coefficient [m^2 s^-1]
    :return:
    sh: shear production [m^2 s^-3]
    '''

    # Define shear production term [m^2 s^-3]
    sh = numpy.zeros(nz)
    # Set the shear production value at street level to zero
    sh[0] = 0
    # Calculate shear production (eq. 5.2, term II, Krayenhoff, PhD thesis)
    # shear production [m^2 s^-3] = Km*[(du/dz)^2+(dv/dz)^2]
    for i in range(1,nz-1):
        # Discretize gradient of x and y components of velocities
        dudz1 = (vx[i]-vx[i-1])/dz
        dvdz1 = (vy[i]-vy[i-1])/dz
        dudz2 = (vx[i+1]-vx[i])/dz
        dvdz2 = (vy[i+1]-vy[i])/dz

        cdm = max(0.5*(km[i]+km[i+1]),cdmin)

        dumdz = 0.5*((dudz1**2+dvdz1**2)+(dudz2**2+dvdz2**2))

        sh[i] = cdm*dumdz

    sh[nz-1] = 0

    return sh