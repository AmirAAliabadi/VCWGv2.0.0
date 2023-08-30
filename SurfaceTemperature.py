import os
import numpy
import math
from pprint import pprint
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from WaterFunctionsDef import IntDef
from scipy.optimize import least_squares
import copy

from Soil_Functions import Soil_Calculations

"""
Compute urban surfaces temperature by solving conduction equation
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Originally developed by Bruno Bueno
Last update: May 2021
"""

class Tsurf_Def(object):
    def __init__(self,thicknessLst, materialLst, T_init,name=None):

        if len(thicknessLst) != len(materialLst):
            raise Exception(self.THICKNESSLST_EQ_MATERIALLST_MSG)
        else:
            self._name = name                                             # purely for internal process
            self.layerThickness = thicknessLst                            # vector of layer thicknesses [m]
            self.layerThermalCond = list(map(lambda z: 0, materialLst))   # vector of layer thermal conductivity [W m^-1 K^-1]
            self.layerVolHeat = list(map(lambda z: 0, materialLst))       # vector of layer volumetric heat capacity [J m^-3 K^-1]

            # Create list of layer k and (Cp*density) from materialLst properties
            for i in range(len(materialLst)):
              self.layerThermalCond[i] = materialLst[i].thermalCond
              self.layerVolHeat[i] = materialLst[i].volHeat

            # vector of layer temperatures [K]
            self.layerTemp = [T_init] * len(thicknessLst)
            # external surface temperature [K]
            self.Text = self.layerTemp[0]
            # internal surface temperature [K]
            self.Tint = self.layerTemp[-1]

            self.z_depth = numpy.zeros(len(self.layerThickness) + 1)
            for i in range(1,len(self.layerThickness)+1):
                self.z_depth[i] = self.z_depth[i-1]+self.layerThickness[i-1]


    def Element(self,SWRabs,LWRabs,LE,Hf,dt,Tdp,bc,flux1=None,flux2=0):
        """
        ------
        INPUT:
        SWRabs: Net shortwave radiation [W m^-2]
        LWRabs: Net longwave radiation [W m^-2]
        LE: Latent heat flux [W m^-2]
        Hf: Sensible heat flux [W m^-2]
        dt: Time step [s]
        Tdp: Deep soil temperature [K]
        bc: Boundary condition at the bottom layer. 1: flux, 2:constant value
        flux1: Heat flux at the surface [W m^-2]
        flux2: Heat flux at the last layer [W m^-2]
        -------
        OUTPUT:
        layerTemp: Layer temperatures [K]
        Text: Exterior temperature [K]
        Tint: Last layer temperature [K]
        """

        # If no heat flux at the exterior surface is not given, it should be calculated from energy balance equation [W m^-2]
        if flux1 == None:
            flux1 = SWRabs + LWRabs - LE - Hf

        # Calculate soil layer temperatures
        # Here is assumed that the interior temperature for all surfaces is known
        self.layerTemp = self.Conduction(dt, flux1, bc, Tdp, flux2, self.layerTemp,self.layerVolHeat,self.layerThermalCond,self.layerThickness)
        # external surface temperature [K]
        self.Text = self.layerTemp[0]
        # internal surface temperature [K]
        self.Tint = self.layerTemp[-1]

    def Conduction(self,dt, flx1, bc, temp2, flx2, layerTemp_old, layerVolHeat, layerThermalCond, layerThickness):

        def is_near_zero(num, eps=1e-10):
            return abs(float(num)) < eps

        def invert(nz, A, C):
            X = [0 for i in range(nz)]

            for i in reversed(range(nz - 1)):
                C[i] = C[i] - A[i][2] * C[i + 1] / A[i + 1][1]
                A[i][1] = A[i][1] - A[i][2] * A[i + 1][0] / A[i + 1][1]

            for i in range(1, nz, 1):
                C[i] = C[i] - A[i][0] * C[i - 1] / A[i - 1][1]

            for i in range(nz):
                X[i] = C[i] / A[i][1]

            return X

        t = copy.copy(layerTemp_old)     # vector of layer temperatures [K]
        hc = copy.copy(layerVolHeat)     # vector of layer volumetric heat [J m^-3 K^-1]
        tc = copy.copy(layerThermalCond) # vector of layer thermal conductivities [W m^-1 K^-1]
        d = copy.copy(layerThickness)    # vector of layer thicknesses [m]

        # flx1                      : net heat flux on surface
        # bc                        : boundary condition parameter (1 or 2)
        # temp2                     : deep soil temperature (avg of air temperature)
        # flx2                      : surface flux (sum of absorbed, emitted, etc.)

        fimp = 0.5    # implicit coefficient
        fexp = 0.5    # explicit coefficient
        num = len(t)  # number of layers

        # Mean thermal conductivity over distance between 2 layers [W m^-1 K^-1]
        tcp = [0 for x in range(num)]
        # Thermal capacity times layer depth [J m^-2 K^-1]
        hcp = [0 for x in range(num)]
        # lower, main, and upper diagonals
        za = [[0 for y in range(3)] for x in range(num)]
        # RHS
        zy = [0 for x in range(num)]

        # --------------------------------------------------------------------------
        # Define the column vectors for heat capacity and conductivity
        hcp[0] = hc[0] * d[0]
        for j in range(1, num):
            tcp[j] = 2. / (d[j - 1] / tc[j - 1] + d[j] / tc[j])
            hcp[j] = hc[j] * d[j]

        # --------------------------------------------------------------------------
        # Define the first row of za matrix, and RHS column vector
        za[0][0] = 0.
        za[0][1] = hcp[0] / dt + fimp * tcp[1]
        za[0][2] = -fimp * tcp[1]
        zy[0] = hcp[0] / dt * t[0] - fexp * tcp[1] * (t[0] - t[1]) + flx1

        # --------------------------------------------------------------------------
        # Define other rows
        for j in range(1, num - 1):
            za[j][0] = fimp * (-tcp[j])
            za[j][1] = hcp[j] / dt + fimp * (tcp[j] + tcp[j + 1])
            za[j][2] = fimp * (-tcp[j + 1])
            zy[j] = hcp[j] / dt * t[j] + fexp * \
                    (tcp[j] * t[j - 1] - tcp[j] * t[j] - tcp[j + 1] * t[j] + tcp[j + 1] * t[j + 1])

        # --------------------------------------------------------------------------
        # Boundary conditions
        # heat flux
        if is_near_zero(bc - 1.):
            za[num - 1][0] = fimp * (-tcp[num - 1])
            za[num - 1][1] = hcp[num - 1] / dt + fimp * tcp[num - 1]
            za[num - 1][2] = 0.
            zy[num - 1] = hcp[num - 1] / dt * t[num - 1] + fexp * tcp[num - 1] * (t[num - 2] - t[num - 1]) + flx2
        # deep-temperature
        elif is_near_zero(bc - 2.):
            za[num - 1][0] = 0.
            za[num - 1][1] = 1.
            za[num - 1][2] = 0.
            zy[num - 1] = temp2
        else:
            print('ERROR: check input parameters in the Conduction routine')

        # --------------------------------------------------------------------------

        zx = invert(num, za, zy)
        # return zx as 1d vector of temperature layers
        return zx



