import numpy
import math
import copy

'''
Rural model: The Monin-Obukhov Similarity Theory (MOST) is used to solve for the vertical profile of potential temperature
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: May 2021
'''

class RSMDef(object):

    def __init__(self, RSMParam,z,nz,dz,T_init,P_init,S_init):

        self.z0 = RSMParam.z0overh_MOST*RSMParam.h_obs       # Rural Aerodynamic Roughness Length [m]
        self.z_T = RSMParam.zToverz0_MOST*self.z0            # Rural thermodynamic length scale [m]
        self.disp = RSMParam.dispoverh_MOST * RSMParam.h_obs # Rural displacement length [m]
        self.RSMParam = RSMParam                             # Constant parameters
        self.z = z                                           # Grid points [m]
        self.nz = nz                                         # Number of grid points in vertical column
        self.dz = dz                                         # Grid resolution [m]

        # Initialize potential temperature profile in the rural area [K]
        self.T_rural = [T_init for x in range(self.nz)]
        # Initialize specific humidity profile in the rural area [g kg^-1]
        self.q_rural = [0 for x in range(self.nz)]
        # Initialize pressure profile in the rural area [Pa]
        self.presProf = [P_init for x in range(self.nz)]
        # Initialize real temperature profile in the rural area [K]
        self.tempRealProf = [T_init for x in range(self.nz)]
        # Initialize density profile at the center of layers in the rural area [kg m^-3]
        self.densityProfC = [None for x in range(self.nz)]
        # Initialize wind speed
        self.S_rural = [S_init for x in range(self.nz)]
        self.u_Zp = S_init
        self.u_veg = S_init

        # Calculate pressure profile [Pa]
        for iz in range(1, self.nz):
            self.presProf[iz] = (self.presProf[iz - 1] ** (287 / RSMParam.cp) -
                                9.81 / RSMParam.cp * (P_init ** (287 / RSMParam.cp)) *
                                 (1. / self.T_rural[iz] + 1. / self.T_rural[iz - 1]) * 0.5 * self.dz) ** \
                                (1. / (287 / RSMParam.cp))

        # Calculate real temperature profile [K]
        for iz in range(self.nz):
            self.tempRealProf[iz] = self.T_rural[iz] * (self.presProf[iz] / P_init) ** (RSMParam.r / RSMParam.cp)

        # Calculate density profile [kg m^-3]
        for iz in range(self.nz):
            self.densityProfC[iz] = self.presProf[iz] / RSMParam.r / self.tempRealProf[iz]

    def MOST(self,EBRural,Text_Rural,MeteoData):
        """
        ------
        INPUT:
        EBRural: Energy balance terms at the surface
        Text_Rural: Surface temperature [K]
        MeteoData: Forcing variables
        -------
        OUTPUT:
        presProf: Pressure profile [Pa]
        tempRealProf: Real temperature profile [K]
        densityProfC: Density profile [kg m^-3]
        u_star: Friction velocity [m s^-1]
        T_rural: Potential temperature profile [K]
        q_rural: Specific humidity profile [kg kg^-1]
        """

        # Calculate pressure profile [Pa]
        for iz in reversed(list(range(self.nz))[1:]):
            self.presProf[iz - 1] = (math.pow(self.presProf[iz], self.RSMParam.r/self.RSMParam.cp) +
                                     9.81/self.RSMParam.cp*(math.pow(MeteoData.Pre, self.RSMParam.r/self.RSMParam.cp)) *
                                     (1./self.T_rural[iz] + 1./self.T_rural[iz-1])*0.5*self.dz) ** (1./(self.RSMParam.r/self.RSMParam.cp))

        # Calculate the real temperature profile [K]
        for iz in range(self.nz):
            self.tempRealProf[iz] = self.T_rural[iz] * (self.presProf[iz] / MeteoData.Pre) ** (self.RSMParam.r / self.RSMParam.cp)

        # Calculate the density profile [kg m^-3]
        for iz in range(self.nz):
            self.densityProfC[iz] = self.presProf[iz] / self.RSMParam.r / self.tempRealProf[iz]

        # Temperature at the lower bound of integral
        # Option 1 (theta_lb): extrapolation
        Slope_T = (MeteoData.Tatm - Text_Rural) / self.z_T
        theta_lb = Text_Rural + Slope_T * self.z_T
        # Option 2 (theta_lb): forcing temperature
        # theta_lb = MeteoData.Tatm
        
        # Number of iteration used to determine friction velocity
        N_iter = 5

        self.ur_wind = copy.copy(MeteoData.Uatm)
        # Check minimum wind speed
        if self.ur_wind < self.RSMParam.WindMin_MOST:
            self.ur_wind = self.RSMParam.WindMin_MOST

        g = 9.81
        # Density at the reference level [kg m^-2]
        rho_0 = MeteoData.Pre / (287 * theta_lb)

        # Calculate turbulent sensible heat flux [K m s^-1]
        self.wt = EBRural.EnergyFlux.HfluxRural / (rho_0 * self.RSMParam.cp)
        # Calculate turbulent latent heat flux [kg kg^-1 m s^-1]
        self.wq = EBRural.EnergyFlux.LEfluxRural / (rho_0 * self.RSMParam.lv)

        # Calculate friction velocity iteratively considering roughness and stability effect
        # Friction velocity at the reference level [m s^-1]
        u_star_init = self.RSMParam.vk * self.ur_wind / math.log((self.RSMParam.h_wind-self.disp) / self.z0)
        # Obukhov length at the reference level [m]
        L_init = -(theta_lb * u_star_init ** 3) / (self.RSMParam.vk * g * self.wt)

        #Iterate until we converge at a stability-corrected friction velocity
        #Iterate maximum number of times
        for j in range(1, N_iter):
            #Stable
            # Solve equation using Businger et al. 1971 and Dyer 1970
            if (self.RSMParam.h_wind-self.disp) / L_init > self.RSMParam.ZL_Pos_cutoff:
                # Option 1 (stable u_star): u_star based on roughness length and stability
                # self.u_star = (self.ur_wind * self.RSMParam.vk)/(math.log((self.RSMParam.h_wind-self.disp) / self.z0)+(4.7*(self.RSMParam.h_wind-self.disp))/L_init-4.7*self.z0/L_init)
                # Option 2 (stable u_star): fraction of mean wind
                self.u_star = 0.07*self.ur_wind

            #Unstable
            elif (self.RSMParam.h_wind-self.disp) / L_init < self.RSMParam.ZL_Neg_cutoff:
                
                # Solve equation using parameterization in Paulson 1970
                alfa2 = (1-16*(self.RSMParam.h_wind-self.disp)/L_init)**(0.25)
                alfa1 = (1 - 16 * (self.z0) / L_init) ** (0.25)
                ksi2 = -2*math.log((1+alfa2)/2)-math.log((1+alfa2**2)/2)+2*math.atan(alfa2)-math.pi/2
                ksi1 = -2*math.log((1+alfa1)/2)-math.log((1+alfa1**2)/2)+2*math.atan(alfa1)-math.pi/2
                self.u_star = self.ur_wind*self.RSMParam.vk/(math.log((self.RSMParam.h_wind-self.disp)/self.z0)+ksi2-ksi1)
                
            #Neutral
            else:
                self.u_star = self.RSMParam.vk * self.ur_wind / math.log((self.RSMParam.h_wind-self.disp) / self.z0)

            if self.u_star < self.RSMParam.u_star_min_MOST:
                self.u_star = self.RSMParam.u_star_min_MOST
            self.L = -(theta_lb * self.u_star ** 3) / (self.RSMParam.vk * g * self.wt)
            L_init = self.L


        if self.L > 0 and self.L < self.RSMParam.L_Pos_min:
            self.L = self.RSMParam.L_Pos_min
        if self.L > 0 and self.L > self.RSMParam.L_Pos_max:
                self.L = self.RSMParam.L_Pos_max
        if self.L < 0 and self.L > self.RSMParam.L_Neg_max:
            self.L = self.RSMParam.L_Neg_max
        if self.L < 0 and self.L < self.RSMParam.L_Neg_min:
            self.L = self.RSMParam.L_Neg_min

        # Surface layer temperature scale
        self.theta_sl = -self.wt / self.u_star
        # Surface layer humidity scale
        self.q_sl = -self.wq / self.u_star

        for j in range(0, self.nz):
            # Stable
            if (self.z[j]-self.disp) / self.L > self.RSMParam.ZL_Pos_cutoff:

                # Businger et al. 1971 and Dyer 1970
                self.T_rural[j] = (self.theta_sl / self.RSMParam.vk) * (numpy.log((self.z[j]-self.disp) / self.z_T) +
                                                                        5*((self.z[j]-self.disp) / self.L)- 5*(self.z_T/self.L)) + theta_lb
                # Specific humidity equation
                # Calculate forcing specific humidity
                # Closed formula
                self.q_rural[j] = (self.q_sl / self.RSMParam.vk) * (numpy.log((self.z[j]-self.disp)/self.z_T) +
                                                                     5*((self.z[j]-self.disp)/self.L) - 5*(self.z_T / self.L)) + MeteoData.q_atm
    
            # Unstable
            elif self.z[1] / self.L < self.RSMParam.ZL_Neg_cutoff:
                
                # Calculate alpha for heat according to Paulson 1970 / Garratt 1994
                alphaHMO2 = (1 - 16 * (self.z[j]-self.disp) / self.L) ** 0.5
                alphaHMO1 = (1 - 16 * self.z_T / self.L) ** 0.5
                # Calculate Psi for heat according to Paulson 1970 / Garratt 1994
                PsiHMO2 = 2 * numpy.log((1 + alphaHMO2) / 2)
                PsiHMO1 = 2 * numpy.log((1 + alphaHMO1) / 2)
                # Calculate wind profile according to Paulson 1970
                self.T_rural[j] = (self.theta_sl / self.RSMParam.vk) * (numpy.log((self.z[j]-self.disp)/self.z_T) - PsiHMO2 + PsiHMO1) + theta_lb

                # Calculate alpha for heat according to Paulson 1970 / Garratt 1994
                alphaHMO2 = (1 - 16 * (self.z[j]-self.disp) / self.L) ** 0.5
                alphaHMO1 = (1 - 16 * self.z_T / self.L) ** 0.5
                # Calculate Psi for heat according to Paulson 1970 / Garratt 1994
                PsiHMO2 = 2 * numpy.log((1 + alphaHMO2) / 2)
                PsiHMO1 = 2 * numpy.log((1 + alphaHMO1) / 2)
                # Calculate wind profile according to Paulson 1970
                self.q_rural[j] = (self.q_sl / self.RSMParam.vk) * (numpy.log((self.z[j]-self.disp) / self.z_T) - PsiHMO2 + PsiHMO1) + MeteoData.q_atm

            # Neutral
            else:
                self.T_rural[j] = theta_lb

                self.q_rural[j] = MeteoData.q_atm

