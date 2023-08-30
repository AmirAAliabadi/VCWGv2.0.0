import numpy
import math
from Radiation_Functions import RadiationFunctions
import copy
from psychrometrics import psychrometrics

"""
Surface energy balance model in the rural area 
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: May 2021
"""

class EnergyBalanceRural_Def(object):

    def __init__(self,Td_initial,CiCO2LeafGroundVegSun_init,CiCO2LeafGroundVegShd_init):

        class EnergyFlux_Def():
            pass
        self.EnergyFlux = EnergyFlux_Def()
        self.EnergyFlux.HfluxRural = 0
        self.EnergyFlux.LEfluxRural = 0
        self.EnergyFlux.SWRabsRural = 0
        self.EnergyFlux.LWRabsRural = 0
        self.EnergyFlux.SWRabsImpRural = 0
        self.EnergyFlux.LWRabsImpRural = 0
        self.EnergyFlux.SWRabsBareRural = 0
        self.EnergyFlux.LWRabsBareRural = 0
        self.EnergyFlux.SWRabsVegRural = 0
        self.EnergyFlux.LWRabsVegRural = 0
        self.EnergyFlux.SWRinRural = 0
        self.EnergyFlux.LWRinRural = 0
        self.EnergyFlux.SWRoutRural = 0
        self.EnergyFlux.LWRoutRural = 0
        self.EnergyFlux.GfluxRural = 0
        self.EnergyFlux.GfluxGroundImp = 0
        self.EnergyFlux.GfluxGroundBare = 0
        self.EnergyFlux.GfluxGroundVeg = 0
        self.EnergyFlux.HfluxGroundImp = 0
        self.EnergyFlux.HfluxGroundBare = 0
        self.EnergyFlux.HfluxGroundVeg = 0
        self.EnergyFlux.LEfluxGroundImp = 0
        self.EnergyFlux.LEfluxGroundBarePond = 0
        self.EnergyFlux.LEfluxGroundBareSoil = 0
        self.EnergyFlux.LEfluxGroundVegInt = 0
        self.EnergyFlux.LEfluxGroundVegPond = 0
        self.EnergyFlux.LEfluxGroundVegSoil = 0
        self.EnergyFlux.LTEfluxGroundVeg = 0
        
        class Eflux_Def():
            pass
        self.Eflux = Eflux_Def
        self.Eflux.EfluxGroundImp = 0
        self.Eflux.EfluxGroundBarePond = 0
        self.Eflux.EfluxGroundBareSoil = 0
        self.Eflux.EfluxGroundBare = 0
        self.Eflux.EfluxGroundVegInt = 0
        self.Eflux.EfluxGroundVegPond = 0
        self.Eflux.EfluxGroundVegSoil = 0
        self.Eflux.TEfluxGroundVeg = 0
        self.Eflux.EfluxGroundVeg = 0
        
        class CiCO2_Def():
            pass
        self.CiCO2 = CiCO2_Def()
        self.CiCO2.CiCO2LeafGroundVegSun = CiCO2LeafGroundVegSun_init
        self.CiCO2.CiCO2LeafGroundVegShd = CiCO2LeafGroundVegShd_init

        class Src_Def():
            pass
        self.Src = Src_Def()
        self.Src.thb = 0
        self.Src.qhb = 0

        self.Td = Td_initial


    def EBSolver_Rural(self,MeteoData,RSMParam,Text,SunPosition,simTime,ParCalculation,RSM,ParThermalGround):
        """
        ------
        INPUT:
        MeteoData: Forcing variables
        RSMParam: Rural model parameters
        Text: Exterior surface temperature [K]
        SunPosition: Sun angles
        simTime: Simulation time parameters
        ParCalculation: General calculation parameters
        RSM: Calculated rural model variables
        -------
        OUTPUT:
        Eflux: Surface heat fluxes (Radiation, sensible, latent, and ground heat fluxes) [W m^-2]
        CiCO2: Leaf Interior  CO2 concentration [umolCO2 mol^-1]
        Src: Sink/source terms at the surface
        Td: Deep soil temperature [K]
        """

        # Calculate radiation [W m^-2]
        RadFun = RadiationFunctions()
        SWR_Rural, LWR_Rural = RadFun.TotalSWR_LWR_Rrual(RSMParam, Text, MeteoData, SunPosition,simTime)

        # Air density [kg m^-3]
        rho_atm = copy.copy(RSM.densityProfC[0])
        # Specific heat air  [J kg^-1 K^-1]
        cp_atm = 1005 + (((RSM.T_rural[0] - 273.15) + 23.15) ** 2) / 3364
        # Latent heat vaporization/condensation [J kg^-1]
        L_heat = 1000 * (2501.3 - 2.361 * (RSM.T_rural[0] - 273.15))

        # Calculate sensible, latent and ground heat fluxes [W m^-2]
        # Option 1 (EB_RuralModel_name): Using "Louis" formulation to calculate sensible heat flux and using Bowen ratio
        # to calculate latent heat flux. Then, ground heat flux is calculated using energy balance at the surface (Louis 1979).
        if RSMParam.EB_RuralModel_name == 'Louis':

            # Calculate sensible heat flux caused by vegetation [W m^-2]
            # Biogenic sensible heat released from plants can be ignored (Plants are cold-blooded)
            vegSens = 0

            # Parameterization using Louis, 1979
            U_nearground = MeteoData.Uatm*math.log(RSMParam.h_temp/(RSMParam.h_obs*RSMParam.z0overh_MOST))/\
                           math.log(RSMParam.h_wind/(RSMParam.h_obs*RSMParam.z0overh_MOST))
            # Calculate total sensible heat flux [W m^-2]
            sens = vegSens + rho_atm*cp_atm*self.Louis_SensHeatFlux(RSMParam.z0_Louis, U_nearground, Text, MeteoData.Tatm,
                                                                    RSMParam.MinWind_rural, RSMParam.h_temp)
            # Calculate latent heat flux [W m^-2]
            # Switch control on the presence or absence of evaporation in the simulation
            lat = sens / RSMParam.Bowen

            # Calculate ground heat flux using energy balance equation
            self.EnergyFlux.GfluxRural = SWR_Rural.SWRabsRural + LWR_Rural.LWRabsRural - sens - lat
            self.EnergyFlux.HfluxRural = sens
            self.EnergyFlux.LEfluxRural = lat

            # Latent heat vaporization/condensation [J kg^-1]
            self.Src.thb = sens / (rho_atm * cp_atm)
            self.Src.qhb = lat / (rho_atm * L_heat)

            # Radiation fluxes [W m^-2]
            self.EnergyFlux.SWRabsRural = copy.copy(SWR_Rural.SWRabsRural)
            self.EnergyFlux.LWRabsRural = copy.copy(LWR_Rural.LWRabsRural)
            self.EnergyFlux.SWRinRural = copy.copy(SWR_Rural.SWRinRural)
            self.EnergyFlux.LWRinRural = copy.copy(LWR_Rural.LWRinRural)
            self.EnergyFlux.SWRoutRural = copy.copy(SWR_Rural.SWRoutRural)
            self.EnergyFlux.LWRoutRural = copy.copy(LWR_Rural.LWRoutRural)

        # Option 2 (EB_RuralModel_name): Using "Penman_Monteith" formulation to calculate latent heat flux and consider
        # fraction of total radiation as ground heat flux. Then, sensible heat flux is calculated using energy balance at the surface (Allen et al. 1998).
        elif RSMParam.EB_RuralModel_name == 'Penman_Monteith':

            # Calculate ground heat flux using Penman-Monteith
            if SWR_Rural.SWRabsRural + LWR_Rural.LWRabsRural > 0:
                MultiplierDay = 0.1
                self.EnergyFlux.GfluxRural = MultiplierDay * (SWR_Rural.SWRabsRural+LWR_Rural.LWRabsRural)
            else:
                MultiplierNight = 0.5
                self.EnergyFlux.GfluxRural = MultiplierNight * (SWR_Rural.SWRabsRural+LWR_Rural.LWRabsRural)

            # Calculate latent heat flux [W m^-2]
            # Using Penman-Monteith formulation
            lat = self.Penman_Monteith(MeteoData.Tatm-273.15,MeteoData.Pre,SWR_Rural.SWRabsRural,LWR_Rural.LWRabsRural,
                                           self.EnergyFlux.GfluxRural,ParCalculation.Elevation,MeteoData.Uatm,RSMParam.h_wind,MeteoData.q_atm)

            # Calculate sensible heat flux using energy balance equation
            sens = SWR_Rural.SWRabsRural + LWR_Rural.LWRabsRural - self.EnergyFlux.GfluxRural - lat
            self.EnergyFlux.HfluxRural = sens
            self.EnergyFlux.LEfluxRural = lat

            # Sink/source term in temperature equation [K m s^-1]
            self.Src.thb = sens / (rho_atm * cp_atm)
            # Sink/source term in temperature equation [kg kg^-1 m s^-1]
            self.Src.qhb = lat / (rho_atm * L_heat)

            # Radiation fluxes [W m^-2]
            self.EnergyFlux.SWRabsRural = copy.copy(SWR_Rural.SWRabsRural)
            self.EnergyFlux.LWRabsRural = copy.copy(LWR_Rural.LWRabsRural)
            self.EnergyFlux.SWRinRural = copy.copy(SWR_Rural.SWRinRural)
            self.EnergyFlux.LWRinRural = copy.copy(LWR_Rural.LWRinRural)
            self.EnergyFlux.SWRoutRural = copy.copy(SWR_Rural.SWRoutRural)
            self.EnergyFlux.LWRoutRural = copy.copy(LWR_Rural.LWRoutRural)

        # Calculate deep soil temperature
        if ParThermalGround.Tdeep_ctrl == 'Force_Restore':
            # Option 1 (Tdeep_ctrl): Using force-restore method
            self.Td = self.Soil_Heat(ParCalculation.dts,Text,numpy.NaN,self.Td,numpy.NaN)[1]
        elif ParThermalGround.Tdeep_ctrl == 'Climate_Data':
            # Option 2 (Tdeep_ctrl): Using climate data
            self.Td = copy.copy(MeteoData.TdeepSoil)

    def Louis_SensHeatFlux(self,z0, Wind, Ts, Ta, windMin, h_temp):
        """
        -----
        INPUT:
        z0: Aerodynamic roughness length [m]
        Wind: Wind speed at the height of measured temperature [m s^-1]
        Ts: Surface temperature [K]
        Ta: Air temperature [K]
        windMin: Minimum wind speed
        h_temp: The height of measured temperature [m]
        -------
        OUTPUT:
        KinematicSensHeatFlux: Kinematic sens heat flux [K m s^-1]
        """

        zz = h_temp
        Ck = 0.4

        Utot = max(Wind, windMin)

        # Compute bulk Richardson number using near surface temperatures
        Ri = 2 * 9.81 * zz * (Ta - Ts) / ((Ta + Ts) * (Utot ** 2))
        # Calculation from Louis, 1979 (eq. 11 and 12)
        b = 9.4
        cm = 7.4
        ch = 5.3
        R = 0.74
        a = Ck / math.log(zz / z0)

        if Ri > 0:
            fh = 1 / ((1 + 0.5 * b * Ri) ** 2)
        else:
            c = b * cm * a * a * (zz / z0) ** 0.5
            c = c * ch / cm
            fh = 1 - b * Ri / (1 + c * (-Ri) ** 0.5)

        KinematicSensHeatFlux = -(a ** 2) * Utot * (Ta - Ts) * fh / R

        return KinematicSensHeatFlux

    def Soil_Heat(self,dt,Ts,Tstm1,Tdptm1,CTt):
        """
        ------
        INPUT:
        dt: Time step [s]
        Ts: Exterior surface temperature [K]
        Tstm1: Exterior surface temperature from two steps back [K]. Note: it is not used, because ground heat flux is calculated somewhere else.
        Tdptm1: Deep soil temperature from previous step [K]
        CTt: Total thermal capacity [K m^2 J^-1]
        -------
        OUTPUT:
        G: Ground heat flux [W m^-2]. Note: it is not used, because ground heat flux is calculated somewhere else.
        Tdp: Deep soil temperature [K]
        """

        # Time canstant [s]
        tau = 86400

        # Temperature Variation [C]
        dTs = Ts - Tstm1

        # Depth Temperature [C]
        Tdp = (1 / (1 + dt / tau)) * (Tdptm1 + (dt / tau) * Ts)

        # Soil Heat Flux [W m^-2]
        G = (1 / CTt) * (2 * numpy.pi * (Ts - Tdp) / tau + dTs / dt)

        return G,Tdp

    def Penman_Monteith(self,T,Pre,SWRabs,LWRabs,G,El,uz,zw,q_atm):
        """
        ------
        INPUT:
        T: Air temperature [C]
        Pre: Air pressure [Pa]
        SWRabs: Net shortwave radiation at the surface [W m^-2]
        LWRabs: Net longwave radiation at the surface [W m^-2]
        G: Ground heat flux [W m^-2]
        El: Elevation above mean sea level [m]
        uz: Forcing wind speed [m s^-1]
        zw: Height of measured wind speed [m]
        q_atm: Air specific humidity [kg kg^-1]
        -------
        OUTPUT:
        LET_os: Latent heat flux [W m^-2]
        """

        # Wind speed at 2 m
        u2 = uz*(4.87/(math.log(67.8*zw-5.42)))

        # Calculate dew point temperature at 2 m [C]
        Td = psychrometrics(T+273.15,q_atm,Pre)[0]

        # Saturation vapor pressure at the mean hourly air temperature [kPa]
        es = 0.6108*numpy.exp(17.27*T/(T+237.3))
        # Actual vapor pressure or saturation vapor pressure at the mean dew point temperature [kPa]
        ea = 0.6108*numpy.exp(17.27*Td/(Td+237.3))

        # Net radiation [MJ m^-2 day^-1]
        Rn = (SWRabs+LWRabs)*3600*1e-6*24

        # Barometric pressure as a function of elevation [kPa]
        Bp = 101.3*((293-0.0065*El)/293)**5.26

        # Latent heat of vaporization [MJ kg^-1]
        Lambda = 2.45

        # Psychrometric constant [kPa C^-1]
        gamma = 0.00163*Bp/Lambda

        # Slope of the saturation vapor pressure curve at mean air temperature [kPa C^-1]
        Delta = (4099*es)/(T+237.3)**2

        # Ground heat flux [MJ m^-2 day^-1]
        G_new = G*1e-6*24*3600

        # Penman evapotranspiration [mm day^-1]
        ET_os = (0.408*Delta*(Rn-G_new) + gamma*(37/(T+273.15))*u2*(es-ea)) / (Delta+gamma*(1+0.34*u2))

        # Density of water [kg m^-3]
        rhow = 1000

        # Latent heat flux [W m^-2]
        LET_os = (Lambda*1e6)*(rhow/(1000*3600*24))*ET_os

        return LET_os

