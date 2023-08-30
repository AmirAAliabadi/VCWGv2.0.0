import os
import numpy
import math
from pprint import pprint
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from Soil_Functions import Soil_Calculations
from scipy.integrate import odeint
from Resistance_Functions import Ressitance_Calculations
from Radiation_Functions import RadiationFunctions
from scipy.optimize import least_squares
import copy

"""
Compute energy balance at the roof
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: January 2021
"""


class EnergyBalanceRoof_Def(object):

    def __init__(self,CiCO2_Sun_init,CiCO2_Shade_init):

        class SWR_Def():
            pass
        self.SWR = SWR_Def()
        self.SWR.SWRabsRoofImp = 0
        self.SWR.SWRoutRoofImp = 0
        self.SWR.SWRinRoofImp = 0
        self.SWR.SWREBRoofImp = 0
        self.SWR.SWRabsRoofVeg = 0
        self.SWR.SWRoutRoofVeg = 0
        self.SWR.SWRinRoofVeg = 0
        self.SWR.SWREBRoofVeg = 0
        self.SWR.SWRabsTotalRoof = 0
        self.SWR.SWRoutTotalRoof = 0
        self.SWR.SWRinTotalRoof = 0
        self.SWR.SWREBTotalRoof = 0

        class LWR_Def():
            pass
        self.LWR = LWR_Def()
        self.LWR.LWRabsRoofImp = 0
        self.LWR.LWRoutRoofImp = 0
        self.LWR.LWRinRoofImp = 0
        self.LWR.LWREBRoofImp = 0
        self.LWR.LWRabsRoofVeg = 0
        self.LWR.LWRoutRoofVeg = 0
        self.LWR.LWRinRoofVeg = 0
        self.LWR.LWREBRoofVeg = 0
        self.LWR.LWRabsTotalRoof = 0
        self.LWR.LWRoutTotalRoof = 0
        self.LWR.LWRinTotalRoof = 0
        self.LWR.LWREBTotalRoof = 0

        class Hflux_Def():
            pass
        self.Hflux = Hflux_Def()
        self.Hflux.HfluxRoofVeg = 0
        self.Hflux.HfluxRoofImp = 0
        self.Hflux.HfluxRoof = 0

        class LEflux_Def():
            pass
        self.LEflux = LEflux_Def()
        self.LEflux.LEfluxRoofImp = 0
        self.LEflux.LEfluxRoofVegInt = 0
        self.LEflux.LEfluxRoofVegPond = 0
        self.LEflux.LEfluxRoofVegSoil = 0
        self.LEflux.LTEfluxRoofVeg = 0
        self.LEflux.LEfluxRoofVeg = 0
        self.LEflux.LEfluxRoof = 0

        class Eflux_Def():
            pass
        self.Eflux = Eflux_Def()
        self.Eflux.EfluxRoofImp = 0
        self.Eflux.EfluxRoofVegInt = 0
        self.Eflux.EfluxRoofVegPond = 0
        self.Eflux.EfluxRoofVegSoil = 0
        self.Eflux.TEfluxRoofVeg = 0
        self.Eflux.EfluxRoofVeg = 0
        self.Eflux.EfluxRoof = 0

        class Gflux_Def():
            pass
        self.Gflux = Gflux_Def()
        self.Gflux.GfluxRoofImp = 0
        self.Gflux.GfluxRoofVeg = 0
        self.Gflux.GfluxRoof = 0

        class Res_Def():
            pass
        self.Res = Res_Def()
        self.Res.raRooftoAtm = 0
        self.Res.rb_LRoof = 0
        self.Res.rap_LRoof = 0
        self.Res.r_soilRoof = 0
        self.Res.rs_sunRoof = 0
        self.Res.rs_shdRoof = 0

        class CiCO2_Def():
            pass
        self.CiCO2 = CiCO2_Def()
        self.CiCO2.CiCO2LeafRoofVegSun = CiCO2_Sun_init
        self.CiCO2.CiCO2LeafRoofVegShd = CiCO2_Shade_init

        class Src_Def():
            pass
        self.Src = Src_Def()
        self.Src.thb = None
        self.Src.qhb = None

    def EBSolver_Roof(self, TemperatureR,MeteoData,Int,ExWater,Vwater,Owater,SoilPotW,Geometry_m,FractionsRoof,ParSoilRoof,
                      PropOpticalRoof,ParVegRoof,ParCalculation,VerticalProfUrban):
        """
        ------
        INPUT:
        TemperatureR: Temperature of the exterior surfaces [K]
        MeteoData: Forcing variables
        Int: Intercepted water
        ExWater: Maximum extractable water [mm m^2 m^-2 ground s^-1]
        Vwater: Volume of water in each soil layer per unit area [mm]
        Owater: Water content (soil moisture) in each layer [-]
        SoilPotW: Soil Water Potential for Second Layer of Vegetation [MPa]
        Geometry_m: Geometric parameters
        FractionsRoof: Fractions of roof covered by vegetation and impervious surface
        ParSoilRoof: Soil parameters
        PropOpticalRoof: Optical properties of the roof (albedo and emissivity)
        ParVegRoof: Roof vegetation parameters
        ParCalculation: General calculation parameters
        VerticalProfUrban: Vertical profiles of variables obtained from 1-D model
        -------
        OUTPUT:
        SWR: Shortwave radiation
             SWRabsRoofImp: Net shortwave radiation at the impervious surface [W m^-2]
             SWRoutRoofImp: Outgoing shortwave radiation from the the impervious surface [W m^-2]
             SWRinRoofImp: Incoming shortwave radiation to the the impervious surface [W m^-2]
             SWREBRoofImp: Shortwave radiation energy balance at the impervious surface [W m^-2]
             SWRabsRoofVeg: Net shortwave radiation at the vegetated surface [W m^-2]
             SWRoutRoofVeg: Outgoing shortwave radiation from the the vegetated surface [W m^-2]
             SWRinRoofVeg: Incoming shortwave radiation to the the vegetated surface [W m^-2]
             SWREBRoofVeg: Shortwave radiation energy balance at the vegetated surface [W m^-2]
             SWRabsTotalRoof: Total net shortwave radiation [W m^-2]
             SWRoutTotalRoof: Total outgoing shortwave radiation [W m^-2]
             SWRinTotalRoof: Total incoming shortwave radiation [W m^-2]
             SWREBTotalRoof: Total shortwave radiation energy balance [W m^-2]
        LWR: Longwave radiation
             LWRabsRoofImp: Net longwave radiation at the impervious surface [W m^-2]
             LWRoutRoofImp: Outgoing longwave radiation from the the impervious surface [W m^-2]
             LWRinRoofImp: Incoming longwave radiation to the the impervious surface [W m^-2]
             LWREBRoofImp: Longwave radiation energy balance at the impervious surface [W m^-2]
             LWRabsRoofVeg: Net longwave radiation at the vegetated surface [W m^-2]
             LWRoutRoofVeg: Outgoing longwave radiation from the the impervious surface [W m^-2]
             LWRinRoofVeg: Incoming longwave radiation to the the vegetated surface [W m^-2]
             LWREBRoofVeg: Longwave radiation energy balance at the vegetated surface [W m^-2]
             LWRabsTotalRoof: Total net longwave radiation [W m^-2]
             LWRoutTotalRoof: Total outgoing longwave radiation [W m^-2]
             LWRinTotalRoof: Total incoming longwave radiation [W m^-2]
             LWREBTotalRoof: Total longwave radiation energy balance [W m^-2]
        Hflux: Sensible heat fluxes
               HfluxRoofVeg: Sensible heat flux from vegetated surface [W m^-2]
               HfluxRoofImp: Sensible heat flux from impervious surface [W m^-2]
               HfluxRoof: Total sensible heat flux [W m^-2]
        LEflux: Latent heat fluxes
                LEfluxRoofImp: Latent heat flux from impervious surface [W m^-2]
                LEfluxRoofVegInt: Latent heat flux from intercepted water on the vegetation [W m^-2]
                LEfluxRoofVegPond: Latent heat flux from surface under vegetation [W m^-2]
                LEfluxRoofVegSoil: Latent heat flux from soil [W m^-2]
                LTEfluxRoofVeg: Latent heat flux from transpiration [W m^-2]
                LEfluxRoofVeg: Total latent heat flux from vegetated surface and soil [W m^-2]
                LEfluxRoof: Total latent heat flux [W m^-2]
        Eflux: Evaporative fluxes
               EfluxRoofImp: Evaporative flux from impervious surface [kg m^-2 s^-1]
               EfluxRoofVegInt: Evaporative flux from intercepted water on the vegetation [kg m^-2 s^-1]
               EfluxRoofVegPond: Evaporative flux from surface under vegetation [kg m^-2 s^-1]
               EfluxRoofVegSoil: Evaporative flux from soil [kg m^-2 s^-1]
               TEfluxRoofVeg: Transpiration flux [kg m^-2 s^-1]
               EfluxRoofVeg: Total evaporative flux from vegetated surface [kg m^-2 s^-1]
               EfluxRoof: Total evaporative flux [kg m^-2 s^-1]
        Gflux: Ground (conductive) heat fluxes
               GfluxRoofImp: Conductive heat flux from impervious surface [W m^-2]
               GfluxRoofVeg: Conductive heat flux from vegetated surface [W m^-2]
               GfluxRoof: Total conductive heat flux [W m^-2]
        CiCO2: Leaf Interior CO2 concentration
               CiCO2LeafRoofVegSun: Sunlit leaf Interior CO2 concentration [umolCO2 mol^-1]
               CiCO2LeafRoofVegShd: Shaded leaf Interior CO2 concentration [umolCO2 mol^-1]
        Src: Sink/Source terms in 1-D model
             thb: Kinematic sensible heat flux [K m s^-1]
             qhb: Evaporative flux [kg kg^-1 m s^-1]
        Res: Resistances
             raRooftoAtm: Aerodynamic resistance [s m^-1]
             rb_LRoof: Leaf boundary layer resistance [s m^-1]
             r_soilRoof: Soil resistance [s m^-1]
             rs_sunRoof: Stomatal resistence of the sunlit vegetation [s m^-1]
             rs_shdRoof: Stomatal resistence of the shaded vegetation [s m^-1]
        """

        # Shortwave radiation
        # Absorbed direct shortwave radiation by vegetated roof [W m^-2]
        SWRabs_dir_veg = (1 - PropOpticalRoof.aveg) * MeteoData.SW_dir
        # Absorbed diffuse shortwave radiation by vegetated roof [W m^-2]
        SWRabs_diff_veg = (1 - PropOpticalRoof.aveg) * MeteoData.SW_diff
        # Absorbed total shortwave radiation by vegetated roof [W m^-2]
        SWR_abs_roofveg = (1 - PropOpticalRoof.aveg) * (MeteoData.SW_dir + MeteoData.SW_diff)
        # Absorbed total shortwave radiation by impervious roof [W m^-2]
        SWR_abs_roofimp = (1 - PropOpticalRoof.aimp) * (MeteoData.SW_dir + MeteoData.SW_diff)
        # Absorbed total shortwave radiation by the roof [W m^-2]
        SWR_abs_totalroof = SWR_abs_roofimp * FractionsRoof.fimp + SWR_abs_roofveg * FractionsRoof.fveg

        # Outgoing total shortwave radiation from the vegetated roof [W m^-2]
        SWR_out_roofveg = PropOpticalRoof.aveg * (MeteoData.SW_dir + MeteoData.SW_diff)
        # Outgoing total shortwave radiation from the impervious roof [W m^-2]
        SWR_out_roofimp = PropOpticalRoof.aimp * (MeteoData.SW_dir + MeteoData.SW_diff)
        # Outgoing total shortwave radiation from the total roof [W m^-2]
        SWR_out_totalroof = SWR_out_roofimp * FractionsRoof.fimp + SWR_out_roofveg * FractionsRoof.fveg

        # Incoming total shortwave radiation to the impervious roof [W m^-2]
        SWR_in_roofimp = (MeteoData.SW_dir + MeteoData.SW_diff)
        # Incoming total shortwave radiation to the vegetated roof [W m^-2]
        SWR_in_roofveg = (MeteoData.SW_dir + MeteoData.SW_diff)
        # Incoming total shortwave radiation to the total roof [W m^-2]
        SWR_in_totalroof = (MeteoData.SW_dir + MeteoData.SW_diff)

        # Verify conservation of shortwave radiation fluxes (must be zero)
        SWR_EB_roofimp = SWR_in_roofimp - SWR_out_roofimp - SWR_abs_roofimp
        SWR_EB_roofveg = SWR_in_roofveg - SWR_out_roofveg - SWR_abs_roofveg
        SWR_EB_totalroof = SWR_in_totalroof - SWR_out_totalroof - SWR_abs_totalroof

        # Longwave radiation
        # Boltzman constant [J K^-1]
        bolzm = 5.67 * 10 ** (-8)
        # Total absorbed longwave radiation by vegetated roof [W m^-2]
        LWR_abs_roofveg = PropOpticalRoof.eveg*(MeteoData.LWR - bolzm*(TemperatureR[1])**4)
        # Total absorbed longwave radiation by impervious roof [W m^-2]
        LWR_abs_roofimp = PropOpticalRoof.eimp*(MeteoData.LWR - bolzm*(TemperatureR[0])**4)
        # Total absorbed longwave radiation by the roof [W m^-2]
        LWR_abs_totalroof = LWR_abs_roofimp * FractionsRoof.fimp + LWR_abs_roofveg * FractionsRoof.fveg

        # Total outgoing longwave radiation by vegetated roof [W m^-2]
        LWR_out_roofveg = PropOpticalRoof.eveg*bolzm*(TemperatureR[1])**4 + (1-PropOpticalRoof.eveg)*MeteoData.LWR
        # Total outgoing longwave radiation by impervious roof [W m^-2]
        LWR_out_roofimp = PropOpticalRoof.eimp*bolzm*(TemperatureR[0])**4 + (1-PropOpticalRoof.eimp)*MeteoData.LWR
        # Total outgoing longwave radiation from the roof [W m^-2]
        LWR_out_totalroof = LWR_out_roofimp * FractionsRoof.fimp + LWR_out_roofveg * FractionsRoof.fveg

        # Total incoming longwave radiation to impervious roof [W m^-2]
        LWR_in_roofimp = MeteoData.LWR
        # Total incoming longwave radiation to vegetated roof [W m^-2]
        LWR_in_roofveg = MeteoData.LWR
        # Total incoming longwave radiation to the roof [W m^-2]
        LWR_in_totalroof = MeteoData.LWR

        # Verify conservation of longwave radiation fluxes (must be zero)
        LWR_EB_roofimp = LWR_in_roofimp - LWR_out_roofimp - LWR_abs_roofimp
        LWR_EB_roofveg = LWR_in_roofveg - LWR_out_roofveg - LWR_abs_roofveg
        LWR_EB_totalroof = LWR_in_totalroof - LWR_out_totalroof - LWR_abs_totalroof

        # Calculate sensible and latent heat
        Hroof_imp, Hroof_veg, Eroof_imp, Eroof_veg, Eroof_ground, Eroof_soil, TEroof_veg, LEroof_imp, LEroof_veg, LEroof_ground, \
        LEroof_soil, LTEroof_veg, Ci_sun_roof, Ci_shd_roof, ra, rb_L, r_soil, rs_sun, rs_shd, thb_roof,qhb_roof = \
            self.HeatFlux_roof_1D(TemperatureR, MeteoData, ParVegRoof, FractionsRoof, Geometry_m,ParSoilRoof,
                                  ParCalculation, SoilPotW, Owater, Vwater, ExWater, Int, self.CiCO2,SWRabs_dir_veg,
                                  SWRabs_diff_veg, VerticalProfUrban)
        # Calculate total sensible heat flux [W m^-2]
        HfluxRoof = Hroof_imp * FractionsRoof.fimp + Hroof_veg * FractionsRoof.fveg

        # Total latent heat flux from vegetation [W m^-2]
        LEfluxRoofVeg = LEroof_veg + LEroof_ground + LEroof_soil + LTEroof_veg
        # Total latent heat flux [W m^-2]
        LEfluxRoof = LEroof_imp * FractionsRoof.fimp + LEfluxRoofVeg * FractionsRoof.fveg
        # Total evaporation flux from vegetation [kg m^-2 s^-1]
        EfluxRoofVeg = Eroof_veg + Eroof_ground + Eroof_soil + TEroof_veg
        # Total evaporation flux [kg m^-2 s^-1]
        EfluxRoof = Eroof_imp * FractionsRoof.fimp + EfluxRoofVeg * FractionsRoof.fveg

        # Save sink/source terms, which will be used in 1-D model for specific humidity (evaporative flux) [kg kg^-1 m s^-1]
        # and for temperature (sensible kinematic heat flux) [K m s^-1]
        self.Src.thb = thb_roof
        self.Src.qhb = qhb_roof
        # Radiation fluxes [W m^-2]
        self.SWR.SWRabsRoofImp = SWR_abs_roofimp
        self.SWR.SWRoutRoofImp = SWR_out_roofimp
        self.SWR.SWRinRoofImp = SWR_in_roofimp
        self.SWR.SWREBRoofImp = SWR_EB_roofimp
        self.LWR.LWRabsRoofImp = LWR_abs_roofimp
        self.LWR.LWRoutRoofImp = LWR_out_roofimp
        self.LWR.LWRinRoofImp = LWR_in_roofimp
        self.LWR.LWREBRoofImp = LWR_EB_roofimp
        self.SWR.SWRabsRoofVeg = SWR_abs_roofveg
        self.SWR.SWRoutRoofVeg = SWR_out_roofveg
        self.SWR.SWRinRoofVeg = SWR_in_roofveg
        self.SWR.SWREBRoofVeg = SWR_EB_roofveg
        self.LWR.LWRabsRoofVeg = LWR_abs_roofveg
        self.LWR.LWRoutRoofVeg = LWR_out_roofveg
        self.LWR.LWRinRoofVeg = LWR_in_roofveg
        self.LWR.LWREBRoofVeg = LWR_EB_roofveg
        self.SWR.SWRabsTotalRoof = SWR_abs_totalroof
        self.SWR.SWRoutTotalRoof = SWR_out_totalroof
        self.SWR.SWRinTotalRoof = SWR_in_totalroof
        self.SWR.SWREBTotalRoof = SWR_EB_totalroof
        self.LWR.LWRabsTotalRoof = LWR_abs_totalroof
        self.LWR.LWRoutTotalRoof = LWR_out_totalroof
        self.LWR.LWRinTotalRoof = LWR_in_totalroof
        self.LWR.LWREBTotalRoof = LWR_EB_totalroof
        # Sensible heat fluxes [W m^-2]
        self.Hflux.HfluxRoofVeg = Hroof_veg
        self.Hflux.HfluxRoofImp = Hroof_imp
        self.Hflux.HfluxRoof = HfluxRoof
        # Latent heat fluxes [W m^-2]
        self.LEflux.LEfluxRoofImp = LEroof_imp
        self.LEflux.LEfluxRoofVegInt = LEroof_veg
        self.LEflux.LEfluxRoofVegPond = LEroof_ground
        self.LEflux.LEfluxRoofVegSoil = LEroof_soil
        self.LEflux.LTEfluxRoofVeg = LTEroof_veg
        self.LEflux.LEfluxRoofVeg = LEfluxRoofVeg
        self.LEflux.LEfluxRoof = LEfluxRoof
        # Evaporative fluxes [kg m^-2 s^-1]
        self.Eflux.EfluxRoofImp = Eroof_imp
        self.Eflux.EfluxRoofVegInt = Eroof_veg
        self.Eflux.EfluxRoofVegPond = Eroof_ground
        self.Eflux.EfluxRoofVegSoil = Eroof_soil
        self.Eflux.TEfluxRoofVeg = TEroof_veg
        self.Eflux.EfluxRoofVeg = EfluxRoofVeg
        self.Eflux.EfluxRoof = EfluxRoof
        # Ground heat fluxes [W m^-2]
        self.Gflux.GfluxRoofImp = SWR_abs_roofimp + LWR_abs_roofimp - Hroof_imp - LEroof_imp
        self.Gflux.GfluxRoofVeg = SWR_abs_roofveg + LWR_abs_roofveg - Hroof_veg - LEroof_veg - LEroof_ground - LEroof_soil - LTEroof_veg
        self.Gflux.GfluxRoof = FractionsRoof.fimp * self.Gflux.GfluxRoofImp + FractionsRoof.fveg * self.Gflux.GfluxRoofVeg
        # Resistances
        self.Res.raRooftoAtm = ra
        self.Res.rb_LRoof = rb_L
        self.Res.r_soilRoof = r_soil
        self.Res.rs_sunRoof = rs_sun
        self.Res.rs_shdRoof = rs_shd
        # Leaf Interior CO2 mixing ratio [umolCO2 mol^-1]
        self.CiCO2.CiCO2LeafRoofVegSun = Ci_sun_roof
        self.CiCO2.CiCO2LeafRoofVegShd = Ci_shd_roof

    def HeatFlux_roof_1D(self,TemperatureR,MeteoData,ParVegRoof,FractionsRoof,Gemeotry_m,ParSoilRoof,ParCalculation,
                         SoilPotW,Owater,Vwater,ExWater,Int,CiCO2Leaf,SWRabs_dir,SWRabs_diff,VerticalProfUrban):
        """
        ------
        INPUT:
        TemperatureR: Temperature of the exterior surfaces [K]
        MeteoData: Forcing variables
        ParVegRoof: Roof vegetation parameters
        FractionsRoof: Fractions of roof covered by vegetation and impervious surface
        Gemeotry_m: Geometric parameters
        ParSoilRoof: Roof soil parameters
        ParCalculation: General calculation parameters
        SoilPotW: Soil Water Potential for Second Layer of Vegetation [MPa]
        Owater: Water content (soil moisture) in each layer [-]
        Vwater: Volume of water in each soil layer per unit area [mm]
        ExWater: Maximum extractable water [mm m^2 m^-2 ground s^-1]
        Int: Intercepted water [mm]
        CiCO2Leaf: Leaf Interior CO2 concentration [umolCO2 mol^-1]
        SWRabs_dir: Absorbed direct shortwave radiation by vegetated roof [W m^-2]
        SWRabs_diff: Absorbed diffuse shortwave radiation by vegetated roof [W m^-2]
        VerticalProfUrban: Vertical profiles of variables obtained from 1-D model
        -------
        OUTPUT:
        Hroof_imp: Sensible heat flux from impervious surface [W m^-2]
        Hroof_veg: Sensible heat flux from vegetated surface [W m^-2]
        Eroof_imp: Evaporative flux from impervious roof [kg m^-2 s^-1]
        Eroof_veg: Evaporative flux from intercepted water by vegetation [kg m^-2 s^-1]
        Eroof_ground: Evaporative flux from runon water on ground under vegetation [kg m^-2 s^-1]
        Eroof_soil: Evaporative flux from soil under vegetation (first layer) [kg m^-2 s^-1]
        TEroof_veg: Transpiration flux [kg m^-2 s^-1]
        LEroof_imp: Latent heat flux from impervious surface [W m^-2]
        LEroof_veg: Latent heat flux from intercepted water on the vegetation [W m^-2]
        LEroof_ground: Latent heat flux from surface under vegetation [W m^-2]
        LEroof_soil: Latent heat flux from soil [W m^-2]
        LTEroof_veg: Latent heat flux from transpiration [W m^-2]
        Ci_sun_roof: Sunlit leaf Interior CO2 concentration [umolCO2 mol^-1]
        Ci_shd_roof: Shaded leaf Interior CO2 concentration [umolCO2 mol^-1]
        ra: Aerodynamic resistance [s m^-1]
        rb: Leaf boundary layer resistance [s m^-1]
        r_soil: Soil resistance [s m^-1]
        rs_sun: Stomatal resistence of the sunlit vegetation [s m^-1]
        rs_shd: Stomatal resistence of the shaded vegetation [s m^-1]
        thb_roof: Evaporative flux [kg kg^-1 m s^-1]
        qhb_roof: Kinematic sensible heat flux [K m s^-1]
        """

        Troof_imp = TemperatureR[0]              # Impervious surface temperature [K]
        Troof_veg = TemperatureR[1]              # Vegetated surface temperature [K]
        T_above_canyon = copy.copy(VerticalProfUrban.th[Gemeotry_m.nz_u])  # Air temperature just above the roof level [K]
        Pre = copy.copy(VerticalProfUrban.presProf[Gemeotry_m.nz_u])       # Air pressure just above the roof level [Pa]
        q_above_canyon = copy.copy(VerticalProfUrban.qn[Gemeotry_m.nz_u])  # Specific humidity just above the roof level [K]
        ea = (q_above_canyon*Pre)/(0.621945 + q_above_canyon)              # Vapor pressure [Pa]
        esat_Tatm = 611*numpy.exp(17.27*(T_above_canyon-273.16) / (237.3+(T_above_canyon-273.16))) # vapor pressure saturation [Pa]
        LAI_roof = copy.copy(ParVegRoof.LAI)
        hc_roof = copy.copy(ParVegRoof.hc)
        d_leaf_roof = copy.copy(ParVegRoof.d_leaf)
        SPAR = copy.copy(ParSoilRoof.SPAR)
        Kbot = copy.copy(ParSoilRoof.Kbot)
        Zs = copy.copy(ParSoilRoof.Zs)

        ResistanceCal = Ressitance_Calculations()
        SoilCal = Soil_Calculations()

        # specific heat air  [J kg^-1 K^-1]
        cp_atm = 1005 + (((T_above_canyon - 273.15) + 23.15) ** 2) / 3364
        # dry air density at atmosphere [kg m^-3]
        rho_atm = copy.copy(VerticalProfUrban.rho[Gemeotry_m.nz_u])
        # Latent heat of vaporization/condensation [J kg^-1]
        L_heat = 1000 * (2501.3 - 2.361 * (T_above_canyon - 273.15))

        # vapor pressure saturation at Troof_imp [Pa]
        esat_T_rimp = 611 * numpy.exp(17.27 * (Troof_imp - 273.16) / (237.3 + (Troof_imp - 273.16)))
        # Saturated specific humidity at Troof_imp [kg kg^-1]
        qsat_T_rimp = (0.622 * esat_T_rimp) / (Pre - 0.378 * esat_T_rimp)
        # vapor pressure saturation at Troof_veg [Pa]
        esat_T_rveg = 611 * numpy.exp(17.27 * (Troof_veg - 273.16) / (237.3 + (Troof_veg - 273.16)))
        # Saturated specific humidity at Troof_veg [kg kg^-1]
        qsat_T_rveg = (0.622 * esat_T_rveg) / (Pre - 0.378 * esat_T_rveg)
        # Average temperature of the roof [K]
        Troof = FractionsRoof.fveg*Troof_veg + FractionsRoof.fimp*Troof_imp

        # Ds = Vapor Pressure Deficit [Pa]
        Ds_atm = esat_Tatm - ea

        # Partitioning of radiation into sunlit and shaded area
        Fsun = (1.0 - numpy.exp(-ParVegRoof.Kopt * (LAI_roof))) / (ParVegRoof.Kopt * (LAI_roof))
        if Fsun < 0.01:
            Fsun = 0
        if Fsun > 1:
            Fsun = 1
        Fshd = 1 - Fsun

        # absorbed direct and diffuse shortwave radiation of the sunlit surface [W m^-2]
        PAR_sun = SWRabs_dir + Fsun * SWRabs_diff
        # absorbed direct and diffuse shortwave radiation of the shaded surface  [W m^-2]
        PAR_shd = Fshd * SWRabs_diff

        # Fraction of vegetation covered by intercepted water (Deardorff(1978))
        # Maximum interception capacity [mm]
        In_max_veg = ParSoilRoof.Sp_In * (LAI_roof + ParVegRoof.SAI)
        # Fraction of vegetation covered by intercepted water
        dw_veg_roof = min(1, (Int.IntRoofVegPlant / In_max_veg) ** (2 / 3))

        # Calculation of resistances
        if FractionsRoof.fveg == 0:
            hc_roof = 0
            LAI_roof = 0
            d_leaf_roof = 0

        if FractionsRoof.fimp == 0:
            Croof = 0
        else:
            Croof = 1

        # Roughness height and length
        zom, zoh, _zom_ground_, _zoh_ground_, disp_h, _zom_H_, _zom_L_, _zoh_H_, _zoh_L_, _d_H_, _d_L_, _zom_other_ = \
            ResistanceCal.Urban_roughness(hc_roof, 0, 0, 0, Croof)

        # Wind profile
        # To use aerodynamic resistance we need to extract at the elevation of typically 2[m] away from the surface
        u_ref, u_Hveg = ResistanceCal.WindProfile_Roof(Gemeotry_m.Height_canyon, hc_roof, VerticalProfUrban,Gemeotry_m)

        # Calculate aerodynamic resistance [s m^-1]
        ra = ResistanceCal.Roof_Aerodynamic_Resistance_1D(VerticalProfUrban,Gemeotry_m,zom,Troof)

        # Calculate soil resistance
        _Zs_, _dz_, _ms_, Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, _EvL_Zs_, _Inf_Zs_, _RfH_Zs_, Rf_Zs, _Zinf_, _Kbot_, \
        _Slo_pot_, _Dz_, _aR_, _aTop_, _rsd_, _lan_dry_, _lan_s_, _cv_s_ = \
            SoilCal.Soil_Parameters_Total(ParSoilRoof.Pcla,ParSoilRoof.Psan,ParSoilRoof.Porg,ParSoilRoof.Kfc,ParSoilRoof.Phy,SPAR,Kbot,
                                          ParVegRoof.CASE_ROOT,ParVegRoof.CASE_ROOT,0,ParVegRoof.ZR95,0,ParVegRoof.ZR50,0,ParVegRoof.ZRmax,Zs)

        r_soil, _b_soil_, alp_soil = \
            ResistanceCal.Soil_Resistance(Troof_veg - 273.15, Pre / 100, u_ref, ea, Int.IntRoofVegGround, Owater[0], Ks_Zs[0],
                                          Osat[0], Ohy[0], L[0], Pe[0], O33[0], alpVG[0], nVG[0], SPAR)

        # sensible heat from impervious area [W m^-2]
        Hroof_imp = cp_atm * rho_atm * (Troof_imp - T_above_canyon) / ra
        # Potential evaporation from runon water on impervious area [kg m^-2 s^-1]
        Eroof_imp_pot = rho_atm * (qsat_T_rimp - q_above_canyon) / ra
        # Potential evaporation from first soil layer [kg m^-2 s^-1]
        Eroof_soil_pot = rho_atm * alp_soil * (qsat_T_rveg - q_above_canyon) / (ra + r_soil)

        if (LAI_roof > 0):

            # Leaf boundary resistance
            rb_L = ResistanceCal.Leaf_Boundary_Resistance(u_Hveg, Troof_veg - 273.15, T_above_canyon - 273.15, hc_roof,d_leaf_roof, LAI_roof,
                                                          MeteoData.Zatm, disp_h, zom)

            rs_sun, rs_shd, Ci_sun, Ci_shd, _An_, _Rdark_, _Lpho_, _SIF_, _DCi_ = \
                ResistanceCal.Canopy_Resistance_An_Evolution(PAR_sun,PAR_shd,LAI_roof,ParVegRoof.Kopt, ParVegRoof.Knit, Fsun,Fshd,
                                                             CiCO2Leaf.CiCO2LeafRoofVegSun,CiCO2Leaf.CiCO2LeafRoofVegShd,MeteoData.Catm_CO2,
                                                             ra,rb_L,Troof_veg-273.15,T_above_canyon-273.15,Pre/100,Ds_atm,SoilPotW,ParVegRoof.Psi_sto_50,
                                                             ParVegRoof.Psi_sto_00,ParVegRoof.CT,ParVegRoof.Vmax,ParVegRoof.DSE,ParVegRoof.Ha,ParVegRoof.FI,
                                                             MeteoData.Catm_O2,ParVegRoof.Do,ParVegRoof.a1,ParVegRoof.go,ParVegRoof.e_rel,ParVegRoof.e_relN,
                                                             ParVegRoof.gmes,ParVegRoof.rjv,ParVegRoof.mSl,ParVegRoof.Sl)

            rb = ResistanceCal.Leaf_Boundary_Resistance(u_Hveg, Troof_veg - 273.15, T_above_canyon - 273.15, hc_roof,d_leaf_roof,LAI_roof,
                                                        MeteoData.Zatm, disp_h, zom)

            # sensible heat from vegetated area [W m^-2]
            Hroof_veg = cp_atm * rho_atm * (Troof_veg - T_above_canyon) / (rb / (2 * (LAI_roof + ParVegRoof.SAI)) + ra)

            # Potential evaporation from intercepted water on vegetation [kg m^-2 s^-1]
            Eroof_veg_pot = rho_atm * (qsat_T_rveg - q_above_canyon) / (rb / ((LAI_roof + ParVegRoof.SAI) * dw_veg_roof) + ra)

            TEroof_veg_sun_pot = rho_atm * (qsat_T_rveg - q_above_canyon) / (rb / ((LAI_roof) * (1 - dw_veg_roof)) + ra + rs_sun /
                                                                ((LAI_roof) * Fsun * (1 - dw_veg_roof)))
            TEroof_veg_shd_pot = rho_atm * (qsat_T_rveg - q_above_canyon) / (rb / ((LAI_roof) * (1 - dw_veg_roof)) + ra + rs_shd /
                                                                ((LAI_roof) * Fshd * (1 - dw_veg_roof)))
            # Total transpiration from vegetation  [kg m^-2 s^-1]
            TEroof_veg_pot = TEroof_veg_sun_pot + TEroof_veg_shd_pot

        else:
            rs_sun = numpy.NaN
            rs_shd = numpy.NaN
            rb = numpy.NaN
            Ci_sun = 0
            Ci_shd = 0

            Hroof_veg = 0
            Eroof_veg_pot = 0
            TEroof_veg_sun_pot = 0
            TEroof_veg_shd_pot = 0
            TEroof_veg_pot = 0

        # Ci = Leaf Interior  CO2 concentration [umolCO2 mol^-1]
        Ci_sun_roof = Ci_sun
        Ci_shd_roof = Ci_shd

        # Condition that evapotranspiration does not exceed available water
        # Real max water evaporation from interception [kg m^-2 s^-1]
        Eroof_imp = min(Eroof_imp_pot, (Int.IntRoofImp / (1000 * ParCalculation.dts) * ParCalculation.rhow))
        Eroof_veg = min(Eroof_veg_pot, (Int.IntRoofVegPlant / (1000 * ParCalculation.dts) * ParCalculation.rhow))
        Eroof_ground = min(Eroof_soil_pot, (Int.IntRoofVegGround / (1000 * ParCalculation.dts) * ParCalculation.rhow))
        Eroof_soil_pot = Eroof_soil_pot - Eroof_ground

        # Water available for Transpiration and Evaporation in a given time step
        # Water flux in each soil layer [kg m^-2 s^-1]
        Eavail_tm1 = [(Vwater[i] / ParCalculation.dts) * (ParCalculation.rhow/1000) for i in range(len(Vwater))]
        Eroof_soil = min(Eroof_soil_pot, Eavail_tm1[0])
        Eavail_tm1[0] = Eavail_tm1[0] - Eroof_soil

        # Water flux for plant respiration in each soil layer [kg m^-2 s^-1]
        Eavail_plant_tm1 = min([Eavail_tm1[i] * Rf_Zs[i] for i in range(len(Eavail_tm1))], [ExWater[i] * (ParCalculation.rhow / 1000) for i in range(len(ExWater))])
        # Cumulated water flux for plant transpiration [kg m^-2 s^-1]
        Eavail_plant_tm1 = sum(Eavail_plant_tm1)
        TEroof_veg = min(TEroof_veg_pot, Eavail_plant_tm1)

        if FractionsRoof.fveg == 0:
            Hroof_veg = 0
            Eroof_veg = 0
            Eroof_ground = 0
            Eroof_soil = 0
            TEroof_veg = 0
        if FractionsRoof.fimp == 0:
            Hroof_imp = 0
            Eroof_imp = 0

        # Latent heat from runon water on impervious area [W m^-2]
        LEroof_imp = L_heat * Eroof_imp
        # Latent heat from intercepted water on vegetation [W m^-2]
        LEroof_veg = L_heat * Eroof_veg
        # Latent heat from runon water on ground under vegetation [W m^-2]
        LEroof_ground = L_heat * Eroof_ground
        # Latent heat from first soil layer [W m^-2]
        LEroof_soil = L_heat * Eroof_soil
        # Latent heat from vegetation  [W m^-2]
        LTEroof_veg = L_heat * TEroof_veg

        class thb_roof_Def():
            pass
        thb_roof = thb_roof_Def()
        # Calculate kinematic sensible heat flux to be used as source/sink term in 1-D model temperature equation[K m s^-1]
        thb_roof.roof_imp = Hroof_imp/(cp_atm*rho_atm)
        thb_roof.roof_veg = Hroof_veg/(cp_atm*rho_atm)

        class qhb_roof_Def():
            pass
        qhb_roof = qhb_roof_Def()
        # Calculate evaporative flux to be used as source/sink term in 1-D model specific humidity equation [kg kg^-1 m s^-1]
        qhb_roof.roof_imp = LEroof_imp/(L_heat*rho_atm)
        qhb_roof.roof_veg = (LEroof_veg+LEroof_ground+LEroof_soil+LTEroof_veg)/(L_heat*rho_atm)

        return Hroof_imp,Hroof_veg,Eroof_imp,Eroof_veg,Eroof_ground,Eroof_soil,TEroof_veg,LEroof_imp,LEroof_veg,LEroof_ground,\
               LEroof_soil,LTEroof_veg,Ci_sun_roof,Ci_shd_roof,ra,rb,r_soil,rs_sun,rs_shd,thb_roof,qhb_roof
