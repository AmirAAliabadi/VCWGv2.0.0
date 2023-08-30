import numpy
from scipy.optimize import fsolve
from UrbanSurfaceHeatFlux import Surface_HeatFlux
from Radiation_Functions import RadiationFunctions
from scipy.optimize import least_squares
import copy

"""
Compute energy balance at the surfaces in the canyon
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: May 2021
"""

class EnergyBalanceCanyon_Def():

    def __init__(self,CiCO2LeafTreeSun_init,CiCO2LeafTreeShd_init,CiCO2LeafGroundVegSun_init,CiCO2LeafGroundVegShd_init,
                 Tdp_gimp_init,Tdp_gbare_init,Tdp_gveg_init,Ttree_init):

        class Src_Def():
            pass
        self.Src = Src_Def()
        self.Src.thb = None
        self.Src.qhb = None
        self.Src.tvb = None

        class Gflux_Def():
            pass
        self.Gflux = Gflux_Def()
        self.Gflux.GfluxGroundImp = 0
        self.Gflux.GfluxGroundBare = 0
        self.Gflux.GfluxGroundVeg = 0
        self.Gflux.GfluxGround = 0
        self.Gflux.GfluxWallSun = 0
        self.Gflux.GfluxWallShade = 0

        class Hflux_Def():
            pass
        self.Hflux = Hflux_Def()
        self.Hflux.HfluxGroundImp = 0
        self.Hflux.HfluxGroundBare = 0
        self.Hflux.HfluxGroundVeg = 0
        self.Hflux.HfluxWallSun = 0
        self.Hflux.HfluxWallShade = 0
        self.Hflux.HfluxTree = 0
        self.Hflux.HfluxGround = 0
        self.Hflux.HfluxCanyon = 0

        class LEflux_Def():
            pass
        self.LEflux = LEflux_Def()
        self.LEflux.LEfluxGroundImp = 0
        self.LEflux.LEfluxGroundBare = 0
        self.LEflux.LEfluxGroundBarePond = 0
        self.LEflux.LEfluxGroundBareSoil = 0
        self.LEflux.LEfluxGroundVeg = 0
        self.LEflux.LEfluxGroundVegInt = 0
        self.LEflux.LEfluxGroundVegPond = 0
        self.LEflux.LEfluxGroundVegSoil = 0
        self.LEflux.LTEfluxGroundVeg = 0
        self.LEflux.LEfluxGround = 0
        self.LEflux.LEfluxWallSun = 0
        self.LEflux.LEfluxWallShade = 0
        self.LEflux.LEfluxTreeInt = 0
        self.LEflux.LTEfluxTree = 0
        self.LEflux.LEfluxTree = 0
        self.LEflux.LEfluxCanyon = 0

        class Eflux_Def():
            pass
        self.Eflux = Eflux_Def()
        self.Eflux.EfluxGroundImp = 0
        self.Eflux.EfluxGroundBarePond = 0
        self.Eflux.EfluxGroundBareSoil = 0
        self.Eflux.EfluxGroundBare = 0
        self.Eflux.EfluxGroundVegInt = 0
        self.Eflux.EfluxGroundVegPond = 0
        self.Eflux.EfluxGroundVegSoil = 0
        self.Eflux.TEfluxGroundVeg = 0
        self.Eflux.EfluxGroundVeg = 0
        self.Eflux.EfluxTreeInt = 0
        self.Eflux.TEfluxTree = 0
        self.Eflux.EfluxGround = 0
        self.Eflux.EfluxTree = 0
        self.Eflux.EfluxWallSun = 0
        self.Eflux.EfluxWallShade = 0

        class Tdepth_Def():
            pass
        self.Tdepth = Tdepth_Def()
        self.Tdepth.TDampGroundImp = Tdp_gimp_init
        self.Tdepth.TDampGroundBare = Tdp_gbare_init
        self.Tdepth.TDampGroundVeg = Tdp_gveg_init

        self.Ttree = [Ttree_init]

        class SWR_Def():
            pass
        self.SWR = SWR_Def()
        self.SWR.SWRin = None
        self.SWR.SWRout = None
        self.SWR.SWRabs = None
        self.SWR.SWRabsDir = None
        self.SWR.SWRabsDiff = None
        self.SWR.SWREB = None
        self.SWR.SWRabsTotalUrban = None
        self.SWR.SWRinTotalUrban = None
        self.SWR.SWRoutTotalUrban = None
        self.SWR.SWREBTotalUrban = None

        class LWR_Def():
            pass
        self.LWR = LWR_Def()
        self.LWR.LWRin = 0
        self.LWR.LWRout = 0
        self.LWR.LWRabs = 0
        self.LWR.LWREB = 0
        self.LWR.LWRabsTotalUrban = 0
        self.LWR.LWRinTotalUrban = 0
        self.LWR.LWRoutTotalUrban = 0
        self.LWR.LWREBTotalUrban = 0

        class CiCO2_Def():
            pass
        self.CiCO2 = CiCO2_Def()
        self.CiCO2.CiCO2LeafTreeSun = CiCO2LeafTreeSun_init
        self.CiCO2.CiCO2LeafTreeShd = CiCO2LeafTreeShd_init
        self.CiCO2.CiCO2LeafGroundVegSun = CiCO2LeafGroundVegSun_init
        self.CiCO2.CiCO2LeafGroundVegShd = CiCO2LeafGroundVegShd_init

        class Res_Def():
            pass
        self.Res = Res_Def()
        self.Res.rap_can = 0
        self.Res.rap_Htree_In = 0
        self.Res.rb_HGround = 0
        self.Res.rb_LGround = 0
        self.Res.r_soilGroundbare = 0
        self.Res.r_soilGroundveg = 0
        self.Res.rs_sunGround = 0
        self.Res.rs_shdGround = 0
        self.Res.rs_sunTree = 0
        self.Res.rs_shdTree = 0

        class OtherParam_Def():
            pass
        self.OtherParam = OtherParam_Def()
        self.OtherParam.Ri_nearGround = 0
        self.OtherParam.Utotal_nearGround = 0
        self.OtherParam.T_nearGround = 0
        self.OtherParam.Tground = 0
        self.OtherParam.alp_soil_bare = 0
        self.OtherParam.alp_soil_veg = 0
        self.OtherParam.Fsun_L = 0
        self.OtherParam.Fshd_L = 0
        self.OtherParam.dw_L = 0


    def EBSolver_Canyon(self,TemperatureC,TempVec,MeteoData,Int,ExWater,Vwater,Owater,SoilPotW,ViewFactor,Geometry_m,ParTree,
                        geometry,FractionsGround,ParSoilGround,ParInterceptionTree,PropOpticalGround,PropOpticalWall,
                        PropOpticalTree,ParThermalGround,ParVegGround,ParVegTree,SunPosition,Anthropogenic,
                        ParCalculation,VerticalProfUrban,ColParam):
        """
        ------
        INPUT:
        TemperatureC: Temperature of the exterior surfaces [K]
        TempVec: Temperature of the exterior surfaces from two steps back [K]
        MeteoData: Forcing variables
        Int: Water interception [mm]
        ExWater: Extractable water [mm m^2 m^-2 ground s^-1]
        Vwater: [mm]
        Owater: Water Content [-]
        SoilPotW: [MPa]
        ViewFactor: View factors
        Geometry_m: Geometric parameters
        ParTree: Trees parameter
        geometry: Normalized geometric parameters
        FractionsGround: Fractions of ground covered by vegetation, impervious, and bare surface
        ParSoilGround: Soil parameters
        ParInterceptionTree: specific water retained by the tree [mm m^2 VEG area m^-2 plant area]
        PropOpticalGround: Optical properties of the ground (albedo and emissivity)
        PropOpticalWall: Optical properties of the wall (albedo and emissivity)
        PropOpticalTree: Optical properties of the tree (albedo and emissivity)
        ParThermalGround: Thermal properties of the ground (Thermal conductivity[W m^-1 K^-1], Volumetric heat capacity [J m^-3 K^-1])
        ParVegGround: Ground vegetation parameters
        ParVegTree: Trees parameter
        SunPosition: sun angles
        Anthropogenic: Anthropogenic heat
        ParCalculation: General calculation parameters
        VerticalProfUrban: Vertical profile of variables obtained from 1-D model
        ColParam: 1-D model parameters
        -------
        OUTPUT:
        Src: Sink/source terms in column model
             thb: Kinematic sensible heat flux from ground [K m s^-1]
             tvb: Kinematic sensible heat flux from wall [K m s^-1]
             qhb: Evaporative flux [kg kg^-1 m s^-1]
        Gflux: Ground (Conductive) heat fluxes
               GfluxGroundImp: Ground heat flux from impervious surface [W m^-2]
               GfluxGroundBare: Ground heat flux from bare surface [W m^-2]
               GfluxGroundVeg: Ground heat flux from vegetated surface [W m^-2]
               GfluxGround: Total ground heat flux [W m^-2]
               GfluxWallSun: Conductive heat flux from sunlit wall [W m^-2]
               GfluxWallShade: Conductive heat flux from shaded wall [W m^-2]
        Hflux: Sensible heat fluxes
               HfluxGroundImp: Sensible heat flux from impervious surface [W m^-2]
               HfluxGroundBare: Sensible heat flux from bare surface [W m^-2]
               HfluxGroundVeg: Sensible heat flux from vegetated surface [W m^-2]
               HfluxWallSun: Sensible heat flux from sunlit wall [W m^-2]
               HfluxWallShade: Sensible heat flux from shaded wall [W m^-2]
               HfluxTree: Sensible heat flux from tree [W m^-2]
               HfluxGround: Total sensible heat flux from ground [W m^-2]
               HfluxCanyon: Total sensible heat flux in the canyon [W m^-2]
        LEflux: Latent heat fluxes
                LEfluxGroundImp: Latent heat flux from impervious surface [W m^-2]
                LEfluxGroundBare: Total latent heat flux from bare ground [W m^-2]
                LEfluxGroundBarePond: Latent heat flux from bare surface [W m^-2]
                LEfluxGroundBareSoil: Latent heat flux from soil under bare surface [W m^-2]
                LEfluxGroundVeg: Total latent heat flux from vegetated ground [W m^-2]
                LEfluxGroundVegInt: Latent heat flux from intercepted water on the low vegetation [W m^-2]
                LEfluxGroundVegPond: Latent heat flux from surface under low vegetation [W m^-2]
                LEfluxGroundVegSoil: Latent heat flux from soil under vegetated surface [W m^-2]
                LTEfluxGroundVeg: Latent heat flux from transpiration of low vegetation [W m^-2]
                LEfluxGround: Total latent heat flux from ground [W m^-2]
                LEfluxWallSun: Latent heat flux from sunlit wall [W m^-2]
                LEfluxWallShade: Latent heat flux from shaded wall [W m^-2]
                LEfluxTreeInt: Latent heat flux from intercepted water on trees [W m^-2]
                LTEfluxTree: Latent heat flux from transpiration of trees [W m^-2]
                LEfluxTree: Total latent heat flux from trees [W m^-2]
                LEfluxCanyon: Total latent heat flux in the canyon [W m^-2]
        Eflux: Evaporative fluxes
               EfluxGroundImp: Evaporation from intercepted water on impervious surface [kg m^-2 s^-1]
               EfluxGroundBarePond: Evaporation from intercepted water on bare ground [kg m^-2 s^-1]
               EfluxGroundBareSoil: Evaporation from soil under bare ground [kg m^-2 s^-1]
               EfluxGroundBare: Total evaporative flux from bare ground [kg m^-2 s^-1]
               EfluxGroundVegInt: Evaporation from intercepted water on low vegetation [kg m^-2 s^-1]
               EfluxGroundVegPond: Evaporation from intercepted water on ground under vegetation [kg m^-2 s^-1]
               EfluxGroundVegSoil: Evaporation from soil under vegetated ground [kg m^-2 s^-1]
               TEfluxGroundVeg: Transpiration from low vegetation [kg m^-2 s^-1]
               EfluxGroundVeg: Total evaporative flux from vegetated ground [kg m^-2 s^-1]
               EfluxTreeInt: Evaporation from intercepted water on trees [kg m^-2 s^-1]
               TEfluxTree: Transpiration from trees [kg m^-2 s^-1]
               EfluxGround: Total evaporative flux from ground [kg m^-2 s^-1]
               EfluxTree: Total evaporative flux from tree [kg m^-2 s^-1]
               EfluxWallSun: Evaporative flux from sunlit wall [kg m^-2 s^-1]
               EfluxWallShade: Evaporative flux from shaded wall [kg m^-2 s^-1]
        Tdepth: Deep soil temperature
                TDampGroundImp: Deep soil temperature of impervious ground [K]
                TDampGroundBare: Deep soil temperature of bare ground [K]
                TDampGroundVeg: Deep soil temperature of vegetated ground [K]
        SWR: Shortwave radiations
             SWRin: Incoming shortwave radiation [W m^-2]
             SWRout: Outgoing shortwave radiation [W m^-2]
             SWRabs: Net shortwave radiation [W m^-2]
             SWRabsDir: Absorbed direct shortwave radiation [W m^-2]
             SWRabsDiff: Absorbed diffuse shortwave radiation [W m^-2]
             SWREB: Energy balance for shortwave radiation [W m^-2]
             SWRabsTotalUrban: Total absorbed shortwave radiation in the urban [W m^-2]
             SWRinTotalUrban: Total incoming shortwave radiation to the urban [W m^-2]
             SWRoutTotalUrban: Total outgoing shortwave radiation from the urban [W m^-2]
             SWREBTotalUrban: Total energy balance for shortwave radiation [W m^-2]
        LWR: Longwave radiations
             LWRin: Incoming longwave radiation [W m^-2]
             LWRout: Outgoing longwave radiation [W m^-2]
             LWRabs: Net longwave radiation [W m^-2]
             LWREB: Energy balance for longwave radiation [W m^-2]
             LWRabsTotalUrban: Total absorbed longwave radiation in the urban [W m^-2]
             LWRinTotalUrban: Total incoming longwave radiation to the urban [W m^-2]
             LWRoutTotalUrban: Total outgoing longwave radiation from the urban [W m^-2]
             LWREBTotalUrban: Total energy balance for longwave radiation [W m^-2]
        CiCO2: Leaf Interior CO2 concentration
               CiCO2LeafTreeSun: Leaf Interior CO2 concentration of sunlit tree [umolCO2 mol^-1]
               CiCO2LeafTreeShd: Leaf Interior CO2 concentration shaded tree [umolCO2 mol^-1]
               CiCO2LeafGroundVegSun: Leaf Interior CO2 concentration of sunlit low vegetation [umolCO2 mol^-1]
               CiCO2LeafGroundVegShd: Leaf Interior CO2 concentration of shaded low vegetation [umolCO2 mol^-1]
        Res: Resistances
             rap_can: Aerodynamic resistance near the ground [s m^-1]
             rap_Htree_In: Aerodynamic resistance between trees and canyon air [s m^-1]
             rb_HGround: Leaf boundary layer resistance [s m^-1]
             rb_LGround: Leaf boundary layer resistance for low vegetation [s m^-1]
             r_soilGroundbare: Bare soil resistance [s m^-1]
             r_soilGroundveg: Vegetated soil resistance [s m^-1]
             rs_sunGround: Stomatal resistance of sunlit low vegetation [s m^-1]
             rs_shdGround: Stomatal resistance of shaded low vegetation [s m^-1]
             rs_sunTree: Stomatal resistance of sunlit trees [s m^-1]
             rs_shdTree: Stomatal resistance of shaded trees [s m^-1]
        OtherParam: Other parameters
                    Ri_nearGround: Richardson number near ground
                    Utotal_nearGround:
                    T_nearGround:
                    Tground:
                    alp_soil_bare: Relative humidity of bare soil [0-1]
                    alp_soil_veg: Relative humidity of vegetated soil [0-1]
                    Fsun_L: Fraction of sunlit for low vegetation [-]
                    Fshd_L: Fraction of shaded for low vegetation [-]
                    dw_L: Fraction of ground vegetation covered by intercepted water [-]
        """


        RadiationCal = RadiationFunctions()
        SurfaceHeatFlux = Surface_HeatFlux()

        # Shortwave radiation
        SWRin_t, SWRout_t, SWRabs_t, SWRabsDir_t, SWRabsDiff_t, SWREB_t = \
            RadiationCal.TotalSWRabsorbed(geometry, FractionsGround, ParTree, PropOpticalGround, PropOpticalWall,
                                          PropOpticalTree, ParVegTree, MeteoData, SunPosition, ViewFactor)

        # Tree absorbed: conversion from spherical to horizontal projected area
        # In the 2-D radiation model, radiation fluxes are incident on straight lines representing flat surface area
        # To convert these to fluxes on circular lines representing spherical surface area a multiplier should be used
        Line2Circle_Multiplier = numpy.pi
        SWRabs_t.SWRabsTree = SWRabs_t.SWRabsTree * Line2Circle_Multiplier
        SWRabsDir_t.SWRabsTree = SWRabsDir_t.SWRabsTree * Line2Circle_Multiplier
        SWRabsDiff_t.SWRabsTree = SWRabsDiff_t.SWRabsTree * Line2Circle_Multiplier

        # Longwave radiation
        LWRin_t, LWRout_t, LWRabs_t, LWREB_t = \
            RadiationCal.TotalLWRabsorbed(TemperatureC, geometry, MeteoData, FractionsGround, PropOpticalGround,
                                          PropOpticalWall, PropOpticalTree, ParTree, ViewFactor)
        # Tree absorbed: conversion from sphere to horizontal projected area
        LWRabs_t.LWRabsTree = LWRabs_t.LWRabsTree * Line2Circle_Multiplier

        # Calculate deep soil temperature
        # Option 1: Using force-restore method
        if ParThermalGround.Tdeep_ctrl == 'Force_Restore':
            # Impervious ground
            _G1GroundImp_, Tdp_ground_imp = \
                SurfaceHeatFlux.ForceRestore_Conductive_Heat_Imp(TemperatureC,self.Tdepth,TempVec,ParCalculation,ParThermalGround,FractionsGround)
            # Bare ground
            type = 0
            _G1GroundBare_, Tdp_ground_bare = \
                SurfaceHeatFlux.ForceRestore_Conductive_Heat_Soil(TemperatureC,self.Tdepth,Owater,TempVec,ParCalculation,
                                                                  ParSoilGround,ParVegGround,ParVegTree,FractionsGround,type)
            # Vegetated ground
            type = 1
            _G1GroundVeg_, Tdp_ground_veg = \
                SurfaceHeatFlux.ForceRestore_Conductive_Heat_Soil(TemperatureC,self.Tdepth,Owater,TempVec,ParCalculation,
                                                                  ParSoilGround,ParVegGround,ParVegTree,FractionsGround,type)
        # Option 2: Using climate data
        elif ParThermalGround.Tdeep_ctrl == 'Climate_Data':
            Tdp_ground_imp = copy.copy(MeteoData.TdeepSoil)
            Tdp_ground_bare = copy.copy(MeteoData.TdeepSoil)
            Tdp_ground_veg = copy.copy(MeteoData.TdeepSoil)

        # Calculate turbulent heat fluxes from ground and trees to the canyon
        HfluxGroundImp, HfluxGroundBare, HfluxGroundVeg, HfluxTree, Eground_imp_pond, Eground_bare_pond, Eground_bare_soil,\
        Eground_veg_int, Eground_veg_pond, Eground_veg_soil, TEground_veg, E_tree_int, TE_tree, Ebare, Eveg, Etree, LEfluxGroundImp,\
        LEfluxGroundBarePond, LEfluxGroundBareSoil, LEfluxGroundVegInt, LEfluxGroundVegPond, LEfluxGroundVegSoil, LTEfluxGroundVeg,\
        LEfluxTreeInt, LTEfluxTree, LEbare, LEveg, LEtree, Ci_sun_tree, Ci_shd_tree, Ci_sun_ground, Ci_shd_ground, rap_can,\
        rap_Htree_In, rb_H, rb_L, r_soil_bare, r_soil_veg, alp_soil_bare, alp_soil_veg, rs_sun_L, rs_shd_L, rs_sun_H,\
        rs_shd_H, Fsun_L, Fshd_L, dw_L = \
            SurfaceHeatFlux.HeatFlux_ground_1D(TemperatureC, MeteoData, Geometry_m, geometry, FractionsGround,ParTree,
                                               ParVegGround, ParVegTree, ParSoilGround, SoilPotW, Owater, Vwater, ExWater,
                                               Int, self.CiCO2, ParInterceptionTree, ParCalculation,SWRabsDir_t.SWRabsTree,
                                               SWRabsDiff_t.SWRabsTree, SWRabsDir_t.SWRabsGroundVeg,SWRabsDiff_t.SWRabsGroundVeg,
                                               VerticalProfUrban,ColParam)

        # Calculate turbulent heat fluxes from sunlit and shaded wall to canyon
        HfluxWallSun, HfluxWallShade, Ewsun, Ewshade, LEwsun, LEwshade = \
            SurfaceHeatFlux.HeatFlux_wall(TemperatureC, Geometry_m,VerticalProfUrban,ParCalculation,ColParam)

        # Calculate total heat flux from ground [W m^-2]
        HfluxGround = HfluxGroundImp*FractionsGround.fimp + HfluxGroundBare*FractionsGround.fbare + HfluxGroundVeg*FractionsGround.fveg

        LEfluxGroundBare = LEfluxGroundBarePond + LEfluxGroundBareSoil
        LEfluxGroundVeg = LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg
        LEfluxGround = LEfluxGroundImp * FractionsGround.fimp + LEfluxGroundBare * FractionsGround.fbare + LEfluxGroundVeg * FractionsGround.fveg
        LEfluxTree = LEfluxTreeInt + LTEfluxTree

        # Boolean operator for presence and absence of trees, vegetated ground, bare ground, and impervious ground
        Ctree = int(ParTree.trees == 1)
        Cveg = int(FractionsGround.fveg > 0)
        Cbare = int(FractionsGround.fbare > 0)
        Cimp = int(FractionsGround.fimp > 0)

        # Calculate total latent heat flux of canyon [W m^-2]
        LEflux_canyon = Cimp*FractionsGround.fimp*LEfluxGroundImp + Cbare*FractionsGround.fbare*(LEfluxGroundBarePond+LEfluxGroundBareSoil) + \
                        Cveg*FractionsGround.fveg*(LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg) + \
                        Ctree*4*geometry.radius_tree*(LEfluxTreeInt + LTEfluxTree)
        # Calculate total sensible heat flux of canyon [W m^-2]
        Hflux_canyon = Anthropogenic.Qf_canyon + Cimp*FractionsGround.fimp*HfluxGroundImp + Cbare*FractionsGround.fbare*HfluxGroundBare +\
                       Cveg*FractionsGround.fveg*HfluxGroundVeg + geometry.hcanyon*HfluxWallSun + geometry.hcanyon*HfluxWallShade +\
                       Ctree*4*geometry.radius_tree*HfluxTree

        # Total evaporative flux from bare soil [kg m^-2 s^-1]
        EfluxGroundBare = Eground_bare_pond + Eground_bare_soil
        # Total evaporative flux from vegetated ground [kg m^-2 s^-1]
        EfluxGroundVeg = Eground_veg_int + Eground_veg_pond + Eground_veg_soil + TEground_veg
        # Total evaporative flux from ground [kg m^-2 s^-1]
        EfluxGround = Eground_imp_pond*FractionsGround.fimp + EfluxGroundBare*FractionsGround.fbare + EfluxGroundVeg*FractionsGround.fveg

        if ParTree.trees == 0:
            SWRabs_t.SWRabsTree = 0
            LWRabs_t.LWRabsTree = 0

        # Calculate trees temperature using least square method
        Ttree = self.FSolver_Tree(TempVec,TemperatureC,MeteoData,Int,ExWater,Vwater,Owater,SoilPotW,self.CiCO2,ViewFactor,Geometry_m,ParTree,
                                  geometry,FractionsGround,ParSoilGround,ParInterceptionTree,PropOpticalGround,PropOpticalWall,
                                  PropOpticalTree,ParVegGround,ParVegTree,SunPosition,ParCalculation,VerticalProfUrban,ColParam)
        self.Ttree = numpy.array(Ttree)

        # Sink/source terms in temperature equation in the column model [K m s^-1]
        self.Src.thb = SurfaceHeatFlux.thb_ground
        self.Src.tvb = SurfaceHeatFlux.tvb_wall
        # Sink/source terms in humidity equation in the column model [kg kg^-1 m s^-1]
        self.Src.qhb = SurfaceHeatFlux.qhb_ground

        # Ground heat fluxes [W m^-2]
        self.Gflux.GfluxGroundImp = SWRabs_t.SWRabsGroundImp + LWRabs_t.LWRabsGroundImp - HfluxGroundImp - LEfluxGroundImp
        self.Gflux.GfluxGroundBare = SWRabs_t.SWRabsGroundBare + LWRabs_t.LWRabsGroundBare - HfluxGroundBare - (LEfluxGroundBarePond + LEfluxGroundBareSoil)
        self.Gflux.GfluxGroundVeg = SWRabs_t.SWRabsGroundVeg + LWRabs_t.LWRabsGroundVeg - HfluxGroundVeg - (LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg)
        self.Gflux.GfluxGround = Cimp*FractionsGround.fimp*self.Gflux.GfluxGroundImp + Cveg*FractionsGround.fveg*self.Gflux.GfluxGroundVeg + Cbare*FractionsGround.fbare*self.Gflux.GfluxGroundBare
        self.Gflux.GfluxWallSun = SWRabs_t.SWRabsWallSun + LWRabs_t.LWRabsWallSun - HfluxWallSun - 0.0
        self.Gflux.GfluxWallShade = SWRabs_t.SWRabsWallShade + LWRabs_t.LWRabsWallShade - HfluxWallShade - 0.0

        # Sensible heat fluxes [W m^-2]
        self.Hflux.HfluxGroundImp = HfluxGroundImp
        self.Hflux.HfluxGroundBare = HfluxGroundBare
        self.Hflux.HfluxGroundVeg = HfluxGroundVeg
        self.Hflux.HfluxWallSun = HfluxWallSun
        self.Hflux.HfluxWallShade = HfluxWallShade
        self.Hflux.HfluxTree =  HfluxTree
        self.Hflux.HfluxGround = HfluxGround
        self.Hflux.HfluxCanyon = Hflux_canyon

        # Latent heat fluxes [W m^-2] (latent heat flux for the walls are not considered)
        self.LEflux.LEfluxGroundImp = LEfluxGroundImp
        self.LEflux.LEfluxGroundBare = LEfluxGroundBarePond + LEfluxGroundBareSoil
        self.LEflux.LEfluxGroundBarePond = LEfluxGroundBarePond
        self.LEflux.LEfluxGroundBareSoil = LEfluxGroundBareSoil
        self.LEflux.LEfluxGroundVeg = LEfluxGroundVegInt + LEfluxGroundVegPond + LEfluxGroundVegSoil + LTEfluxGroundVeg
        self.LEflux.LEfluxGroundVegInt = LEfluxGroundVegInt
        self.LEflux.LEfluxGroundVegPond = LEfluxGroundVegPond
        self.LEflux.LEfluxGroundVegSoil = LEfluxGroundVegSoil
        self.LEflux.LTEfluxGroundVeg = LTEfluxGroundVeg
        self.LEflux.LEfluxGround = LEfluxGround
        self.LEflux.LEfluxWallSun = 0
        self.LEflux.LEfluxWallShade = 0
        self.LEflux.LEfluxTreeInt = LEfluxTreeInt
        self.LEflux.LTEfluxTree = LTEfluxTree
        self.LEflux.LEfluxTree = LEfluxTree
        self.LEflux.LEfluxCanyon = LEflux_canyon

        # Evaporative fluxes [kg m^-2 s^-1]
        self.Eflux.EfluxGroundImp = Eground_imp_pond
        self.Eflux.EfluxGroundBarePond = Eground_bare_pond
        self.Eflux.EfluxGroundBareSoil = Eground_bare_soil
        self.Eflux.EfluxGroundBare = EfluxGroundBare
        self.Eflux.EfluxGroundVegInt = Eground_veg_int
        self.Eflux.EfluxGroundVegPond = Eground_veg_pond
        self.Eflux.EfluxGroundVegSoil = Eground_veg_soil
        self.Eflux.TEfluxGroundVeg = TEground_veg
        self.Eflux.EfluxGroundVeg = EfluxGroundVeg
        self.Eflux.EfluxTreeInt = E_tree_int
        self.Eflux.TEfluxTree = TE_tree
        self.Eflux.EfluxGround = EfluxGround
        self.Eflux.EfluxTree = E_tree_int+TE_tree
        self.Eflux.EfluxWallSun = 0
        self.Eflux.EfluxWallShade = 0

        # Deep soil temperature [K]
        self.Tdepth.TDampGroundImp = Tdp_ground_imp
        self.Tdepth.TDampGroundBare = Tdp_ground_bare
        self.Tdepth.TDampGroundVeg = Tdp_ground_veg

        # Radiation [W m^-2]
        self.SWR.SWRin = SWRin_t
        self.SWR.SWRout = SWRout_t
        self.SWR.SWRabs = SWRabs_t
        self.SWR.SWRabsDir = SWRabsDir_t
        self.SWR.SWRabsDiff = SWRabsDiff_t
        self.SWR.SWREB = SWREB_t
        self.LWR.LWRin = LWRin_t
        self.LWR.LWRout = LWRout_t
        self.LWR.LWRabs = LWRabs_t
        self.LWR.LWREB = LWREB_t

        # Leaf Interior  CO2 concentration [umolCO2 mol^-1]
        self.CiCO2.CiCO2LeafTreeSun = Ci_sun_tree
        self.CiCO2.CiCO2LeafTreeShd = Ci_shd_tree
        self.CiCO2.CiCO2LeafGroundVegSun = Ci_sun_ground
        self.CiCO2.CiCO2LeafGroundVegShd = Ci_shd_ground

        # Resistances
        self.Res.rap_can = rap_can
        self.Res.rap_Htree_In = rap_Htree_In
        self.Res.rb_HGround = rb_H
        self.Res.rb_LGround = rb_L
        self.Res.r_soilGroundbare = r_soil_bare
        self.Res.r_soilGroundveg = r_soil_veg
        self.Res.rs_sunGround = rs_sun_L
        self.Res.rs_shdGround = rs_shd_L
        self.Res.rs_sunTree = rs_sun_H
        self.Res.rs_shdTree = rs_shd_H

        # Other parameters
        self.OtherParam.Ri_nearGround = SurfaceHeatFlux.HfGruond_otherparam.Ri_nearGround
        self.OtherParam.Utotal_nearGround = SurfaceHeatFlux.HfGruond_otherparam.Utotal_nearGround
        self.OtherParam.T_nearGround = SurfaceHeatFlux.HfGruond_otherparam.T_nearGround
        self.OtherParam.Tground = SurfaceHeatFlux.HfGruond_otherparam.Tground
        self.OtherParam.alp_soil_bare = alp_soil_bare
        self.OtherParam.alp_soil_veg = alp_soil_veg
        self.OtherParam.Fsun_L = Fsun_L
        self.OtherParam.Fshd_L = Fshd_L
        self.OtherParam.dw_L = dw_L

    def EBSolver_Tree(self,Ttree,TemperatureC,MeteoData,Int,ExWater,Vwater,Owater,SoilPotW,CiCO2Leaf,ViewFactor,Gemeotry_m,ParTree,
                      geometry,FractionsGround,ParSoilGround,ParInterceptionTree,PropOpticalGround,PropOpticalWall,PropOpticalTree,
                      ParVegGround,ParVegTree,SunPosition,ParCalculation,VerticalProfUrban,ColParam):

        """
        ------
        INPUT:
        Ttree: Trees temperature [K]
        TemperatureC: Temperature of the exterior surfaces [K]
        MeteoData: Forcing variables
        Int: Water interception [mm]
        ExWater: Extractable water [mm m^2 m^-2 ground s^-1]
        Vwater: [mm]
        Owater: Water Content [-]
        SoilPotW: [MPa]
        CiCO2Leaf: Leaf Interior  CO2 concentration [umolCO2 mol^-1]
        ViewFactor: View factors
        Gemeotry_m: Geometric parameters
        ParTree: Trees parameter
        geometry: Normalized geometric parameters
        FractionsGround: Fractions of ground covered by vegetation, impervious, and bare surface
        ParSoilGround: Soil parameters
        ParInterceptionTree: Specific water retained by the tree [mm m^2 VEG area m^-2 plant area]
        PropOpticalGround: Optical properties of the ground (albedo and emissivity)
        PropOpticalWall: Optical properties of the wall (albedo and emissivity)
        PropOpticalTree: Optical properties of the tree (albedo and emissivity)
        ParVegGround: Ground vegetation parameters
        ParVegTree: Trees parameter
        SunPosition: Sun angles
        ParCalculation: General calculation parameters
        VerticalProfUrban: Vertical profile of variables obtained from 1-D model
        ColParam: 1-D model parameters
        -------
        OUTPUT:
        YTree: Energy balance
        EBTree: sensible and latent
        """

        # Update trees temperature in the temperature vector
        TemperatureC_update = [TemperatureC[0],TemperatureC[1],TemperatureC[2],TemperatureC[3],TemperatureC[4],Ttree[0]]


        RadiationCal = RadiationFunctions()
        SurfaceHeatFlux = Surface_HeatFlux()

        # Shortwave radiation
        SWRin_t, SWRout_t, SWRabs_t, SWRabsDir_t, SWRabsDiff_t, SWREB_t = \
            RadiationCal.TotalSWRabsorbed(geometry, FractionsGround, ParTree, PropOpticalGround, PropOpticalWall,
                                          PropOpticalTree, ParVegTree, MeteoData, SunPosition, ViewFactor)

        # Tree absorbed: conversion from sphere to horizontal projected area
        Line2Circle_Multiplier = numpy.pi
        SWRabs_t.SWRabsTree = SWRabs_t.SWRabsTree * Line2Circle_Multiplier
        SWRabsDir_t.SWRabsTree = SWRabsDir_t.SWRabsTree * Line2Circle_Multiplier
        SWRabsDiff_t.SWRabsTree = SWRabsDiff_t.SWRabsTree * Line2Circle_Multiplier

        # Longwave radiation
        LWRin_t, LWRout_t, LWRabs_t, LWREB_t = \
            RadiationCal.TotalLWRabsorbed(TemperatureC_update, geometry, MeteoData, FractionsGround, PropOpticalGround,
                                          PropOpticalWall, PropOpticalTree, ParTree, ViewFactor)
        # Tree absorbed: conversion from sphere to horizontal projected area
        LWRabs_t.LWRabsTree = LWRabs_t.LWRabsTree * Line2Circle_Multiplier

        # Turbulent heat fluxes from ground and trees to canyon
        HfluxGroundImp, HfluxGroundBare, HfluxGroundVeg, HfluxTree, Eground_imp_pond, Eground_bare_pond, Eground_bare_soil,\
        Eground_veg_int, Eground_veg_pond, Eground_veg_soil, TEground_veg, E_tree_int, TE_tree, Ebare, Eveg, Etree, LEfluxGroundImp,\
        LEfluxGroundBarePond, LEfluxGroundBareSoil, LEfluxGroundVegInt, LEfluxGroundVegPond, LEfluxGroundVegSoil, LTEfluxGroundVeg,\
        LEfluxTreeInt, LTEfluxTree, LEbare, LEveg, LEtree, Ci_sun_tree, Ci_shd_tree, Ci_sun_ground, Ci_shd_ground, rap_can,\
        rap_Htree_In, rb_H, rb_L, r_soil_bare, r_soil_veg, alp_soil_bare, alp_soil_veg, rs_sun_L, rs_shd_L, rs_sun_H,\
        rs_shd_H, Fsun_L, Fshd_L, dw_L = \
            SurfaceHeatFlux.HeatFlux_ground_1D(TemperatureC_update, MeteoData, Gemeotry_m, geometry, FractionsGround,ParTree,
                                               ParVegGround, ParVegTree, ParSoilGround, SoilPotW, Owater, Vwater, ExWater,
                                               Int, CiCO2Leaf, ParInterceptionTree, ParCalculation,SWRabsDir_t.SWRabsTree,
                                               SWRabsDiff_t.SWRabsTree, SWRabsDir_t.SWRabsGroundVeg,SWRabsDiff_t.SWRabsGroundVeg,
                                               VerticalProfUrban,ColParam)

        if ParTree.trees == 0:
            SWRabs_t.SWRabsTree = 0
            LWRabs_t.LWRabsTree = 0

        # Energy balance for trees
        if ParTree.trees > 0:
            # Net energy flux to leaves [W m^-2]
            YTree = SWRabs_t.SWRabsTree + LWRabs_t.LWRabsTree - HfluxTree - LEfluxTreeInt - LTEfluxTree
        else:
            # When there are no trees we do not need to iteratively solve for tree temperature
            YTree = 0

        class EBTree_Def():
            pass
        self.EBTree = EBTree_Def()
        self.EBTree.HfluxTree = HfluxTree
        self.EBTree.LEfluxTreeInt = LEfluxTreeInt
        self.EBTree.LTEfluxTree = LTEfluxTree
        self.EBTree.LEfluxTree = LEfluxTreeInt + LTEfluxTree
        self.EBTree.EfluxTreeInt = E_tree_int
        self.EBTree.TEfluxTree = TE_tree
        self.EBTree.EfluxTree = E_tree_int + TE_tree
        self.EBTree.CiCO2LeafTreeSun = Ci_sun_tree
        self.EBTree.CiCO2LeafTreeShd = Ci_shd_tree
        self.EBTree.rap_Htree_In = rap_Htree_In
        self.EBTree.rb_HGround = rb_H
        self.EBTree.rs_sunTree = rs_sun_H
        self.EBTree.rs_shdTree = rs_shd_H

        return YTree

    def FSolver_Tree(self,TempVec,TemperatureC,MeteoData,Int,ExWater,Vwater,Owater,SoilPotW,CiCO2Leaf,ViewFactor,Gemeotry_m,
                     ParTree,geometry,FractionsGround,ParSoilGround,ParInterceptionTree,PropOpticalGround,PropOpticalWall,
                     PropOpticalTree,ParVegGround,ParVegTree,SunPosition,ParCalculation,VerticalProfUrban,ColParam):

        # Use temperature from previous time step as a starting point
        TTree = copy.copy(TemperatureC[5])

        # For least squares fitting we need a lower bound and an upper bound temperature
        lb = 243
        ub = 373

        Result = least_squares(self.EBSolver_Tree, TTree, bounds=(lb,ub),
                               args=(TemperatureC,MeteoData,Int,ExWater,Vwater,Owater,SoilPotW,CiCO2Leaf,ViewFactor,Gemeotry_m,
                                     ParTree,geometry,FractionsGround,ParSoilGround,ParInterceptionTree,PropOpticalGround,
                                     PropOpticalWall,PropOpticalTree,ParVegGround,ParVegTree,SunPosition,ParCalculation,
                                     VerticalProfUrban,ColParam))
        T = Result.x

        # If solver failed, retry with different starting value
        if Result.status < 1:
            TTree = copy.copy(MeteoData.Tatm[itt])

            Result = least_squares(self.EBSolver_Tree, TTree, bounds=(lb, ub),
                                   args=(TemperatureC, MeteoData, Int, ExWater, Vwater, Owater, SoilPotW, CiCO2Leaf,
                                         ViewFactor, Gemeotry_m,ParTree, geometry, FractionsGround, ParSoilGround, ParInterceptionTree,
                                         PropOpticalGround,PropOpticalWall, PropOpticalTree, ParVegGround, ParVegTree, SunPosition,
                                         ParCalculation,VerticalProfUrban,ColParam))
            T = Result.x

        for i in range(3):
            if Result.status < 1:
                TTree = TempVec.TTree + i

                Result = least_squares(self.EBSolver_Tree, TTree, bounds=(lb, ub),
                                       args=(TemperatureC, MeteoData, Int, ExWater, Vwater, Owater, SoilPotW, CiCO2Leaf,
                                             ViewFactor, Gemeotry_m, ParTree, geometry, FractionsGround, ParSoilGround,
                                             ParInterceptionTree,PropOpticalGround, PropOpticalWall, PropOpticalTree, ParVegGround,
                                             ParVegTree, SunPosition,ParCalculation, VerticalProfUrban,ColParam))
                T = Result.x

        for i in range(3):
            if Result.status < 1:
                TTree = TempVec.TTree - i

                Result = least_squares(self.EBSolver_Tree, TTree, bounds=(lb, ub),
                                       args=(TemperatureC, MeteoData, Int, ExWater, Vwater, Owater, SoilPotW, CiCO2Leaf,
                                             ViewFactor, Gemeotry_m, ParTree, geometry, FractionsGround, ParSoilGround,
                                             ParInterceptionTree,PropOpticalGround, PropOpticalWall, PropOpticalTree, ParVegGround,
                                             ParVegTree, SunPosition,ParCalculation, VerticalProfUrban,ColParam))
                T = Result.x

        return T
