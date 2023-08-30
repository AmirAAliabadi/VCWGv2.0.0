import os
import numpy
import math
from scipy.interpolate import interp1d
from Soil_Functions import Soil_Calculations
from Resistance_Functions import Ressitance_Calculations
from SurfaceHeatFluxDef import ForceRestoreConductiveHeatImp,ForceRestoreConductiveHeatSoil
import copy

'''
Surface Energy Balance Model (SEBM): Calculate turbulent heat fluxes at the surface of urban elements
Developed by Mohsen Moradi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: June 2021
'''

class Surface_HeatFlux(object):

    def ForceRestore_Conductive_Heat_Imp(self,TemperatureC,TempDamp,TempVec,ParCalculation,ParThermalGround,FractionsGround):

        """
        ------
        INPUT:
        TemperatureC: Temperature of the exterior surfaces [K]
        TempDamp: Deep soil temperature from previous step [K]
        TempVec: Temperature of the exterior surfaces from two steps back [K]
        ParCalculation: General calculation parameters
        ParThermalGround: Thermal properties of the ground (Thermal conductivity[W m^-1 K^-1], Volumetric heat capacity [J m^-3 K^-1])
        FractionsGround: Fractions of ground covered by vegetated, bare,and impervious surface
        -------
        OUTPUT:
        G: Ground heat flux [W m^-2] (Note: the force-restore method is only used to calculate deep soil temperature.)
        Tdp: Deep soil temperature [K]
        """

        Ts = TemperatureC[0]

        if FractionsGround.fimp > 0:
            Cimp = 1
        else:
            Cimp = 0

        # Time constant [s]
        tau = 86400
        # Total Thermal Capacity Soil [K m^2 J^-1]
        CTt = 2 * (numpy.sqrt(numpy.pi / (ParThermalGround.lan_dry_imp * ParThermalGround.cv_s_imp * tau)))

        SoilHeat = self.Soil_Heat(ParCalculation.dts,Ts-273.15,TempVec.TGroundImp-273.15,TempDamp.TDampGroundImp-273.15,CTt)
        Tdp = SoilHeat[1]
        G = SoilHeat[0]*Cimp

        self.CondHeatImp = ForceRestoreConductiveHeatImp()
        self.CondHeatImp.G = G*Cimp
        self.CondHeatImp.Tdp = Tdp + 273.15

        return G,Tdp + 273.15

    def Impervious_Conductive_Heat(self,TemperatureC,TempVec,Anthropogenic,ParThermalWall,WallLayers,ParCalculation,type):
        # Sun
        if type == 1:
            Ts = TemperatureC[3]
            Tb = Anthropogenic.Tb
            Tint = TemperatureC[6]
            Tint_tm1 = TempVec.TWallIntSun
            lan_dry1 = ParThermalWall.lan_dry
            lan_dry2 = ParThermalWall.lan_dry
            dz1 = WallLayers.dz1_wall
            dz2 = WallLayers.dz2_wall
            cv_s1 = ParThermalWall.cv_s
            cv_s2 = ParThermalWall.cv_s
            dts = ParCalculation.dts
        # Shade
        elif type == 0:
            Ts = TemperatureC[4]
            Tb = Anthropogenic.Tb
            Tint = TemperatureC[7]
            Tint_tm1 = TempVec.TWallIntShade
            lan_dry1 = ParThermalWall.lan_dry
            lan_dry2 = ParThermalWall.lan_dry
            dz1 = WallLayers.dz1_wall
            dz2 = WallLayers.dz2_wall
            cv_s1 = ParThermalWall.cv_s
            cv_s2 = ParThermalWall.cv_s
            dts = ParCalculation.dts
        else:
            print('please, enter sun or shade for sunlit or shaded wall')

        # Soil Heat Flux [W m^-2]
        G1 = lan_dry1 * (Ts - Tint) / dz1
        G2 = lan_dry2 * (Tint - Tb) / dz2
        dS = (cv_s1 + cv_s2) / 2 * (dz1 + dz2) / dts * (Tint - Tint_tm1)

        return G1,G2,dS

    def Impervious_Conductive_HeatRoof(self,TemperatureR,TempVec,Anthropogenic,ParThermalRoof,ParSoilRoof,ParCalculation):

        Ts = TemperatureR[0]
        Tint = TemperatureR[2]
        Tb = Anthropogenic.Tb
        Tint_tm1 = TempVec.TRoofIntImp
        lan_dry1 = ParThermalRoof.lan_dry_imp
        lan_dry2 = ParThermalRoof.lan_dry_imp
        dz1 = ParSoilRoof.dz1
        dz2 = ParSoilRoof.dz2
        cv_s1 = ParThermalRoof.cv_s_imp
        cv_s2 = ParThermalRoof.cv_s_imp
        dts = ParCalculation.dts

        G1 = lan_dry1 * (Ts - Tint) / dz1
        G2 = lan_dry2 * (Tint - Tb) / dz2
        dS = (cv_s1 + cv_s2) / 2 * (dz1 + dz2) / dts * (Tint - Tint_tm1)

        return G1,G2,dS

    def Soil_Heat(self,dt,Ts,Tstm1,Tdptm1,CTt):

        """
        ------
        INPUT:
        dt: Time step [s]
        Ts: Surface temperature [C]
        Tstm1: Surface temperature from two steps back [C]
        Tdptm1: Deep soil temperature from previous step [C]
        CTt: Total thermal capacity [K m^2 J^-1]
        -------
        OUTPUT:
        G: Ground heat flux [W m^-2] (Note: the force-restore method is only used to calculate deep soil temperature.)
        Tdp: Deep soil temperature [C]
        """

        # Time constant [s]
        tau = 86400

        # Temperature Variation [C]
        dTs = Ts - Tstm1

        # Depth Temperature [C]
        Tdp = (1 / (1 + dt / tau)) * (Tdptm1 + (dt / tau) * Ts)
        # Soil Heat Flux [W m^-2]
        G = (1 / CTt) * (2 * numpy.pi * (Ts - Tdp) / tau + dTs / dt)

        return G,Tdp

    def Soil_Conductive_Heat(self,TemperatureR,TempVec,Anthropogenic,Owater,ParVegRoof,ParSoilRoof,ParThermalRoof,ParCalculation):

        Troof = TemperatureR[1]
        Tint = TemperatureR[3]
        Tint_tm1 = TempVec.TRoofIntImp
        Tb = Anthropogenic.Tb
        Otm1 = Owater.OwRoofSoilVeg
        Rrootl = ParVegRoof.Rrootl
        PsiL50 = ParVegRoof.PsiL50
        PsiX50 = ParVegRoof.PsiX50
        CASE_ROOT = ParVegRoof.CASE_ROOT
        ZR95 = ParVegRoof.ZR95
        ZR50 = ParVegRoof.ZR50
        ZRmax = ParVegRoof.ZRmax
        Pcla = ParSoilRoof.Pcla
        Psan = ParSoilRoof.Psan
        Porg = ParSoilRoof.Porg
        Kfc = ParSoilRoof.Kfc
        Phy = ParSoilRoof.Phy
        SPAR = ParSoilRoof.SPAR
        Kbot = ParSoilRoof.Kbot
        Zs = ParSoilRoof.Zs
        dz1 = ParSoilRoof.dz1
        dz2 = ParSoilRoof.dz2
        cv_s2 = ParThermalRoof.cv_s_imp
        lan_dry2 = ParThermalRoof.lan_dry_imp
        dts = ParCalculation.dts

        SoilCal = Soil_Calculations()

        # Calculation of soil parameters
        _Zs_, dz, ms, Osat, Ohy, _nVG_, _alpVG_, _Ks_Zs_, _L_, _Pe_, _O33_, _SPAR_, _EvL_Zs_, \
        _Inf_Zs_, _RfH_Zs_, _RfL_Zs_, _Zinf_, _Kbot_, _Slo_pot_, _Dz_, _aR_, _aTop_, rsd, lan_dry, lan_s, cv_s = \
            SoilCal.Soil_Parameters_Total(Pcla, Psan, Porg, Kfc, Phy, SPAR, Kbot, CASE_ROOT, CASE_ROOT, 0, ZR95, 0, ZR50,
                                          0, ZRmax, Zs)

        Tdamptm1 = Tb
        # Soil temperature of the layer
        Tdp = [Tdamptm1 for i in range(ms)]
        # Soil temperature of the layer at the previous time step
        Tdptm1 = Tdp

        lanS, cv_soil, _CTt_ = SoilCal.Soil_Thermal_Properties([Tdptm1[i] - 273.15 for i in range(len(Tdptm1))],
                                                                   rsd, lan_dry, lan_s, cv_s, Osat, Ohy, Otm1)

        # Average soil parameters according to soil layer thickness
        dz = dz / sum(dz, 2)
        lanS = numpy.dot(dz,lanS)
        cv_soil = numpy.dot(dz,cv_soil)

        # Average roof parameters according to roof thickness
        dz_roof = [dz1, dz2]
        cv_roof = [cv_soil, cv_s2]
        dz_roof = [dz_roof[i] / sum(dz_roof) for i in range(len(dz_roof))]
        cv_roof = numpy.dot(dz_roof,cv_roof)

        # Computation of heat fluxes
        # Soil Heat Flux [W m^-2]
        G1 = lanS * (Troof - Tint) / dz1
        G2 = lan_dry2 * (Tint - Tb) / dz2
        dS = cv_roof * (dz1 + dz2) / dts * (Tint - Tint_tm1)

        return G1,G2,dS


    def ForceRestore_Conductive_Heat_Soil(self,TemperatureC,TempDamp,Owater,TempVec,ParCalculation,ParSoilGround,ParVegGround,
                                          ParVegTree,FractionsGround,type):
        """
        ------
        INPUT:
        TemperatureC: Temperature of the exterior surfaces [K]
        TempDamp: Deep Soil Temperature from previous step [K]
        Owater: Water Content [-]
        TempVec: Temperature of the exterior surfaces from two steps back [K]
        ParCalculation: General calculation parameters
        ParSoilGround: Soil parameters
        ParVegGround: Ground vegetation parameters
        ParVegTree: Trees parameter
        FractionsGround: Fractions of ground covered by vegetated, bare,and impervious surface
        type: type of the surface (vegetated or bare)
        -------
        OUTPUT:
        G: Ground heat flux [W m^-2] (Note: the force-restore method is only used to calculate deep soil temperature.)
        Tdp: Deep soil temperature [K]
        """

        if type == 0: # bare soil
            Ts = TemperatureC[1]
            Tdptm1 = TempDamp.TDampGroundBare
            Otm1 = Owater.OwGroundSoilBare[0]
            Tstm1 = TempVec.TGroundBare
            Pcla = ParSoilGround.Pcla
            Psan = ParSoilGround.Psan
            Porg = ParSoilGround.Porg
            Kfc = ParSoilGround.Kfc
            Phy = ParSoilGround.Phy
            SPAR = ParSoilGround.SPAR
            Kbot = ParSoilGround.Kbot
            CASE_ROOT_H = ParVegTree.CASE_ROOT
            CASE_ROOT_L = ParVegGround.CASE_ROOT
            ZR95_H = ParVegTree.ZR95
            ZR95_L = ParVegGround.ZR95
            ZR50_H = ParVegTree.ZR50
            ZR50_L = ParVegGround.ZR50
            ZRmax_H = ParVegTree.ZRmax
            ZRmax_L = ParVegGround.ZRmax
            Zs = ParSoilGround.Zs
            if FractionsGround.fbare > 0:
                Csoil = 1
            else:
                Csoil = 0

        elif type == 1: # vegetated ground
            Ts = TemperatureC[2]
            Tdptm1 = TempDamp.TDampGroundVeg
            Otm1 = Owater.OwGroundSoilVeg[0]
            Tstm1 = TempVec.TGroundVeg
            Pcla = ParSoilGround.Pcla
            Psan = ParSoilGround.Psan
            Porg = ParSoilGround.Porg
            Kfc = ParSoilGround.Kfc
            Phy = ParSoilGround.Phy
            SPAR = ParSoilGround.SPAR
            Kbot = ParSoilGround.Kbot
            CASE_ROOT_H = ParVegTree.CASE_ROOT
            CASE_ROOT_L = ParVegGround.CASE_ROOT
            ZR95_H = ParVegTree.ZR95
            ZR95_L = ParVegGround.ZR95
            ZR50_H = ParVegTree.ZR50
            ZR50_L = ParVegGround.ZR50
            ZRmax_H = ParVegTree.ZRmax
            ZRmax_L = ParVegGround.ZRmax
            Zs = ParSoilGround.Zs
            if FractionsGround.fveg > 0:
                Csoil = 1
            else:
                Csoil = 0
        else:
            print('please enter a valid specification of ground type. 0 = bare, 1 = vegetated')


        SoilCal = Soil_Calculations()
        SoilCal.Soil_Parameters_Total(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,
                                      ZRmax_H,ZRmax_L,Zs)

        SoilCal.Soil_Thermal_Properties(Tdptm1-273.15,SoilCal.SoilParamTotal.rsd,SoilCal.SoilParamTotal.lan_dry,SoilCal.SoilParamTotal.lan_s,
                                        SoilCal.SoilParamTotal.cv_s,SoilCal.SoilParamTotal.Osat,SoilCal.SoilParamTotal.Ohy,Otm1)

        SoilHeat = self.Soil_Heat(ParCalculation.dts,Ts-273.15,Tstm1-273.15,Tdptm1-273.15,SoilCal.SoilThProp.CTt)
        Tdp = SoilHeat[1]+273.15
        G = SoilHeat[0]*Csoil

        self.CondHeatSoil = ForceRestoreConductiveHeatSoil()
        self.CondHeatSoil.G = G
        self.CondHeatSoil.Tdp = Tdp

        return G, Tdp

    def HeatFlux_wall(self,TemperatureC,Geometry_m,VerticalProfUrban,ParCalculation,ColParam):

        """
        ------
        INPUT:
        TemperatureC: Vector of canyon surface temperatures [K]
        Geometry_m: Geometric parameters
        VerticalProfUrban: Vertical profile of variables obtained from 1-D model
        ParCalculation: General calculation parameters
        ColParam: Column model parameters
        -------
        OUTPUT:
        Hwsun: Sensible heat flux from sunlit wall [W m^-2]
        Hwshade: Sensible heat flux from shaded wall [W m^-2]
        Ewsun: Evaporative flux from sunlit wall [kg m^-2 s^-1]. Note: It is zero in the current version of VCWG.
        Ewshade: Evaporative flux from shaded wall [kg m^-2 s^-1]. Note: It is zero in the current version of VCWG.
        LEwsun: Latent heat flux from sunlit wall [W m^-2]. Note: It is zero in the current version of VCWG.
        LEwshade: Latent heat flux from shaded wall [W m^-2]. Note: It is zero in the current version of VCWG.
        """

        # Parameter definitions
        Twsun = TemperatureC[3]
        Twshade = TemperatureC[4]
        T_canyon = copy.copy(VerticalProfUrban.th)
        Tavg_canyon = numpy.mean(VerticalProfUrban.th[0:Geometry_m.nz_u])

        # Air specific heat  [J kg^-1 K^-1]
        cp_atm = 1005 + (((Tavg_canyon - 273.15) + 23.15) ** 2) / 3364
        # dry air density [kg m^-3]
        rho_atm = numpy.mean(VerticalProfUrban.rho[:Geometry_m.nz_u])
        # Latent heat vaporization/condensation [J kg^-1]
        L_heat = 1000 * (2501.3 - 2.361 * (Tavg_canyon - 273.15))

        # Call the resistance function
        ResistanceCal = Ressitance_Calculations()

        Hwsun_z = []
        Hwshade_z = []
        for i_z in range(Geometry_m.nz_u):
            # Calculate wall resistance [s m^-1]
            RES_w = ResistanceCal.Wall_Aerodynamic_Resistance(VerticalProfUrban,Geometry_m,ColParam.WindMin_Urban,cp_atm,
                                                              i_z,ParCalculation)
            # Calculate sensible heat flux from sunlit wall [W m^-2]
            Hwsun_z.append(cp_atm * VerticalProfUrban.rho[i_z] * (Twsun-T_canyon[i_z]) / (RES_w))
            # Calculate sensible heat flux from shaded wall [W m^-2]
            Hwshade_z.append(cp_atm * VerticalProfUrban.rho[i_z] * (Twshade-T_canyon[i_z]) / (RES_w))

        # Calculate total sensible heat flux as the area weighted average of sensible heat fluxes from wall layers [W m^-2]
        Hwsun = (Geometry_m.dz/Geometry_m.Height_canyon)*sum(Hwsun_z)
        Hwshade = (Geometry_m.dz/Geometry_m.Height_canyon)*sum(Hwshade_z)

        # Calculate evaporative fluxes [kg m^-2 s^-1]. Note: It is zero in the current version of VCWG.
        Ewsun = 0
        Ewshade = 0
        # Calculate latent heat fluxes [W m^-2]. Note: It is zero in the current version of VCWG.
        LEwsun = L_heat * Ewsun
        LEwshade = L_heat * Ewshade

        class tvb_wall_Def():
            pass
        self.tvb_wall = tvb_wall_Def()
        # Sink/source terms in 1-D model [K m s^-1]
        self.tvb_wall.Hwsun_z = [Hwsun_z[i]/(cp_atm*rho_atm) for i in range(len(Hwsun_z))]
        self.tvb_wall.Hwshade_z = [Hwshade_z[i]/(cp_atm*rho_atm) for i in range(len(Hwshade_z))]
        self.tvb_wall.Hwsun = Hwsun/(cp_atm*rho_atm)
        self.tvb_wall.Hwshade = Hwshade/(cp_atm*rho_atm)
        # Total sensible heat fluxes [W m^-2]
        self.tvb_wall.Hfluxwsun = Hwsun
        self.tvb_wall.Hfluxwshade = Hwshade

        return Hwsun,Hwshade,Ewsun,Ewshade,LEwsun,LEwshade

    def HeatFlux_ground_1D(self,TemperatureC,MeteoData,Gemeotry_m,geometry,FractionsGround,ParTree,ParVegGround,ParVegTree,
                           ParSoilGround,SoilPotW,Owater,Vwater,ExWater,Int,CiCO2Leaf,ParInterceptionTree,ParCalculation,
                           SWRdir_abs_tree,SWRdiff_abs_tree,SWRdir_abs_groundveg,SWRdiff_abs_groundveg,VerticalProfUrban,ColParam):
        """
        ------
        INPUT:
        TemperatureC: Temperature of the exterior surfaces [K]
        MeteoData: Forcing variables
        Gemeotry_m: Geometric parameters
        geometry: Normalized geometric parameters
        FractionsGround: Fractions of ground covered by vegetated, bare, and impervious surface
        ParTree: Trees parameter
        ParVegGround: Ground vegetation parameters
        ParVegTree: Trees parameter
        ParSoilGround: Soil parameters
        SoilPotW: Soil water potential [MPa]
        Owater: Water Content [-]
        Vwater: Volume of water [mm]
        ExWater: extractable water [mm m^2 m^-2 ground s^-1]
        Int: Water interception [mm]
        CiCO2Leaf: Leaf Interior  CO2 mixing ratio [umolCO2 mol^-1]
        ParInterceptionTree: specific water retained by the tree [mm m^2 VEG area m^-2 plant area]
        ParCalculation: General calculation parameters
        SWRdir_abs_tree: Direct shortwave radiation absorbed by trees [W m^-2]
        SWRdiff_abs_tree: Diffusive shortwave radiation absorbed by trees [W m^-2]
        SWRdir_abs_groundveg: Direct shortwave radiation absorbed by ground vegetation [W m^-2]
        SWRdiff_abs_groundveg: Diffusive shortwave radiation absorbed by ground vegetation [W m^-2]
        VerticalProfUrban: Vertical profile of variables obtained from 1-D model
        ColParam: 1-D model parameters
        -------
        OUTPUT:
        Himp: Sensible heat flux from impervious surface [W m^-2]
        Hbare: Sensible heat flux from bare surface [W m^-2]
        Hveg: Sensible heat flux from vegetated surface [W m^-2]
        Htree: Sensible heat flux from tree [W m^-2]
        Eimp: Evaporation from intercepted water on impervious surface [kg m^-2 s^-1]
        Ebare_pond: Evaporation from intercepted water on bare ground [kg m^-2 s^-1]
        Ebare_soil: Evaporation from soil under bare ground [kg m^-2 s^-1]
        Eveg_int: Evaporation from intercepted water on low vegetation [kg m^-2 s^-1]
        Eveg_pond: Evaporation from intercepted water on ground under vegetation [kg m^-2 s^-1]
        Eveg_soil: Evaporation from soil under vegetated ground [kg m^-2 s^-1]
        TEveg: Transpiration from low vegetation [kg m^-2 s^-1]
        Etree_int: Evaporation from intercepted water on trees [kg m^-2 s^-1]
        TEtree: Transpiration from trees [kg m^-2 s^-1]
        Ebare: Evaporation from bare surface [kg m^-2 s^-1]
        Eveg: Evaporation from vegetated surface [kg m^-2 s^-1]
        Etree: Evaporation from tree [kg m^-2 s^-1]
        LEimp: Latent heat flux from impervious surface [W m^-2]
        LEbare_pond: Latent heat flux from runon on bare surface [W m^-2]
        LEbare_soil: Latent heat flux from soil of bare surface [W m^-2]
        LEveg_int: Latent heat flux from intercepted water on low vegetation [W m^-2]
        LEveg_pond: Latent heat flux from vegetated surface [W m^-2]
        LEveg_soil: Latent heat flux from soil underneath the low vegetation [W m^-2]
        LTEveg: Latent heat of transpiration from vegetated surface [W m^-2]
        LEtree_int: Latent heat flux from intercepted water on tree [W m^-2]
        LTEtree: Latent heat of transpiration from tree [W m^-2]
        LEbare: Latent heat flux from bare surface [W m^-2]
        LEveg: Latent heat flux from vegetated surface [W m^-2]
        LEtree: Latent heat flux from tree [W m^-2]
        Ci_sun_H: Leaf Interior CO2 mixing ratio of sunlit tree [umolCO2 mol^-1]
        Ci_shd_H: Leaf Interior CO2 mixing ratio shaded tree [umolCO2 mol^-1]
        Ci_sun_L: Leaf Interior CO2 mixing ratio of sunlit low vegetation [umolCO2 mol^-1]
        Ci_shd_L: Leaf Interior CO2 mixing ratio of shaded low vegetation [umolCO2 mol^-1]
        rap_can: Aerodynamic resistance near the ground [s m^-1]
        rap_Htree_In: Aerodynamic resistance between trees and canyon air [s m^-1]
        rb_H: Leaf boundary layer resistance for trees [s m^-1]
        rb_L: Leaf boundary layer resistance for low vegetation [s m^-1]
        r_soil_bare: Bare soil resistance [s m^-1]
        r_soil_veg: Vegetated soil resistance [s m^-1]
        alp_soil_bare: Relative humidity of bare soil [0-1]
        alp_soil_veg: Relative humidity of vegetated soil [0-1]
        rs_sun_L: Stomatal resistance of sunlit low vegetation [s m^-1]
        rs_shd_L: Stomatal resistance of shaded low vegetation [s m^-1]
        rs_sun_H: Stomatal resistance of sunlit trees [s m^-1]
        rs_shd_H: Stomatal resistance of shaded trees [s m^-1]
        Fsun_L: Fraction of sunlit for low vegetation [-]
        Fshd_L: Fraction of shaded for low vegetation [-]
        dw_L: Fraction of ground vegetation covered by intercepted water [-]
        """

        # Parameter specification
        Timp = TemperatureC[0]
        Tbare = TemperatureC[1]
        Tveg = TemperatureC[2]
        Ttree = TemperatureC[5]
        Tcanyon = copy.copy(VerticalProfUrban.th[0])
        qcanyon = copy.copy(VerticalProfUrban.qn[0])
        WindSpeed_top = numpy.sqrt(VerticalProfUrban.vx[-1]**2 + VerticalProfUrban.vy[-1]**2)
        WindSpeed_bottom = numpy.sqrt(VerticalProfUrban.vx[0]**2 + VerticalProfUrban.vy[0]**2)
        q_intp = interp1d(Gemeotry_m.z[:-1], VerticalProfUrban.qn)
        vx_intp = interp1d(Gemeotry_m.z[:-1], VerticalProfUrban.vx)
        vy_intp = interp1d(Gemeotry_m.z[:-1], VerticalProfUrban.vy)
        th_intp = interp1d(Gemeotry_m.z[:-1], VerticalProfUrban.th)
        Pre_intp = interp1d(Gemeotry_m.z[:-1], VerticalProfUrban.presProf)

        rad_tree = copy.copy(geometry.radius_tree)
        SPARTREE = ParVegTree.SPARTREE
        Pcla = ParSoilGround.Pcla
        Psan = ParSoilGround.Psan
        Porg = ParSoilGround.Porg
        Kfc = ParSoilGround.Kfc
        Phy = ParSoilGround.Phy
        SPAR = ParSoilGround.SPAR
        Kbot = ParSoilGround.Kbot
        CASE_ROOT_H = ParVegTree.CASE_ROOT
        CASE_ROOT_L = ParVegGround.CASE_ROOT
        ZR95_H = ParVegTree.ZR95
        ZR95_L = ParVegGround.ZR95
        ZR50_H = ParVegTree.ZR50
        ZR50_L = ParVegGround.ZR50
        ZRmax_H = ParVegTree.ZRmax
        ZRmax_L = ParVegGround.ZRmax
        Zs = ParSoilGround.Zs
        Psi_L_tm1 = SoilPotW.SoilPotWGroundVeg_L
        Psi_H_tm1 = SoilPotW.SoilPotWGroundTot_H


        ResistanceCal = Ressitance_Calculations()
        SoilCal = Soil_Calculations()

        # Boolean operator for presence and absence of trees, vegetated ground, bare ground, and impervious ground
        Ctree = int(ParTree.trees == 1)
        Cveg = int(FractionsGround.fveg > 0)
        Cbare = int(FractionsGround.fbare > 0)
        Cimp = int(FractionsGround.fimp > 0)

        if numpy.isnan(rad_tree):
            rad_tree = 0

        # Specific heat air  [J kg^-1 K^-1]
        cp_atm = 1005 + (((Tcanyon - 273.15) + 23.15) ** 2) / 3364
        # dry air density at atmosphere [kg m^-3]
        rho_atm = copy.copy(VerticalProfUrban.rho[0])
        # Latent heat vaporization/condensation [J kg^-1]
        L_heat = 1000 * (2501.3 - 2.361 * (Tcanyon - 273.15))

        # Average temperature of the ground [K]
        Tground = FractionsGround.fveg*Tveg + FractionsGround.fbare*Tbare + FractionsGround.fimp*Timp

        # vapor pressure saturation at Tground_imp[Pa]
        esat_T_imp = 611 * numpy.exp(17.27 * (Timp - 273.16) / (237.3 + (Timp - 273.16)))
        # Saturated specific humidity at Tground_imp [kg kg^-1]
        qsat_T_imp = round((0.622 * esat_T_imp) / (MeteoData.Pre - 0.378 * esat_T_imp),4)
        # vapor pressure saturation at Tground_bare [Pa]
        esat_T_bare = 611 * numpy.exp(17.27 * (Tbare - 273.16) / (237.3 + (Tbare - 273.16)))
        # Saturated specific humidity at Tground_bare [kg kg^-1]
        qsat_T_bare = round((0.622 * esat_T_bare) / (MeteoData.Pre - 0.378 * esat_T_bare),4)
        # vapor pressure saturation at Tground_veg [Pa]
        esat_T_veg = 611 * numpy.exp(17.27 * (Tveg - 273.16) / (237.3 + (Tveg - 273.16)))
        # Saturated specific humidity at Tground_veg [kg kg^-1]
        qsat_T_veg = round((0.622 * esat_T_veg) / (MeteoData.Pre - 0.378 * esat_T_veg),4)
        # vapor pressure saturation at Ttree [Pa]
        esat_T_tree = 611 * numpy.exp(17.27 * (Ttree - 273.16) / (237.3 + (Ttree - 273.16)))
        # Saturated specific humidity at Ttree [kg kg^-1]
        qsat_T_tree = round((0.622 * esat_T_tree) / (MeteoData.Pre - 0.378 * esat_T_tree),4)
        # vapor pressure saturation at T_canyon [Pa]
        esat_T_canyon = 611 * numpy.exp(17.27 * (Tcanyon - 273.16) / (237.3 + (Tcanyon - 273.16)))
        # Vapor pressure at T_canyon [Pa]
        e_T_canyon = qcanyon * MeteoData.Pre / (0.622 + 0.378 * qcanyon)
        # Vapor pressure at T_tree [Pa]
        q_tree = q_intp(Gemeotry_m.Height_tree)
        e_T_tree = q_tree * MeteoData.Pre / (0.622 + 0.378 * q_tree)

        # Parameters for stomata resistance
        # Leaf Interior CO2 mixing ratio [umolCO2 mol^-1]
        Citm1_sun_H = Ctree * CiCO2Leaf.CiCO2LeafTreeSun
        Citm1_shd_H = Ctree * CiCO2Leaf.CiCO2LeafTreeShd
        Citm1_sun_L = Cveg * CiCO2Leaf.CiCO2LeafGroundVegSun
        Citm1_shd_L = Cveg * CiCO2Leaf.CiCO2LeafGroundVegShd
        # Vapor Pressure Deficit [Pa]
        Ds_canyon = esat_T_canyon - e_T_canyon
        Ds_tree = esat_T_tree - e_T_tree


        # Partitioning of radiation into sunlit and shaded area
        # Fraction of sunlit/shaded for high vegetation
        Fsun_H = Ctree * ((1.0 - numpy.exp(-ParVegTree.Kopt * (ParVegTree.LAI))) / (ParVegTree.Kopt * (ParVegTree.LAI)))
        if Fsun_H < 0.01:
            Fsun_H = 0
        if Fsun_H > 1:
            Fsun_H = 1
        Fshd_H = Ctree * (1 - Fsun_H)
        # Absorbed direct and diffuse shortwave radiation of the sunlit surface [W m^-2]
        PAR_sun_H = Ctree * (SWRdir_abs_tree + Fsun_H*SWRdiff_abs_tree)
        # Absorbed direct and diffuse shortwave radiation of the shaded surface [W m^-2]
        PAR_shd_H = Ctree * (Fshd_H*SWRdiff_abs_tree)

        # Fraction of sunlit/shaded for low vegetation
        Fsun_L = Cveg * ((1.0 - numpy.exp(-ParVegGround.Kopt * (ParVegGround.LAI))) / (ParVegGround.Kopt * (ParVegGround.LAI)))
        if Fsun_L < 0.01:
            Fsun_L = 0
        if Fsun_L > 1:
            Fsun_L = 1
        Fshd_L = Cveg * (1 - Fsun_L)
        # Absorbed direct and diffuse shortwave radiation of the sunlit surface [W m^-2]
        PAR_sun_L = Cveg * (SWRdir_abs_groundveg + Fsun_L*SWRdiff_abs_groundveg)
        # Absorbed direct and diffuse shortwave radiation of the shaded surface  [W m^-2]
        PAR_shd_L = Cveg * (Fshd_L*SWRdiff_abs_groundveg)

        # Calculate fraction of vegetation covered by intercepted water (Deardorff 1978)
        # Maximum possible interception on high vegetation [mm]
        In_max_H = Ctree * (ParInterceptionTree.Sp_In * (ParVegTree.LAI + ParVegTree.SAI))
        # Fraction of tree covered by intercepted water [-]
        dw_H = Ctree * (min(1, (Int.IntTree / In_max_H)**(2/3)))
        # Maximum possible interception on low vegetation [mm]
        In_max_L = Cveg * (ParSoilGround.Sp_In * (ParVegGround.LAI + ParVegGround.SAI))
        # Fraction of ground vegetation covered by intercepted water [-]
        dw_L = Cveg * (min(1, (Int.IntGroundVegPlant / In_max_L)**(2/3)))

        # Calculate ground aerodynamic roughness length
        _zom_, _zoh_, zom_ground, _zoh_ground_, _disp_h_, _zom_H_, _zom_L_, _zoh_H_, _zoh_L_, _d_H_, _d_L_, _zom_other_ = \
            ResistanceCal.Urban_roughness(Ctree*Gemeotry_m.Height_tree, Cveg*ParVegGround.hc, Cbare, Cimp, 0)

        # Calculate wind speed at the level of high vegetation [m s^-1]
        if Ctree > 0:
            u_tree = numpy.sqrt(vx_intp(Gemeotry_m.Height_tree)**2 + vy_intp(Gemeotry_m.Height_tree)**2)
        else:
            u_tree = 0.1

        # Calculate wind speed at the level of low vegetation [m s^-1]
        if Cveg > 0:
            u_Lveg = WindSpeed_bottom
        else:
            u_Lveg = 0.1

        # Calculate displacement height and roughness length of the canyon
        dcan, zomcan, _u_Hcan_, _u_Zref_und_, _w_Zref_und_, _alpha_ = \
            ResistanceCal.WindProfile_Canyon(Gemeotry_m.Height_canyon,Gemeotry_m.Height_tree,Gemeotry_m.Radius_tree,Gemeotry_m.Width_canyon,
                                             Gemeotry_m.Width_roof,ParVegTree.Kopt,ParVegTree.LAI,MeteoData.Zatm,WindSpeed_top,
                                             1.5,ParTree.trees,1.5, zom_ground)


        # Calculate aerodynamic resistance (Louis 1979)
        rap_can, rap_Htree_In, _u_Hcan_, alpha, Ri_nearGround = \
            ResistanceCal.Ground_Aerodynamic_Resistance_1D(WindSpeed_top,MeteoData.Zatm,VerticalProfUrban,Gemeotry_m,Tcanyon,Tground,
                                                           Gemeotry_m.Height_canyon,dcan,zomcan,zom_ground,Gemeotry_m.Height_tree,
                                                           Gemeotry_m.Radius_tree,ColParam)

        # Calculate stomatal and leaf boundary resistances of tree
        if Ctree == 1 and ParVegTree.LAI > 0:

            # Calculate temperature and pressure at the height of tree [K]
            Tcanyon_tree = th_intp(Gemeotry_m.Height_tree)
            Pre_tree = Pre_intp(Gemeotry_m.Height_tree)

            # Leaf boundary resistance [s m^-1]
            rb_H = ResistanceCal.Leaf_BR(u_tree, Ttree - 273.15, Tcanyon_tree - 273.15, ParVegTree.d_leaf, alpha)
            # Stomatal resistance [s m^-1]
            rs_sun_H, rs_shd_H, Ci_sun_H, Ci_shd_H, _An_, _Rdark_, _Lpho_, _SIF_, _DCi_ = \
                ResistanceCal.Canopy_Resistance_An_Evolution(PAR_sun_H, PAR_shd_H, ParVegTree.LAI, ParVegTree.Kopt, ParVegTree.Knit,
                                                             Fsun_H, Fshd_H, Citm1_sun_H, Citm1_shd_H, MeteoData.Catm_CO2, rap_Htree_In,
                                                             rb_H, Ttree - 273.15, Pre_tree/100, Ds_tree, Psi_H_tm1,
                                                             ParVegTree.Psi_sto_50, ParVegTree.Psi_sto_00, ParVegTree.CT, ParVegTree.Vmax,
                                                             ParVegTree.DSE, ParVegTree.Ha, ParVegTree.FI, MeteoData.Catm_O2, ParVegTree.Do,
                                                             ParVegTree.a1,ParVegTree.go, ParVegTree.e_rel, ParVegTree.e_relN,
                                                             ParVegTree.gmes, ParVegTree.rjv)

        else:
            rb_H = numpy.inf
            rs_sun_H = numpy.inf
            rs_shd_H = numpy.inf
            Ci_sun_H = 0
            Ci_shd_H = 0
            Tcanyon_tree = Ttree

        # Calculate stomatal and leaf boundary resistances of ground vegetation
        if Cveg == 1 and ParVegGround.LAI > 0:

            # Leaf boundary resistance [s m^-1]
            rb_L = ResistanceCal.Leaf_BR(u_Lveg, Tveg - 273.15, Tcanyon - 273.15, ParVegGround.d_leaf, alpha)
            # Stomatal resistance [s m^-1]
            rs_sun_L, rs_shd_L, Ci_sun_L, Ci_shd_L, _An_, _Rdark_, _Lpho_, _SIF_, _DCi_ = \
                ResistanceCal.Canopy_Resistance_An_Evolution(PAR_sun_L, PAR_shd_L, ParVegGround.LAI, ParVegGround.Kopt, ParVegGround.Knit,
                                                             Fsun_L,Fshd_L, Citm1_sun_L, Citm1_shd_L, MeteoData.Catm_CO2, rap_can, rb_L,
                                                             Tveg - 273.15, MeteoData.Pre/100, Ds_canyon, Psi_L_tm1,
                                                             ParVegGround.Psi_sto_50, ParVegGround.Psi_sto_00, ParVegGround.CT,
                                                             ParVegGround.Vmax, ParVegGround.DSE, ParVegGround.Ha, ParVegGround.FI,
                                                             MeteoData.Catm_O2, ParVegGround.Do, ParVegGround.a1, ParVegGround.go,
                                                             ParVegGround.e_rel, ParVegGround.e_relN, ParVegGround.gmes, ParVegGround.rjv)
        else:
            rb_L = numpy.inf
            rs_sun_L = numpy.inf
            rs_shd_L = numpy.inf
            Ci_sun_L = 0
            Ci_shd_L = 0


        numpy.nan, numpy.nan, numpy.nan, Osat, Ohy, nVG, alpVG, Ks_Zs, L, Pe, O33, SPAR, numpy.nan, numpy.nan, RfH_Zs, \
        RfL_Zs, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan, numpy.nan = \
            SoilCal.Soil_Parameters_Total(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,
                                          ZR50_L,ZRmax_H,ZRmax_L,Zs)

        # Calculate soil resistance
        r_soil_bare, numpy.nan, alp_soil_bare = \
            ResistanceCal.Soil_Resistance(Tbare - 273.15, VerticalProfUrban.presProf[0]/100, WindSpeed_bottom, e_T_canyon,
                                          Int.IntGroundBare, Owater.OwGroundSoilBare[0],Ks_Zs[0], Osat[0], Ohy[0], L[0], Pe[0],
                                          O33[0], alpVG[0], nVG[0], SPAR)

        r_soil_veg, numpy.nan, alp_soil_veg = \
            ResistanceCal.Soil_Resistance(Tveg - 273.15, VerticalProfUrban.presProf[0]/100, WindSpeed_bottom, e_T_canyon,
                                          Int.IntGroundBare, Owater.OwGroundSoilVeg[0], Ks_Zs[0], Osat[0], Ohy[0], L[0], Pe[0],
                                          O33[0], alpVG[0], nVG[0], SPAR)

        # Calculate sensible heat fluxes [W m^-2]
        Himp = Cimp * (cp_atm*rho_atm * (Timp-Tcanyon) / rap_can)
        Hbare = Cbare * (cp_atm*rho_atm * (Tbare-Tcanyon) / rap_can)
        Hveg = Cveg * (cp_atm*rho_atm * (Tveg-Tcanyon) / (rb_L / (2*(ParVegGround.LAI+ParVegGround.SAI))+rap_can))
        Htree = Ctree * (cp_atm*rho_atm * (Ttree-Tcanyon_tree) / (rb_H/(2*(ParVegTree.LAI+ParVegTree.SAI))+rap_Htree_In))


        ###############
        leaf_dim = 0.72 * 0.05
        gHa = 1.4 * 0.135 * numpy.sqrt(u_tree / leaf_dim)
        cp_mol = 29.3
        Htree = 2*gHa*cp_mol*(Ttree-Tcanyon_tree)
        if Htree != 0:
            rap_Htree_In = (cp_atm*rho_atm * (Ttree-Tcanyon_tree))/Htree
        ###############

        # Potential evaporation from impervious surface [kg m^-2 s^-1]
        Eimp_pot = Cimp * (rho_atm * (qsat_T_imp-qcanyon)/ rap_can)
        # Potential evaporation from bare surface [kg m^-2 s^-1]
        Ebare_soil_pot = Cbare * (rho_atm * (alp_soil_bare*qsat_T_bare-qcanyon) / (rap_can+r_soil_bare))

        if Cveg == 0:
            # Potential evaporation from intercepted water on vegetated ground [kg m^-2 s^-1]
            Eveg_int_pot = 0
            # Potential evaporation from soil under vegetated ground [kg m^-2 s^-1]
            Eveg_soil_pot = 0
            # Potential transpiration from sunlit vegetated ground [kg m^-2 s^-1]
            TEveg_sun_pot = 0
            # Potential transpiration from shaded vegetated ground [kg m^-2 s^-1]
            TEveg_shd_pot = 0
        else:
            # Potential evaporation from intercepted water on vegetated ground [kg m^-2 s^-1]
            Eveg_int_pot = Cveg * (rho_atm*(qsat_T_veg-qcanyon)/ (rb_L/((ParVegGround.LAI+ParVegGround.SAI)*dw_L)+rap_can))
            # Potential evaporation from soil under vegetated ground [kg m^-2 s^-1]
            Eveg_soil_pot = Cveg * (rho_atm*(alp_soil_veg*qsat_T_veg-qcanyon) / (rap_can+r_soil_veg))
            # Potential transpiration from sunlit vegetated ground [kg m^-2 s^-1]
            TEveg_sun_pot = Cveg * (rho_atm*(qsat_T_veg-qcanyon)/ (rb_L/((ParVegGround.LAI)*Fsun_L*(1-dw_L))+rap_can +
                                                                   rs_sun_L / ((ParVegGround.LAI)*Fsun_L*(1-dw_L))))
            # Potential transpiration from shaded vegetated ground [kg m^-2 s^-1]
            TEveg_shd_pot = Cveg * (rho_atm * (qsat_T_veg - qcanyon)/ (rb_L / ((ParVegGround.LAI) * Fshd_L * (1 - dw_L)) + rap_can +
                                                                   rs_shd_L / ((ParVegGround.LAI) * Fshd_L * (1 - dw_L))))
        # Total potential transpiration from shaded and sunlit vegetated ground [kg m^-2 s^-1]
        TEveg_pot = Cveg * (TEveg_sun_pot + TEveg_shd_pot)

        if Ctree == 0:
            # Potential evaporation from intercepted water on tree [kg m^-2 s^-1]
            Etree_int_pot = 0
            # Potential transpiration from sunlit tree [kg m^-2 s^-1]
            TEtree_sun_pot = 0
            # Potential transpiration from shaded tree [kg m^-2 s^-1]
            TEtree_shd_pot = 0
        else:
            # Potential evaporation from intercepted water on tree [kg m^-2 s^-1]
            Etree_int_pot = Ctree * (rho_atm * (qsat_T_tree - q_tree)/ (rb_H / ((ParVegTree.LAI + ParVegTree.SAI) * dw_H) + rap_Htree_In))
            # Potential transpiration from sunlit tree [kg m^-2 s^-1]
            TEtree_sun_pot = Ctree * (rho_atm * (qsat_T_tree - q_tree)/ (rb_H / ((ParVegTree.LAI) * Fsun_H * (1 - dw_H)) +
                                                                      rap_Htree_In + rs_sun_H / ((ParVegTree.LAI) * Fsun_H * (1 - dw_H))))
            # Potential transpiration from shaded tree [kg m^-2 s^-1]
            TEtree_shd_pot = Ctree * (rho_atm * (qsat_T_tree - q_tree)/ (rb_H / ((ParVegTree.LAI) * Fshd_H * (1 - dw_H)) +
                                                                      rap_Htree_In + rs_shd_H / ((ParVegTree.LAI) * Fsun_H * (1 - dw_H))))
        # Total potential transpiration from shaded and sunlit tree [kg m^-2 s^-1]
        TEtree_pot = Ctree * (TEtree_sun_pot + TEtree_shd_pot)

        # Condition that evapotranspiration does not exceed available water
        # Water limitations of interception and ponding
        # Real max water evaporation from interception [kg m^-2 s^-1]
        Eimp = min(Eimp_pot, (Int.IntGroundImp / (1000 * ParCalculation.dts) * ParCalculation.rhow))

        # Real max water evaporation from interception [kg m^-2 s^-1]
        Ebare_pond = min(Ebare_soil_pot, (Int.IntGroundBare / (1000 * ParCalculation.dts) * ParCalculation.rhow))
        # The potential evaporation from bare surface is sum of evaporation from ponding + evaporation from soil under bare surface
        # Evaporation from soil under bare surface [kg m^-2 s^-1]
        Ebare_soil_pot = Ebare_soil_pot - Ebare_pond

        # Real max water evaporation from interception [kg m^-2 s^-1]
        Eveg_int = min(Eveg_int_pot, (Int.IntGroundVegPlant / (1000 * ParCalculation.dts) * ParCalculation.rhow))

        # Real max water evaporation from interception [kg m^-2 s^-1]
        Eveg_pond = min(Eveg_soil_pot, (Int.IntGroundVegGround / (1000 * ParCalculation.dts) * ParCalculation.rhow))
        # The potential evaporation from vegetated surface is sum of evaporation from ponding + evaporation from soil under vegetated surface
        # Evaporation from soil under vegetated surface [kg m^-2 s^-1]
        Eveg_soil_pot = Eveg_soil_pot - Eveg_pond

        # Real max water evaporation from interception [kg m^-2 s^-1]
        Etree_int = min(Etree_int_pot, (Int.IntTree / (1000 * ParCalculation.dts) * ParCalculation.rhow))

        # Water limitation to soil evaporation
        # Water mass flux in each soil layer [kg m^-2 s^-1]
        VavailImp_tm1 = [(Vwater.VGroundSoilImp[i] / ParCalculation.dts) * (ParCalculation.rhow / 1000) for i in range(len(Vwater.VGroundSoilImp))]
        # Water mass flux in each soil layer [kg m^-2 s^-1]
        VavailBare_tm1 = [(Vwater.VGroundSoilBare[i] / ParCalculation.dts) * (ParCalculation.rhow / 1000) for i in range(len(Vwater.VGroundSoilBare))]
        # Water mass flux in each soil layer [kg m^-2 s^-1]
        VavailVeg_tm1 = [(Vwater.VGroundSoilVeg[i] / ParCalculation.dts) * (ParCalculation.rhow / 1000) for i in range(len(Vwater.VGroundSoilVeg))]

        # Real evaporation from soil under bare surface [kg m^-2 s^-1]
        Ebare_soil = min(Ebare_soil_pot, VavailBare_tm1[0])
        # Real evaporation from soil under vegetated surface [kg m^-2 s^-1]
        Eveg_soil = min(Eveg_soil_pot, VavailVeg_tm1[0])
        # Available water in the first soil layer of bare and vegetated ground after soil evaporation
        VavailBare_tm1[0] = VavailBare_tm1[0] - Ebare_soil
        VavailVeg_tm1[0] = VavailVeg_tm1[0] - Eveg_soil

        # Water limitation to ground vegetation transpiration
        # Root fraction in the soil layers [-]
        RfH_Zs_Imp = numpy.zeros(len(Zs) - 1)
        RfH_Zs_Imp[:] = numpy.NaN

        _Zs_, _dz_, _ms_, _Osat_, _Ohy_, _nVG_, _alpVG_, _Ks_Zs_, _L_, _Pe_, _O33_,_SPAR_, _EvL_Zs_, _Inf_Zs_, RfH_Zs_ImpL2,\
        _RfL_Zs_, _Zinf_, _Kbot_, _Slo_pot_, _Dz_, _aR_, _aTop_, _rsd_, _lan_dry_, _lan_s_, _cv_s_ = \
            SoilCal.Soil_Parameters_Total(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,
                                          ZR50_L,ZRmax_H,ZRmax_L,Zs[2:])
        RfH_Zs_Imp[0:2] = [0,0]
        RfH_Zs_Imp[2:] = RfH_Zs_ImpL2

        # How much water is available per crown area
        if SPARTREE == 1:
            # Tree roots can access all water in the soil (imp, bare, veg)
            Ccrown = [Cveg*FractionsGround.fveg, Cbare*FractionsGround.fbare, Cimp*FractionsGround.fimp, Ctree*4*rad_tree]

            if Ccrown[0] == 0:
                Vavail_Veg_tm1_L = [numpy.NaN for i in range(len(VavailVeg_tm1))]
            else:
                Vavail_Veg_tm1_L = [Ccrown[0] / (Ccrown[0] + Ccrown[0] * Ccrown[3])* VavailVeg_tm1[i] for i in range(len(VavailVeg_tm1))]

            if Ccrown[0] * Ccrown[3] == 0 and (Ccrown[0] + Ccrown[0] * Ccrown[3]) == 0:
                Vavail_Veg_tm1_H = [numpy.NaN for i in range(len(VavailVeg_tm1))]
            else:
                Vavail_Veg_tm1_H = [Ccrown[0] * Ccrown[3] / (Ccrown[0] + Ccrown[0] * Ccrown[3])* VavailVeg_tm1[i] for i in range(len(VavailVeg_tm1))]

            if Ccrown[1] == 0 or Ccrown[3] == 0:
                Vavail_Bare_tm1_H = [numpy.NaN for i in range(len(VavailBare_tm1))]
            else:
                Vavail_Bare_tm1_H = [Ccrown[1] * Ccrown[3] / (Ccrown[1] * Ccrown[3])* VavailBare_tm1[i] for i in range(len(VavailBare_tm1))]

            if Ccrown[2] == 0 or Ccrown[3] == 0:
                Vavail_Imp_tm1_H = [numpy.NaN for i in range(len(VavailImp_tm1))]
            else:
                Vavail_Imp_tm1_H = [Ccrown[2] * Ccrown[3] / (Ccrown[2] * Ccrown[3])* VavailImp_tm1[i] for i in range(len(VavailImp_tm1))]

        # If the tree crown is smaller than the combined vegetated and bare fraction, then the trees only transpire from
        # these fractions. Otherwise, they also transpire from the impervious ground fraction.
        else:
            if (4*rad_tree) <= (FractionsGround.fveg+FractionsGround.fbare):
                Ccrown = [Cveg*FractionsGround.fveg, Cbare*FractionsGround.fbare, Cimp*FractionsGround.fimp,
                          Ctree*4*rad_tree/(FractionsGround.fveg+FractionsGround.fbare)]

                if Ccrown[0] == 0:
                    Vavail_Veg_tm1_L = [numpy.NaN for i in range(len(VavailVeg_tm1))]
                else:
                    Vavail_Veg_tm1_L = [Ccrown[0] / (Ccrown[0] + Ccrown[0]* Ccrown[3])* VavailVeg_tm1[i] for i in range(len(VavailVeg_tm1))]

                if Ccrown[0] * Ccrown[3] == 0 and (Ccrown[0] + Ccrown[0] * Ccrown[3]) == 0:
                    Vavail_Veg_tm1_H = [numpy.NaN for i in range(len(VavailVeg_tm1))]
                else:
                    Vavail_Veg_tm1_H = [Ccrown[0] * Ccrown[3] / (Ccrown[0] + Ccrown[0] * Ccrown[3]) * VavailVeg_tm1[i] for i in range(len(VavailVeg_tm1))]

                if Ccrown[1] == 0 or Ccrown[3] == 0:
                    Vavail_Bare_tm1_H = [numpy.NaN for i in range(len(VavailBare_tm1))]
                else:
                    Vavail_Bare_tm1_H = [Ccrown[1] * Ccrown[3] / (Ccrown[1] * Ccrown[3]) * VavailBare_tm1[i] for i in range(len(VavailBare_tm1))]

                Vavail_Imp_tm1_H = [0. * VavailImp_tm1[i] for i in range(len(VavailImp_tm1))]

            elif (4*rad_tree)>(FractionsGround.fveg+FractionsGround.fbare):
                Ccrown = [Cveg*FractionsGround.fveg, Cbare*FractionsGround.fbare, Cimp*FractionsGround.fimp, Ctree*1,
                          Ctree*((4*rad_tree)-(FractionsGround.fveg+FractionsGround.fbare))/FractionsGround.fimp]

                if Ccrown[0] == 0 and (Ccrown[0] + Ccrown[0] * Ccrown[3]) == 0:
                    Vavail_Veg_tm1_L = [numpy.NaN for i in range(len(VavailVeg_tm1))]
                else:
                    Vavail_Veg_tm1_L = [Ccrown[0] / (Ccrown[0] + Ccrown[0] * Ccrown[3]) * VavailVeg_tm1[i] for i in range(len(VavailVeg_tm1))]

                if Ccrown[0] * Ccrown[3] == 0 and (Ccrown[0] + Ccrown[0] * Ccrown[3]) == 0:
                    Vavail_Veg_tm1_H = [numpy.NaN for i in range(len(VavailVeg_tm1))]
                else:
                    Vavail_Veg_tm1_H = [Ccrown[0] * Ccrown[3] / (Ccrown[0] + Ccrown[0] * Ccrown[3])* VavailVeg_tm1[i] for i in range(len(VavailVeg_tm1))]

                if Ccrown[1] == 0 or Ccrown[3] == 0:
                    Vavail_Bare_tm1_H = [numpy.NaN for i in range(len(VavailBare_tm1))]
                else:
                    Vavail_Bare_tm1_H = [Ccrown[1] * Ccrown[3] / (Ccrown[1] * Ccrown[3]) * VavailBare_tm1[i] for i in range(len(VavailBare_tm1))]

                if Ccrown[2] == 0 or Ccrown[4] == 0:
                    Vavail_Imp_tm1_H = [numpy.NaN for i in range(len(VavailImp_tm1))]
                else:
                    Vavail_Imp_tm1_H = [Ccrown[2] * Ccrown[4] / (Ccrown[2] * Ccrown[4]) * VavailImp_tm1[i] for i in range(len(VavailImp_tm1))]

        # Minimum available and extractable water for plants [kg m^-2 s^-1]
        Vavail_Veg_tm1_L = sum(min([Vavail_Veg_tm1_L[i]*RfL_Zs[i] for i in range(len(Vavail_Veg_tm1_L))],
                                   [ExWater.ExWaterGroundVeg_L[i]*(ParCalculation.rhow / 1000) for i in range(len(ExWater.ExWaterGroundVeg_L))]))
        Vavail_Veg_tm1_H = sum(min([Vavail_Veg_tm1_H[i]*RfH_Zs[i] for i in range(len(Vavail_Veg_tm1_H))],
                                   [ExWater.ExWaterGroundVeg_H[i]*(ParCalculation.rhow / 1000) for i in range(len(ExWater.ExWaterGroundVeg_H))]))
        Vavail_Bare_tm1_H = sum(min([Vavail_Bare_tm1_H[i]*RfH_Zs[i] for i in range(len(Vavail_Bare_tm1_H))],
                                    [ExWater.ExWaterGroundBare_H[i]*(ParCalculation.rhow / 1000) for i in range(len(ExWater.ExWaterGroundBare_H))]))
        Vavail_Imp_tm1_H = sum(min([Vavail_Imp_tm1_H[i]*RfH_Zs_Imp[i] for i in range(len(Vavail_Imp_tm1_H))],
                                   [ExWater.ExWaterGroundImp_H[i]*(ParCalculation.rhow / 1000) for i in range(len(ExWater.ExWaterGroundImp_H))]))

        # Water limitation for tree and ground vegetation transpiration [kg m^-2 s^-1]
        Vavail_plant_tm1_H = (Ccrown[0] * Vavail_Veg_tm1_H + Ccrown[1] * Vavail_Bare_tm1_H + Ccrown[2] * Vavail_Imp_tm1_H) / Ccrown[3]
        Vavail_plant_tm1_L = Vavail_Veg_tm1_L
        TEveg = min(TEveg_pot, Vavail_plant_tm1_L)
        TEtree = min(TEtree_pot, Vavail_plant_tm1_H)

        # Evapotranspiration [kg m^-2 s^-1]
        Ebare = Cbare * (Ebare_pond + Ebare_soil)
        Eveg = Cveg * (Eveg_int + Eveg_pond + Eveg_soil + TEveg)
        Etree = Ctree * (Etree_int + TEtree)

        # Latent heat [W m^-2]
        LEimp = Cimp * (L_heat * Eimp)
        LEbare_pond = Cbare * (L_heat * Ebare_pond)
        LEbare_soil = Cbare * (L_heat * Ebare_soil)
        LEbare = Cbare * (LEbare_pond + LEbare_soil)
        LEveg_int = Cveg * (L_heat * Eveg_int)
        LEveg_pond = Cveg * (L_heat * Eveg_pond)
        LEveg_soil = Cveg * (L_heat * Eveg_soil)
        LTEveg = Cveg * (L_heat * TEveg)
        LEveg = Cveg * (LEveg_int + LEveg_pond + LEveg_soil + LTEveg)
        LEtree_int = Ctree * (L_heat * Etree_int)
        LTEtree = Ctree * (L_heat * TEtree)
        LEtree = Ctree * (LEtree_int + LTEtree)

        class thb_ground_Def():
            pass
        self.thb_ground = thb_ground_Def()
        self.thb_ground.ground_imp = Himp/(cp_atm*rho_atm)
        self.thb_ground.ground_bare = Hbare/(cp_atm*rho_atm)
        self.thb_ground.ground_veg = Hveg/(cp_atm*rho_atm)
        self.thb_ground.tree = Htree/(cp_atm*rho_atm)
        self.thb_ground.ustar = ResistanceCal.Ustar_Atm

        class qhb_ground_Def():
            pass
        self.qhb_ground = qhb_ground_Def()
        self.qhb_ground.ground_imp = LEimp/(L_heat*rho_atm)
        self.qhb_ground.ground_bare = (LEbare_pond+LEbare_soil)/(L_heat*rho_atm)
        self.qhb_ground.ground_veg = (LEveg_int+LEveg_pond+LEveg_soil+LTEveg)/(L_heat*rho_atm)
        self.qhb_ground.tree = LEtree/(L_heat*rho_atm)

        class HfGruond_otherparam_Def():
            pass
        self.HfGruond_otherparam = HfGruond_otherparam_Def()
        self.HfGruond_otherparam.Ri_nearGround = 0
        self.HfGruond_otherparam.Utotal_nearGround = 0
        self.HfGruond_otherparam.T_nearGround = Tcanyon
        self.HfGruond_otherparam.Tground = Tground


        return Himp,Hbare,Hveg,Htree,Eimp,Ebare_pond,Ebare_soil,Eveg_int,Eveg_pond,Eveg_soil,TEveg,Etree_int,TEtree, Ebare,\
               Eveg,Etree, LEimp,LEbare_pond,LEbare_soil,LEveg_int, LEveg_pond,LEveg_soil,LTEveg,LEtree_int,LTEtree, LEbare,\
               LEveg,LEtree, Ci_sun_H,Ci_shd_H,Ci_sun_L,Ci_shd_L, rap_can,rap_Htree_In,rb_H,rb_L, r_soil_bare,r_soil_veg,\
               alp_soil_bare,alp_soil_veg, rs_sun_L,rs_shd_L,rs_sun_H,rs_shd_H,Fsun_L,Fshd_L,dw_L



