import os
import numpy
import math
from Water_Functions import Water_Calculations
import copy

"""
Compute water balance at the roof
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Originally developed by Naika Meili
Last update: May 2021
"""

class WaterBalanceRoof_Def(object):

    def __init__(self,Int_init,TE_init,Owater_init,E_init,ParSoilRoof,ExWater_init,V_init,SoilPotW_init):

        class Int_Def():
            pass
        self.Int = Int_Def()
        self.Int.IntRoofImp = Int_init
        self.Int.IntRoofVegPlant = Int_init
        self.Int.IntRoofVegGround = Int_init
        self.Int.IntRooftot = Int_init

        self.Owater_OwRoofSoilVeg = [Owater_init for i in range(ParSoilRoof.ms)]
        self.TEfluxRoofVeg = TE_init
        self.EfluxRoofVegSoil = E_init
        self.ExWaterRoof_L = [ExWater_init for i in range(ParSoilRoof.ms)]
        self.VRoofSoil = V_init
        self.SoilPotWRoof_L = SoilPotW_init

        class Runoff_Def():
            pass
        self.Runoff = Runoff_Def()
        self.Runoff.QRoofImp = 0
        self.Runoff.QRoofVegDrip = 0
        self.Runoff.QRoofVegPond = 0
        self.Runoff.QRoofVegSoil = 0

        class Leakage_Def():
            pass
        self.Leakage = Leakage_Def()
        self.Leakage.LkRoofImp = 0
        self.Leakage.LkRoofVeg = 0
        self.Leakage.LkRoof = 0

        self.RunonRoofTot = 0
        self.RunoffRoofTot = 0
        self.Anthp = 0

        class dInt_dt_Def():
            pass
        self.dInt_dt = dInt_dt_Def()
        self.dInt_dt.dInt_dtRoofImp = 0
        self.dInt_dt.dInt_dtRoofVegPlant = 0
        self.dInt_dt.dInt_dtRoofVegGround = 0
        self.dInt_dt.dInt_dtRooftot = 0

        class WB_Def():
            pass
        self.WB = WB_Def()
        self.WB.In_imp = 0
        self.WB.In_veg = 0
        self.WB.In_ground = 0
        self.WB.soil = 0
        self.WB.imp_tot = 0
        self.WB.veg_tot = 0
        self.WB.Total = 0


    def WBSolver_Roof(self,Eroof_imp,Eroof_veg,Eroof_ground,EfluxRoofVegSoil_EB,TEfluxRoofVeg_EB,MeteoData,FractionsRoof,
                      ParSoilRoof,ParCalculation,ParVegRoof,Anthropogenic):
        """
        ------
        INPUT:
        Eroof_imp: Evaporative flux from impervious roof [kg m^-2 s^-1]
        Eroof_veg: Evaporative flux from intercepted water by vegetation [kg m^-2 s^-1]
        Eroof_ground: Evaporative flux from ponding water on ground under vegetation [kg m^-2 s^-1]
        EfluxRoofVegSoil_EB: Evaporative flux from soil under vegetation (first layer) [kg m^-2 s^-1]
        TEfluxRoofVeg_EB: Transpiration flux [kg m^-2 s^-1]
        MeteoData: Forcing variables
        FractionsRoof:
        ParSoilRoof: Roof soil parameters
        ParCalculation: General calculation parameters
        ParVegRoof: Roof vegetation parameters
        Anthropogenic:
        -------
        OUTPUT:
        Runoff: Runoff and water table rise
                QRoofImp: Runoff from impervious surface [mm s^-1]
                QRoofVegDrip: Runoff from vegetation (Summation of dripping, saturation excess, and throughfall) [mm s^-1]
                QRoofVegPond: Runoff from ground under vegetation [mm s^-1]
                QRoofVegSoil: Water table rise that reaches the surface level (Dunne Runoff) [mm]
        Leakage: Leakage (infiltration) from the surfaces
                 LkRoofImp: leakage(infiltration) of impervious surface [mm s^-1]
                 LkRoofVeg: leakage at the bedrock [mm s^-1]
                 LkRoof: total leakage [mm s^-1]
        RunonRoofTot: Total runon [mm s^-1]
        RunoffRoofTot: Total runoff [mm s^-1]
        Anthp: Anthropogenic water [mm s^-1]
        dInt_dt: Change in intercepted water (storage)
                 dInt_dtRoofImp: Water storage at the impervious surface [mm s^-1]
                 dInt_dtRoofVegPlant: Water storage at the vegetation [mm s^-1]
                 dInt_dtRoofVegGround: Water storage at the ground under vegetation [mm s^-1]
                 dInt_dtRooftot: Total water storage [mm s^-1]
        dVRoofSoil_dt: Change in the soil water content (storage term) [mm s^-1]
        fRoofVeg: Infiltration rate [mm s^-1]
        VRoofSoil: Volume of water in each soil layer per unit area [mm]
        SoilPotWRoof_L: Soil Water Potential for Second Layer of Vegetation [MPa]
        ExWaterRoof_L: Maximum extractable water [mm m^2 m^-2 ground s^-1]
        Int: Intercepted water
             IntRoofImp: Intercepted water at the impervious surface [mm]
             IntRoofVegPlant: Intercepted water at the vegetation [mm]
             IntRoofVegGround: Intercepted water at the surface under vegetation [mm]
             IntRooftot: Total Intercepted water [mm]
        WB: Water balance (it should be close to zero)
            In_imp: Water balance for the impervious surface [mm s^-1]
            In_veg: Water balance for the vegetation [mm s^-1]
            In_ground: Water balance for the ground under vegetation [mm s^-1]
            soil: Water balance for the soil column [mm s^-1]
            imp_tot: Water balance for the impervious surface [mm s^-1]
            veg_tot: Water balance for the vegetation and underneath ground [mm s^-1]
            Total: Total water balance [mm s^-1]
        Owater_OwRoofSoilVeg: Water content (soil moisture) in each layer [-]
        """

        Water_Cal = Water_Calculations()

        Runon_tm1 = copy.copy(self.RunonRoofTot)

        # Calculate Water fluxes
        # Impervious ground
        q_runoff_imp, self.Int.IntRoofImp, dIn_imp_dt, Lk_imp, WBalance_In_imp = \
            Water_Cal.Water_Impervious(MeteoData.Rain, Runon_tm1,Eroof_imp, self.Int.IntRoofImp,ParCalculation.dts,ParCalculation.rhow,
                                       ParSoilRoof.In_max_imp,ParSoilRoof.Kimp)
        # Vegetation
        q_runoff_veg, self.Int.IntRoofVegPlant, dIn_veg_dt, WBalance_In_veg = \
            Water_Cal.Water_Vegetation(MeteoData.Rain,Eroof_veg,self.Int.IntRoofVegPlant,ParSoilRoof.Sp_In,ParVegRoof.LAI,
                                       ParVegRoof.SAI,ParCalculation.rhow,ParCalculation.dts)
        # Ground under vegetation
        q_runoff_ground, self.Int.IntRoofVegGround, dIn_ground_dt, f_ground, WBalance_In_ground = \
            Water_Cal.Water_Ground(q_runoff_veg + Anthropogenic.Waterf_roof/3600,Runon_tm1,Eroof_ground,self.Owater_OwRoofSoilVeg,
                                   self.Int.IntRoofVegGround,ParSoilRoof.In_max_ground,ParSoilRoof.Pcla,ParSoilRoof.Psan,
                                   ParSoilRoof.Porg,ParSoilRoof.Kfc,ParSoilRoof.Phy,ParSoilRoof.SPAR,ParSoilRoof.Kbot,
                                   ParVegRoof.CASE_ROOT,ParVegRoof.CASE_ROOT,0,ParVegRoof.ZR95,0,ParVegRoof.ZR50,0,ParVegRoof.ZRmax,
                                   ParSoilRoof.Zs, ParCalculation.dts, ParCalculation.rhow)
        # Soil column
        V, self.Owater_OwRoofSoilVeg, self.OSwRoofSoil, Lk, _Psi_s_H_, Psi_s, _Exwat_H_, Exwat, Rd, self.TEfluxRoofVeg, _TE_H_, \
        self.EfluxRoofVegSoil, dV_dt, WBalance_soil, Psi_soil, Ko = \
            Water_Cal.Water_Soil(self.Owater_OwRoofSoilVeg, f_ground, 0, TEfluxRoofVeg_EB, EfluxRoofVegSoil_EB,
                                 [0.*ParSoilRoof.Zs[i] for i in range(len(ParSoilRoof.Zs))],ParCalculation.dts,ParSoilRoof.Pcla,
                                 ParSoilRoof.Psan,ParSoilRoof.Porg,ParSoilRoof.Kfc,ParSoilRoof.Phy,ParSoilRoof.SPAR,ParSoilRoof.Kbot,
                                 ParVegRoof.CASE_ROOT,ParVegRoof.CASE_ROOT,0,ParVegRoof.ZR95,0,ParVegRoof.ZR50,0,ParVegRoof.ZRmax,0,
                                 ParVegRoof.Rrootl,0,ParVegRoof.PsiL50,0,ParVegRoof.PsiX50,ParSoilRoof.Zs,ParCalculation.rhow)

        # Calculate total runoff as the fraction of total excess water that leaves the system [mm s^-1]
        self.RunoffRoofTot = FractionsRoof.Per_runoff * (FractionsRoof.fimp*q_runoff_imp + FractionsRoof.fveg*(q_runoff_ground + Rd/ParCalculation.dts))
        # Calculate total runon as the fraction of total excess water that leaves the system [mm s^-1]
        self.RunonRoofTot = (1-FractionsRoof.Per_runoff) * (FractionsRoof.fimp*q_runoff_imp +
                                                            FractionsRoof.fveg*(q_runoff_ground + Rd/ParCalculation.dts))

        # Water balance for the impervious surface [mm s^-1]
        WBalance_imp_tot = MeteoData.Rain + Runon_tm1 - Eroof_imp*1000/ParCalculation.rhow - q_runoff_imp - Lk_imp - dIn_imp_dt
        # Water balance for the vegetated surface [mm s^-1]
        WBalance_veg_tot = MeteoData.Rain + Runon_tm1 + Anthropogenic.Waterf_roof/3600 - \
                           (Eroof_veg + Eroof_ground + self.EfluxRoofVegSoil + self.TEfluxRoofVeg)*1000/ParCalculation.rhow - \
                           Lk - q_runoff_ground - Rd/ParCalculation.dts - dIn_veg_dt - dIn_ground_dt - dV_dt

        # Calculate total water evaporated [mm s^-1]
        E_tot = (FractionsRoof.fimp*Eroof_imp +
                 FractionsRoof.fveg*(Eroof_veg + Eroof_ground + self.EfluxRoofVegSoil + self.TEfluxRoofVeg))*1000/ParCalculation.rhow
        # Calculate total leakage [mm s^-1]
        Leak_tot = FractionsRoof.fimp*Lk_imp + FractionsRoof.fveg*Lk
        # Calculate total storage [mm s^-1]
        Storage_tot = FractionsRoof.fimp*dIn_imp_dt + FractionsRoof.fveg*(dIn_veg_dt + dIn_ground_dt + dV_dt)

        # Calculate total water balance [mm s^-1]
        WBalance_tot = MeteoData.Rain + Runon_tm1 + FractionsRoof.fveg*Anthropogenic.Waterf_roof/3600 - E_tot - Leak_tot - \
                       self.RunoffRoofTot - self.RunonRoofTot - Storage_tot


        self.Runoff.QRoofImp = q_runoff_imp
        self.Runoff.QRoofVegDrip = q_runoff_veg
        self.Runoff.QRoofVegPond = q_runoff_ground
        self.Runoff.QRoofVegSoil = Rd

        self.Leakage.LkRoofImp = Lk_imp
        self.Leakage.LkRoofVeg = Lk
        self.Leakage.LkRoof = Leak_tot

        self.Anthp = Anthropogenic.Waterf_roof/3600

        self.dInt_dt.dInt_dtRoofImp = dIn_imp_dt
        self.dInt_dt.dInt_dtRoofVegPlant = dIn_veg_dt
        self.dInt_dt.dInt_dtRoofVegGround = dIn_ground_dt
        self.dInt_dt.dInt_dtRooftot = dIn_imp_dt*FractionsRoof.fimp + (dIn_veg_dt+dIn_ground_dt)*FractionsRoof.fveg

        self.dVRoofSoil_dt = dV_dt
        self.fRoofVeg = f_ground
        self.VRoofSoil = V
        self.SoilPotWRoof_L = Psi_s
        self.ExWaterRoof_L = Exwat

        self.WB.In_imp = WBalance_In_imp
        self.WB.In_veg = WBalance_In_veg
        self.WB.In_ground = WBalance_In_ground
        self.WB.soil = WBalance_soil
        self.WB.imp_tot = WBalance_imp_tot
        self.WB.veg_tot = WBalance_veg_tot
        self.WB.Total = WBalance_tot

        self.Int.IntRooftot = self.Int.IntRoofImp*FractionsRoof.fimp + (self.Int.IntRoofVegGround+self.Int.IntRoofVegPlant)*FractionsRoof.fveg