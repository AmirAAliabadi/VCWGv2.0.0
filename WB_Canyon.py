import numpy
from Soil_Functions import Soil_Calculations
from scipy.integrate import odeint
from Water_Functions import Water_Calculations
import copy

"""
Compute water balance in the canyon
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Originally developed by Naika Meili
Last update: May 2021
"""

class WaterBalanceCanyon_Def(object):

    def __init__(self,Owater_init,Int_init,Runon_init,ParSoilGround,Vwater_init,ExWater_init):

        class Int_Def():
            pass
        self.Int = Int_Def()
        self.Int.IntGroundImp = Int_init
        self.Int.IntGroundBare = Int_init
        self.Int.IntGroundVegPlant = Int_init
        self.Int.IntGroundVegGround = Int_init
        self.Int.IntTree = Int_init

        class Owater_Def():
            pass
        self.Owater = Owater_Def()
        self.Owater.OwGroundSoilImp = [Owater_init for i in range(ParSoilGround.ms)]
        self.Owater.OwGroundSoilBare = [Owater_init for i in range(ParSoilGround.ms)]
        self.Owater.OwGroundSoilVeg = [Owater_init for i in range(ParSoilGround.ms)]
        self.Owater.OwGroundSoilTot = [Owater_init for i in range(ParSoilGround.ms)]

        class Runon_Def():
            pass
        self.Runon = Runon_Def()
        self.Runon.RunonGroundTot = Runon_init
        self.Runon.RunonUrban = Runon_init

        class Qinlat_Def():
            pass
        self.Qinlat = Qinlat_Def()
        self.Qinlat.Qin_imp = [0.0 for i in range(ParSoilGround.ms)]
        self.Qinlat.Qin_bare = [0.0 for i in range(ParSoilGround.ms)]
        self.Qinlat.Qin_veg = [0.0 for i in range(ParSoilGround.ms)]
        self.Qinlat.Qin_bare2imp = [0.0 for i in range(ParSoilGround.ms)]
        self.Qinlat.Qin_veg2imp = [0.0 for i in range(ParSoilGround.ms)]
        self.Qinlat.Qin_veg2bare = [0.0 for i in range(ParSoilGround.ms)]
        self.Qinlat.Qin_imp2bare = [0.0 for i in range(ParSoilGround.ms)]
        self.Qinlat.Qin_bare2veg = [0.0 for i in range(ParSoilGround.ms)]
        self.Qinlat.Qin_imp2veg = [0.0 for i in range(ParSoilGround.ms)]

        class Runoff_Def():
            pass
        self.Runoff = Runoff_Def()
        self.Runoff.QGroundImp = 0.0
        self.Runoff.QGroundBarePond = 0.0
        self.Runoff.QGroundBareSoil = 0.0
        self.Runoff.QTree = 0.0
        self.Runoff.QGroundVegDrip = 0.0
        self.Runoff.QGroundVegPond = 0.0
        self.Runoff.QGroundVegSoil = 0.0
        self.Runoff.RunoffGroundTot = Runon_init
        self.Runoff.RunoffUrban = Runon_init

        class dInt_dt_Def():
            pass
        self.dInt_dt = dInt_dt_Def()
        self.dInt_dt.dInt_dtGroundImp = 0.0
        self.dInt_dt.dInt_dtGroundBare = 0.0
        self.dInt_dt.dInt_dtGroundVegPlant = 0.0
        self.dInt_dt.dInt_dtGroundVegGround = 0.0
        self.dInt_dt.dInt_dtTree = 0.0

        class Infiltration_Def():
            pass
        self.Infiltration = Infiltration_Def()
        self.Infiltration.fGroundBare = 0.0
        self.Infiltration.fGroundVeg = 0.0
        self.Infiltration.fGroundImp = 0.0

        class Vwater_Def():
            pass
        self.Vwater = Vwater_Def()
        self.Vwater.VGroundSoilImp = Vwater_init
        self.Vwater.VGroundSoilBare = Vwater_init
        self.Vwater.VGroundSoilVeg = Vwater_init
        self.Vwater.VGroundSoilTot = Vwater_init

        class Leakage_Def():
            pass
        self.Leakage = Leakage_Def()
        self.Leakage.LkGroundImp = 0
        self.Leakage.LkGroundBare = 0
        self.Leakage.LkGroundVeg = 0
        self.Leakage.LkGround = 0

        class ExWater_Def():
            pass
        self.ExWater = ExWater_Def()
        self.ExWater.ExWaterGroundImp_H = [ExWater_init for i in range(ParSoilGround.ms)]
        self.ExWater.ExWaterGroundImp_L = [ExWater_init for i in range(ParSoilGround.ms)]
        self.ExWater.ExWaterGroundBare_H = [ExWater_init for i in range(ParSoilGround.ms)]
        self.ExWater.ExWaterGroundBare_L = [ExWater_init for i in range(ParSoilGround.ms)]
        self.ExWater.ExWaterGroundVeg_H = [ExWater_init for i in range(ParSoilGround.ms)]
        self.ExWater.ExWaterGroundVeg_L = [ExWater_init for i in range(ParSoilGround.ms)]
        self.ExWater.ExWaterGroundTot_H = [ExWater_init for i in range(ParSoilGround.ms)]
        self.ExWater.ExWaterGroundTot_L = [ExWater_init for i in range(ParSoilGround.ms)]

        class SoilPotW_Def():
            pass
        self.SoilPotW = SoilPotW_Def()
        self.SoilPotW.SoilPotWGroundImp_H = 0.0
        self.SoilPotW.SoilPotWGroundBare_H = 0.0
        self.SoilPotW.SoilPotWGroundVeg_H = 0.0
        self.SoilPotW.SoilPotWGroundTot_H = 0.0
        self.SoilPotW.SoilPotWGroundImp_L = 0.0
        self.SoilPotW.SoilPotWGroundBare_L = 0.0
        self.SoilPotW.SoilPotWGroundVeg_L = 0.0
        self.SoilPotW.SoilPotWGroundTot_L = 0.0

        class dVwater_dt_Def():
            pass
        self.dVwater_dt = dVwater_dt_Def()
        self.dVwater_dt.dVGroundSoilImp_dt = 0
        self.dVwater_dt.dVGroundSoilBare_dt = 0
        self.dVwater_dt.dVGroundSoilVeg_dt = 0
        self.dVwater_dt.dVGroundSoilTot_dt = 0

        class WBIndv_Def:
            pass
        self.WBIndv = WBIndv_Def()
        self.WBIndv.WB_In_tree = 0
        self.WBIndv.WB_In_gveg = 0
        self.WBIndv.WB_In_gimp = 0
        self.WBIndv.WB_In_gbare = 0
        self.WBIndv.WB_Pond_gveg = 0
        self.WBIndv.WB_Soil_gimp = 0
        self.WBIndv.WB_Soil_gbare = 0
        self.WBIndv.WB_Soil_gveg = 0

        class WBTot_Def():
            pass
        self.WBTot = WBTot_Def()
        self.WBTot.WBSurf_tree = 0
        self.WBTot.WBSurf_imp = 0
        self.WBTot.WBSurf_bare = 0
        self.WBTot.WBSurf_veg = 0
        self.WBTot.WBSoil_imp = 0
        self.WBTot.WBSoil_bare = 0
        self.WBTot.WBSoil_veg = 0
        self.WBTot.WBImp_tot = 0
        self.WBTot.WBBare_tot = 0
        self.WBTot.WBVeg_tot = 0
        self.WBTot.WBCanyon_flux = 0
        self.WBTot.WBTree_level = 0
        self.WBTot.WBGround_level = 0
        self.WBTot.WBSoil_level = 0
        self.WBTot.WBCanyon_level = 0

        class Rd_Def():
            pass
        self.Rd = Rd_Def()
        self.Rd.Rd_gbare = 0
        self.Rd.Rd_gveg = 0
        self.Rd.Rd_gimp = 0

        class TE_Def():
            pass
        self.TE = TE_Def()
        self.TE.TEgveg_imp = 0
        self.TE.TEtree_imp = 0
        self.TE.TEgveg_bare = 0
        self.TE.TEtree_bare = 0
        self.TE.TEgveg_veg = 0
        self.TE.TEtree_veg = 0

        self.RainGround = 0
        self.Anth_gbare = 0
        self.Anth_gveg = 0


    def WBSolver_Canyon(self,MeteoData,Etree_In,Egveg_In,Egimp_Pond,Egbare_Pond,Egveg_Pond,Egbare_soil,Egveg_soil,TEgveg,
                        TEtree,ParSoilGround,ParInterceptionTree,ParCalculation,ParVegGround,ParVegTree,FractionsGround,
                        geometry,ParTree,Gemeotry_m,Anthropogenic):
        """
        ------
        INPUT:
        MeteoData: Forcing variables
        Etree_In: Evaporation from intercepted water on trees [kg m^-2 s^-1]
        Egveg_In: Evaporation from intercepted water on low vegetation [kg m^-2 s^-1]
        Egimp_Pond: Evaporation from intercepted water on impervious surface [kg m^-2 s^-1]
        Egbare_Pond: Evaporation from intercepted water on bare ground [kg m^-2 s^-1]
        Egveg_Pond: Evaporation from intercepted water on ground under vegetation [kg m^-2 s^-1]
        Egbare_soil: Evaporation from soil under bare ground [kg m^-2 s^-1]
        Egveg_soil: Evaporation from soil under vegetated ground [kg m^-2 s^-1]
        TEgveg: Transpiration from low vegetation [kg m^-2 s^-1]
        TEtree: Transpiration from trees [kg m^-2 s^-1]
        ParSoilGround: Ground soil parameters
        ParInterceptionTree: specific water retained by the tree [mm m^2 VEG area m^-2 plant area]
        ParCalculation: General calculation parameters
        ParVegGround: Ground vegetation parameters
        ParVegTree: Trees parameter
        FractionsGround: Fractions of ground covered by vegetation,bare soil, and impervious surface
        geometry: Normalized geometric parameters
        ParTree: Trees parameter
        Gemeotry_m: Geometric parameters
        Anthropogenic: Anthropogenic water
        -------
        OUTPUT:
        Int: Intercepted water
             IntGroundImp: Intercepted water on the impervious surface [mm]
             IntGroundBare: Intercepted water on the bare surface [mm]
             IntGroundVegPlant: Intercepted water on the low vegetation [mm]
             IntGroundVegGround: Intercepted water on the vegetated surface [mm]
             IntTree: Intercepted water on the tree [mm]
        Owater: Water content (soil moisture) in each layer
                OwGroundSoilImp: Water content in the soil column layers under impervious surface [-]
                OwGroundSoilBare: Water content in the soil column layers under bare surface [-]
                OwGroundSoilVeg: Water content in the soil column layers under vegetated surface [-]
                OwGroundSoilTot: Total water content in the soil column layers [-]
        Runon: Runon at the ground surfaces
               RunonGroundTot: Total runon at the ground [mm s^-1]
               RunonUrban: Total runon in the urban [mm s^-1]. Note: it is calculated in the VCWG_Hydrology.py
        Qinlat: Lateral water flux between soil columns
                Qin_imp: Net change in water of each layer in impervious ground soil column [mm]
                Qin_bare: Net change in water of each layer in bare ground soil column [mm]
                Qin_veg: Net change in water of each layer in vegetated ground soil column [mm]
                Qin_bare2imp: Water flux from bare ground soil column to impervious ground soil column [mm]
                Qin_veg2imp: Water flux from vegetated ground soil column to impervious ground soil column [mm]
                Qin_veg2bare: Water flux from vegetated ground soil column to bare ground soil column [mm]
                Qin_imp2bare: Water flux from impervious ground soil column to bare ground soil column [mm]
                Qin_bare2veg: Water flux from bare ground soil column to vegetated ground soil column [mm]
                Qin_imp2veg: Water flux from impervious ground soil column to vegetated ground soil column [mm]
        Runoff: Runoff from different surfaces
                QGroundImp: Runoff from impervious surface [mm s^-1]
                QGroundBarePond: Runoff from bare surface [mm s^-1]
                QGroundBareSoil: Water table rise that reaches the bare surface level (Dunne Runoff) [mm]
                QTree: Runoff from trees (dripping, saturation excess, and throughfall) [mm s^-1]
                QGroundVegDrip: Runoff from low vegetation (dripping, saturation excess, and throughfall) [mm s^-1]
                QGroundVegPond: Runoff from surface under low vegetation [mm s^-1]
                QGroundVegSoil: Water table rise that reaches the vegetated surface level (Dunne Runoff) [mm]
                RunoffGroundTot: Total runoff at the ground [mm s^-1]
                RunoffUrban: Total runoff in the urban [mm s^-1]. Note: it is calculated in the VCWG_Hydrology.py
        dInt_dt: Storage term in water balance equation
                dInt_dtGroundImp: Change in intercepted water on the impervious surface [mm s^-1]
                dInt_dtGroundBare: Change in intercepted water on the bare surface [mm s^-1]
                dInt_dtGroundVegPlant: Change in intercepted water on the low vegetation [mm s^-1]
                dInt_dtGroundVegGround: Change in intercepted water on the surface under low vegetation [mm s^-1]
                dInt_dtTree: Change in intercepted water on the trees [mm s^-1]
        Infiltration: Infiltration rate
                      fGroundBare: Infiltration to the bare ground [mm s^-1]
                      fGroundVeg: Infiltration to the vegetated ground [mm s^-1]
                      fGroundImp: Infiltration to the impervious ground [mm s^-1]
        Vwater: Volume of water
                VGroundSoilImp: Volume of water in each layer of impervious ground soil column per unit area [mm]
                VGroundSoilBare: Volume of water in each layer of bare ground soil column per unit area [mm]
                VGroundSoilVeg: Volume of water in each layer of vegetated ground soil column per unit area [mm]
                VGroundSoilTot: Total volume of water in each layer soil column per unit area [mm]
        Leakage: Leakage at bedrock
                 LkGroundImp: Leakage at the bedrock of impervious ground column [mm s^-1]
                 LkGroundBare: Leakage at the bedrock of bare ground column [mm s^-1]
                 LkGroundVeg: Leakage at the bedrock of vegetated ground column [mm s^-1]
                 LkGround: Total leakage [mm s^-1]
        ExWater: Maximum extractable water by soil columns
                 ExWaterGroundImp_H: Maximum extractable water by tree from impervious ground soil column [mm m^2 m^-2 ground s^-1]
                 ExWaterGroundImp_L: Maximum extractable water by low vegetation from impervious ground soil column [mm m^2 m^-2 ground s^-1]
                 ExWaterGroundBare_H: Maximum extractable water by tree from bare ground soil column [mm m^2 m^-2 ground s^-1]
                 ExWaterGroundBare_L: Maximum extractable water by low vegetation from bare ground soil column [mm m^2 m^-2 ground s^-1]
                 ExWaterGroundVeg_H: Maximum extractable water by tree from vegetated ground soil column [mm m^2 m^-2 ground s^-1]
                 ExWaterGroundVeg_L: Maximum extractable water by low vegetation from vegetated ground soil column [mm m^2 m^-2 ground s^-1]
                 ExWaterGroundTot_H: Total maximum extractable water by tree from ground soil column [mm m^2 m^-2 ground s^-1]
                 ExWaterGroundTot_L: Total maximum extractable water by low vegetation from ground soil column [mm m^2 m^-2 ground s^-1]
        SoilPotW: Soil water potential
                  SoilPotWGroundImp_H: Soil water potential for first layer of vegetation in the impervious ground soil column [MPa]
                  SoilPotWGroundBare_H: Soil water potential for first layer of vegetation in the bare ground soil column [MPa]
                  SoilPotWGroundVeg_H: Soil water potential for first layer of vegetation in the vegetated ground soil column [MPa]
                  SoilPotWGroundTot_H: Total Soil water potential for first layer of vegetation in the soil column [MPa]
                  SoilPotWGroundImp_L: Soil water potential for second layer of vegetation in the impervious ground soil column [MPa]
                  SoilPotWGroundBare_L: Soil water potential for second layer of vegetation in the bare ground soil column [MPa]
                  SoilPotWGroundVeg_L: Soil water potential for second layer of vegetation in the vegetated ground soil column [MPa]
                  SoilPotWGroundTot_L: Total Soil water potential for second layer of vegetation in the soil column [MPa]
        dVwater_dt: Total water change in each soil column
                    dVGroundSoilImp_dt: Total water change in soil column of impervious ground [mm s^-1]
                    dVGroundSoilBare_dt: Total water change in soil column of bare ground [mm s^-1]
                    dVGroundSoilVeg_dt: Total water change in soil column of vegetated ground [mm s^-1]
                    dVGroundSoilTot_dt: Total water change in soil columns [mm s^-1]
        WBIndv: Water balance equation
                WB_In_tree: Water balance for trees [mm s^-1]
                WB_In_gveg: Water balance for low vegetation [mm s^-1]
                WB_In_gimp: Water balance at the impervious surface [mm s^-1]
                WB_In_gbare: Water balance at the bare surface [mm s^-1]
                WB_Pond_gveg: Water balance at the surface under vegetation [mm s^-1]
                WB_Soil_gimp: Water balance for the impervious ground soil column [mm s^-1]
                WB_Soil_gbare: Water balance for the bare ground soil column [mm s^-1]
                WB_Soil_gveg: Water balance for the vegetated ground soil column [mm s^-1]
        WBTot: Water balance at different surfaces
               WBSurf_tree: Water balance for tree [mm s^-1]
               WBSurf_imp: Water balance at the impervious surface [mm s^-1]
               WBSurf_bare: Water balance at the bare surface [mm s^-1]
               WBSurf_veg: Water balance at the vegetated surface [mm s^-1]
               WBSoil_imp: Water balance for the impervious ground soil column [mm s^-1]
               WBSoil_bare: Water balance for the bare ground soil column [mm s^-1]
               WBSoil_veg: Water balance for the vegetated ground soil column [mm s^-1]
               WBImp_tot: Total water balance for the impervious ground (surface and soil column) [mm s^-1]
               WBBare_tot: Total water balance for the bare ground (surface and soil column) [mm s^-1]
               WBVeg_tot: Total water balance for the vegetated ground (surface and soil column) [mm s^-1]
               WBCanyon_flux: Water balance for the canyon [mm s^-1]
               WBTree_level: Water balance at the tree level [mm s^-1]
               WBGround_level: Water balance at the ground level (impervious, bare and vegetated ground) [mm s^-1]
               WBSoil_level: Water balance at the soil level (impervious, bare and vegetated column) [mm s^-1]
               WBCanyon_level: Water balance at the canyon level [mm s^-1]
        Rd: Water table rise
            Rd_gbare: Water table rise that reaches the bare surface level (Dunne Runoff) [mm]
            Rd_gveg: Water table rise that reaches the vegetated surface level (Dunne Runoff) [mm]
            Rd_gimp: Water table rise that reaches the impervious surface level (Dunne Runoff) [mm]
        TE:
            TEgveg_imp: Transpiration from low vegetation in the impervious ground column [kg m^-2 s^-1]
            TEtree_imp: Transpiration from tree in the impervious ground column [kg m^-2 s^-1]
            TEgveg_bare: Transpiration from low vegetation in the bare ground column [kg m^-2 s^-1]
            TEtree_bare: Transpiration from tree in the bare ground column [kg m^-2 s^-1]
            TEgveg_veg: Transpiration from low vegetation in the vegetated ground column [kg m^-2 s^-1]
            TEtree_veg: Transpiration from tree in the vegetated ground column [kg m^-2 s^-1]
        RainGround: Water received by any ground fraction including rain and dripping from trees [mm s^-1]
        Anth_gbare: Anthropogenic water at the bare surface [mm s^-1]
        Anth_gveg: Anthropogenic water at the vegetated surface [mm s^-1]
        """

        Water_Cal = Water_Calculations()

        # Re-define input parameters which are overwritten in this function
        Egbare_soil_local = copy.copy(Egbare_soil)
        Egveg_soil_local = copy.copy(Egveg_soil)
        TEgveg_local = copy.copy(TEgveg)

        # Boolean operator for presence and absence of trees, vegetated ground, bare ground, and impervious ground
        Ctree = int(ParTree.trees == 1)
        Cveg = int(FractionsGround.fveg > 0)
        Cbare = int(FractionsGround.fbare > 0)
        Cimp = int(FractionsGround.fimp > 0)

        #------------------------
        # Vegetation interception
        #------------------------
        # Trees water intercepted
        Water_Cal.Water_Vegetation(MeteoData.Rain,Etree_In,self.Int.IntTree,ParInterceptionTree.Sp_In,ParVegTree.LAI,
                                   ParVegTree.SAI,ParCalculation.rhow,ParCalculation.dts)
        q_tree_dwn = Water_Cal.WaterVegetation.q_runoff_veg * Ctree     # [mm s^-1]
        In_tree = Water_Cal.WaterVegetation.In_veg * Ctree              # [mm]
        dIn_tree_dt = Water_Cal.WaterVegetation.dIn_veg_dt * Ctree      # [mm s^-1]
        WB_In_tree = Water_Cal.WaterVegetation.WBalance_In_veg * Ctree  # [mm s^-1]

        # Water received by any ground fraction including rain and dripping from trees [mm s^-1]
        Rain_ground = 4*geometry.radius_tree*Ctree*q_tree_dwn + (1-4*geometry.radius_tree*Ctree)*MeteoData.Rain

        # Vegetated ground water intercepted
        Water_Cal.Water_Vegetation(Rain_ground,Egveg_In,self.Int.IntGroundVegPlant,ParSoilGround.Sp_In,ParVegGround.LAI,
                                   ParVegGround.SAI,ParCalculation.rhow,ParCalculation.dts)
        q_gveg_dwn = Water_Cal.WaterVegetation.q_runoff_veg * Cveg    # [mm s^-1]
        In_gveg = Water_Cal.WaterVegetation.In_veg * Cveg             # [mm]
        dIn_gveg_dt = Water_Cal.WaterVegetation.dIn_veg_dt * Cveg     # [mm s^-1]
        WB_In_gveg = Water_Cal.WaterVegetation.WBalance_In_veg * Cveg # [mm s^-1]

        # -----------------------
        # Water runon (ponding) on ground
        # -----------------------
        # Impervious surface
        Water_Cal.Water_Impervious(Rain_ground,self.Runon.RunonGroundTot,Egimp_Pond,self.Int.IntGroundImp,ParCalculation.dts,
                                   ParCalculation.rhow,ParSoilGround.In_max_imp,ParSoilGround.Kimp)
        q_gimp_runoff = Water_Cal.WaterImpervious.q_runoff_imp * Cimp  # [mm s^-1]
        In_gimp = Water_Cal.WaterImpervious.In_imp * Cimp              # [mm]
        dIn_gimp_dt = Water_Cal.WaterImpervious.dIn_imp_dt * Cimp      # [mm s^-1]
        f_inf_gimp = Water_Cal.WaterImpervious.Lk_imp * Cimp           # [mm s^-1]
        WB_In_gimp = Water_Cal.WaterImpervious.WBalance_In_imp * Cimp  # [mm s^-1]

        # Bare ground
        Water_Cal.Water_Ground(Rain_ground+Anthropogenic.Waterf_canyonBare/3600,self.Runon.RunonGroundTot,Egbare_Pond,
                               self.Owater.OwGroundSoilBare,self.Int.IntGroundBare,ParSoilGround.In_max_bare,ParSoilGround.Pcla,
                               ParSoilGround.Psan,ParSoilGround.Porg,ParSoilGround.Kfc,ParSoilGround.Phy,ParSoilGround.SPAR,
                               ParSoilGround.Kbot,ParVegTree.CASE_ROOT,ParVegGround.CASE_ROOT,ParVegTree.ZR95,ParVegGround.ZR95*numpy.NaN,
                               ParVegTree.ZR50,ParVegGround.ZR50*numpy.NaN,ParVegTree.ZRmax,ParVegGround.ZRmax*numpy.NaN,
                               ParSoilGround.Zs,ParCalculation.dts,ParCalculation.rhow)
        q_gbare_runoff = Water_Cal.WaterGround.q_runoff_ground * Cbare  # [mm s^-1]
        In_gbare = Water_Cal.WaterGround.In_ground * Cbare              # [mm]
        dIn_gbare_dt = Water_Cal.WaterGround.dIn_ground_dt * Cbare      # [mm s^-1]
        f_inf_gbare = Water_Cal.WaterGround.f_ground * Cbare            # [mm s^-1]
        WB_In_gbare = Water_Cal.WaterGround.WBalance_In_ground * Cbare  # [mm s^-1]

        # Ground under vegetation
        Water_Cal.Water_Ground(q_gveg_dwn+Anthropogenic.Waterf_canyonVeg/3600,self.Runon.RunonGroundTot,Egveg_Pond,
                               self.Owater.OwGroundSoilVeg,self.Int.IntGroundVegGround,ParSoilGround.In_max_underveg,
                               ParSoilGround.Pcla,ParSoilGround.Psan,ParSoilGround.Porg,ParSoilGround.Kfc,ParSoilGround.Phy,
                               ParSoilGround.SPAR,ParSoilGround.Kbot,ParVegTree.CASE_ROOT,ParVegGround.CASE_ROOT,ParVegTree.ZR95,
                               ParVegGround.ZR95,ParVegTree.ZR50,ParVegGround.ZR50,ParVegTree.ZRmax,ParVegGround.ZRmax,
                               ParSoilGround.Zs,ParCalculation.dts,ParCalculation.rhow)
        q_gveg_runoff = Water_Cal.WaterGround.q_runoff_ground * Cveg    # [mm s^-1]
        In_gveg_pond = Water_Cal.WaterGround.In_ground * Cveg           # [mm]
        dIn_gveg_pond_dt = Water_Cal.WaterGround.dIn_ground_dt * Cveg   # [mm s^-1]
        f_inf_gveg = Water_Cal.WaterGround.f_ground * Cveg              # [mm s^-1]
        WB_Pond_gveg = Water_Cal.WaterGround.WBalance_In_ground * Cveg  # [mm s^-1]

        #--------------------------------------------------------------------------
        # Water distribution in different soil columns (impervious, bare, vegetated)
        #--------------------------------------------------------------------------
        # Transpiration assumptions: the tree can extract water from all soil layers.
        # The low vegetation can only take water from the vegetated soil layer
        # TEtree: Tree evaporation per m^2 tree crown area

        # Tree roots can access all water in the soil (imp, bare, veg)
        if ParVegTree.SPARTREE == 1:
            TEtree0_imp = TEtree * Ctree * (4*geometry.radius_tree) * Cimp   # [kg m^-2 s^-1]
            TEtree0_bare = TEtree * Ctree * (4*geometry.radius_tree) * Cbare # [kg m^-2 s^-1]
            TEtree0_veg = TEtree * Ctree * (4*geometry.radius_tree) * Cveg   # [kg m^-2 s^-1]

        # If the tree crown is smaller than the combined vegetated and bare fraction, then the trees only transpire
        # from these fractions. Otherwise, they also transpire from the impervious ground fraction.
        elif ParVegTree.SPARTREE == 2:
            if (4*geometry.radius_tree) <= (FractionsGround.fveg+FractionsGround.fbare):
                TEtree0_imp = 0 * Ctree * Cimp   # [kg m^-2 s^-1]
                TEtree0_bare = TEtree * (4*geometry.radius_tree) / (FractionsGround.fveg + FractionsGround.fbare) * Ctree * Cbare # [kg m^-2 s^-1]
                TEtree0_veg = TEtree * (4*geometry.radius_tree) / (FractionsGround.fveg + FractionsGround.fbare) * Ctree * Cveg   # [kg m^-2 s^-1]
            elif (4*geometry.radius_tree) > (FractionsGround.fveg+FractionsGround.fbare):
                TEtree0_imp = ((4*geometry.radius_tree) - (FractionsGround.fveg+FractionsGround.fbare))*TEtree/FractionsGround.fimp *Ctree*Cimp # [kg m^-2 s^-1]
                TEtree0_bare = TEtree * Ctree * Cbare # [kg m^-2 s^-1]
                TEtree0_veg = TEtree * Ctree * Cveg   # [kg m^-2 s^-1]

        TEgveg_local = TEgveg_local * Cveg             # [kg m^-2 s^-1]
        Egbare_soil_local = Egbare_soil_local * Cbare  # [kg m^-2 s^-1]
        Egveg_soil_local = Egveg_soil_local * Cveg     # [kg m^-2 s^-1]

        # Impervious column
        Water_Cal.Water_Soil(self.Owater.OwGroundSoilImp[2:],f_inf_gimp,TEtree0_imp,0,0,[0 for i in range(len(self.Qinlat.Qin_imp[2:]))],
                             ParCalculation.dts,ParSoilGround.Pcla,ParSoilGround.Psan,ParSoilGround.Porg,ParSoilGround.Kfc,
                             ParSoilGround.Phy,ParSoilGround.SPAR,ParSoilGround.Kbot,ParVegTree.CASE_ROOT,ParVegGround.CASE_ROOT,
                             ParVegTree.ZR95,ParVegGround.ZR95,ParVegTree.ZR50,ParVegGround.ZR50,ParVegTree.ZRmax,ParVegGround.ZRmax,
                             ParVegTree.Rrootl,ParVegGround.Rrootl,ParVegTree.PsiL50,ParVegGround.PsiL50,ParVegTree.PsiX50,
                             ParVegGround.PsiX50,[ParSoilGround.Zs[i]-ParSoilGround.Zs[2] for i in range(2,len(ParSoilGround.Zs))],
                             ParCalculation.rhow)
        V_gimp1 = numpy.zeros(len(Water_Cal.WaterSoil.V)+2)
        V_gimp1[0:2] = [numpy.NaN,numpy.NaN]
        V_gimp1[2:] = Water_Cal.WaterSoil.V               # [mm]
        Lk_gimp1 = Water_Cal.WaterSoil.Lk                 # [mm s^-1]
        self.Rd.Rd_gimp = Water_Cal.WaterSoil.Rd          # [mm]
        self.TE.TEgveg_imp = Water_Cal.WaterSoil.TE_L     # [kg m^-2 s^-1]
        self.TE.TEtree_imp = Water_Cal.WaterSoil.TE_H     # [kg m^-2 s^-1]
        self.Egimp_soil1 = Water_Cal.WaterSoil.E_soil     # [kg m^-2 s^-1]
        WB_Soil_gimp1 = Water_Cal.WaterSoil.WBalance_soil # [mm s^-1]

        # Bare ground column
        Water_Cal.Water_Soil(self.Owater.OwGroundSoilBare,f_inf_gbare,TEtree0_bare,0,Egbare_soil_local,
                             [0 for i in range(len(self.Qinlat.Qin_bare))],ParCalculation.dts,ParSoilGround.Pcla,
                             ParSoilGround.Psan,ParSoilGround.Porg,ParSoilGround.Kfc,ParSoilGround.Phy,ParSoilGround.SPAR,
                             ParSoilGround.Kbot,ParVegTree.CASE_ROOT,ParVegGround.CASE_ROOT,ParVegTree.ZR95,ParVegGround.ZR95,
                             ParVegTree.ZR50,ParVegGround.ZR50,ParVegTree.ZRmax,ParVegGround.ZRmax,ParVegTree.Rrootl,
                             ParVegGround.Rrootl,ParVegTree.PsiL50,ParVegGround.PsiL50,ParVegTree.PsiX50,ParVegGround.PsiX50,
                             ParSoilGround.Zs,ParCalculation.rhow)
        V_gbare1 = Water_Cal.WaterSoil.V                   # [mm]
        Lk_gbare1 = Water_Cal.WaterSoil.Lk                 # [mm s^-1]
        self.Rd.Rd_gbare = Water_Cal.WaterSoil.Rd          # [mm]
        self.TE.TEgveg_bare = Water_Cal.WaterSoil.TE_L     # [kg m^-2 s^-1]
        self.TE.TEtree_bare = Water_Cal.WaterSoil.TE_H     # [kg m^-2 s^-1]
        self.Egbare_Soil1 = Water_Cal.WaterSoil.E_soil     # [kg m^-2 s^-1]
        WB_Soil_gbare1 = Water_Cal.WaterSoil.WBalance_soil # [mm s^-1]

        # Vegetated ground column
        Water_Cal.Water_Soil(self.Owater.OwGroundSoilVeg,f_inf_gveg,TEtree0_veg,TEgveg_local,Egveg_soil_local,
                             [0 for i in range(len(self.Qinlat.Qin_veg))],ParCalculation.dts,ParSoilGround.Pcla,
                             ParSoilGround.Psan,ParSoilGround.Porg,ParSoilGround.Kfc,ParSoilGround.Phy,ParSoilGround.SPAR,
                             ParSoilGround.Kbot,ParVegTree.CASE_ROOT,ParVegGround.CASE_ROOT,ParVegTree.ZR95,ParVegGround.ZR95,
                             ParVegTree.ZR50,ParVegGround.ZR50,ParVegTree.ZRmax,ParVegGround.ZRmax,ParVegTree.Rrootl,
                             ParVegGround.Rrootl,ParVegTree.PsiL50,ParVegGround.PsiL50,ParVegTree.PsiX50,ParVegGround.PsiX50,
                             ParSoilGround.Zs,ParCalculation.rhow)
        V_gveg1 = Water_Cal.WaterSoil.V                    # [mm]
        Lk_gveg1 = Water_Cal.WaterSoil.Lk                  # [mm s^-1]
        self.Rd.Rd_gveg = Water_Cal.WaterSoil.Rd           # [mm]
        self.TE.TEgveg_veg = Water_Cal.WaterSoil.TE_L      # [kg m^-2 s^-1]
        self.TE.TEtree_veg = Water_Cal.WaterSoil.TE_H      # [kg m^-2 s^-1]
        self.Egveg_Soil1 = Water_Cal.WaterSoil.E_soil      # [kg m^-2 s^-1]
        WB_Soil_gveg1 = Water_Cal.WaterSoil.WBalance_soil  # [mm s^-1]

        # Lateral Richards equation
        T_SPAN = [0,ParCalculation.dts]
        # Thickness of the layers [mm]
        dz = numpy.diff(ParSoilGround.Zs)

        SoilCal = Soil_Calculations()
        SoilCal.Soil_Parameters_Total(ParSoilGround.Pcla,ParSoilGround.Psan,ParSoilGround.Porg,ParSoilGround.Kfc,ParSoilGround.Phy,
                                      ParSoilGround.SPAR,ParSoilGround.Kbot,ParVegTree.CASE_ROOT,ParVegGround.CASE_ROOT,
                                      ParVegTree.ZR95,ParVegGround.ZR95,ParVegTree.ZR50,ParVegGround.ZR50,ParVegTree.ZRmax,
                                      ParVegGround.ZRmax,ParSoilGround.Zs)
        Osat = SoilCal.SoilParamTotal.Osat
        Ohy = SoilCal.SoilParamTotal.Ohy
        nVG = SoilCal.SoilParamTotal.nVG
        alpVG = SoilCal.SoilParamTotal.alpVG
        Ks_Zs = SoilCal.SoilParamTotal.Ks_Zs
        L = SoilCal.SoilParamTotal.L
        Pe = SoilCal.SoilParamTotal.Pe
        O33 = SoilCal.SoilParamTotal.O33

        V_gimp2 = numpy.zeros(len(dz))
        V_gbare2 = numpy.zeros(len(dz))
        V_gveg2 = numpy.zeros(len(dz))

        # Solve the ordinary differential equation (ODE) for each soil layer, describing the change in soil moisture over time
        # The first two soil layers are only calculated as the exchange between two soil columns as the impervious soil
        # column is assumed to be impervious there
        for i in range(0,2):
            Vlat1 = [V_gbare1[i],V_gveg1[i]]

            Vout = odeint(SoilCal.Soil_Moistures_Rich_Comp_Lat2, Vlat1, T_SPAN,
                          args=(dz[i],ParSoilGround.SPAR,Ks_Zs[i],Osat[i],Ohy[i],L[i],Pe[i],O33[i],alpVG[i],nVG[i],Cbare,
                                Cveg,FractionsGround.fbare,FractionsGround.fveg,Gemeotry_m.Width_canyon))
            V_gimp2[i] = numpy.NaN   # [mm]
            V_gbare2[i] = Vout[-1,0] # [mm]
            V_gveg2[i] = Vout[-1,1]  # [mm]

        # Solve the ODE for third layer to the end
        for i in range(2,len(dz)):
            Vlat1 = [V_gimp1[i], V_gbare1[i], V_gveg1[i]]

            Vout = odeint(SoilCal.Soil_Moistures_Rich_Comp_Lat3, Vlat1, T_SPAN,
                          args=(dz[i],ParSoilGround.SPAR,Ks_Zs[i],Osat[i],Ohy[i],L[i],Pe[i],O33[i],alpVG[i],nVG[i],Cimp,
                                Cbare,Cveg,FractionsGround.fimp,FractionsGround.fbare,FractionsGround.fveg,Gemeotry_m.Width_canyon))
            V_gimp2[i] = Vout[-1, 0]  # [mm]
            V_gbare2[i] = Vout[-1, 1] # [mm]
            V_gveg2[i] = Vout[-1, 2]  # [mm]


        # Back compute lateral differentials
        dVin_imp = [V_gimp2[i] - V_gimp1[i] for i in range(0,len(V_gimp2))]      # [mm]
        dVin_bare = [V_gbare2[i] - V_gbare1[i] for i in range(0,len(V_gbare2))]  # [mm]
        dVin_veg = [V_gveg2[i] - V_gveg1[i] for i in range(0,len(V_gveg2))]      # [mm]

        Qin_bare2imp = numpy.zeros(len(dz))
        Qin_bare2veg = numpy.zeros(len(dz))
        Qin_imp2bare = numpy.zeros(len(dz))
        Qin_imp2veg = numpy.zeros(len(dz))
        Qin_veg2imp = numpy.zeros(len(dz))
        Qin_veg2bare = numpy.zeros(len(dz))
        Qin_bare2imp[0:2] = [numpy.NaN,numpy.NaN]
        Qin_imp2bare[0:2] = [numpy.NaN,numpy.NaN]
        Qin_imp2veg[0:2] = [numpy.NaN,numpy.NaN]
        Qin_veg2imp[0:2] = [numpy.NaN,numpy.NaN]

        # Impervious surface
        SoilCal.Soil_Moisture_Conductivity_Update(V_gimp2[2:],ParSoilGround.Pcla,ParSoilGround.Psan,ParSoilGround.Porg,
                                                  ParSoilGround.Kfc,ParSoilGround.Phy,ParSoilGround.SPAR,ParSoilGround.Kbot,
                                                  ParVegTree.CASE_ROOT,ParVegGround.CASE_ROOT,ParVegTree.ZR95,ParVegGround.ZR95,
                                                  ParVegTree.ZR50,ParVegGround.ZR50,ParVegTree.ZRmax,ParVegGround.ZRmax,
                                                  [ParSoilGround.Zs[i]-ParSoilGround.Zs[2] for i in range(2,len(ParSoilGround.Zs))],
                                                  ParVegTree.Rrootl,ParVegGround.Rrootl,ParVegTree.PsiL50,ParVegGround.PsiL50,
                                                  ParVegTree.PsiX50,ParVegGround.PsiX50)
        V_gimp2 = SoilCal.SoilMoistCond.V                # [mm]
        O_gimp2 = SoilCal.SoilMoistCond.O                # [-]
        OS_gimp2 = SoilCal.SoilMoistCond.OS              # [-]
        Psi_Soil_gimp2 = SoilCal.SoilMoistCond.Psi_soil  # [mm]
        Psi_s_H_gimp2 = SoilCal.SoilMoistCond.Psi_s_H    # [MPa]
        Psi_s_L_gimp2 = SoilCal.SoilMoistCond.Psi_s_L    # [MPa]
        Exwat_H_gimp2 = SoilCal.SoilMoistCond.Exwat_H    # [mm m^2 m^-2 ground s^-1]
        Exwat_L_gimp2 = SoilCal.SoilMoistCond.Exwat_L    # [mm m^2 m^-2 ground s^-1]
        Kf_gimp2 = SoilCal.SoilMoistCond.Ko              # [mm s^-1]

        V_gimp2_bckup = V_gimp2
        V_gimp2 = numpy.zeros(len(V_gimp2_bckup)+2)
        V_gimp2[0:2] = [numpy.NaN,numpy.NaN]
        V_gimp2[2:] = V_gimp2_bckup

        O_gimp2_bckup = O_gimp2
        O_gimp2 = numpy.zeros(len(O_gimp2_bckup)+2)
        O_gimp2[0:2] = [numpy.NaN,numpy.NaN]
        O_gimp2[2:] = O_gimp2_bckup

        Exwat_H_gimp2_bckup = Exwat_H_gimp2
        Exwat_H_gimp2 = numpy.zeros(numpy.size(Exwat_H_gimp2_bckup)+2)
        Exwat_H_gimp2[0:2] = [numpy.NaN, numpy.NaN]
        Exwat_H_gimp2[2:] = Exwat_H_gimp2_bckup

        Exwat_L_gimp2_bckup = Exwat_L_gimp2
        Exwat_L_gimp2 = numpy.zeros(numpy.size(Exwat_L_gimp2)+2)
        Exwat_L_gimp2[0:2] = [numpy.NaN, numpy.NaN]
        Exwat_L_gimp2[2:] = Exwat_L_gimp2_bckup

        Psi_Soil_gimp2_bckup = Psi_Soil_gimp2
        Psi_Soil_gimp2 = numpy.zeros(len(Psi_Soil_gimp2)+2)
        Psi_Soil_gimp2[0:2] = [numpy.NaN, numpy.NaN]
        Psi_Soil_gimp2[2:] = Psi_Soil_gimp2_bckup

        Kf_gimp2_bckup = Kf_gimp2
        Kf_gimp2 = numpy.zeros(len(Kf_gimp2)+2)
        Kf_gimp2[0:2] = [numpy.NaN, numpy.NaN]
        Kf_gimp2[2:] = Kf_gimp2_bckup

        # Bare ground
        SoilCal.Soil_Moisture_Conductivity_Update(V_gbare2,ParSoilGround.Pcla,ParSoilGround.Psan,ParSoilGround.Porg,
                                                  ParSoilGround.Kfc,ParSoilGround.Phy,ParSoilGround.SPAR,ParSoilGround.Kbot,
                                                  ParVegTree.CASE_ROOT,ParVegGround.CASE_ROOT,ParVegTree.ZR95,ParVegGround.ZR95,
                                                  ParVegTree.ZR50,ParVegGround.ZR50,ParVegTree.ZRmax,ParVegGround.ZRmax,
                                                  ParSoilGround.Zs,ParVegTree.Rrootl,ParVegGround.Rrootl,ParVegTree.PsiL50,
                                                  ParVegGround.PsiL50,ParVegTree.PsiX50,ParVegGround.PsiX50)
        V_gbare2 = SoilCal.SoilMoistCond.V               # [mm]
        O_gbare2 = SoilCal.SoilMoistCond.O               # [-]
        OS_gbare2 = SoilCal.SoilMoistCond.OS             # [-]
        Psi_soil_gbare2 = SoilCal.SoilMoistCond.Psi_soil # [mm]
        Psi_s_H_gbare2 = SoilCal.SoilMoistCond.Psi_s_H   # [MPa]
        Psi_s_L_gbare2 = SoilCal.SoilMoistCond.Psi_s_L   # [MPa]
        Exwat_H_gbare2 = SoilCal.SoilMoistCond.Exwat_H   # [mm m^2 m^-2 ground s^-1]
        Exwat_L_gbare2 = SoilCal.SoilMoistCond.Exwat_L   # [mm m^2 m^-2 ground s^-1]
        Kf_gbare2 = SoilCal.SoilMoistCond.Ko             # [mm s^-1]

        # Vegetated ground
        SoilCal.Soil_Moisture_Conductivity_Update(V_gveg2,ParSoilGround.Pcla,ParSoilGround.Psan,ParSoilGround.Porg,
                                                  ParSoilGround.Kfc,ParSoilGround.Phy,ParSoilGround.SPAR,ParSoilGround.Kbot,
                                                  ParVegTree.CASE_ROOT,ParVegGround.CASE_ROOT,ParVegTree.ZR95,ParVegGround.ZR95,
                                                  ParVegTree.ZR50,ParVegGround.ZR50,ParVegTree.ZRmax,ParVegGround.ZRmax,
                                                  ParSoilGround.Zs,ParVegTree.Rrootl,ParVegGround.Rrootl,ParVegTree.PsiL50,
                                                  ParVegGround.PsiL50,ParVegTree.PsiX50,ParVegGround.PsiX50)
        V_gveg2 = SoilCal.SoilMoistCond.V                # [mm]
        O_gveg2 = SoilCal.SoilMoistCond.O                # [-]
        OS_gveg2 = SoilCal.SoilMoistCond.OS              # [-]
        Psi_soil_gveg2 = SoilCal.SoilMoistCond.Psi_soil  # [mm]
        Psi_s_H_gveg2 = SoilCal.SoilMoistCond.Psi_s_H    # [MPa]
        Psi_s_L_gveg2 = SoilCal.SoilMoistCond.Psi_s_L    # [MPa]
        Exwat_H_gveg2 = SoilCal.SoilMoistCond.Exwat_H    # [mm m^2 m^-2 ground s^-1]
        Exwat_L_gveg2 = SoilCal.SoilMoistCond.Exwat_L    # [mm m^2 m^-2 ground s^-1]
        Kf_gveg2 = SoilCal.SoilMoistCond.Ko              # [mm s^-1]

        # Change in water column
        # Water volume in each layer [mm]
        Vtm1_imp = [(self.Owater.OwGroundSoilImp[i] - Ohy[i]) * dz[i] for i in range(0,len(dz))]
        Vtm1_imp[0] = numpy.NaN
        Vtm1_imp[1] = numpy.NaN
        # Water volume in each layer [mm]
        Vtm1_bare = [(self.Owater.OwGroundSoilBare[i] - Ohy[i]) * dz[i] for i in range(0,len(dz))]
        # Water volume in each layer [mm]
        Vtm1_veg = [(self.Owater.OwGroundSoilVeg[i] - Ohy[i]) * dz[i] for i in range(0,len(dz))]

        # Calculate total change in water in each soil column
        dV_dt_gimpTot = (numpy.nansum(V_gimp2) - numpy.nansum(Vtm1_imp))/ParCalculation.dts      # [mm s^-1]
        dV_dt_gbareTot = (numpy.nansum(V_gbare2) - numpy.nansum(Vtm1_bare))/ParCalculation.dts   # [mm s^-1]
        dV_dt_gvegTot = (numpy.nansum(V_gveg2) - numpy.nansum(Vtm1_veg))/ParCalculation.dts      # [mm s^-1]

        #-----------------------------------------------------------------
        # Surface water balance fluxes: Tree, impervious, bare, vegetated ground
        #-----------------------------------------------------------------
        # Water balance flux [mm s^-1]: Tree
        WBsurf_tree = Ctree*MeteoData.Rain - Etree_In*1000/ParCalculation.rhow - q_tree_dwn - dIn_tree_dt
        # Water balance flux [mm s^-1]: Impervious surface
        WBsurf_imp = Cimp*Rain_ground + Cimp*self.Runon.RunonGroundTot - Egimp_Pond*1000/ParCalculation.rhow - q_gimp_runoff - f_inf_gimp - dIn_gimp_dt
        # Water balance flux [mm s^-1]: Bare ground
        WBsurf_bare = Cbare*Rain_ground + Cbare*self.Runon.RunonGroundTot + Anthropogenic.Waterf_canyonBare/3600 - Egbare_Pond*1000/ParCalculation.rhow - \
                      f_inf_gbare - q_gbare_runoff - dIn_gbare_dt
        # Water balance flux [mm s^-1]: Vegetated ground (including both low vegetation and ground underneath the vegetation)
        WBsurf_veg = Cveg*Rain_ground + Cveg*self.Runon.RunonGroundTot + Anthropogenic.Waterf_canyonVeg/3600 - (Egveg_In+Egveg_Pond)*1000/ParCalculation.rhow - \
                     f_inf_gveg - q_gveg_runoff - dIn_gveg_dt - dIn_gveg_pond_dt

        #-------------------------------------------------------
        # Soil water balance fluxes: impervious, bare, vegetated ground
        #-------------------------------------------------------
        # Water balance flux [mm s^-1]: Soil column under impervious surface
        WBsoil_imp = f_inf_gimp + numpy.nansum(dVin_imp)/ParCalculation.dts - (self.TE.TEgveg_imp+self.TE.TEtree_imp+self.Egimp_soil1)*1000/ParCalculation.rhow -\
                     Lk_gimp1 - self.Rd.Rd_gimp/ParCalculation.dts - dV_dt_gimpTot
        # Water balance flux [mm s^-1]: Soil column under bare ground
        WBsoil_bare = f_inf_gbare + numpy.nansum(dVin_bare)/ParCalculation.dts - \
                      (self.TE.TEgveg_bare+self.TE.TEtree_bare+self.Egbare_Soil1)*1000/ParCalculation.rhow - \
                      Lk_gbare1 - self.Rd.Rd_gbare/ParCalculation.dts - dV_dt_gbareTot
        # Water balance flux [mm s^-1]: Soil column under vegetated ground
        WBsoil_veg = f_inf_gveg + numpy.nansum(dVin_veg)/ParCalculation.dts - (self.TE.TEgveg_veg+self.TE.TEtree_veg+self.Egveg_Soil1)*1000/ParCalculation.rhow - \
                     Lk_gveg1 - self.Rd.Rd_gveg/ParCalculation.dts - dV_dt_gvegTot

        #---------------------------------------------------------
        # Ground water balance fluxes: impervious, bare, vegetated ground
        #---------------------------------------------------------
        # Water balance flux [mm s^-1]: Impervious ground including surface and soil column
        WBimp_tot = Cimp*Rain_ground + Cimp*self.Runon.RunonGroundTot + numpy.nansum(dVin_imp)/ParCalculation.dts - \
                    (Egimp_Pond+self.TE.TEgveg_imp+self.TE.TEtree_imp+self.Egimp_soil1)*1000/ParCalculation.rhow - \
                    q_gimp_runoff - dIn_gimp_dt - Lk_gimp1 - self.Rd.Rd_gimp/ParCalculation.dts - dV_dt_gimpTot
        # Water balance flux [mm s^-1]: Bare ground including surface and soil column
        WBbare_tot = Cbare*Rain_ground + Cbare*self.Runon.RunonGroundTot + numpy.nansum(dVin_bare)/ParCalculation.dts + \
                     Anthropogenic.Waterf_canyonBare/3600 - \
                     (Egbare_Pond+self.TE.TEgveg_bare+self.TE.TEtree_bare+self.Egbare_Soil1)*1000/ParCalculation.rhow - \
                     q_gbare_runoff - dIn_gbare_dt - Lk_gbare1 - self.Rd.Rd_gbare/ParCalculation.dts - dV_dt_gbareTot
        # Water balance flux [mm s^-1]: Vegetated ground including low vegetation, surface and soil column
        WBveg_tot = Cveg*Rain_ground + Cveg*self.Runon.RunonGroundTot + numpy.nansum(dVin_veg)/ParCalculation.dts + \
                    Anthropogenic.Waterf_canyonVeg/3600 - \
                    (Egveg_In+Egveg_Pond+self.TE.TEgveg_veg+self.TE.TEtree_veg+self.Egveg_Soil1)*1000/ParCalculation.rhow - \
                    q_gveg_runoff - dIn_gveg_dt - dIn_gveg_pond_dt - Lk_gveg1 - self.Rd.Rd_gveg/ParCalculation.dts - dV_dt_gvegTot

        # Calculate total runoff [mm s^-1]
        Runoff = FractionsGround.Per_runoff * (FractionsGround.fimp*(q_gimp_runoff + self.Rd.Rd_gimp/ParCalculation.dts) +
                                               FractionsGround.fbare*(q_gbare_runoff+self.Rd.Rd_gbare/ParCalculation.dts) +
                                               FractionsGround.fveg*(q_gveg_runoff+self.Rd.Rd_gveg/ParCalculation.dts))
        # Calculate total runon [mm s^-1]
        Runon_local = (1-FractionsGround.Per_runoff) * (FractionsGround.fimp*(q_gimp_runoff+self.Rd.Rd_gimp/ParCalculation.dts) +
                                                        FractionsGround.fbare*(q_gbare_runoff+self.Rd.Rd_gbare/ParCalculation.dts) +
                                                        FractionsGround.fveg*(q_gveg_runoff+self.Rd.Rd_gveg/ParCalculation.dts))
        # Calculate total water flux in the canyon (ground and trees)  [mm s^-1]
        Etot = (FractionsGround.fimp*(Egimp_Pond+self.TE.TEgveg_imp+self.TE.TEtree_imp+self.Egimp_soil1) +
                FractionsGround.fbare *(Egbare_Pond+self.Egbare_Soil1+self.TE.TEgveg_bare+self.TE.TEtree_bare) +
                FractionsGround.fveg*(Egveg_In+Egveg_Pond+self.Egveg_Soil1+self.TE.TEgveg_veg+self.TE.TEtree_veg) +
                Ctree*(4*geometry.radius_tree)*(Etree_In)) * 1000/ParCalculation.rhow
        # Calculate total deep leakage [mm s^-1]
        DeepGLk = (FractionsGround.fimp*Lk_gimp1 + FractionsGround.fbare*Lk_gbare1 + FractionsGround.fveg*Lk_gveg1)
        # Calculate total storage in the canyon [mm s^-1]
        StorageTot = FractionsGround.fimp*(dIn_gimp_dt+dV_dt_gimpTot) + FractionsGround.fbare*(dIn_gbare_dt+dV_dt_gbareTot) + \
                     FractionsGround.fveg*(dIn_gveg_dt+dIn_gveg_pond_dt+dV_dt_gvegTot) + (4*geometry.radius_tree)*dIn_tree_dt
        # Water balance for the canyon [mm s^-1]
        WBcanyon_flux = MeteoData.Rain + self.Runon.RunonGroundTot + FractionsGround.fveg*Anthropogenic.Waterf_canyonVeg/3600 + \
                        FractionsGround.fbare*Anthropogenic.Waterf_canyonBare/3600 - Etot - DeepGLk - Runoff - Runon_local - StorageTot

        #--------------------
        # Total water balance fluxes
        #--------------------
        # Water balance flux [mm s^-1]: Tree level
        WBtree_level = (4*geometry.radius_tree) * (Ctree*MeteoData.Rain - Etree_In*1000/ParCalculation.rhow - q_tree_dwn - dIn_tree_dt)
        # Water balance flux [mm s^-1]: Ground level including impervious, bare and vegetated ground
        WBground_level = FractionsGround.fimp*(Rain_ground + self.Runon.RunonGroundTot - dIn_gimp_dt - Egimp_Pond*1000/ParCalculation.rhow -
                                               q_gimp_runoff - self.Rd.Rd_gimp/ParCalculation.dts - f_inf_gimp) + \
                         FractionsGround.fbare*(Rain_ground + self.Runon.RunonGroundTot + Anthropogenic.Waterf_canyonBare/3600 -
                                                dIn_gbare_dt - Egbare_Pond*1000/ParCalculation.rhow - q_gbare_runoff -
                                                self.Rd.Rd_gbare/ParCalculation.dts - f_inf_gbare) + \
                         FractionsGround.fveg*(Rain_ground + self.Runon.RunonGroundTot + Anthropogenic.Waterf_canyonVeg/3600 -
                                               dIn_gveg_dt - dIn_gveg_pond_dt - (Egveg_In+Egveg_Pond)*1000/ParCalculation.rhow -
                                               q_gveg_runoff - self.Rd.Rd_gveg/ParCalculation.dts - f_inf_gveg)
        # Water balance flux [mm s^-1]: Soil level
        WBsoil_level = FractionsGround.fimp*(f_inf_gimp + (-self.TE.TEgveg_imp-self.TE.TEtree_imp-self.Egimp_soil1)*1000/ParCalculation.rhow -
                                             Lk_gimp1 - dV_dt_gimpTot) + \
                       FractionsGround.fbare*(f_inf_gbare + (-self.TE.TEgveg_bare-self.TE.TEtree_bare-self.Egbare_Soil1)*1000/ParCalculation.rhow -
                                              Lk_gbare1 - dV_dt_gbareTot) + \
                       FractionsGround.fveg*(f_inf_gveg + (- self.TE.TEgveg_veg-self.TE.TEtree_veg-self.Egveg_Soil1)*1000/ParCalculation.rhow -
                                             Lk_gveg1 - dV_dt_gvegTot)
        # Water balance flux [mm s^-1]: Canyon level
        WBcanyon_level = MeteoData.Rain + self.Runon.RunonGroundTot + (4*geometry.radius_tree)*(-Etree_In*1000/ParCalculation.rhow-dIn_tree_dt) + \
                         FractionsGround.fimp*(-dIn_gimp_dt-Egimp_Pond*1000/ParCalculation.rhow-q_gimp_runoff-self.Rd.Rd_gimp/ParCalculation.dts -
                                               (self.TE.TEgveg_imp + self.TE.TEtree_imp + self.Egimp_soil1)*1000/ParCalculation.rhow -
                                               Lk_gimp1 - dV_dt_gimpTot)+ \
                         FractionsGround.fbare*(Anthropogenic.Waterf_canyonBare/3600 - dIn_gbare_dt - Egbare_Pond*1000/ParCalculation.rhow -
                                                q_gbare_runoff - self.Rd.Rd_gbare -
                                                (self.TE.TEgveg_bare + self.TE.TEtree_bare + self.Egbare_Soil1)*1000/ParCalculation.rhow -
                                                Lk_gbare1 - dV_dt_gbareTot) + \
                         FractionsGround.fveg*(Anthropogenic.Waterf_canyonVeg/3600 - dIn_gveg_dt - dIn_gveg_pond_dt -
                                               (Egveg_In + Egveg_Pond)*1000/ParCalculation.rhow - q_gveg_runoff - self.Rd.Rd_gveg/ParCalculation.dts -
                                               (self.TE.TEgveg_veg + self.TE.TEtree_veg + self.Egveg_Soil1)*1000/ParCalculation.rhow -
                                               Lk_gveg1 - dV_dt_gvegTot)

        # Replace the nan values with finite numbers
        V_gimp = copy.copy(V_gimp2)
        V_gimp[numpy.isnan(V_gimp)] = 0

        O_gimp = copy.copy(O_gimp2)
        O_gimp[numpy.isnan(O_gimp)] = Ohy[numpy.isnan(O_gimp)]

        Exwat_H_gimp = copy.copy(Exwat_H_gimp2)
        Exwat_H_gimp[numpy.isnan(Exwat_H_gimp)] = 0

        Exwat_L_gimp = copy.copy(Exwat_L_gimp2)
        Exwat_L_gimp[numpy.isnan(Exwat_L_gimp)] = 0

        # Average properties and rescaling
        V = numpy.nansum([[FractionsGround.fimp*V_gimp[i],FractionsGround.fbare*V_gbare2[i],FractionsGround.fveg*V_gveg2[i]]
                          for i in range(0,len(V_gbare2))],axis=1)
        O = numpy.nansum([[FractionsGround.fimp*O_gimp[i],FractionsGround.fbare*O_gbare2[i],FractionsGround.fveg*O_gveg2[i]]
                          for i in range(0,len(O_gbare2))],axis=1)
        OS = numpy.nansum([[FractionsGround.fimp*OS_gimp2,FractionsGround.fbare*OS_gbare2,FractionsGround.fveg*OS_gveg2]],axis=1)
        Lk = numpy.nansum([[FractionsGround.fimp*Lk_gimp1,FractionsGround.fbare*Lk_gbare1,FractionsGround.fveg*Lk_gveg1]],axis=1)
        Rd = numpy.nansum([[FractionsGround.fimp*self.Rd.Rd_gimp,FractionsGround.fbare*self.Rd.Rd_gbare,FractionsGround.fveg*self.Rd.Rd_gveg]],axis=1)
        dV_dt = numpy.nansum([[FractionsGround.fimp*dV_dt_gimpTot,FractionsGround.fbare*dV_dt_gbareTot,FractionsGround.fveg*dV_dt_gvegTot]],axis=1)

        V[0] = V[0]/(FractionsGround.fbare+FractionsGround.fveg)
        V[1] = V[1]/(FractionsGround.fbare+FractionsGround.fveg)
        V[numpy.isnan(V)] = 0

        O[0] = O[0]/(FractionsGround.fbare+FractionsGround.fveg)
        O[1] = O[1]/(FractionsGround.fbare+FractionsGround.fveg)
        O[numpy.isnan(O)] = Ohy[numpy.isnan(O)]

        Exwat_H = numpy.nansum([[FractionsGround.fimp*Exwat_H_gimp[i],FractionsGround.fbare*Exwat_H_gbare2[i],FractionsGround.fveg*Exwat_H_gveg2[i]]
                                for i in range(0,numpy.size(Exwat_H_gveg2))], axis=1)

        if ParVegTree.SPARTREE == 1:
            # Tree roots can access all water in the soil (imp, bare, veg)
            self.TE.TEtree_imp = self.TE.TEtree_imp / (4 * geometry.radius_tree) * Cimp
            self.TE.TEtree_bare = self.TE.TEtree_bare / (4 * geometry.radius_tree) * Cbare
            self.TE.TEtree_veg = self.TE.TEtree_veg / (4 * geometry.radius_tree) * Cveg

        else:
            # If the tree crown is smaller than the combined vegetated and bare fraction, then only the trees transpire
            #  from these fractions. Otherwise, there is also transpiration from roots in the impervious ground fraction.
            if 4*geometry.radius_tree <= (FractionsGround.fveg+FractionsGround.fbare):
                self.TE.TEtree_imp = 0
                self.TE.TEtree_bare = self.TE.TEtree_bare / (4 * geometry.radius_tree) * Cbare
                self.TE.TEtree_veg = self.TE.TEtree_veg / (4 * geometry.radius_tree) * Cveg
            elif 4*geometry.radius_tree > (FractionsGround.fveg+FractionsGround.fbare):
                self.TE.TEtree_imp = self.TE.TEtree_imp / ((4*geometry.radius_tree) - (FractionsGround.fveg+FractionsGround.fbare)) \
                                     * FractionsGround.fimp*Cimp
                self.TE.TEtree_bare = self.TE.TEtree_bare * Cbare
                self.TE.TEtree_veg = self.TE.TEtree_veg * Cveg

        TEtree_tot = (FractionsGround.fimp*self.TE.TEtree_imp + FractionsGround.fbare*self.TE.TEtree_bare + FractionsGround.fveg*self.TE.TEtree_veg)

        self.Runoff.QGroundImp = q_gimp_runoff
        self.Runoff.QGroundBarePond = q_gbare_runoff
        self.Runoff.QGroundBareSoil = self.Rd.Rd_gbare
        self.Runoff.QTree = q_tree_dwn
        self.Runoff.QGroundVegDrip = q_gveg_dwn
        self.Runoff.QGroundVegPond = q_gveg_runoff
        self.Runoff.QGroundVegSoil = self.Rd.Rd_gveg
        self.Runoff.RunoffGroundTot = Runoff

        self.dInt_dt.dInt_dtGroundImp = dIn_gimp_dt
        self.dInt_dt.dInt_dtGroundBare = dIn_gbare_dt
        self.dInt_dt.dInt_dtGroundVegPlant = dIn_gveg_dt
        self.dInt_dt.dInt_dtGroundVegGround = dIn_gveg_pond_dt
        self.dInt_dt.dInt_dtTree = dIn_tree_dt

        self.Infiltration.fGroundBare = f_inf_gbare
        self.Infiltration.fGroundVeg = f_inf_gveg
        self.Infiltration.fGroundImp = f_inf_gimp

        self.Vwater.VGroundSoilImp = V_gimp
        self.Vwater.VGroundSoilBare = V_gbare2
        self.Vwater.VGroundSoilVeg = V_gveg2
        self.Vwater.VGroundSoilTot = V

        self.Leakage.LkGroundImp = Lk_gimp1
        self.Leakage.LkGroundBare = Lk_gbare1
        self.Leakage.LkGroundVeg = Lk_gveg1
        self.Leakage.LkGround = Lk

        self.ExWater.ExWaterGroundImp_H = Exwat_H_gimp
        self.ExWater.ExWaterGroundImp_L = Exwat_L_gimp
        self.ExWater.ExWaterGroundBare_H = Exwat_H_gbare2
        self.ExWater.ExWaterGroundBare_L = Exwat_L_gbare2
        self.ExWater.ExWaterGroundVeg_H = Exwat_H_gveg2
        self.ExWater.ExWaterGroundVeg_L = Exwat_L_gveg2
        self.ExWater.ExWaterGroundTot_H = Exwat_H
        self.ExWater.ExWaterGroundTot_L = Exwat_L_gveg2

        self.SoilPotW.SoilPotWGroundImp_H = Psi_s_H_gimp2
        self.SoilPotW.SoilPotWGroundBare_H = Psi_s_H_gbare2
        self.SoilPotW.SoilPotWGroundVeg_H = Psi_s_H_gveg2
        self.SoilPotW.SoilPotWGroundTot_H = Psi_s_H_gveg2
        self.SoilPotW.SoilPotWGroundImp_L = Psi_s_L_gimp2
        self.SoilPotW.SoilPotWGroundBare_L = Psi_s_L_gbare2
        self.SoilPotW.SoilPotWGroundVeg_L = Psi_s_L_gveg2
        self.SoilPotW.SoilPotWGroundTot_L = Psi_s_L_gveg2

        self.dVwater_dt.dVGroundSoilImp_dt = dV_dt_gimpTot
        self.dVwater_dt.dVGroundSoilBare_dt = dV_dt_gbareTot
        self.dVwater_dt.dVGroundSoilVeg_dt = dV_dt_gvegTot
        self.dVwater_dt.dVGroundSoilTot_dt = dV_dt

        self.WBIndv.WB_In_tree = WB_In_tree
        self.WBIndv.WB_In_gveg = WB_In_gveg
        self.WBIndv.WB_In_gimp = WB_In_gimp
        self.WBIndv.WB_In_gbare = WB_In_gbare
        self.WBIndv.WB_Pond_gveg = WB_Pond_gveg
        self.WBIndv.WB_Soil_gimp = WB_Soil_gimp1
        self.WBIndv.WB_Soil_gbare = WB_Soil_gbare1
        self.WBIndv.WB_Soil_gveg = WB_Soil_gveg1

        self.WBTot.WBSurf_tree = WBsurf_tree
        self.WBTot.WBSurf_imp = WBsurf_imp
        self.WBTot.WBSurf_bare = WBsurf_bare
        self.WBTot.WBSurf_veg = WBsurf_veg
        self.WBTot.WBSoil_imp = WBsoil_imp
        self.WBTot.WBSoil_bare = WBsoil_bare
        self.WBTot.WBSoil_veg = WBsoil_veg
        self.WBTot.WBImp_tot = WBimp_tot
        self.WBTot.WBBare_tot = WBbare_tot
        self.WBTot.WBVeg_tot = WBveg_tot
        self.WBTot.WBCanyon_flux = WBcanyon_flux
        self.WBTot.WBTree_level = WBtree_level
        self.WBTot.WBGround_level = WBground_level
        self.WBTot.WBSoil_level = WBsoil_level
        self.WBTot.WBCanyon_level = WBcanyon_level

        self.Runon.RunonGroundTot = Runon_local

        self.Owater.OwGroundSoilImp = O_gimp
        self.Owater.OwGroundSoilBare = O_gbare2
        self.Owater.OwGroundSoilVeg = O_gveg2
        self.Owater.OwGroundSoilTot = O

        self.Int.IntGroundImp = In_gimp
        self.Int.IntGroundBare = In_gbare
        self.Int.IntGroundVegPlant = In_gveg
        self.Int.IntGroundVegGround = In_gveg_pond
        self.Int.IntTree = In_tree

        self.Qinlat.Qin_imp = dVin_imp
        self.Qinlat.Qin_bare = dVin_bare
        self.Qinlat.Qin_veg = dVin_veg
        self.Qinlat.Qin_bare2imp = Qin_bare2imp
        self.Qinlat.Qin_veg2imp = Qin_veg2imp
        self.Qinlat.Qin_veg2bare = Qin_veg2bare
        self.Qinlat.Qin_imp2bare = Qin_imp2bare
        self.Qinlat.Qin_bare2veg = Qin_bare2veg
        self.Qinlat.Qin_imp2veg = Qin_imp2veg

        self.RainGround = Rain_ground

        self.EfluxCanyon = Etot

        self.WaterStorageCanyon = StorageTot

        self.Anth_gbare = Anthropogenic.Waterf_canyonBare/3600
        self.Anth_gveg = Anthropogenic.Waterf_canyonVeg/3600