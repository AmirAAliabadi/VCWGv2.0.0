import os
import numpy
import math
from pprint import pprint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from Soil_Functions import Soil_Calculations
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from WaterFunctionsDef import WaterCanyonDef,WaterVegetationDef,WaterImperviousDef,WaterGroundDef,WaterSoilDef
import copy

'''
Water Functions:
Developed by Mohsen Moradi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Originally developed by Naika Meili
Last update: February 2021
'''

class Water_Calculations(object):

    def Water_Ground(self,q_runon_veg,Runon_tm1,E_ground,Otm1,In_ground_tm1,In_max_ground,Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,
                     CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs,dts,rhow):
        """
        ------
        INPUT:
        q_runon_veg: Runoff from vegetation plus anthropogenic water (Summation of dripping, saturation excess, and throughfall) [mm s^-1]
        Runon_tm1: Runon [mm s^-1]
        E_ground: Evaporative flux from runon water on ground under vegetation [kg m^-2 s^-1]
        Otm1: Soil moisture/water content in the different soil layers [-]
        In_ground_tm1: Intercepted water [mm]
        In_max_ground: Maximum interception capacity of ground under vegetation [mm]
        Pcla: Fraction of clay in the soil [-]
        Psan: Fraction of sand in the soil [-]
        Porg: Fraction of organic material in the soil [-]
        Kfc: Conductivity at field capacity [mm h^-1]
        Phy: Suction at the residual/hygroscopic water content [kPa]
        SPAR: Soil parameter type
        Kbot: Conductivity at the bedrock layer [mm h^-1]
        CASE_ROOT_H: Type of Root Profile of high vegetation
        CASE_ROOT_L: Type of Root Profile of low vegetation
        ZR95_H: Root depth 95 percentile of high vegetation [mm]
        ZR95_L: Root depth 95 percentile of low vegetation [mm]
        ZR50_H: Root depth 50 percentile of high vegetation [mm]
        ZR50_L: Root depth 50 percentile of low vegetation [mm]
        ZRmax_H: Maximum Root depth of high vegetation [mm]
        ZRmax_L: Maximum Root depth of low vegetation [mm]
        Zs: Soil layer discretization [mm]
        dts: Time step [s]
        rhow: Density of water [kg m^-3]
        -------
        OUTPUT:
        q_runoff_ground: Runoff [mm s^-1]
        In_ground: Intercepted water [mm]
        dIn_ground_dt: Change in intercepted water (storage term) [mm s^-1]
        f_ground: Infiltration rate [mm s^-1]
        WBalance_In_ground: Volume balance [mm s^-1]
        """

        # Re-define input parameter which is overwritten in this function
        q_runon_veg_local = copy.copy(q_runon_veg)

        # Calculate soil parameters depending on soil composition
        SoilCal = Soil_Calculations()
        SoilCal.Soil_Parameters_Total(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,
                                      ZRmax_H,ZRmax_L,Zs)
        Osat = SoilCal.SoilParamTotal.Osat   # [-]
        Ohy = SoilCal.SoilParamTotal.Ohy     # [-]
        nVG = SoilCal.SoilParamTotal.nVG     # [mm^-1]
        alpVG = SoilCal.SoilParamTotal.alpVG # [mm^-1]
        Ks_Zs = SoilCal.SoilParamTotal.Ks_Zs # [mm s^-1]
        L = SoilCal.SoilParamTotal.L         # [-]
        Pe = SoilCal.SoilParamTotal.Pe       # [kPa]
        O33 = SoilCal.SoilParamTotal.O33     # [-]
        Zinf =  SoilCal.SoilParamTotal.Zinf  # [mm]

        # Calculate input fluxes into the system ((Runoff from vegetation + anthropogenic water) + Runon(Ponding)) [mm s^-1]
        q_runon_veg_local = q_runon_veg_local+Runon_tm1
        # Calculate intercepted water for the new time step (water incoming to first soil layer) [mm]
        In_ground = In_ground_tm1 + q_runon_veg_local*dts - E_ground*dts*1000/rhow
        # Total water flux incoming to Soil Layer [mm s^-1]
        WIS = In_ground/dts
        # Interception/ponding from the previous time step as an indicator if there is ponding infiltration [mm s^-1]
        In_flux = In_ground_tm1/dts

        # Infiltration rate [mm s^-1]
        SoilCal.Infiltration_2(Osat[0],Ohy[0],L[0],alpVG[0],nVG[0],Pe[0],Ks_Zs[0],O33[0],SPAR,Otm1[0],Zinf,WIS,1,In_flux,dts)
        f_ground = SoilCal.Infiltration2.f

        # Update intercepted water [mm]
        In_ground = In_ground_tm1 + q_runon_veg_local*dts - E_ground*dts*1000/rhow - f_ground*dts

        # Calculate runoff by checking if the intercepted water exceed its maximum [mm s^-1]
        q_runoff_ground = ((In_ground - In_max_ground)/dts) * int(In_ground > In_max_ground)

        # Update interception by checking if the intercepted water exceed its maximum [mm]
        In_ground = In_ground - (In_ground-In_max_ground)*int(In_ground > In_max_ground)

        # Calculate change in intercepted water (storage term) [mm s^-1]
        dIn_ground_dt = (In_ground-In_ground_tm1)/dts

        # Volume balance check [mm s^-1]
        WBalance_In_ground = q_runon_veg_local - E_ground*1000/rhow - f_ground - q_runoff_ground - dIn_ground_dt

        self.WaterGround = WaterGroundDef()
        self.WaterGround.q_runoff_ground = q_runoff_ground
        self.WaterGround.In_ground = In_ground
        self.WaterGround.dIn_ground_dt = dIn_ground_dt
        self.WaterGround.f_ground = f_ground
        self.WaterGround.WBalance_In_ground = WBalance_In_ground

        return q_runoff_ground,In_ground,dIn_ground_dt,f_ground,WBalance_In_ground

    def Water_Impervious(self,Rain,Runon_tm1,E_imp,In_imp_tm1,dts,rhow,In_max_imp,K_imp):

        """
        ------
        INPUT:
        Rain: Precipitation [mm s^-1]
        Runon_tm1: Runon from previous time step [mm s^-1]
        E_imp: Evaporation from impervious surfaces [kg m^-2 s^-1]
        In_imp_tm1: Interception from previous time step [mm]
        dts: Time step [s]
        rhow: Density of water [kg m^-3]
        In_max_imp: Maximum interception capacity of urban area [mm]
        K_imp: Hydraulic conductivity of impervious area [mm h^-1]
        ------
        OUTPUT:
        q_runoff_imp: Runoff [mm s^-1]
        In_imp: Intercepted water [mm]
        dIn_imp_dt: Storage term in water balance equation (change in intercepted water) [mm s^-1]
        Lk_imp: Leakage (infiltration) from impervious area [mm s^-1]
        WBalance_In_imp: Volume balance [mm s^-1]
        """

        # Re-define input parameters which are overwritten in this function
        Rain_local = copy.copy(Rain)
        # Calculate input fluxes into the system (Rain + Runon(Ponding)) [mm s^-1]
        Rain_local = Rain_local+Runon_tm1
        # Calculate intercepted water for the new time step [mm]
        In_imp = In_imp_tm1 + Rain_local*dts - E_imp*dts*1000/rhow

        # Leakage (infiltration) from impervious area [mm s^-1]
        Lk_imp = 0
        # Update intercepted water by removing leakage [mm]
        In_imp = In_imp_tm1 + Rain_local*dts - E_imp*dts*1000/rhow - Lk_imp*dts

        # Calculate runoff by checking if the intercepted water exceed its maximum [mm s^-1]
        q_runoff_imp = ( (In_imp -In_max_imp)/dts ) * int(In_imp > In_max_imp)

        # Update interception by checking if the intercepted water exceed its maximum [mm]
        In_imp = In_imp - q_runoff_imp*dts

        # Calculate change in intercepted water [mm s^-1]
        dIn_imp_dt = (In_imp-In_imp_tm1)/dts

        # Water Balance check [mm s^-1]
        WBalance_In_imp = Rain_local - (E_imp*1000/rhow) - q_runoff_imp - Lk_imp - dIn_imp_dt

        self.WaterImpervious = WaterImperviousDef()
        self.WaterImpervious.q_runoff_imp = q_runoff_imp
        self.WaterImpervious.In_imp = In_imp
        self.WaterImpervious.dIn_imp_dt = dIn_imp_dt
        self.WaterImpervious.Lk_imp = Lk_imp
        self.WaterImpervious.WBalance_In_imp = WBalance_In_imp

        return q_runoff_imp,In_imp,dIn_imp_dt,Lk_imp,WBalance_In_imp

    def Water_Soil(self,Otm1,f,TE_H,TE_L,E_soil,Qlat_in,dts,Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,CASE_ROOT_H,CASE_ROOT_L,ZR95_H,
                   ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L,Zs,rhow):
        """
        ------
        INPUT:
        Otm1: Soil moisture in the different soil layers [-]
        f: Infiltration rate into first soil layer [mm s^-1]
        TE_H: Transpiration from high vegetation [kg m^-2 s^-1]
        TE_L: Transpiration from low vegetation [kg m^-2 s^-1]
        E_soil: Evaporative flux from soil under vegetation (first layer) [kg m^-2 s^-1]
        Qlat_in: Lateral water flux [mm]
        dts: Time step [s]
        Pcla: Fraction of clay in the soil [-]
        Psan: Fraction of sand in the soil [-]
        Porg: Fraction of organic material in the soil [-]
        Kfc: Conductivity at field capacity [mm h^-1]
        Phy: Suction at the residual/hygroscopic water content [kPa]
        SPAR: Soil parameter type
        Kbot: Conductivity at the bedrock layer [mm h^-1]
        CASE_ROOT_H: Type of Root Profile of high vegetation
        CASE_ROOT_L: Type of Root Profile of low vegetation
        ZR95_H: Root depth 95 percentile of high vegetation [mm]
        ZR95_L: Root depth 95 percentile of low vegetation [mm]
        ZR50_H: Root depth 50 percentile of high vegetation [mm]
        ZR50_L: Root depth 50 percentile of low vegetation [mm]
        ZRmax_H: Maximum Root depth of high vegetation [mm]
        ZRmax_L: Maximum Root depth of low vegetation [mm]
        Rrootl_H: Root length index of high vegetation [m root m^-2 PFT]
        Rrootl_L: Root length index of low vegetation [m root m^-2 PFT]
        PsiL50_H: Water Potential at 50% loss conductivity of high vegetation [MPa]
        PsiL50_L: Water Potential at 50% loss conductivity of low vegetation [MPa]
        PsiX50_H: Water potential at 50% of xylem hydraulic conductivity and limit for water extraction from soil of high vegetation [MPa]
        PsiX50_L: Water potential at 50% of xylem hydraulic conductivity and limit for water extraction from soil of low vegetation [MPa]
        Zs: Soil layer discretization [mm]
        rhow: Density of water [kg m^-3]
        -------
        OUTPUT:
        V: Volume of water in each soil layer per unit area [mm]. Never includes the residual water content Ohy
        O: Soil moisture in each soil layer [-]. Always includes the residual water content Ohy
        OS: Water Content for evaporation [-]
        Lk: Leakage at bedrock [mm s^-1]
        Psi_s_H: Soil water potential for first layer of vegetation [MPa]
        Psi_s_L: Soil Water Potential  for Second layer of Vegetation [MPa]
        Exwat_H: Maximum extractable water [mm m^2 m^-2 ground s^-1]
        Exwat_L: Maximum extractable water [mm m^2 m^-2 ground s^-1]
        Rd: Water table rise that reaches the surface level (Dunne Runoff) [mm]
        TE_L_local: Transpiration from low vegetation [kg m^-2 s^-1]
        TE_H_local: Transpiration from high vegetation [kg m^-2 s^-1]
        E_soil_local: Evaporative flux from soil under vegetation (first layer) [kg m^-2 s^-1]
        dV_dt: Change in the soil water content (storage term) [mm s^-1]
        WBalance_soil: Soil water balance [mm s^-1]
        Psi_soil: Tension at O [mm]
        Ko: Hydraulic conductivity at O [mm s^-1]
        """

        # Re-define input parameters which are overwritten in this function
        TE_H_local = copy.copy(TE_H)
        TE_L_local = copy.copy(TE_L)
        E_soil_local = copy.copy(E_soil)
        Qlat_in_local = copy.copy(Qlat_in)

        # Calculate soil parameters depending on soil composition
        SoilCal = Soil_Calculations()
        SoilCal.Soil_Parameters_Total(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,
                                      ZRmax_H,ZRmax_L,Zs)
        dz = SoilCal.SoilParamTotal.dz          # [mm]
        ms = SoilCal.SoilParamTotal.ms          # [-]
        Osat = SoilCal.SoilParamTotal.Osat      # [-]
        Ohy = SoilCal.SoilParamTotal.Ohy        # [-]
        nVG = SoilCal.SoilParamTotal.nVG        # [mm^-1]
        alpVG = SoilCal.SoilParamTotal.alpVG    # [mm^-1]
        Ks_Zs = SoilCal.SoilParamTotal.Ks_Zs    # [mm s^-1]
        L = SoilCal.SoilParamTotal.L            # [-]
        Pe = SoilCal.SoilParamTotal.Pe          # [kPa]
        O33 = SoilCal.SoilParamTotal.O33        # [-]
        EvL_Zs = SoilCal.SoilParamTotal.EvL_Zs  # [-]
        Inf_Zs = SoilCal.SoilParamTotal.Inf_Zs  # [-]
        RfH_Zs = SoilCal.SoilParamTotal.RfH_Zs  # [%]
        RfL_Zs = SoilCal.SoilParamTotal.RfL_Zs  # [%]
        Slo_pot = SoilCal.SoilParamTotal.Slo_pot
        Dz = SoilCal.SoilParamTotal.Dz          # [mm]
        aR = SoilCal.SoilParamTotal.aR
        aTop = SoilCal.SoilParamTotal.aTop      # [mm]

        E_soil_local = E_soil_local*1000/rhow # [mm s^-1]
        TE_L_local = TE_L_local*1000/rhow     # [mm s^-1]
        TE_H_local = TE_H_local*1000/rhow     # [mm s^-1]

        # Distributed Sink: Evaporation from bare soil and Transpiration
        # Evaporation from bare soil * ratio of soil depth [mm s^-1]
        E_soil_dis = E_soil_local * EvL_Zs
        # Root water uptake from different soil layers of high vegetation [mm s^-1]
        TE_dis_H = TE_H_local * RfH_Zs
        # Root water uptake from different soil layers of low vegetation [mm s^-1]
        TE_dis_L = TE_L_local * RfL_Zs

        # Calculate leakage at bedrock
        SoilCal.Leakage_Bottom(Otm1,Ks_Zs,Osat,Ohy,L,nVG,SoilCal.SoilParamTotal.Kbot,ms,SPAR)
        # Leakage at bedrock [mm s^-1]
        Lk = SoilCal.LeakageBottom.Lk
        # Lateral water flux [mm s^-1]
        Qlat_in_local = [Qlat_in_local[i]/dts for i in range(0,len(Qlat_in_local))]

        # Solving Richards equation and calculating the change in water volume
        # Initial value for water volume in each layer [mm]
        V0 = [(Otm1[i] - Ohy[i])*dz[i] for i in range(0,len(dz))]
        T_SPAN = [0,dts]
        ISeep = numpy.ones(ms)
        # Water content in each soil layer [mm]. V never includes the residual water content Ohy
        Vout = odeint(SoilCal.Soil_Moistures_Rich_Comp,V0,T_SPAN,args=(Lk,f,E_soil_dis,TE_dis_H,TE_dis_L,Qlat_in_local,Slo_pot,
                                                                       ISeep,SPAR,Osat,Ohy,O33,dz,Ks_Zs,Dz,ms,L,Pe,aR,aTop,
                                                                       alpVG,nVG,Zs,1,0,0))
        V = [Vout[-1,iV] for iV in range(len(V0))]
        if numpy.isnan(sum(V)):
            print('NaN values in the Volumes')

        # Update soil moisture content in different soil layers
        SoilCal.Soil_Water_MultiLayer(V,Zs,dz,ms,Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,SPAR,EvL_Zs,Inf_Zs,RfH_Zs,RfL_Zs,
                                      Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L)

        O = SoilCal.SoilWaterMultLay.O             # [-]
        OS = SoilCal.SoilWaterMultLay.OS           # [-]
        Psi_s_H = SoilCal.SoilWaterMultLay.Psi_s_H # [MPa]
        Psi_s_L = SoilCal.SoilWaterMultLay.Psi_s_L # [MPa]
        Exwat_H = SoilCal.SoilWaterMultLay.Exwat_H # [mm m^2 m^-2 ground s^-1]
        Exwat_L = SoilCal.SoilWaterMultLay.Exwat_L # [mm m^2 m^-2 ground s^-1]
        Rd = SoilCal.SoilWaterMultLay.Rd           # [mm]
        WTR = SoilCal.SoilWaterMultLay.WTR         # [mm]

        Ko = numpy.zeros(len(O))
        Psi_soil = numpy.zeros(len(O))
        for i in range(0,len(O)):
            SoilCal.Conductivity_Suction(SPAR,Ks_Zs[i],Osat[i],Ohy[i],L[i],Pe[i],O33[i],alpVG[i],nVG[i],O[i])
            Ko[i] = SoilCal.CondSuc.Ko       # [mm s^-1]
            Psi_soil[i] = SoilCal.CondSuc.Po # [mm]

        # Volume Correction for Rd and WTR
        V[0] = V[0]+WTR[1]-Rd
        for i in range(1,len(V)-1):
            V[i] = V[i]+(WTR[i+1]-WTR[i])
        V[-1] = V[-1]-WTR[-1]

        # Volume Compensation - Correct negative volumes
        if sum(i < 0 for i in V):
            SoilCal.Volume_Correction(V,EvL_Zs,RfH_Zs,RfL_Zs,dts*E_soil_local,dts*TE_dis_H,dts*TE_dis_L,dts*Lk)
            V = SoilCal.VolCorr.V             # [mm]
            TE_dis_H = SoilCal.VolCorr.T_H    # [mm]
            TE_dis_L = SoilCal.VolCorr.T_L    # [mm]
            E_soil_local = SoilCal.VolCorr.EG # [mm]
            Lk = SoilCal.VolCorr.Lk           # [mm]

            TE_H_local = sum(TE_dis_H)/dts  # [mm s^-1]
            TE_L_local = sum(TE_dis_L)/dts  # [mm s^-1]
            E_soil_local = E_soil_local/dts # [mm s^-1]
            Lk = Lk/dts                     # [mm s^-1]



        # Change in the soil water content (storage term) [mm s^-1]
        dV_dt = (sum(V)-sum(V0))/dts
        # Volume Balance check [mm s^-1]
        WBalance_soil = f + numpy.nansum(Qlat_in_local) - Lk - TE_H_local - TE_L_local - E_soil_local - Rd/dts - dV_dt

        E_soil_local = E_soil_local / (1000/rhow) # [kg m^-2 s^-1]
        TE_L_local = TE_L_local / (1000/rhow)     # [kg m^-2 s^-1]
        TE_H_local = TE_H_local / (1000/rhow)     # [kg m^-2 s^-1]


        self.WaterSoil = WaterSoilDef()
        self.WaterSoil.V = V                         # Volume of water in each soil layer per unit area [mm]. Never includes the residual water content Ohy
        self.WaterSoil.O = O                         # Soil moisture in each soil layer [-]. Always includes the residual water content Ohy
        self.WaterSoil.OS = OS                       # [-]
        self.WaterSoil.Lk = Lk                       # Leakage at bedrock [mm s^-1]
        self.WaterSoil.Psi_s_H = Psi_s_H             # [MPa]
        self.WaterSoil.Psi_s_L = Psi_s_L             # [MPa]
        self.WaterSoil.Exwat_H = Exwat_H             # Max extractable water for Vegetation [mm m^2 m^-2 ground s^-1]
        self.WaterSoil.Exwat_L = Exwat_L             # Max extractable water for Vegetation [mm m^2 m^-2 ground s^-1]
        self.WaterSoil.Rd = Rd                       # saturation excess runoff / Dunne Runoff [mm]
        self.WaterSoil.TE_L = TE_L_local             # [kg m^-2 s^-1]
        self.WaterSoil.TE_H = TE_H_local             # [kg m^-2 s^-1]
        self.WaterSoil.E_soil = E_soil_local         # [kg m^-2 s^-1]
        self.WaterSoil.dV_dt = dV_dt                 # [mm s^-1]
        self.WaterSoil.WBalance_soil = WBalance_soil # [mm s^-1]
        self.WaterSoil.Psi_soil = Psi_soil           # [mm]
        self.WaterSoil.Ko = Ko                       # [mm s^-1]

        return V,O,OS,Lk,Psi_s_H,Psi_s_L,Exwat_H,Exwat_L,Rd,TE_L_local,TE_H_local,E_soil_local,dV_dt,WBalance_soil,Psi_soil,Ko

    def Water_Vegetation(self,Rain,E_veg,In_veg_tm1,Sp_In,LAI,SAI,rhow,dts):

        """
        ------
        INPUT:
        Rain: Precipitation [mm s^-1]
        E_veg: Evaporative flux from intercepted water by vegetation [kg m^-2 s^-1]
        In_veg_tm1: Interception from previous time step [mm]
        Sp_In: Specific water retained by a vegetated surface [mm m^2 VEG area m^-2 plant area]
        LAI: Leaf area index [m^2 m^-2]
        SAI: Stem area index [m^2 m^-2]
        rhow: Density f water [kg m^-3]
        dts: Time step [s]
        -------
        OUTPUT:
        q_runoff_veg: Runoff [mm s^-1]
        In_veg: Intercepted water [mm]
        dIn_veg_dt: Change in interception (storage term in water balance equation) [mm s^-1]
        WBalance_In_veg: Volume balance [mm s^-1]
        """

        # Throughfall through foliage and dripping parameters
        PAI = LAI+SAI
        # Model constant (Ramirez and Senarath, 2000)
        Kthroughfall = 0.75
        # Projected leaf area fraction onto the ground [-]
        Cfol = 1 - numpy.exp(-Kthroughfall * PAI)
        # Exponential decay parameter [mm^-1]
        gc = 3.7
        # Drainage rate coefficient [mm s^-1]
        Kc = (0.001*60)/3600

        # Maximum interception capacity of vegetation [mm]
        In_max_veg = Sp_In * (LAI + SAI)

        # One vegetation layer
        # Precipitation onto the canopy foliage [mm s^-1]
        Rain_fol = Rain * Cfol
        # Throughfall [mm s^-1]
        Rain_tf = Rain * (1-Cfol)

        # Calculate intercepted water for the new time step [mm]
        In_veg = In_veg_tm1 + Rain_fol*dts - E_veg*dts*1000/rhow
        # Calculate storage (saturation) excess [mm]
        SE_veg = (In_veg - In_max_veg) * int(In_veg > In_max_veg)

        # Update intercepted water [mm]
        In_veg = In_veg - SE_veg

        # Calculate dripping [mm s^-1]
        Dr_veg = Kc*numpy.exp(gc*(In_veg-In_max_veg)) * int(In_veg > 0)

        # Update intercepted water [mm]
        In_veg = In_veg - Dr_veg*dts

        # Update dripping [mm s^-1]
        Dr_veg = Dr_veg + (In_veg/dts)*int(In_veg < 0)

        # Update intercepted water [mm]
        In_veg = In_veg * int(In_veg > 0)

        # Summation of dripping, saturation excess, and throughfall [mm s^-1]
        q_runoff_veg = Dr_veg + SE_veg/dts + Rain_tf

        # Calculate change in interception (storage term in water balance equation) [mm s^-1]
        dIn_veg_dt = (In_veg - In_veg_tm1)/dts

        # Volume Balance check [mm s^-1]
        WBalance_In_veg = Rain - (E_veg*1000/rhow) - q_runoff_veg - dIn_veg_dt

        self.WaterVegetation = WaterVegetationDef()
        self.WaterVegetation.q_runoff_veg = q_runoff_veg
        self.WaterVegetation.In_veg = In_veg
        self.WaterVegetation.dIn_veg_dt = dIn_veg_dt
        self.WaterVegetation.WBalance_In_veg = WBalance_In_veg

        return q_runoff_veg,In_veg,dIn_veg_dt,WBalance_In_veg
