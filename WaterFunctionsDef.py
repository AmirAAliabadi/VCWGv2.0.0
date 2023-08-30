import os
import numpy
import math
from pprint import pprint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.integrate import odeint


class WaterCanyonDef(object):

    def __init__(self,q_tree_dwn=None,In_tree=None,dIn_tree_dt=None,q_gveg_dwn=None,In_gveg=None,dIn_gveg_dt=None,
                 q_gimp_runoff=None,In_gimp=None,dIn_gimp_dt=None,f_inf_gimp=None,q_gbare_runoff=None,In_gbare=None,
                 dIn_gbare_dt=None,f_inf_gbare=None,q_gveg_runoff=None,In_gveg_pond=None,dIn_gveg_pond_dt=None,
                 f_inf_gveg=None,V_gimp=None,O_gimp=None,OS_gimp=None,Lk_gimp=None,Psi_s_H_gimp=None,Psi_s_L_gimp=None,
                 Exwat_H_gimp=None,Exwat_L_gimp=None,Rd_gimp=None,TEgveg_imp=None,TEtree_imp=None,Egimp_soil=None,
                 dV_dt_gimp=None,Psi_soil_gimp=None,Kf_gimp=None,V_gbare=None,O_gbare=None,OS_gbare=None,Lk_gbare=None,
                 Psi_s_H_gbare=None,Psi_s_L_gbare=None,Exwat_H_gbare=None,Exwat_L_gbare=None,Rd_gbare=None,TEgveg_bare=None,
                 TEtree_bare=None,Egbare_Soil=None,dV_dt_gbare=None,Psi_soil_gbare=None,Kf_gbare=None,V_gveg=None,
                 O_gveg=None,OS_gveg=None,Lk_gveg=None,Psi_s_H_gveg=None,Psi_s_L_gveg=None,Exwat_H_gveg=None,Exwat_L_gveg=None,
                 Rd_gveg=None,TEgveg_veg=None,TEtree_veg=None,Egveg_Soil=None,dV_dt_gveg=None,Psi_soil_gveg=None,Kf_gveg=None,
                 Qin_imp=None,Qin_bare=None,Qin_veg=None,Qin_bare2imp=None,Qin_bare2veg=None,Qin_imp2bare=None,Qin_imp2veg=None,
                 Qin_veg2imp=None,Qin_veg2bare=None,V=None,O=None,OS=None,Lk=None,Rd=None,dV_dt=None,Psi_s_L=None,Exwat_L=None,
                 TEgveg_tot=None,Psi_s_H_tot=None,Exwat_H=None,TEtree_tot=None,EB_TEtree=None,EB_TEgveg=None,WBIndv=None,
                 WBTot=None,Runoff=None,Runon=None,Etot=None,DeepGLk=None,StorageTot=None):
        self.q_tree_dwn = q_tree_dwn
        self.In_tree = In_tree
        self.dIn_tree_dt = dIn_tree_dt
        self.q_gveg_dwn = q_gveg_dwn
        self.In_gveg = In_gveg
        self.dIn_gveg_dt = dIn_gveg_dt
        self.q_gimp_runoff = q_gimp_runoff
        self.In_gimp = In_gimp
        self.dIn_gimp_dt = dIn_gimp_dt
        self.f_inf_gimp = f_inf_gimp
        self.q_gbare_runoff = q_gbare_runoff
        self.In_gbare = In_gbare
        self.dIn_gbare_dt = dIn_gbare_dt
        self.f_inf_gbare = f_inf_gbare
        self.q_gveg_runoff = q_gveg_runoff
        self.In_gveg_pond = In_gveg_pond
        self.dIn_gveg_pond_dt = dIn_gveg_pond_dt
        self.f_inf_gveg = f_inf_gveg
        self.V_gimp = V_gimp
        self.O_gimp = O_gimp
        self.OS_gimp = OS_gimp
        self.Lk_gimp = Lk_gimp
        self.Psi_s_H_gimp = Psi_s_H_gimp
        self.Psi_s_L_gimp = Psi_s_L_gimp
        self.Exwat_H_gimp = Exwat_H_gimp
        self.Exwat_L_gimp = Exwat_L_gimp
        self.Rd_gimp = Rd_gimp
        self.TEgveg_imp = TEgveg_imp
        self.TEtree_imp = TEtree_imp
        self.Egimp_soil = Egimp_soil
        self.dV_dt_gimp = dV_dt_gimp
        self.Psi_soil_gimp = Psi_soil_gimp
        self.Kf_gimp = Kf_gimp
        self.V_gbare = V_gbare
        self.O_gbare = O_gbare
        self.OS_gbare = OS_gbare
        self.Lk_gbare = Lk_gbare
        self.Psi_s_H_gbare = Psi_s_H_gbare
        self.Psi_s_L_gbare = Psi_s_L_gbare
        self.Exwat_H_gbare = Exwat_H_gbare
        self.Exwat_L_gbare = Exwat_L_gbare
        self.Rd_gbare = Rd_gbare
        self.TEgveg_bare = TEgveg_bare
        self.TEtree_bare = TEtree_bare
        self.Egbare_Soil = Egbare_Soil
        self.dV_dt_gbare = dV_dt_gbare
        self.Psi_soil_gbare = Psi_soil_gbare
        self.Kf_gbare = Kf_gbare
        self.V_gveg = V_gveg
        self.O_gveg = O_gveg
        self.OS_gveg = OS_gveg
        self.Lk_gveg = Lk_gveg
        self.Psi_s_H_gveg = Psi_s_H_gveg
        self.Psi_s_L_gveg = Psi_s_L_gveg
        self.Exwat_H_gveg = Exwat_H_gveg
        self.Exwat_L_gveg = Exwat_L_gveg
        self.Rd_gveg = Rd_gveg
        self.TEgveg_veg = TEgveg_veg
        self.TEtree_veg = TEtree_veg
        self.Egveg_Soil = Egveg_Soil
        self.dV_dt_gveg = dV_dt_gveg
        self.Psi_soil_gveg = Psi_soil_gveg
        self.Kf_gveg = Kf_gveg
        self.Qin_imp = Qin_imp
        self.Qin_bare = Qin_bare
        self.Qin_veg = Qin_veg
        self.Qin_bare2imp = Qin_bare2imp
        self.Qin_bare2veg = Qin_bare2veg
        self.Qin_imp2bare = Qin_imp2bare
        self.Qin_imp2veg = Qin_imp2veg
        self.Qin_veg2imp = Qin_veg2imp
        self.Qin_veg2bare = Qin_veg2bare
        self.V = V
        self.O = O
        self.OS = OS
        self.Lk = Lk
        self.Rd = Rd
        self.dV_dt = dV_dt
        self.Psi_s_L = Psi_s_L
        self.Exwat_L = Exwat_L
        self.TEgveg_tot = TEgveg_tot
        self.Psi_s_H_tot = Psi_s_H_tot
        self.Exwat_H = Exwat_H
        self.TEtree_tot = TEtree_tot
        self.EB_TEtree = EB_TEtree
        self.EB_TEgveg = EB_TEgveg
        self.WBIndv = WBIndv
        self.WBTot = WBTot
        self.Runoff = Runoff
        self.Runon = Runon
        self.Etot = Etot
        self.DeepGLk = DeepGLk
        self.StorageTot = StorageTot

class IntDef(object):

    def __init__(self,IntGroundImp=None,IntGroundBare=None,IntGroundVegPlant=None,IntGroundVegGround=None,IntTree=None):
        self.IntGroundImp = IntGroundImp
        self.IntGroundBare = IntGroundBare
        self.IntGroundVegPlant = IntGroundVegPlant
        self.IntGroundVegGround = IntGroundVegGround
        self.IntTree = IntTree

class RunonDef(object):

    def __init__(self, RunonGroundTot=None,RunoffGroundTot=None):
        self.RunonGroundTot = RunonGroundTot
        self.RunoffGroundTot = RunoffGroundTot

class QinlatDef(object):

    def __init__(self,Qin_bare2imp=None,Qin_veg2imp=None,Qin_veg2bare=None,Qin_imp2bare=None,Qin_bare2veg=None,
                 Qin_imp2veg=None,Qin_imp=None,Qin_bare=None,Qin_veg=None):
        self.Qin_bare2imp = Qin_bare2imp
        self.Qin_veg2imp = Qin_veg2imp
        self.Qin_veg2bare = Qin_veg2bare
        self.Qin_imp2bare = Qin_imp2bare
        self.Qin_bare2veg = Qin_bare2veg
        self.Qin_imp2veg = Qin_imp2veg
        self.Qin_imp = Qin_imp
        self.Qin_bare = Qin_bare
        self.Qin_veg = Qin_veg

class WaterVegetationDef(object):

    def __init__(self,q_runon_veg=None,In_veg=None,dIn_veg_dt=None,WBalance_In_veg=None):
        self.q_runon_veg = q_runon_veg
        self.In_veg = In_veg
        self.dIn_veg_dt = dIn_veg_dt
        self.WBalance_In_veg = WBalance_In_veg

class WaterImperviousDef(object):

    def __init__(self,q_runon_imp=None,In_imp=None,dIn_imp_dt=None,Lk_imp=None,WBalance_In_imp=None):
        self.q_runon_imp = q_runon_imp
        self.In_imp = In_imp
        self.dIn_imp_dt = dIn_imp_dt
        self.Lk_imp = Lk_imp
        self.WBalance_In_imp = WBalance_In_imp

class WaterGroundDef(object):

    def __init__(self,q_runon_ground=None,In_ground=None,dIn_ground_dt=None,f_ground=None,WBalance_In_ground=None):
        self.q_runon_ground = q_runon_ground
        self.In_ground = In_ground
        self.dIn_ground_dt = dIn_ground_dt
        self.f_ground = f_ground
        self.WBalance_In_ground = WBalance_In_ground


class WaterSoilDef(object):

    def __init__(self,V=None,O=None,OS=None,Lk=None,Psi_s_H=None,Psi_s_L=None,Exwat_H=None,Exwat_L=None,Rd=None,TE_L=None,
                 TE_H=None,E_soil=None,dV_dt=None,WBalance_soil=None,Psi_soil=None,Ko=None):
        self.V = V
        self.O = O
        self.OS = OS
        self.Lk = Lk
        self.Psi_s_H = Psi_s_H
        self.Psi_s_L = Psi_s_L
        self.Exwat_H = Exwat_H
        self.Exwat_L = Exwat_L
        self.Rd = Rd
        self.TE_L = TE_L
        self.TE_H = TE_H
        self.E_soil = E_soil
        self.dV_dt = dV_dt
        self.WBalance_soil = WBalance_soil
        self.Psi_soil = Psi_soil
        self.Ko = Ko

class SoilPotWDef(object):

    def __init__(self,SoilPotWRoofVeg_H=None,SoilPotWGroundImp_H=None,SoilPotWGroundBare_H=None,SoilPotWGroundVeg_H=None,
                 SoilPotWGroundTot_H=None,SoilPotWRoofVeg_L=None,SoilPotWGroundImp_L=None,SoilPotWGroundBare_L=None,
                 SoilPotWGroundVeg_L=None,SoilPotWGroundTot_L=None):
        self.SoilPotWRoofVeg_H = SoilPotWRoofVeg_H
        self.SoilPotWGroundImp_H = SoilPotWGroundImp_H
        self.SoilPotWGroundBare_H = SoilPotWGroundBare_H
        self.SoilPotWGroundVeg_H = SoilPotWGroundVeg_H
        self.SoilPotWGroundTot_H = SoilPotWGroundTot_H
        self.SoilPotWRoofVeg_L = SoilPotWRoofVeg_L
        self.SoilPotWGroundImp_L = SoilPotWGroundImp_L
        self.SoilPotWGroundBare_L = SoilPotWGroundBare_L
        self.SoilPotWGroundVeg_L = SoilPotWGroundVeg_L
        self.SoilPotWGroundTot_L = SoilPotWGroundTot_L

class ExWaterDef(object):

    def __init__(self,ExWaterGroundImp_H=None,ExWaterGroundBare_H=None,ExWaterGroundVeg_H=None,ExWaterGroundVeg_L=None):
        self.ExWaterGroundImp_H = ExWaterGroundImp_H
        self.ExWaterGroundBare_H = ExWaterGroundBare_H
        self.ExWaterGroundVeg_H = ExWaterGroundVeg_H
        self.ExWaterGroundVeg_L = ExWaterGroundVeg_L
























