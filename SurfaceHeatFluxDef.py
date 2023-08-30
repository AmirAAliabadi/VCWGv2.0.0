import os
import numpy
import math
from pprint import pprint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.integrate import odeint


class ForceRestoreConductiveHeatImp(object):

    def __init__(self,G=None,Tdp=None):
        self.G = G
        self.Tdp = Tdp

class ForceRestoreConductiveHeatSoil(object):

    def __init__(self,G=None,Tdp=None):
        self.G = G
        self.Tdp = Tdp

class AdjustTground(object):

    def __init__(self,T_imp_adj=None,T_soil_adj=None,T_veg_adj=None,LEveg_int=None,LE_imp=None,LE_bare_pond=None,
                 LE_veg_pond=None,LE_bare_soil=None,LE_veg_soil=None,LTE_veg=None,H_imp=None,H_bare_soil=None,H_veg=None):
        self.T_imp_adj = T_imp_adj
        self.T_soil_adj = T_soil_adj
        self.T_veg_adj = T_veg_adj
        self.LEveg_int = LEveg_int
        self.LE_imp = LE_imp
        self.LE_bare_pond = LE_bare_pond
        self.LE_veg_pond = LE_veg_pond
        self.LE_bare_soil = LE_bare_soil
        self.LE_veg_soil = LE_veg_soil
        self.LTE_veg = LTE_veg
        self.H_imp = H_imp
        self.H_bare_soil = H_bare_soil
        self.H_veg = H_veg

class SensibleFluxDef(object):

    def __init__(self,H_ground_imp=None,H_ground_soil=None,H_ground_veg=None):
        self.H_ground_imp = H_ground_imp
        self.H_ground_soil = H_ground_soil
        self.H_ground_veg = H_ground_veg

class LatentFluxDef(object):

    def __init__(self,Etree=None,LEfluxGroundImp=None,LEfluxGroundBarePond=None,LEfluxGroundBareSoil=None,
                 LEfluxGroundVegInt=None,LEfluxGroundVegPond=None,LEfluxGroundVegSoil=None,LTEfluxGroundVeg=None,
                 LEfluxTreeInt=None,LTEfluxTree=None,LEbare=None,LEveg=None,LEtree=None,Ci_sun_tree=None,Ci_shd_tree=None,
                 Ci_sun_ground=None,Ci_shd_ground=None,EfluxGroundImp=None,EfluxGroundBarePond=None,EfluxGroundBareSoil=None,
                 EfluxGroundBare=None,EfluxGroundVegInt=None,EfluxGroundVegPond=None,EfluxGroundVegSoil=None,TEfluxGroundVeg=None,
                 EfluxGroundVeg=None,EfluxGround=None):
        self.EfluxGroundImp = EfluxGroundImp
        self.EfluxGroundBarePond = EfluxGroundBarePond
        self.EfluxGroundBareSoil = EfluxGroundBareSoil
        self.EfluxGroundBare = EfluxGroundBare
        self.EfluxGroundVegInt = EfluxGroundVegInt
        self.EfluxGroundVegPond = EfluxGroundVegPond
        self.EfluxGroundVegSoil = EfluxGroundVegSoil
        self.TEfluxGroundVeg = TEfluxGroundVeg
        self.EfluxGroundVeg = EfluxGroundVeg
        self.EfluxGround = EfluxGround
        self.Etree = Etree
        self.LEfluxGroundImp = LEfluxGroundImp
        self.LEfluxGroundBarePond = LEfluxGroundBarePond
        self.LEfluxGroundBareSoil = LEfluxGroundBareSoil
        self.LEfluxGroundVegInt = LEfluxGroundVegInt
        self.LEfluxGroundVegPond = LEfluxGroundVegPond
        self.LEfluxGroundVegSoil = LEfluxGroundVegSoil
        self.LTEfluxGroundVeg = LTEfluxGroundVeg
        self.LEfluxTreeInt = LEfluxTreeInt
        self.LTEfluxTree = LTEfluxTree
        self.LEbare = LEbare
        self.LEveg = LEveg
        self.LEtree = LEtree
        self.Ci_sun_tree = Ci_sun_tree
        self.Ci_shd_tree = Ci_shd_tree
        self.Ci_sun_ground = Ci_sun_ground
        self.Ci_shd_ground = Ci_shd_ground



