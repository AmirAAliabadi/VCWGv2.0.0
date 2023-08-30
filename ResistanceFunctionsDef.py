import os
import numpy
import math
from pprint import pprint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.integrate import odeint

"""
Define parameters in Resistance function
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: March 2021
"""

class ResistanceAero(object):

    def __init__(self,r_ground=None,r_ground_b=None,r_roof=None,r_roof_b=None,r_wall=None):
        self.r_ground = r_ground
        self.r_ground_b = r_ground_b
        self.r_roof = r_roof
        self.r_roof_b = r_roof_b
        self.r_wall = r_wall

class ResistanceSoil(object):

    def __init__(self,r_soil=None,b_soil=None,alp_soil=None):
        self.r_soil = r_soil
        self.b_soil = b_soil
        self.alp_soil = alp_soil

class CanopyResistanceAnEvolution(object):

    def __init__(self,rs_sun=None,rs_shd=None,Ci_sun=None,Ci_shd=None,An=None,Rdark=None,Lpho=None,SIF=None,DCi=None):
        self.rs_sun = rs_sun
        self.rs_shd = rs_shd
        self.Ci_sun = Ci_sun
        self.Ci_shd = Ci_shd
        self.An = An
        self.Rdark = Rdark
        self.Lpho = Lpho
        self.SIF = SIF
        self.DCi = DCi

class PhotosynthesisBiochemical(object):

    def __init__(self,CcF=None,An=None,rs=None,Rdark=None,F755nm=None,GAM=None,gsCO2=None):
        self.CcF = CcF
        self.An = An
        self.rs = rs
        self.Rdark = Rdark
        self.F755nm = F755nm
        self.GAM = GAM
        self.gsCO2 = gsCO2

class CiCO2LeafDef(object):

    def __init__(self,CiCO2LeafRoofVegSun=None,CiCO2LeafRoofVegShd=None,CiCO2LeafGroundVegSun=None,CiCO2LeafGroundVegShd=None,
                 CiCO2LeafTreeSun=None,CiCO2LeafTreeShd=None):
        self.CiCO2LeafRoofVegSun = CiCO2LeafRoofVegSun
        self.CiCO2LeafRoofVegShd = CiCO2LeafRoofVegShd
        self.CiCO2LeafGroundVegSun = CiCO2LeafGroundVegSun
        self.CiCO2LeafGroundVegShd = CiCO2LeafGroundVegShd
        self.CiCO2LeafTreeSun = CiCO2LeafTreeSun
        self.CiCO2LeafTreeShd = CiCO2LeafTreeShd




