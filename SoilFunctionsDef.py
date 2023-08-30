import os
import numpy
import math
from pprint import pprint
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy.integrate import odeint

"""
Define soil parameters
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: June 2020
"""

class SoilParametersTotal(object):

    def __init__(self,Zs=None,dz=None,ms=None, Osat=None, Ohy=None, nVG=None, alpVG=None, Ks_Zs=None, L=None, Pe=None,
                 O33=None, SPAR=None, EvL_Zs=None, Inf_Zs=None, RfH_Zs=None, RfL_Zs=None, Zinf=None, Kbot=None, Slo_pot=None,
                 Dz=None, aR=None, aTop=None,rsd=None, lan_dry=None, lan_s=None, cv_s=None):
        self.Zs = Zs
        self.dz = dz
        self.ms = ms
        self.Osat = Osat
        self.Ohy = Ohy
        self.nVG = nVG
        self.alpVG = alpVG
        self.Ks_Zs = Ks_Zs
        self.L = L
        self.Pe = Pe
        self.O33 = O33
        self.SPAR = SPAR
        self.EvL_Zs = EvL_Zs
        self.Inf_Zs = Inf_Zs
        self.RfH_Zs = RfH_Zs
        self.RfL_Zs = RfL_Zs
        self.Zinf = Zinf
        self.Kbot = Kbot
        self.Slo_pot = Slo_pot
        self.Dz = Dz
        self.aR = aR
        self.aTop = aTop
        self.rsd = rsd
        self.lan_dry = lan_dry
        self.lan_s = lan_s
        self.cv_s = cv_s

class ConductivitySuction(object):

    def __init__(self,Ko=None,Po=None):
        self.Ko = Ko
        self.Po = Po

class SoilThermalProperties(object):

    def __init__(self,lanS=None,cv_Soil=None,CTt=None):
        self.lanS = lanS
        self.cv_Soil = cv_Soil
        self.CTt = CTt

class SoilParameters(object):

    def __init__(self,Osat=None,L=None,Pe=None,Ks=None,O33=None,rsd=None,lan_dry=None,lan_s=None,cv_s=None,K_usle=None):
        self.Osat = Osat
        self.L = L
        self.Pe = Pe
        self.Ks = Ks
        self.O33 = O33
        self.rsd = rsd
        self.lan_dry = lan_dry
        self.lan_s = lan_s
        self.cv_s = cv_s
        self.K_usle = K_usle

class SoilMoistureConductivityUpdate(object):

    def __init__(self,V=None,O=None,OS=None,Psi_soil=None,Psi_s_H=None,Psi_s_L=None,Exwat_H=None,Exwat_L=None,Ko=None):
        self.V = V
        self.O = O
        self.OS = OS
        self.Psi_soil = Psi_soil
        self.Psi_s_H = Psi_s_H
        self.Psi_s_L = Psi_s_L
        self.Exwat_H = Exwat_H
        self.Exwat_L = Exwat_L
        self.Ko = Ko

class OwaterDef(object):

    def __init__(self,OwGroundSoilImp=None,OwGroundSoilBare=None,OwGroundSoilVeg=None):
        self.OwGroundSoilImp = OwGroundSoilImp
        self.OwGroundSoilBare = OwGroundSoilBare
        self.OwGroundSoilVeg = OwGroundSoilVeg

class VwaterDef(object):

    def __init__(self,VGroundSoilImp=None,VGroundSoilBare=None,VGroundSoilVeg=None,VGroundSoilTot=None):
        self.VGroundSoilImp = VGroundSoilImp
        self.VGroundSoilBare = VGroundSoilBare
        self.VGroundSoilVeg = VGroundSoilVeg
        self.VGroundSoilTot = VGroundSoilTot

class ParSoilGroundDef(object):

    def __init__(self,Zs=None,ms=None,In_max_imp=None,In_max_underveg=None,In_max_bare=None,Sp_In=None,Kimp=None,Kfc=None,
                 Phy=None,SPAR=None,Kbot=None,Pcla=None,Psan=None,Porg=None,O33=None,dz=None):
        self.Zs = Zs
        self.ms = ms
        self.In_max_imp = In_max_imp
        self.In_max_underveg = In_max_underveg
        self.In_max_bare = In_max_bare
        self.Sp_In = Sp_In
        self.Kimp = Kimp
        self.Kfc = Kfc
        self.Phy = Phy
        self.SPAR = SPAR
        self.Kbot = Kbot
        self.Pcla = Pcla
        self.Psan = Psan
        self.Porg = Porg
        self.O33 = O33
        self.dz = dz

class ParSoilRoofDef(object):
    def __init__(self):
        self.Zs = [0,10,20,50,100]
        self.ms = len(self.Zs)-1
        self.dz1 = 0.1
        self.dz2 = 0.1
        self.In_max_imp = 0.25
        self.In_max_ground = 10
        self.Sp_In = 0.2
        self.Kimp = 0
        self.Kfc = 0.2
        self.Phy = 10000
        self.SPAR = 2
        self.Kbot = None
        self.Pcla = 0.2
        self.Psan = 0.4
        self.Porg = 0.025
        self.O33 = None
        self.dz = numpy.diff(self.Zs)


class Infiltration2Def(object):

    def __init__(self,f=None,fpot=None):
        self.f = f
        self.fpot = fpot

class LeakageBottomDef(object):

    def __init__(self,Lk=None):
        self.Lk = Lk

class SoilWaterMultiLayerDef(object):

    def __init__(self,O=None,ZWT=None,OF=None,OSs=None,Psi_s_H=None,Psi_s_L=None,gsr_H=None,gsr_L=None,Exwat_H=None,
                 Exwat_L=None,Rd=None,WTR=None,POT=None,OH=None,OL=None):
        self.O = O
        self.ZWT = ZWT
        self.OF = OF
        self.OSs = OSs
        self.Psi_s_H = Psi_s_H
        self.Psi_s_L = Psi_s_L
        self.gsr_H = gsr_H
        self.gsr_L = gsr_L
        self.Exwat_H = Exwat_H
        self.Exwat_L = Exwat_L
        self.Rd = Rd
        self.WTR = WTR
        self.POT = POT
        self.OH = OH
        self.OL = OL

class VolumeCorrectionDef(object):

    def __init__(self,V=None,T_H=None,T_L=None,EG=None,Lk=None):
        self.V = V
        self.T_H = T_H
        self.T_L = T_L
        self.EG = EG
        self.Lk = Lk














