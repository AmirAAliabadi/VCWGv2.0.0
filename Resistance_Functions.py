import os
import numpy
import math
from scipy.interpolate import interp1d
from scipy.optimize import fsolve
from ResistanceFunctionsDef import ResistanceSoil,CanopyResistanceAnEvolution,PhotosynthesisBiochemical
from Soil_Functions import Soil_Calculations
import copy

'''
Resistance Functions:
Developed by Mohsen Moradi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: June 2021
'''

class Ressitance_Calculations(object):

    def Canopy_Resistance_An_Evolution(self,PAR_sun,PAR_shd,LAI,Kopt,Knit,Fsun,Fshd,Citm1_sun,Citm1_shd,Ca,ra,rb,Ts,Pre,
                                       Ds,Psi_L,Psi_sto_50,Psi_sto_99,CT,Vmax,DS,Ha,FI,Oa,Do,a1,go,e_rel,e_relN,gmes,rjv):
        """
        ------
        INPUT:
        PAR_sun: Absorbed direct and diffuse shortwave radiation of the sunlit surface [W m^-2]
        PAR_shd: Absorbed direct and diffuse shortwave radiation of the shaded surface [W m^-2]
        LAI: Leaf area index [-]
        Kopt: optical depth of direct beam perunit plant area [-]
        Knit: Canopy nitrogen decay coefficient [-]
        Fsun: Partitioning of radiation into sunlit area
        Fshd: Partitioning of radiation into shaded area
        Citm1_sun: Leaf Interior CO2 mixing ratio [umolCO2 mol^-1]
        Citm1_shd: Leaf Interior CO2 mixing ratio [umolCO2 mol^-1]
        Ca: Atmospheric CO2 mixing ratio [umolCO2 mol^-1]
        ra: Aerodynamic resistance [s m^-1]
        rb: Leaf boundary resistance [s m^-1]
        Ts: Temperature of the vegetation [C]
        Pre: Pressure at the height of trees [mbar]
        Ds: Vapor Pressure Deficit [Pa]
        Psi_L: Soil water potential for first/second layer of vegetation [MPa]
        Psi_sto_50: Water Potential at 50% loss conductivity [MPa]
        Psi_sto_99: Water Potential at PLCs loss conductivity [MPa]
        CT: Photosyntesis Typology for Plants, Photosynthetic pathway
        Vmax: Maximum Rubisco Capacity [umolCO2 s^-1 m^-2]
        DS: Plant Dependent, Activation Energy in Photosynthesis for Rubisco Capacity [kJ mol^-1]
        Ha: Plant Dependent, Activation energy. [kJ mol^-1 K^-1]
        FI: Intrinsic quantum Efficiency [umolCO2 umolPhotons^-1]
        Oa: Intercellular Partial Pressure Oxygen [ppm] - [umolO2 mol^-1]
        Do: Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis [Pa]
        a1: Empirical parameter connecting stomatal aperture and net assimilation [-]
        go: minimal Stomatal Conductance [mol s^-1 m^-2]
        e_rel: Relative Efficiency of the photosynthesis apparatus due to Age/Day-length [-]
        e_relN: Relative efficiency of the photosynthesis apparatus due to N limitations [-]
        gmes: Mesophyll conductance, not used [mol CO2 s^-1 m^-2]
        rjv: Scaling factor between Jmax and Vmax
        -------
        OUTPUT:
        rs_sun: stomatal resistence [s m^-1]
        rs_shd: stomatal resistence [s m^-1]
        Ci_sun: Leaf Interior CO2 mixing ratio of sunlit [umolCO2 mol^-1]
        Ci_shd: Leaf Interior CO2 mixing ratio of shaded [umolCO2 mol^-1]
        An: Net Assimiltation Rate [umolCO2 s^-1 m^-2 ]
        Rdark: Surface Leaf Concentration [umolCO2 s^-1 m^-2 ]
        """

        # Re-define input parameters which are overwritten in this function
        PAR_sun_local = copy.copy(PAR_sun)
        PAR_shd_local = copy.copy(PAR_shd)
        Citm1_sun_local = copy.copy(Citm1_sun)
        Citm1_shd_local = copy.copy(Citm1_shd)
        Vmax_local = copy.copy(Vmax)

        if Citm1_sun_local < 200:
            Citm1_sun_local = 200
        if Citm1_shd_local < 200:
            Citm1_shd_local = 200

        # ANSW_SCA is assumed to be one
        # Relative efficiency for age
        Vmax_local = Vmax_local*e_rel*e_relN
        # Scaling from leaf to canopy
        # To be recomputed for Vmax only for LAI and with Kopt to avoid issue with SAI LAIdead
        FsunV = (1 - numpy.exp(-Kopt * (LAI))) / (Kopt * (LAI))
        if FsunV < 0.01:
            FsunV = 0
        if FsunV > 1:
            FsunV = 1
        FshdV = 1 - FsunV

        # Two big leaves with Kn
        Can_sun = (1 - numpy.exp(-(Kopt + Knit) * LAI)) / (Kopt + Knit)
        Can_shd = (1 - numpy.exp(-Knit * LAI)) / Knit - (1 - numpy.exp(-(Kopt + Knit) * LAI)) / (Kopt + Knit)
        # Two big leaves with Kn
        Vmax_sun = Vmax_local * Can_sun / (LAI * FsunV)
        Vmax_shd = Vmax_local * Can_shd / (LAI * FshdV)

        if FsunV == 0:
            Vmax_sun = 0

        # minimum canopy conductance
        go_sun = go
        # Canopy Boundary layer resistance
        rb_sun = rb
        # minimum canopy conductance
        go_shd = go
        # Canopy Boundary layer resistance
        rb_shd = rb

        gmes_sun = gmes
        gmes_shd = gmes

        PAR_sun_local = PAR_sun_local / (LAI * Fsun)
        PAR_shd_local = PAR_shd_local / (LAI * Fshd)

        # sunlit fraction
        if Fsun > 0:
            Ci_sun = fsolve(self.CO2_Concentration,Citm1_sun_local,args=(PAR_sun_local,Ca,ra,rb_sun,Ts,Pre,Ds,Psi_L,Psi_sto_50,
                            Psi_sto_99,CT,Vmax_sun,DS,Ha,FI,Oa,Do,a1,go_sun,gmes_sun,rjv),xtol=1)

            if Ci_sun.size == 1:
                Ci_sun = Ci_sun[0]
            else:
                print('Size of Ci_shd is greater than 1')

            self.Photosynthesis_Biochemical(Ci_sun,PAR_sun_local,Ca,ra,rb_sun,Ts,Pre,Ds,Psi_L,Psi_sto_50,Psi_sto_99,CT,
                                            Vmax_sun,DS,Ha,FI,Oa,Do,a1,go_sun,gmes_sun,rjv)
            CiF_sun = self.PhotoBiochem.CcF
            An_sun = self.PhotoBiochem.An
            rc_sun = self.PhotoBiochem.rs
            Rdark_sun = self.PhotoBiochem.Rdark
            SIF_sun = self.PhotoBiochem.F755nm
        else:
            Ci_sun = 0
            CiF_sun = 0
            An_sun = 0
            Rdark_sun = 0
            rc_sun = numpy.inf
            SIF_sun = 0

        # Shadowed fraction
        if Fshd > 0:
            Ci_shd = fsolve(self.CO2_Concentration,Citm1_shd_local,args=(PAR_shd_local,Ca,ra,rb_shd,Ts,Pre,Ds,Psi_L,Psi_sto_50,
                                                                         Psi_sto_99,CT,Vmax_shd,DS,Ha,FI,Oa,Do,a1,go_shd,gmes_shd,rjv),xtol=1)
            if Ci_shd.size == 1:
                Ci_shd = Ci_shd[0]
            else:
                print('Size of Ci_shd is greater than 1')

            self.Photosynthesis_Biochemical(Ci_shd,PAR_shd_local,Ca,ra,rb_shd,Ts,Pre,Ds,Psi_L,Psi_sto_50,Psi_sto_99,CT,
                                            Vmax_shd,DS,Ha,FI,Oa,Do,a1,go_shd,gmes_shd,rjv)
            CiF_shd = self.PhotoBiochem.CcF
            An_shd = self.PhotoBiochem.An
            rc_shd = self.PhotoBiochem.rs
            Rdark_shd = self.PhotoBiochem.Rdark
            SIF_shd = self.PhotoBiochem.F755nm
        else:
            Ci_shd = 0
            CiF_shd = 0
            An_shd = 0
            Rdark_shd = 0
            rc_shd = numpy.inf
            SIF_shd = 0

        DCi_sun = Ci_sun - CiF_sun
        DCi_shd = Ci_shd - CiF_shd
        DCi = (DCi_sun + DCi_shd) / 2

        An = An_sun * (LAI * Fsun) + An_shd * (LAI * Fshd)
        Rdark = Rdark_sun * (LAI * Fsun) + Rdark_shd * (LAI * Fshd)
        SIF = SIF_sun * (LAI * Fsun) + SIF_shd * (LAI * Fshd)

        # stomatal resistence [s m^-1]
        rs_sun = rc_sun
        rs_shd = rc_shd

        lanp = 0.469 # [J umol^-1 CO2]
        Lpho = (An + Rdark) * lanp # [W m^-2]

        self.CanopyResEvl = CanopyResistanceAnEvolution()
        self.CanopyResEvl.rs_sun = rs_sun
        self.CanopyResEvl.rs_shd = rs_shd
        self.CanopyResEvl.Ci_sun = Ci_sun
        self.CanopyResEvl.Ci_shd = Ci_shd
        self.CanopyResEvl.An = An
        self.CanopyResEvl.Rdark = Rdark
        self.CanopyResEvl.Lpho = Lpho
        self.CanopyResEvl.SIF = SIF
        self.CanopyResEvl.DCi = DCi

        return rs_sun,rs_shd,Ci_sun,Ci_shd,An,Rdark,Lpho,SIF,DCi

    def CO2_Concentration(self,Cc,IPAR,Csl,ra,rb,Ts,Pre,Ds,Psi_L,Psi_sto_50,Psi_sto_99,CT,Vmax,DS,Ha,FI,Oa,Do,a1,go,gmes,rjv):

        self.Photosynthesis_Biochemical(Cc,IPAR,Csl,ra,rb,Ts,Pre,Ds,Psi_L,Psi_sto_50,Psi_sto_99,CT,Vmax,DS,Ha,FI,Oa,Do,a1,go,gmes,rjv)

        CcF = self.PhotoBiochem.CcF
        self.DCi = Cc - CcF

        return self.DCi

    def Photosynthesis_Biochemical(self,Cc,IPAR,Csl,ra,rb,Ts,Pre,Ds,Psi_L,Psi_sto_50,Psi_sto_00,CT,Vmax,DS,Ha,FI,Oa,Do,
                                   a1,go,gmes,rjv):

        # Re-define input parameters which are overwritten in this function
        Cc_local = copy.copy(Cc)
        IPAR_local = copy.copy(IPAR)
        Csl_local = copy.copy(Csl)
        ra_local = copy.copy(ra)
        rb_local = copy.copy(rb)
        Pre_local = copy.copy(Pre)
        DS_local = copy.copy(DS)
        Ha_local = copy.copy(Ha)
        Oa_local = copy.copy(Oa)
        go_local = copy.copy(go)

        Ta = Ts
        Pre0 = 101325 # [Pa]
        Tf = 273.15   # [K]

        # Conversion factors
        Pre_local = Pre_local * 100     # [Pa]
        IPAR_local = IPAR_local * 4.57  # [umolPhotons s^-1 m^-2]
        ra_local = ra_local * (0.0224 * (Ta + 273.15) * Pre0 / (Tf * Pre_local)) * 10 ** (-6) # [m^2 s umolH2O^-1]
        rb_local = rb_local * (0.0224 * (Ta + 273.15) * Pre0 / (Tf * Pre_local)) * 10 ** (-6) # [m^2 s umolH2O^-1]
        Cc_local = Cc_local * 10 ** (-6) * Pre_local   # Partial Pressure [Pa * molCO2 molAIR^-1]
        Oa_local = Oa_local * 10 ** (-6) * Pre_local   # [Pa]
        Csl_local = Csl_local * 10 ** (-6) * Pre_local # Leaf surface CO2 concentration [Pa]

        # Mesophyl Conductance [s m^2 umolCO2^-1]
        rmes = 1 / (1e+6 * gmes)
        go_local = go_local * 10 ** 6 # [umolCO2 s^-1 m^-2]

        # Temperature dependence
        # Maximum Rubisco Capacity Vm
        Ts_k = Ts + 273.15 # [K]
        # Reference Temperature [K]
        Tref = 25 + 273.15
        # Gas Constant [kJ K^-1 mol^-1]
        R = 0.008314

        # ANS_TEMP is assumed t be one
        # Deactivation Energy [kJ mol^-1]
        Hd = 200

        # Mix of Temperature Function and High Temperature Inhibition
        kT = numpy.exp(Ha_local * (Ts_k - Tref)/(Tref * R * Ts_k)) * (1+numpy.exp((Tref*DS_local-Hd) / (Tref*R))) / \
             (1+numpy.exp((Ts_k * DS_local - Hd)/(Ts_k * R)))
        Vm = Vmax * kT # [umolCO2 s^-1 m^-2]

        # Maximum Electron Transport Rate Jm
        # Deactivation Energy [kJ mol^-1]
        Hd = 200
        # Activation Energy [kJ mol^-1]
        Ha_local = 50
        # entropy factor [kJ mol^-1 K^-1]
        DS_local = 0.646

        kT = numpy.exp(Ha_local * (Ts_k - Tref) / (Tref * R * Ts_k)) * (1 + numpy.exp((Tref * DS_local - Hd) / (Tref * R))) / \
             (1 + numpy.exp((Ts_k * DS_local - Hd) / (Ts_k * R)))
        # [umol electrons s^-1 m^-2]
        Jmax = Vmax * rjv
        Jm = Jmax * kT # [umol electrons s^-1 m^-2]

        # Triose Phosphate Utilization
        # Activation Energy [kJ mol^-1]
        Ha_local = 53.1
        # entropy factor [kJ mol^-1 K^-1]
        DS_local = 0.490
        # Deactivation Energy [kJ mol^-1]
        Hd = 150.65

        TPU25 = 0.1182 * Vmax # [umolCO2 s^-1 m^-2]
        kT = numpy.exp(Ha_local * (Ts_k - Tref) / (Tref * R * Ts_k)) * (1 + numpy.exp((Tref * DS_local - Hd) / (Tref * R))) / \
             (1 + numpy.exp((Ts_k * DS_local - Hd) / (Ts_k * R)))
        TPU = TPU25 * kT # [umolCO2 s^-1 m^-2]

        if CT == 4:
            s1 = 0.3  # [1 K^-1]
            s3 = 0.2  # [1 K^-1] 0.3 (Cox 2001)
            Tup = 40  # [C]
            Tlow = 15 # [C]

            # Temperature Function 1 for Maximum Rubisco Capacity
            f1T = 1 / (1 + numpy.exp(s1 * (Ts - Tup)))
            # Temperature Function 2 for Maximum Rubisco Capacity
            f2T = 1 / (1 + numpy.exp(s3 * (Tlow - Ts)))
            fT = 2**(0.1 * (Ts - 25))
            Vm = Vmax * fT * f1T * f2T   # [umolCO2 s^-1 m^-2]

            ke25 = 20000 * Vmax
            ke = ke25 * fT

        # CO2 concentration point
        # ANSG is assumed to be 2
        # Activation Energy [kJ mol^-1]
        Ha_local = 37.83
        kT = numpy.exp(Ha_local * (Ts_k - Tref) / (Tref * R * Ts_k))
        GAM25 = 42.75 # [umol mol^-1]
        GAM25 = GAM25 * 10 ** (-6) * Pre_local # [Pa]
        # Michaelis - Menten Constant for CO2 [Pa]
        GAM = GAM25 * kT

        if CT == 3:
            # Michaelis-Menten Constants for CO2 and O2
            # Activation Energy [kJ mol^-1]
            Ha_local = 79.43
            Kc25 = 404.9 # [umol mol^-1]
            Kc25 = Kc25 * 10 ** (-6) * Pre_local # [Pa]
            kT = numpy.exp(Ha_local * (Ts_k - Tref) / (Tref * R * Ts_k))
            Kc = Kc25 * kT

            # Activation Energy [kJ mol^-1]
            Ha_local = 36.38
            Ko25 = 278.4 # [umol mol^-1]
            Ko25 = Ko25 * 10 ** (-3) * Pre_local # [Pa]
            kT = numpy.exp(Ha_local * (Ts_k - Tref) / (Tref * R * Ts_k))
            # Michaelis-Menten Constant for O2
            Ko = Ko25 * kT

        # Dark Respiration
        if CT == 3:
            Ha_local = 46.39
            DS_local = 0.490
            Hd = 150.65

            Rdark25 = 0.015 * Vmax
            kT = numpy.exp(Ha_local * (Ts_k - Tref) / (Tref * R * Ts_k)) * (1 + numpy.exp((Tref * DS_local - Hd) / (Tref * R))) / \
                 (1 + numpy.exp((Ts_k * DS_local - Hd) / (Ts_k * R)))
            Rdark = Rdark25 * kT

        elif CT == 4:
            fT = 2.0 ** (0.1 * (Ts - 25))
            # Temperature Function 3 for Respiration
            fT3 = 1 / (1 + numpy.exp(1.3 * (Ts - 55)))
            Rdark25 = 0.025 * Vmax
            # Leaf Maintainance Respiration / Dark Respiration [umolCO2 s^-1 m^-2]
            Rdark = Rdark25 * fT * fT3

        # Photosynthesis factors
        # Light Absorbed by Photosystem II in CO2 units [umolCO2 s^-1 m^-2]
        Q = FI * IPAR_local
        d1 = 0.7
        d2 = -(Q + Jm / 4)
        d3 = Q * Jm / 4

        # Electron Transport Rate
        J = min((-d2 + numpy.sqrt(d2 ** 2 - 4 * d1 * d3)) / (2 * d1), (-d2 - numpy.sqrt(d2 ** 2 - 4 * d1 * d3)) / (2 * d1))

        if CT == 3:
            # Gross Assimilation Rate Limited by Rubisco [umolCO2 s^-1 m^-2]
            JC = Vm * (Cc_local - GAM) / (Cc_local + Kc * (1 + Oa_local / Ko))
            # Light Limited
            # Gross Assimilation Rate Limited by Light [umolCO2 s^-1 m^-2]
            JL = J * (Cc_local - GAM) / (Cc_local + 2 * GAM)
            # Capacity of the leaf to export or utilize the products of photosynthesis
            # Gross Assimilation Rate Limited by Export [umolCO2 s^-1 m^-2]
            JE = 3 * TPU
        elif CT == 4:
            # Rubisco Limited
            # Gross Assimilation Rate Limited by Rubisco [umolCO2 s^-1 m^-2]
            JC = Vm
            # Light Limited
            # Gross Assimilation Rate Limited by Light [umolCO2 s^-1 m^-2]
            JL = Q
            # PEP Carboxylase Limited
            JE = ke * Cc_local / Pre_local

        # First Polynomium
        if CT == 3:
            b1 = 0.98
            b2 = -(JC + JL)
            b3 = JC * JL
        elif CT == 4:
            b1 = 0.80
            b2 = -(JC + JL)
            b3 = JC * JL

        # Smoothed Minimum between JC and JE [umolCO2 s^-1 m^-2]
        JP = min((-b2 + numpy.sqrt(b2 ** 2 - 4 * b1 * b3)) / (2 * b1), (-b2 - numpy.sqrt(b2 ** 2 - 4 * b1 * b3)) / (2 * b1))

        # Second Polynomium
        if CT == 3:
            c1 = 0.95
            c2 = -(JP + JE)
            c3 = JP * JE
        elif CT == 4:
            c1 = 0.95
            c2 = -(JP + JE)
            c3 = JP * JE
        # Gross Assimilation Rate Potential [umolCO2 s^-1 m^-2]
        A = min((-c2 + numpy.sqrt(c2 ** 2 - 4 * c1 * c3)) / (2 * c1), (-c2 - numpy.sqrt(c2 ** 2 - 4 * c1 * c3)) / (2 * c1))

        # New Water Stress Function
        Rgsws = 0.02
        p2 = math.log((1 - Rgsws) / Rgsws) / (Psi_sto_00 - Psi_sto_50) # [MPa^-1]
        q2 = -p2 * Psi_sto_50 # [-]
        Rgsw = 1 / (1 + numpy.exp(p2 * Psi_L + q2))
        fO = (1 - Rgsw)
        if fO > 1:
            fO = 1
        if fO < 0:
            fO = 0

        # Solar-induced chlorophyll fluorescence (SIF)
        # Je is the actual electron transport rate calculated from the CO2 exchange data
        if CT == 3:
            Jfe = A * (Cc_local + 2 * GAM) / (Cc_local - GAM)
        elif CT == 4:
            Jfe = A

        fiP0= FI*4 # [umol Electrons umolPhotons^-1]
        fiP = fiP0 * Jfe / Q # [0.4 max - stress decrease ]
        # degree of light saturation
        dls = 1 - fiP / fiP0

        kf = 0.05
        kd = max(0.03 * Ts + 0.0773, 0.087)
        kn = (6.2473 * dls - 0.5944) * dls

        fiF = kf / (kf + kd + kn) * (1 - fiP) # [umol Electrons umolPhotons^-1]
        SIF = IPAR_local * fiF # [umol electrons s^-1 m^-2]

        # k theoretically a function of Vmax and Chlorophyll content
        k = 0.0375 * Vmax + 8.25 # [umol m^-2 s^-1 / W m^-2 sr^-1 um^-1]
        F755nm = SIF / k # [W m^-2 sr^-1 um^-1]

        # Gross Assimilation Rate [umolCO2 s^-1 m^-2]
        A = A * fO
        # Net Assimilation Rate [umolCO2 s^-1 m^-2]
        An = A - Rdark

        # Stomatal Conductance
        gsCO2 = go_local + a1 * An * Pre_local / ((Cc_local - GAM) * (1 + Ds / Do))
        if gsCO2 < go_local:
            gsCO2 = go_local

        # Stomatal resistance or Canopy [s m^2 umolCO2^-1]
        rsCO2 = 1 / gsCO2

        CcF = Csl_local - An * Pre_local * (rsCO2 + rmes + 1.37 * rb_local + ra_local) # [Pa]
        if CcF < 0:
            CcF = 0

        # Stomatal resistance or canopy [s m^2 molH2O^-1]
        rsH20 = (rsCO2 / 1.64) * (10 ** 6)
        # Net Assimilation Rate [umolCO2 s^-1 m^-2]
        An = (Csl_local - CcF) / (Pre_local * (rsCO2 + rmes + 1.37 * rb_local + ra_local))

        CcF = CcF / (Pre_local * 10 ** (-6)) # [umolCO2 molAIR^-1]
        rs = rsH20 * (Tf * Pre_local) / (0.0224 * (Ts + 273.15) * Pre0) # Stomatal resistance or Canopy [s m^-1]

        self.PhotoBiochem = PhotosynthesisBiochemical()
        self.PhotoBiochem.CcF = CcF
        self.PhotoBiochem.An = An
        self.PhotoBiochem.rs = rs
        self.PhotoBiochem.Rdark = Rdark
        self.PhotoBiochem.F755nm = F755nm
        self.PhotoBiochem.GAM = GAM
        self.PhotoBiochem.gsCO2 = gsCO2

    def Soil_Resistance(self,T_soil,Pre,Ws,ea,q_runon,O,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,SPAR):

        """
        ------
        INPUT:
        T_soil: Soil temperature [C]
        Pre: Pressure [mbar]
        Ws: Wind speed [m s^-1]
        ea: Vapor pressure at T_canyon [Pa]
        q_runon: Intercepted water on the surface [mm]
        O: Water Content []
        Ks: Hydraulic conductivity at saturation [mm s^-1]
        Osat: Water content at saturation, saturation moisture 0 kPa  [-]
        Ohy: Hygroscopic Moisture Evaporation cessation []
        L: Slope of logarithmic tension-moisture curve [-]
        Pe: Tension at air antry (bubbling pressure) [kPa]
        O33: Soil water content at -33 [kPa] of water potential
        alpVG: Alpha parameter Van-Genuchten soil water retention curve [mm^-1]
        nVG: n parameter Van-Genuchten soil water retention curve [mm^-1]
        SPAR: Soil parameter type
        -------
        OUTPUT:
        r_soil: Soil resistance [s m^-1]
        b_soil: beta factor [0-1]
        alp_soil: relative humidity [0-1]
        """

        # -------------------------
        # Calculate soil resistance
        # -------------------------
        # Soil Temperature [K]
        Ts_k = T_soil + 273.15
        # Water density [kg m^-3]
        row = 1000
        g = 9.81
        # water vapor gas constant [J kg^-1 K^-1]
        Rd = 461.5
        # vapor molecular diffusivity [m^2 s^-1]
        Da = (2.11 * 1e-5) * (((Ts_k) / 273.15)**1.94) * (Pre * 100 / 101325)
        esat = 611 * numpy.exp(17.27 * T_soil / (237.3 + T_soil)) # [Pa]

        SoilCal = Soil_Calculations()
        SoilCal.Conductivity_Suction(SPAR,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,O)
        Ko = SoilCal.CondSuc.Ko # [mm s^-1]
        Po = SoilCal.CondSuc.Po # [mm]
        if Po < 0:
            Po = 0
        alp_soil = numpy.exp(-Po * g / (1000 * Rd * Ts_k))
        # 40-200 um  Size of the pores --  Particle Size/3 [m]
        Psz = (11.12 * nVG**3.286) * 1e-6
        # Boundary Layer Thickness
        dm = 2.26 * 1e-3 / numpy.sqrt(Ws)
        ### it is only for ANSW = 4
        gammap = (alp_soil * esat - ea) / (row * Rd * Ts_k) # [-]
        if gammap < 0:
            r_soil = 0
        else:
            # Internal soil viscous resistance [s m^-1]
            rsv = gammap / (4 * Ko / (1000 * 3600))
            f_O = (2 / numpy.pi) * (numpy.sqrt(numpy.pi / (4 * O)) - 1) / numpy.sqrt(4 * O) # [-]
            # viscous boundary layer resistance [s m^-1]
            rvbl = (dm + Psz * f_O) / Da
            r_soil = rvbl + rsv

        if O <= Ohy:
            r_soil = numpy.inf
        b_soil = 1

        if q_runon > 0:
            r_soil = 0
            alp_soil = 1
            b_soil = 1

        self.ResSoil = ResistanceSoil()
        self.ResSoil.r_soil = r_soil
        self.ResSoil.alp_soil = alp_soil
        self.ResSoil.b_soil = b_soil

        return r_soil,b_soil,alp_soil

    def Leaf_Boundary_Resistance(self,Ws,Ts,Ta,hc,d_leaf,LAI,zatm,disp_h,zom):

        # Re-define input parameters which are overwritten in this function
        d_leaf_local = copy.copy(d_leaf)

        # Wind speed [m s^-1]
        u = Ws

        d_leaf_local = d_leaf_local / 100  # [m]
        # von Karman constant
        k = 0.4
        # Empirical coefficient [m s^-0.5]
        a = 0.01
        # Zero plane displacement [m]
        d = disp_h
        # Domain height [m]
        z = zatm
        # Hypothesis Logarithmic distribution of wind speed
        # Friction Velocity  [m s^-1]
        us = k * u / math.log((z - d) / zom)
        # Wind Speed top Canopy [m s^-1]
        u_hc = (us / k) * math.log((hc - d) / zom)
        # Attenuation Coefficient
        alpha_den = (z / hc - 1)
        alpha = math.log(u / u_hc) / alpha_den
        alpha = 0.5 * alpha * LAI / 2

        # Expression of Leaf Boundary Layer Resistance
        gb = (2 * a / alpha) * ((u_hc / d_leaf_local) ** 0.5) * (1 - numpy.exp(-alpha / 2)) \
            if d_leaf_local != 0 else (2 * a / alpha) * (numpy.inf ** 0.5) * (1 - numpy.exp(-alpha / 2))

        # Expression for free convection
        # Molecular diffusivity of heat [m^2 s^-1]
        Dh = 1.9e-5
        # Grashof number [-]
        if Ts > Ta:
            Gr = 1.6e8 * (Ts - Ta) * d_leaf_local ** 3
        else:
            Gr = 0
        # The leaf boundary conductance at free convection [m s^-1]
        gb_free = 0.5 * Dh * (Gr ** 0.25) / d_leaf_local if d_leaf_local != 0 else numpy.nan

        gb = gb + gb_free

        # Leaf Boundary Layer Resistance [s m^-1] one-sided for unit leaf
        rb = 1 / gb

        return rb

    def Leaf_BR(self,u_hc,Ts,Ta,d_leaf,alpha):

        """
        ------
        INPUT:
        u_hc: wind speed at the height of trees [m s^-1]
        Ts: Trees temperature [C]
        Ta: Air temperature at the height of trees [C]
        d_leaf: Leaf dimension of trees[cm]
        alpha: Attenuation Coefficient [-]
        -------
        OUTPUT:
        rb: Leaf boundary layer resistance [s m^-1]
        """

        # Re-define input parameters which are overwritten in this function
        d_leaf_local = copy.copy(d_leaf)

        d_leaf_local = d_leaf_local/100 # [m]
        a = 0.01                        # [m s^-0.5] (Chodhury and Monteith 1988)

        # Expression for Leaf Boundary Layer Resistance [m s^-1]
        gb = (2 * a / alpha) * ((u_hc / d_leaf_local)**0.5) * (1 - numpy.exp(-alpha / 2))

        # Expression for free convection  (Leuning 1995, Monteith 1973)
        Dh = 1.9 * 1e-5 # [m^2 s^-1]
        Gr = 1.6 * 1e+8 * (Ts - Ta) * (d_leaf_local ** 3)* (Ts > Ta) # [-]
        gb_free = 0.5 * Dh * Gr ** (0.25) / d_leaf_local             # [m s^-1]
        gb = gb + gb_free

        # Leaf Boundary Layer Resistance [s m^-1] one-sided for unit leaf
        rb = 1 / gb

        return rb

    def Urban_roughness(self,hc_H,hc_L,Csoil,Croad,Croof):

        """
        ------
        INPUT:
        hc_H: Height of high vegetation [m]
        hc_L: Height of low vegetation [m]
        Csoil: boolean operator for presence (1) and absence (0) of soil
        Croad: boolean operator for presence (1) and absence (0) of road
        Croof: boolean operator for presence (1) and absence (0) of roof
        -------
        OUTPUT:
        zom: roughness eddy diffusivities for momentum [m]
        zoh: roughness eddy diffusivities for heat [m]
        disp_h: maximum displacement height [m]
        zom_H: high vegetation roughness momentum [m]
        zom_L: low vegetation roughness momentum [m]
        zoh_H: high vegetation roughness heat [m]
        zoh_L: low vegetation roughness heat [m]
        d_H: displacement height of high vegetation [m]
        d_L: displacement height of low vegetation [m]
        zom_other: roughness momentum for the other urban surfaces [m]
        """

        # bare soil roughness momentum [m]
        if Csoil == 1:
            zom_soil = 0.003
        else:
            zom_soil = 0.0

        # road roughness momentum
        if Croad == 1:
            zom_road = 0.003
        else:
            zom_road = 0.0

        # roof roughness momentum Wang et al. (2013) [m]
        if Croof == 1:
            zom_roof = 0.01
        else:
            zom_roof = 0.0

        # vegetation roughness momentum [m] Brutsaert (1975) high vegetation
        zom_H = 0.123 * hc_H
        # vegetation roughness momentum [m] Brutsaert (1975) low vegetation
        zom_L = 0.123 * hc_L

        zom_other = [zom_soil, zom_road, zom_roof]
        # roughness eddy diffusivities for momentum [m]
        zom_other = max(zom_other)

        # Heat Roughness [m]
        zoh_L = zom_L * 0.1
        zoh_H = zom_H * 0.1
        # roughness  eddy diffusivities for heat  [m] [Brutsaert (1975)]
        zoh_other = 0.1 * zom_other

        # PATCH SCALE ROUGHNESS
        zom = max(max(zom_H, zom_L), zom_other)
        zoh = max(max(zoh_H, zoh_L), zoh_other)
        zom_ground = max(zom_L, zom_other)
        zoh_ground = max(zoh_L, zoh_other)

        # Displacement height
        d_L = 2 / 3 * hc_L
        d_H = 2 / 3 * hc_H
        disp_h = max(d_H, d_L)

        return zom,zoh,zom_ground,zoh_ground,disp_h,zom_H,zom_L,zoh_H,zoh_L,d_H,d_L,zom_other

    def WindProfile_Canyon(self,Hcan,Htree,R_tree,Wcan,Wroof,Kopt,LAI_t,Zatm,WindSpeed_top,Zp,trees,Zref_und,zom_und):

        """
        ------
        INPUT:
        Hcan: canyon height [m]
        Htree: Tree height [m]
        R_tree: Tree radius [m]
        Wcan: Canyon width [m]
        Wroof: Roof width [m]
        Kopt: Optical transmission factor [-]
        LAI_t: Leaf area index of tree [-]
        Zatm: Height of the domain [m]
        WindSpeed_top: wind speed at the top of the domain [m s^-1]
        Zp: Height of interest within the canyon [m]
        trees: Presence of trees [0: No, 1: Yes]
        Zref_und: Refrence height [m]
        zom_und: Aerodynamic roughness length [m]
        -------
        OUTPUT:
        dcan: Urban displacement height including trees [m]
        zomcan:	Urban momentum roughness height including trees [m]
        u_Hcan:	Wind speed at canyon height [m s^-1]
        u_Zpcan: Wind speed within canyon at height Zpcan [m s^-1]
        w_Zpcan: Vertical wind speed within canyon [m s^-1]
        """

        # Re-define input parameters which are overwritten in this function
        Htree_local = copy.copy(Htree)
        R_tree_local = copy.copy(R_tree)
        LAI_t_local = copy.copy(LAI_t)

        if trees == 0:
            Htree_local = 0
            R_tree_local = 0
            LAI_t_local = 0

        # Best fit for staggered arrays , or a = 3.59 for square arrays
        a = 4.43
        # Von Karman constant
        k = 0.4
        # b=1, no incorporation for drag correction factors. Good fit for staggered arrays
        b = 1.0
        # nominal drag for cubical obstacles
        CDb = 1.2

        # Plan area fraction of buildings and vegetation
        Ap_build = Wroof
        Ap_tree = 4 * R_tree_local
        Ap_urb = Wcan + Wroof

        # Frontal area fraction of vegetation and buildings: assumption infinite urban canyon perpendicular to the
        # wind direction (Length of building and plot equals infinity)
        Af_build_s = Hcan
        Af_veg_s = 2 * R_tree_local

        # Tree canopy transmittance (optical = P2D)
        P2D = numpy.exp(-Kopt * LAI_t_local)
        # Guan et al. 2003
        P3D = P2D ** 0.40
        # Guan et al. 2000
        Pv = (-1.251 * P3D ** 2 + + 0.489 * P3D + 0.803) / CDb

        # Plan area fraction of buildings and Calculation of structural parameters and wind profile in the city
        Lp_tot = (Ap_build + (1 - P3D) * Ap_tree) / Ap_urb
        H_tot = (Hcan * Ap_build + (Htree_local + R_tree_local) * (1 - P3D) * Ap_tree) / (Ap_build + (1 - P3D) * Ap_tree)

        # Urban displacement height and roughness length with incorporation of trees (Kent 2017), (MacDonald 1998)
        # displacement height of canyon [m], eq. 23
        dcan = H_tot * (1 + a ** (-Lp_tot) * (Lp_tot - 1))

        Af_build = H_tot / (H_tot - dcan) * Af_build_s
        Af_veg = H_tot / (H_tot - dcan) * Af_veg_s

        zomcan = H_tot * (1 - dcan / H_tot) * numpy.exp( -(1 / k ** 2 * 0.5 * b * CDb * (1 - dcan / H_tot) *
                                                           (Af_build + Pv * Af_veg) / Ap_urb) ** (-0.5))
        zohcan = zomcan / 10

        # Calculation of wind profile above and in the canyon with a logarithmic and exponential wind profile.
        # Friction Velocity Atmosphere [m s^-1]
        us_atm = k * WindSpeed_top / math.log((Zatm - dcan) / zomcan)
        self.Ustar_Atm = us_atm
        # Wind Speed at canyon top [m s^-1]
        u_Hcan = (us_atm / k) * math.log((Hcan - dcan) / zomcan)
        # Attenuation Coefficient Canyon not corrected for presence of trees.
        alpha = math.log(WindSpeed_top / u_Hcan) / (Zatm / Hcan - 1)

        if Zp >= Hcan:
            u_Zp = (us_atm / k) * math.log((Zp - dcan) / zomcan)
            w_Zp = 0
        elif Zp <= Hcan and Zp >= Zref_und:
            u_Zp = u_Hcan * numpy.exp(-alpha * (1 - Zp / Hcan))
            w_Zp = 0
        elif Zp <= Zref_und and Zp >= zom_und:
            uref_und = u_Hcan * numpy.exp(-alpha * (1 - Zref_und / Hcan))
            usref_und = k * uref_und / math.log(Zref_und / zom_und)
            u_Zp = (usref_und / k) * math.log(Zp / zom_und)
            w_Zp = 0
        else:
            u_Zp = 0
            w_Zp = 0
            print('wind speed calculation height higher than reference height or lower than roughness length')

        return  dcan,zomcan,u_Hcan,u_Zp,w_Zp,alpha

    def WindProfile_Roof(self,Hcan,hveg,VerticalProfUrban,Geometry_m):

        # Wind profile from 1-D model
        vx = copy.copy(VerticalProfUrban.vx)
        vy = copy.copy(VerticalProfUrban.vy)
        z_urban = copy.copy(Geometry_m.z[:-1])

        vx_intp = interp1d(z_urban, vx)
        vy_intp = interp1d(z_urban, vy)

        u_Zp = numpy.sqrt(vx_intp(Geometry_m.dz/2+Hcan)**2+vy_intp(Geometry_m.dz/2+Hcan)**2)
        u_Hveg = numpy.sqrt(vx_intp(hveg+Hcan)**2+vy_intp(hveg+Hcan)**2)

        return u_Zp,u_Hveg

    def Wall_Aerodynamic_Resistance(self,VerticalProfUrban,Geometry_m,windMin,Cp,iz_wall,ParCalculation):

        """
        ------
        INPUT:
        VerticalProfUrban: Vertical profile of variables obtained from 1-D model
        Geometry_m: Geometric parameters
        windMin: Minimum wind speed in the urban area [m s^-1]
        Cp: Air specific heat  [J kg^-1 K^-1]
        iz_wall: z index in the urban canyon
        ParCalculation: General calculation parameters
        -------
        OUTPUT:
        RES: Aerodynamic resistance [s m^-1]
        """

        # Air density [kg m^-3]
        rho = copy.copy(VerticalProfUrban.rho[iz_wall])

        vett = numpy.sqrt(VerticalProfUrban.vx[iz_wall] ** 2 + VerticalProfUrban.vy[iz_wall] ** 2)
        vett = max(vett, windMin)

        # Convective heat transfer coefficient [W K^-1 m^-2]
        hc = 5.678 * (1.09 + 0.23 * (vett / 0.3048))
        # Using energy balance for a control volume inside the urban unit, the convective heat transfer coefficient should be limited
        # hc must be less than (rho * cp / dt) * [(1-lambdap) * Hmean / (4 * lambdaf * dz)]
        if hc > ((rho*Cp/ParCalculation.dts) * ((1-Geometry_m.lambdap)*Geometry_m.Height_canyon) / (4*Geometry_m.lambdaf*Geometry_m.dz)):
            hc = (rho*Cp/ParCalculation.dts) * ((1-Geometry_m.lambdap)*Geometry_m.Height_canyon) / (4*Geometry_m.lambdaf*Geometry_m.dz)
        # Term in energy equation [s m^-1]
        RES = ((rho * Cp)/hc)

        return RES

    def Ground_Aerodynamic_Resistance_1D(self,WindSpeed_top,Zatm,VerticalProfUrban,Gemeotry_m,Ta,Ts,hcan,dcan,zomcan,zom_und,Ztree,Rtree,ColParam):
        """
        ------
        INPUT:
        WindSpeed_top: Wind speed at the top of the domain [m s^-1]
        Zatm: Domain height [m]
        VerticalProfUrban: Vertical profile of variables obtained from 1-D model
        Gemeotry_m: Geometric parameters
        Ta: Air temperature near the ground [K]
        Ts: Total ground temperature [K]
        hcan: Canyon height [m]
        dcan: Displacement height of the canyon [m]
        zomcan: Aerodynamic roughness length of the canyon [m]
        zom_und: Aerodynamic roughness length of the ground [m]
        Ztree: Trees height [m]
        Rtree: Trees radius [m]
        ColParam: 1-D model parameters
        -------
        OUTPUT:
        rap_can: Aerodynamic resistance near the ground [s m^-1]
        rap_Ztree_In: Aerodynamic resistance between trees and canyon air [s m^-1]
        u_Hcan: Wind speed at the canyon height [m s^-1]
        alpha: Attenuation Coefficient [-]
        Ri: Bulk Richardson number
        Utot:
        """

        # Interpolate wind speed
        vx = copy.copy(VerticalProfUrban.vx)
        vy = copy.copy(VerticalProfUrban.vy)
        z_urban = copy.copy(Gemeotry_m.z[:-1])

        vx_intp = interp1d(z_urban, vx)
        vy_intp = interp1d(z_urban, vy)

        # Wind speed at the canyon height [m s^-1]
        u_Hcan = numpy.sqrt(vx_intp(hcan)**2+vy_intp(hcan)**2)

        # Make sure that the wind speed at the top of the domain is different from the wind speed at the canyon height
        # (if they are equal, then Attenuation Coefficient will be zero)
        wind_top = WindSpeed_top if WindSpeed_top != u_Hcan else u_Hcan+0.1
        # Attenuation Coefficient Canyon not corrected for presence of trees
        alpha = math.log(wind_top / u_Hcan) / (Zatm / hcan - 1)

        Ck = 0.4

        #---------------------------------------
        # Aerodynamic resistance near the ground
        #---------------------------------------
        zz = Gemeotry_m.dz/2
        Utot = numpy.sqrt(vx[0]**2+vy[0]**2)
        Utot = max(Utot,ColParam.WindMin_Urban)
        # Compute bulk Richardson number
        # Ta and Ts should be in [K]
        Ri = 2 * 9.81 * zz * (Ta - Ts) / ((Ta + Ts) * (Utot ** 2))
        if Ri > 0.16:
            Ri = 0.16

        # Calculation from Louis, 1979 (eq. 11 and 12)
        b = 9.4
        cm = 7.4
        ch = 5.3
        R = 0.74
        a = Ck / math.log(zz / zom_und)
        if Ri > 0:
            fm = 1 / ((1 + 0.5 * b * Ri) ** 2)
            fh = fm
        else:
            c = b * cm * a * a * (zz / zom_und) ** 0.5
            fm = 1 - b * Ri / (1 + c * (-Ri) ** 0.5)
            c = c * ch / cm
            fh = 1 - b * Ri / (1 + c * (-Ri) ** 0.5)

        rap_can = R / ((a ** 2) * Utot * fh)

        # --------------------------------------
        # Aerodynamic resistance above the trees
        # --------------------------------------
        zz = dcan + zomcan
        Utot = numpy.sqrt(vx_intp(zz) ** 2 + vy_intp(zz) ** 2)
        # Calculation from Louis, 1979 (eq. 11 and 12)
        b = 9.4
        cm = 7.4
        ch = 5.3
        R = 0.74
        a = Ck / math.log(zz / zom_und)
        if Ri > 0:
            fm = 1 / ((1 + 0.5 * b * Ri) ** 2)
            fh = fm
        else:
            c = b * cm * a * a * (zz / zom_und) ** 0.5
            fm = 1 - b * Ri / (1 + c * (-Ri) ** 0.5)
            c = c * ch / cm
            fh = 1 - b * Ri / (1 + c * (-Ri) ** 0.5)

        rap_can_AboveTree = R / ((a ** 2) * Utot * fh)

        # ------------------------------------------------
        # Aerodynamic resistance just underneath the trees
        # ------------------------------------------------
        zz = Ztree-Rtree
        Utot = numpy.sqrt(vx_intp(zz) ** 2 + vy_intp(zz) ** 2)
        # Calculation from Louis, 1979 (eq. 11 and 12)
        b = 9.4
        cm = 7.4
        ch = 5.3
        R = 0.74
        a = Ck / math.log(zz / zom_und)
        if Ri > 0:
            fm = 1 / ((1 + 0.5 * b * Ri) ** 2)
            fh = fm
        else:
            c = b * cm * a * a * (zz / zom_und) ** 0.5
            fm = 1 - b * Ri / (1 + c * (-Ri) ** 0.5)
            c = c * ch / cm
            fh = 1 - b * Ri / (1 + c * (-Ri) ** 0.5)

        rap_Ztree = R / ((a ** 2) * Utot * fh)

        # Aerodynamic resistance between trees and canyon air
        rap_Ztree_In = max(rap_can_AboveTree - rap_Ztree, 0)


        return rap_can,rap_Ztree_In,u_Hcan,alpha,Ri

    def Roof_Aerodynamic_Resistance_1D(self,VerticalProfUrban,Geometry_m,z0,Ts):

        vx = copy.copy(VerticalProfUrban.vx)
        vy = copy.copy(VerticalProfUrban.vy)
        th = copy.copy(VerticalProfUrban.th)
        z_urban = copy.copy(Geometry_m.z[:-1])

        vx_intp = interp1d(z_urban, vx)
        vy_intp = interp1d(z_urban, vy)
        th_intp = interp1d(z_urban, th)

        vx_air = vx_intp(Geometry_m.dz/2 + Geometry_m.Height_canyon)
        vy_air = vy_intp(Geometry_m.dz/2 + Geometry_m.Height_canyon)
        th_air = th_intp(Geometry_m.dz/2 + Geometry_m.Height_canyon)

        Ck = 0.4

        zz = Geometry_m.dz/2
        Utot = numpy.sqrt(vx_air**2 + vy_air**2)

        # Compute bulk Richardson number
        # th_air and Ts should be in [K]
        Ri = 2 * 9.81 * zz * (th_air - Ts) / ((th_air + Ts) * (Utot ** 2))

        # Calculation from Louis, 1979 (eq. 11 and 12)
        b = 9.4
        cm = 7.4
        ch = 5.3
        R = 0.74
        a = Ck / math.log(zz / z0)
        if Ri > 0:
            fm = 1 / ((1 + 0.5 * b * Ri) ** 2)
            fh = fm
        else:
            c = b * cm * a * a * (zz / z0) ** 0.5
            fm = 1 - b * Ri / (1 + c * (-Ri) ** 0.5)
            c = c * ch / cm
            fh = 1 - b * Ri / (1 + c * (-Ri) ** 0.5)

        ra = R / ((a ** 2) * Utot * fh)

        return ra























