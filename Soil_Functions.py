import os
import numpy
import math
from SoilFunctionsDef import SoilParametersTotal,ConductivitySuction,SoilThermalProperties,SoilParameters,SoilMoistureConductivityUpdate,Infiltration2Def,LeakageBottomDef,SoilWaterMultiLayerDef,VolumeCorrectionDef
import copy

'''
Soil Functions:
Developed by Mohsen Moradi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: June 2021
'''

class Soil_Calculations(object):

    def Conductivity_Suction(self,SPAR,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,O):

        """
        ------
        INPUT:
        SPAR: Soil parameter type (1:VanGenuchten, 2:Saxton-Rawls)
        Ks: Hydraulic conductivity at saturation for each soil layer [mm s^-1]
        Osat: Water content at saturation, saturation moisture 0 kPa  [-]
        Ohy: Hygroscopic Moisture Evaporation cessation 10000 kPa - 10 MPa - 1000 m [-]
        L: Slope of logarithmic tension-moisture curve [-]
        Pe: Tension at air entry (bubbling pressure) [kPa]
        O33: 33 kPa Moisture [-]
        alpVG: Alpha parameter Van-Genuchten soil water retention curve [mm^-1]
        nVG: n parameter Van-Genuchten soil water retention curve [mm^-1]
        O: Soil moisture/water content in the different soil layers [-]
        -------
        OUTPUT:
        Ko: Hydraulic conductivity at O [mm s^-1]
        Po: Tension at O [mm]
        """

        if SPAR == 1:
            # [-]
            Se = (O - Ohy)/ (Osat - Ohy)
            # [mm]
            mVG = 1 - 1 / nVG
            # Tension at O [mm]
            Po = -(1 / alpVG) * ((Se) ** (-1. / mVG) - 1) ** (1. / nVG)
            # Hydraulic conductivity at O [mm s^-1]
            Ko = Ks * ((Se) ** (0.5)) * (1 - (1 - (Se) ** (1. / mVG)) ** mVG) ** 2
        else:
            # Specific weight water[N m^-3]
            gw = 9810
            B = 1 / L
            # Coefficient of moisture tension
            if O33 == 0:
                A = 0
            else:
                A = numpy.exp(math.log(33) + B * math.log(O33))

            # Equation 17 in Saxton and Rawls, 2006. [mm s^-1]
            Ko = Ks * (O / Osat) ** (3 + (2 / L))

            IFCond = numpy.all(O<O33)
            if IFCond:
                # [kPa]
                Psi = A * (O**(-B))
            else:
                # [kPa]
                Psi = 33 - ((O - O33) * (33 - Pe) / (Osat - O33))

            # Tension at O [mm]
            Po = 1000 * 1000 * Psi / (gw)

        self.CondSuc = ConductivitySuction()
        self.CondSuc.Ko = Ko  # [mm s^-1]
        self.CondSuc.Po = Po  # [mm]

        return Ko,Po

    def Infiltration_2(self,Osat,Ohy,L,alpVG,nVG,Pe,Ks_Zs,O33,SPAR,O,Dz,WIS,cosalp,Pond,dts):

        """
        ------
        INPUT:
        Osat: Water content at saturation, saturation moisture 0 kPa  [-]
        Ohy: Hygroscopic Moisture Evaporation cessation 10000 kPa - 10 MPa - 1000 m [-]
        L: Slope of logarithmic tension-moisture curve [-]
        alpVG: Alpha parameter Van-Genuchten soil water retention curve [mm^-1]
        nVG: n parameter Van-Genuchten soil water retention curve [mm^-1]
        Pe: Tension at air entry (bubbling pressure) [kPa]
        Ks_Zs: Hydraulic conductivity at saturation for each soil layer [mm s^-1]
        O33: 33 kPa Moisture [-]
        SPAR: Soil parameter type (1:VanGenuchten, 2:Saxton-Rawls)
        O: Soil moisture/water content in the different soil layers [-]
        Dz: Distance from surface to half-layer [mm]
        WIS: Water Incoming to Soil Layer [mm s^-1]
        cosalp: model constant
        Pond: Interception/ponding from the previous time step as an indicator if there is ponding infiltration [mm s^-1]
        dts: Time step [s]
        -------
        OUTPUT:
        f: Infiltration rate [mm s^-1]
        fpot: Potential Infiltration rate [mm s^-1]
        """


        if SPAR ==1:
            PAE = 0
        else:
            gw = 9810
            # Suction at air entry (bubbling pressure) [mm]
            PAE = -1000*1000*Pe/(gw)

        # [mm]
        P0 = max(Pond*dts,0.0)


        self.Conductivity_Suction(SPAR,Ks_Zs,Osat,Ohy,L,Pe,O33,alpVG,nVG,O)
        Ko = self.CondSuc.Ko  # [mm s^-1]
        Po = self.CondSuc.Po  # [mm]

        # Water Potential [mm]
        P = -Po
        Khalf = 0.5 * (Ko + Ks_Zs)
        # Potential Infiltration rate [mm s^-1]
        fpot = Khalf*(1*cosalp-(-PAE+(P-P0))/Dz)
        # Infiltration rate [mm s^-1]
        f = max(0, min(fpot, WIS))

        self.Infiltration2 = Infiltration2Def()
        self.Infiltration2.f = f
        self.Infiltration2.fpot = fpot

        return f,fpot

    def Leakage_Bottom(self,O,Ks_Zs,Osat,Ohy,L,nVG,Kbot,ms,SPAR):

        """
        ------
        INPUT:
        O: Soil moisture in the different soil layers [-]
        Ks_Zs: Hydraulic conductivity at saturation for each soil layer [mm s^-1]
        Osat: Water content at saturation, saturation moisture 0 kPa  [-]
        Ohy: Hygroscopic Moisture Evaporation cessation 10000 kPa - 10 MPa - 1000 m [-]
        L: Slope of logarithmic tension-moisture curve
        nVG: n parameter Van-Genuchten soil water retention curve [mm^-1]
        Kbot: Conductivity at the bedrock layer [mm s^-1]
        ms: Number of soil layers [-]
        SPAR: Soil parameter type (1:VanGenuchten, 2:Saxton-Rawls)
        -------
        OUTPUT:
        Lk: Leakage [mm s^-1]
        """

        if numpy.isnan(Kbot):
            if SPAR == 1:
                Se = (O[ms-1] - Ohy[ms-1]) / (Osat[ms-1] - Ohy[ms-1])
                mVG = 1 - 1 / nVG[ms]
                # Estimation Conductivity last layer [mm s^-1]
                Ko = Ks_Zs[ms-1] * ((Se)**(0.5)) * (1 - (1 - (Se)**(1 / mVG))**mVG)**2
            else:
                # Estimation Conductivity last layer [mm s^-1]
                Ko = Ks_Zs[ms-1] * (O[ms-1]/ Osat[ms-1])**(3 + (2 / L[ms-1]))
            Lk = Ko
        else:
            if O[ms-1] > Osat[ms-1]-1e-5:
                Lk = Kbot
            else:
                Lk = 0

        self.LeakageBottom = LeakageBottomDef()
        self.LeakageBottom.Lk = Lk

        return Lk

    def Root_Fraction_General(self,Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L):

        # RfH_Zs : Root Fraction for High Vegetation [1...m] [%]
        # RfL_Zs : Root Fraction for Low  Vegetation [1...m] [%]

        n = len(Zs) - 1
        RfH_Zs = numpy.zeros(n)
        RfL_Zs = numpy.zeros(n)

        if numpy.isnan(ZRmax_H):
            if ZR95_H > Zs[n] or ZR95_L > Zs[n] or ZRmax_L > Zs[n]:
                print('ERROR: LAST LAYER TOO SHALLOW FOR ACCOMODATING ROOTS')
        elif numpy.isnan(ZRmax_L):
            if ZR95_H > Zs[n] or ZR95_L > Zs[n] or ZRmax_H > Zs[n]:
                print('ERROR: LAST LAYER TOO SHALLOW FOR ACCOMODATING ROOTS')
        elif numpy.isnan(ZRmax_H) and numpy.isnan(ZRmax_L):
            if ZR95_H > Zs[n] or ZR95_L > Zs[n]:
                print('ERROR: LAST LAYER TOO SHALLOW FOR ACCOMODATING ROOTS')
        else:
            if ZR95_H > Zs[n] or ZR95_L > Zs[n] or ZRmax_H > Zs[n] or ZRmax_L > Zs[n]:
                print('ERROR: LAST LAYER TOO SHALLOW FOR ACCOMODATING ROOTS')

        if CASE_ROOT ==1: # Exponential Profile
            # Shape of Root Distribution [mm^-1]
            eta_H = 3. / ZR95_H if ZR95_H != 0 else numpy.inf
            eta_L = 3. / ZR95_L if ZR95_L != 0 else numpy.inf

            i = 0
            if ZR95_H != 0:
                while i <= n-1: ### n-1 is right???
                    if ZR95_H > Zs[i+1]:
                        RfH_Zs[i] = numpy.exp(-eta_H * Zs[i]) - numpy.exp(-eta_H * Zs[i + 1])
                    else:
                        RfH_Zs[i] = numpy.exp(-eta_H * Zs[i]) - numpy.exp(-eta_H * ZR95_H)
                        i = n-1
                    i = i+1
            i = 0
            if ZR95_L != 0:
                while i <= n-1:
                    if ZR95_L > Zs[i+1]:
                        RfL_Zs[i] = numpy.exp(-eta_L * Zs[i]) - numpy.exp(-eta_L * Zs[i+1])
                    else:
                        RfL_Zs[i] = numpy.exp(-eta_L * Zs[i]) - numpy.exp(-eta_L * ZR95_L)
                        i = n-1
                    i = i+1

            # Water Content for the Root
            Rto1 = 0.9502
            # Root Proportion in the Layer
            RfH_Zs[:] = RfH_Zs[:] / Rto1
            RfL_Zs[:] = RfL_Zs[:] / Rto1

        # Linear Dose Response
        # Schenk and Jackson 2002, Collins and Bras 2007
        if CASE_ROOT ==2:
            c_H = 2.94 / math.log(ZR50_H / ZR95_H)
            c_L = 2.94 / math.log(ZR50_L / ZR95_L)

            i = 0
            if ZR95_H != 0:
                while i <= n-1:
                    if ZR95_H > Zs[i+1]:
                        RfH_Zs[i] = 1 / (1 + (Zs[i+1] / ZR50_H)**c_H) - 1 / (1+(Zs[i]/ZR50_H)**c_H)
                    else:
                        RfH_Zs[i] = 1 / (1 + (ZR95_H / ZR50_H)**c_H) - 1 / (1+(Zs[i]/ZR50_H)**c_H)
                        i = n-1
                    i = i+1
            i = 0
            if ZR95_L != 0:
                while i <= n-1:
                    if ZR95_L > Zs[i+1]:
                        RfL_Zs[i] = 1 / (1 + (Zs[i+1] / ZR50_L)**c_L) - 1 / (1+(Zs[i]/ZR50_L)**c_L)
                    else:
                        RfL_Zs[i] = 1 / (1 + (ZR95_L / ZR50_L)**c_L) - 1 / (1+(Zs[i]/ZR50_L)**c_L)
                        i = n-1
                    i = i+1
            Rto1 = 0.9498
            # Root Proportion in the Layer
            RfH_Zs[:] = RfH_Zs[:] / Rto1
            RfL_Zs[:] = RfL_Zs[:] / Rto1

        # Constant Profile
        if CASE_ROOT == 3:

            i = 0
            if ZR95_H != 0:
                while i <= n-1:
                    if ZR95_H > Zs[i+1]:
                        RfH_Zs[i] = (Zs[i+1]-Zs[i])/ZR95_H
                    else:
                        RfH_Zs[i] = (ZR95_H - Zs[i]) / ZR95_H
                        i = n-1
                    i = i+1
            i = 1
            if ZR95_L != 0:
                while i <= n-1:
                    if ZR95_L > Zs[i+1]:
                        RfL_Zs[i] = (Zs[i+1]-Zs[i])/ZR95_L
                    else:
                        RfL_Zs[i] = (ZR95_L - Zs[i]) / ZR95_L
                        i = n-1
                    i = i+1

        # Deep (Tap) Root Profile
        if CASE_ROOT == 4:
            c_H = 2.94 / math.log(ZR50_H / ZR95_H)
            c_L = 2.94 / math.log(ZR50_L / ZR95_L)

            i = 0
            if ZR95_H != 0:
                while i <= n-1:
                    if ZR95_H > Zs[i + 1]:
                        RfH_Zs[i] = 1 / (1 + (Zs[i+1] / ZR50_H)**c_H) - 1/(1+(Zs[i]/ZR50_H)**c_H)
                    elif ZR95_H <= Zs[i+1] and ZR95_H > Zs[i]:
                        RfH_Zs[i] = 1 / (1 + (ZR95_H / ZR50_H)**c_H) - 1/(1+(Zs[i]/ZR50_H)**c_H)
                        if ZRmax_H[j] <= Zs[i+1]:
                            RfH_Zs[i] = RfH_Zs[i] + 0.0502*(ZRmax_H- ZR95_H) / (ZRmax_H - ZR95_H)
                            i = n-1
                        else:
                            RfH_Zs[i] = RfH_Zs[i] + 0.0502*(Zs[i+1] - ZR95_H) / (ZRmax_H - ZR95_H)
                    elif ZRmax_H > Zs[i+1]:
                        RfH_Zs[i] = 0.0502 * (Zs[i + 1] - Zs[i]) / (ZRmax_H - ZR95_H)
                    else:
                        RfH_Zs[i] = 0.0502 * (ZRmax_H - Zs[i]) / (ZRmax_H - ZR95_H)
                        i = n-1
                    i = i+1
            i = 1
            if ZR95_L != 0:
                while i <= n-1:
                    if ZR95_L > Zs[i+1]:
                        RfL_Zs[i] =  1/(1 + (Zs[i+1]/ZR50_L)**c_L) -  1/(1 + (Zs[i]/ZR50_L)**c_L)
                    elif ZR95_L <= Zs[i+1] and ZR95_L > Zs[i]:
                        RfL_Zs[i] = 1/(1 + (ZR95_L / ZR50_L)**c_L) - 1/(1 + (Zs[i]/ZR50_L)**c_L)
                        if ZRmax_L <= Zs[i+1]:
                            RfL_Zs[i] = RfL_Zs[i] + 0.0502*(ZRmax_L - ZR95_L) / (ZRmax_L - ZR95_L)
                            i = n-1
                        else:
                            RfL_Zs[i] = RfL_Zs[i] + 0.0502*(Zs[i + 1] - ZR95_L) / (ZRmax_L - ZR95_L)
                    elif ZRmax_L > Zs[i+1]:
                        RfL_Zs[i] = 0.0502 * (Zs[i + 1] - Zs[i]) / (ZRmax_L - ZR95_L)
                    else:
                        RfL_Zs[i] = 0.0502 * (ZRmax_L - Zs[i]) / (ZRmax_L - ZR95_L)
                        i = n-1
                    i = i+1

        return RfH_Zs,RfL_Zs

    def Root_Soil_Conductance(self,Ks,Rl,rcyl,rroot,Zr):

        # Re-define input parameters which are overwritten in this function
        Zr_local = copy.copy(Zr)
        Ks_local = copy.copy(Ks)

        # Rooting depth [m]
        Zr_local = Zr_local / 1000
        # water density [kg m^-3]
        row = 1000
        # gravity acceleration [m s^-2]
        g = 9.81
        # water density [mmolH20 kg^-1]
        rho = 55555

        # [m MPa^-1]
        CF = 10**6 / (row * g)
        # [m s^-1]
        Ks_local = Ks_local / 1000
        OPT = 3
        if OPT == 1:
            # [s^-1]
            gsr = Ks_local * numpy.sqrt(Rl / (2 * Zr_local))
            # [m s^-1 MPa^-1]
            gsr = gsr * CF
            # [mmol H20 m^-2 ground s MPa]
            Ksr = gsr * row * rho

        if OPT == 2:
            # Radial conductivity or root [s^-1]
            Kr = 5 * 1e-08 / CF
            # Radial thickness of the rhizosphere [m]
            Lrs = rcyl
            # [s^-1]
            gsr = numpy.sqrt(Ks_local * Kr / Lrs)
            # [mmol H20 m^-2 ground s MPa]
            Ksr = gsr * CF * row * rho

        if OPT == 3:
            # [s^-1]
            gsr = Ks_local * Rl * (2 * numpy.pi) / (math.log(rcyl / rroot))
            Ksr = gsr * CF * row * rho # [mmol H20 m^-2 ground s MPa]

        return Ksr

    def Soil_Parameters(self,Psan,Pcla,Porg):
        # Note that Pcla+Psan+Porg must be less than 1
        Psil = 1 - Psan - Pcla - Porg
        if Psil < 0:
            print('SOIL PERCENTAGE INPUTS INCONSISTENT')
        # Weight fraction of gravel [g Gravel g^-1 bulk soil]
        Rw = 0
        # Density factor
        DF = 1
        # Fraction of organic material in the soil in percent ### ???
        Porg_perc = Porg * 100
        O1500t = -0.024 * Psan + 0.487 * Pcla + 0.006 * Porg_perc + 0.005 * (Psan * Porg_perc) - 0.013 * (
                Pcla * Porg_perc) + 0.068 * (Psan * Pcla) + 0.031
        O1500 = O1500t + 0.14 * O1500t - 0.02  # 1500 kPa Moisture
        O33t = -0.251 * Psan + 0.195 * Pcla + 0.011 * Porg_perc + 0.006 * (Psan * Porg_perc) - 0.027 * (
                Pcla * Porg_perc) + 0.452 * (Psan * Pcla) + 0.299
        O33 = O33t + (1.283 * O33t ** 2 - 0.374 * O33t - 0.015)  # 33 kPa Moisture
        Os_33t = 0.278 * Psan + 0.034 * Pcla + 0.022 * Porg_perc - 0.018 * (Psan * Porg_perc) - 0.027 * (
                Pcla * Porg_perc) - 0.584 * (Psan * Pcla) + 0.078
        Os_33 = Os_33t + (0.636 * Os_33t - 0.107)  # SAT - 33 kPa Moisture

        # Coefficient of moisture tension
        B = (math.log(1500) - math.log(33)) / (math.log(O33) - math.log(O1500))
        # Coefficient of moisture tension
        A = numpy.exp(math.log(33) + B * math.log(O33))
        # Slope of logarithmic tension-moisture curve [-]
        L = 1 / B

        # Saturation moisture 0 kPa []
        Osat = O33 + Os_33 - 0.097 * Psan + 0.043
        # Normal density dry soil  [kg m^-3]  (Equation 6 in Saxton and Rawls, 2006)
        rsd = (1 - Osat) * 2650

        # Adjust density density dry soil [kg m^-3]
        rsd_df = rsd * DF
        # Saturation moisture 0 kPa [%]
        Osat_df = 1 - (rsd_df / 2650)
        # 33 kPa Moisture [%]
        O33_df = O33 - 0.2 * (Osat - Osat_df)
        # SAT-33 kPa Moisture [%]
        Os_33_df = Osat_df - O33_df

        Osat = Osat_df
        O33 = O33_df
        Os_33 = Os_33_df
        rsd = rsd_df

        # Saturation conductivity [mm s^-1]
        Ks = (1930 * (Osat - O33) ** (3 - L))/3600
        # Tension at air entry, first solution, [kPa]
        Pet = -21.67 * Psan - 27.93 * Pcla - 81.97 * Os_33 + 71.12 * (Psan * Os_33) + 8.29 * (Pcla * Os_33) + 14.05 * (
                Psan * Pcla) + 27.16
        # Tension at air entry (bubbling pressure) [kPa]
        Pe = Pet + (0.02 * Pet ** 2 - 0.113 * Pet - 0.70)
        if Pe < 0.5:
            Pe = 0.5

        # Gravel effects
        # matric soil density / gravel density []
        alpha = rsd / 2650
        # Volume fraction of gravel [g cm^-3]
        Rv = (alpha * Rw) / (1 - Rw * (1 - alpha))
        # Bulk soil density (matric plus gravel), [g cm^-3]
        rhoB = (rsd / 1000) * (1 - Rv) + (Rv * 2.65)
        Osatb = Osat * (1 - Rv)
        O33b = O33 * (1 - Rv)
        # Saturated conductivity bulk soil [mm s^-1]
        Kb = Ks * (1 - Rw) / (1 - Rw * (1 - 3 * alpha / 2))

        rsd = rhoB * 1000  # [kg m^-3]
        Osat = Osatb  # [-]
        O33 = O33b  # [-]
        Ks = Kb  # [mm s^-1]

        # -------------------------------
        # Thermal characteristics of soil
        # -------------------------------
        # Normal density dry soil for thermal properties computing [kg m^-3]
        rsd_s = 2700 * (1 - Osat)
        # Thermal conductivity dry soil [W m^-1 K^-1]
        lan_dry = (0.135 * rsd_s + 64.7) / (2700 - 0.947 * rsd_s)
        # Thermal conductivity soil solid [W m^-1 K^-1]
        lan_s = (8.8 * Psan + 2.92 * Pcla) / (Psan + Pcla)
        # Volumetric heat capacity soil solid [J m^-3 K^-1]
        cv_s = 1e+6 * (2.128 * Psan + 2.385 * Pcla) / (Psan + Pcla)

        # -----------------
        # K USLE Parameter
        # -----------------
        # Organic carbon
        Porg_c = Porg / 1.72
        fsand = (0.2 + 0.3 * numpy.exp(-25.6 * Psan * (1 - Psil)))
        fcli = (Psil / (Pcla + Psil)) ** 0.3
        forg = (1 - (0.25 * Porg_c) / (Porg_c + numpy.exp(3.72 - 2.95 * Porg_c)))
        fhisand = (1 - (0.7 * (1 - Psan)) / ((1 - Psan) + numpy.exp(-5.51 - 22.9 * (1 - Psan))))
        # K_Usle [ton h MJ^-1 mm^-1]
        K_usle = fsand * fcli * forg * fhisand
        # Erosivity factor [kg h J^-1 mm^-1]
        K_usle = K_usle / 1000

        self.SoilParam = SoilParameters()
        self.SoilParam.Osat = Osat
        self.SoilParam.L = L
        self.SoilParam.Pe = Pe
        self.SoilParam.Ks = Ks
        self.SoilParam.O33 = O33
        self.SoilParam.rsd = rsd
        self.SoilParam.lan_dry = lan_dry
        self.SoilParam.lan_s = lan_s
        self.SoilParam.cv_s = cv_s
        self.SoilParam.K_usle = K_usle

        return Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle

    def Soil_ParametersII(self,ms,Osat,L,Pe,Ks,O33,nVG,alpVG,Kfc,Pss,Pwp,Phy,Ohy,SPAR):

        # nVG : n parameter Van-Genuchten soil water retention curve [mm^-1]
        # L : Slope of logaritimc tension-moisture curve
        # Kfc : Conductivity at field capacity [mm h^-1]
        # Ks : saturation conductivity [mm s^-1]
        # Ofc : Field Capacity Moisture [] Kns < 0.2 mm/h
        # Oss : Stomatal closure begin moisture 30 kPa - 0.03 MPa - 3 m []
        # Owp : Stomatal closure end moisture -Wilting point 3000 kPa - 3 MPa -300 m []
        # Ohy : Hygroscopic Moisture Evaporation cessation 10000 kPa - 10 MPa - 1000 m []

        # Re-define input parameters which are overwritten in this function
        Pss_local = copy.copy(Pss)
        Pwp_local = copy.copy(Pwp)
        Ohy_local = copy.copy(Ohy)


        if SPAR == 1:
            mVG = 1 - 1 / nVG # [mm]
            Pss_local = -101.9368 * Pss_local # [mm]
            Se = 1 / ((1 + abs(alpVG * Pss_local)**nVG)**mVG)
            Se[Se > 1] = 1
            O = Ohy_local + (Osat - Ohy_local) * Se
            Oss = O

            Pwp_local = -101.9368 * Pwp_local # [mm]
            Se = 1 / ((1 + abs(alpVG * Pwp_local)**nVG)**mVG)
            Se[Se > 1] = 1
            O = Ohy_local + (Osat - Ohy_local) * Se
            Owp = O

            Ofc = numpy.zeros(ms)
            Ohy_local = numpy.zeros(ms)

        else:
            Ofc = numpy.zeros(ms)
            Oss = numpy.zeros(ms)
            Owp = numpy.zeros(ms)
            Ohy_local = numpy.zeros(ms)
            for i in range(0,ms):
                B = 1/L[i]
                # Coefficient of moisture tension
                A = numpy.exp(math.log(33) + B * math.log(O33[i]))

                if Pss_local < 33:
                    Oss[i] = O33[i] + (33 - Pss_local) * (Osat[i] - O33[i]) / (33 - Pe[i])
                else:
                    Oss[i] = (Pss_local/A)**(-1/B)

                if Pwp_local < 33:
                    Owp[i] = O33[i] + (33 - Pwp_local) * (Osat[i] - O33[i]) / (33 - Pe[i])
                else:
                    Owp[i] = (Pwp_local / A)**(-1 / B)

                if Phy < 33:
                    Ohy_local[i] = O33[i] + (33 - Phy) * (Osat[i] - O33[i]) / (33 - Pe[i])
                else:
                    Ohy_local[i] = (Phy / A)**(-1 / B)

                Ofc[i] = Osat[i] * ((Kfc/3600) / Ks[i])**(1 / (3 + (2 / L[i])))

        return Ofc,Oss,Owp,Ohy_local

    def Soil_Thermal_Properties(self,Tdp,rsd,lan_dry,lan_s,cv_s,Osat,Ohy,O):
        """
        ------
        INPUT:
        Tdp: Dampening temperature [C]
        rsd: Dry soil density [kg m^-3]
        lan_dry: Dry soil thermal conductivity [W m^-1 K^-1]
        lan_s: Solid soil thermal conductivity [W m^-1 K^-1]
        cv_s: Solid soil volumertic heat capacity [J m^-3 K^-1]
        Osat: Water content at saturation [-]
        Ohy: Residual / Hygroscopic / Wilting Point water content [-]
        O: Soil moisture, soil water content at previous time step [-]
        -------
        OUTPUT:
        CTt: Soil total thermal capacity [K m^2 J^-1]
        lanS: Soil thermal conductivity [W m^-1 K^-1]
        cv_Soil: Soil volumetric heat capacity [J m^-3 K^-1]
        """

        # Re-define input parameters which are overwritten in this function
        O_local = copy.copy(O)

        # Water properties
        # Density of water [kg m^-3]
        row = 1000
        # Thermal conductivity of water [W m^-1 K^-1]
        lan_wat = 0.58
        # Thermal conductivity of ice [W m^-1 K^-1]
        lan_ice = 2.29
        # Volumetric heat capacity of water [J m^-3 K^-1]
        cv_w = 4186000

        n = numpy.size(O_local)

        # This if condition is used, because size of O will change
        if n == 1:
            # Thermal conductivity Soil [W m^-1 K^-1]
            lanS = 0
            # Volumetric heat capacity Soil  [J m^-3 K^-1]
            cv_Soil = 0
            # Soil density [kg m^-3]
            rsoil = 0
            # Specific Heat of Soil [J kg^-1 K^-1]
            cs_Soil = 0

            if O_local < Ohy[0]:
                O_local = Ohy[0]
            if O_local > Osat[0]:
                O_local = Osat[0]
            # Frozen layer
            Oice = [Osat[i] if Tdp < 0 else 0 for i in range(len(Osat))]

            # Each soil layer
            if Tdp > 0:
                lan_sat = (lan_wat ** Osat[0]) * (lan_s[0] ** (1 - Osat[0]))
                Ke = math.log((O_local + Oice[0]) / Osat[0]) + 1
                Ke = Ke*int(Ke >= 0)

            else:
                # Liquid water content at saturation
                Oliq = Osat[0] - Oice[0]
                # Saturated Conductivity [W m^-1 K^-1]
                lan_sat = (lan_wat ** Osat[0]) * (lan_s[0] ** (1 - Osat[0])) * (lan_ice ** (Osat[0] - Oliq))
                Ke = (O_local + Oice[0]) / Osat[0]

            if O_local/Osat[0] > 10**(-7):
                 # Thermal conductivity Soil [W m^-1 K^-1]
                lanS = Ke * lan_sat + (1 - Ke) * lan_dry[0]
            else:
                # Thermal conductivity Soil [W m^-1 K^-1]
                lanS = lan_dry[0]

            # Volumetric heat capacity Soil  [J m^-3 K^-1]
            cv_Soil = cv_s[0] * (1 - Osat[0]) + O_local * cv_w
            # Soil Density [kg m^-3]
            rsoil = rsd[0] + (O_local - Ohy[0]) * row
            # Specific Heat Soil [J kg^-1 K^-1]
            cs_Soil = cv_Soil / rsoil

            # time constant [s]
            tau = 86400
            # Total Thermal Capacity Soil [K m^2 J^-1]
            CTt = 2 * (numpy.sqrt(numpy.pi / (lanS * cs_Soil * rsoil * tau)))

        else:
            # Thermal conductivity Soil [W m^-1 K^-1]
            lanS = numpy.zeros(n)
            # Volumetric heat capacity Soil  [J m^-3 K^-1]
            cv_Soil = numpy.zeros(n)
            # Soil density [kg m^-3]
            rsoil = numpy.zeros(n)
            # Specific Heat of Soil [J kg^-1 K^-1]
            cs_Soil = numpy.zeros(n)

            for i in range(0,n):
                if O_local[i] < Ohy[i]:
                    O_local[i] = Ohy[i]
            for i in range(0,n):
                if O_local[i] > Osat[i]:
                    O_local[i] = Osat[i]

            # Frozen layer
            Oice = Osat*int(Tdp < 0)

            # Each soil layer
            for i in range(n):
                if Tdp > 0:
                    # Saturated Conductivity [W m^-1 K^-1]
                    lan_sat = (lan_wat ** Osat[i]) * (lan_s[i] ** (1 - Osat[i]))
                    # Kersten number
                    Ke = math.log((O_local[i] + Oice[i]) / Osat[i]) + 1
                    Ke = Ke * int(Ke >= 0)

                else:
                    # Liquid water content at saturation
                    Oliq = Osat[i] - Oice[i]
                    # Saturated Conductivity [W m^-1 K^-1]
                    lan_sat = (lan_wat ** Osat[i]) * (lan_s[i] ** (1 - Osat[i])) * (lan_ice ** (Osat[i] - Oliq))
                    Ke = (O_local[i] + Oice[i]) / Osat[i]

                if O_local[i]/Osat[i] > 10**(-7):
                    # Thermal conductivity Soil [W m^-1 K^-1]
                    lanS[i] = Ke * lan_sat + (1 - Ke) * lan_dry[i]
                else:
                    # Thermal conductivity Soil [W m^-1 K^-1]
                    lanS[i] = lan_dry[i]

                # Volumetric heat capacity Soil  [J m^-3 K^-1]
                cv_Soil[i] = cv_s[i] * (1 - Osat[i]) + O_local[i] * cv_w
                # Soil Density [kg m^-3]
                rsoil[i] = rsd[i] + (O_local[i] - Ohy[i]) * row
                # Specific Heat Soil [J kg^-1 K^-1]
                cs_Soil[i] = cv_Soil[i] / rsoil[i]

            # time constant [s]
            tau = 86400
            # Total Thermal Capacity Soil [K m^2 J^-1]
            CTt = 2 * (numpy.sqrt(numpy.pi / (lanS[0] * cs_Soil[0] * rsoil[0] * tau)))


        self.SoilThProp = SoilThermalProperties()
        self.SoilThProp.lanS = lanS
        self.SoilThProp.cv_Soil = cv_Soil
        self.SoilThProp.CTt = CTt

        return lanS,cv_Soil,CTt

    def Soil_Water_MultiLayer(self,V,Zs,dz,n,Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,SPAR,EvL_Zs,Inf_Zs,RfH_Zs,RfL_Zs,Rrootl_H,
                              Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L):
        """
        ------
        INPUT:
        V: Water content in each soil layer [mm]
        Zs: Soil layer discretization [mm]
        dz: Thickness of soil layers [mm]
        n: Number of soil layers [-]
        Osat: Water content at saturation, saturation moisture 0 kPa  [-]
        Ohy:  MoistHygroscopicure Evaporation cessation 10000 kPa - 10 MPa - 1000 m [-]
        nVG: n parameter Van-Genuchten soil water retention curve [mm^-1]
        alpVG: Alpha parameter Van-Genuchten soil water retention curve [mm^-1]
        Ks_Zs: Hydraulic conductivity at saturation for each soil layer [mm s^-1]
        L: Slope of logarithmic tension-moisture curve
        Pe: Tension at air entry (bubbling pressure) [kPa]
        O33: 33 kPa Moisture [-]
        SPAR: Soil parameter type (1:VanGenuchten, 2:Saxton-Rawls)
        EvL_Zs: Fraction of evaporation depth in a specific soil layer [-]
        Inf_Zs: Fraction of infiltration depth in a specific soil layer [-]
        RfH_Zs: Root Fraction for High Vegetation [1...m] [%]
        RfL_Zs: Root Fraction for Low Vegetation [1...m] [%]
        Rrootl_H: Root length index of high vegetation [m root m^-2 PFT]
        Rrootl_L: Root length index of low vegetation [m root m^-2 PFT]
        PsiL50_H: Water Potential at 50% loss conductivity of high vegetation [MPa]
        PsiL50_L: Water Potential at 50% loss conductivity of low vegetation [MPa]
        PsiX50_H: Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil of high vegetation [MPa]
        PsiX50_L: Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil of low vegetation [MPa]
        ------
        OUTPUT:
        O: Water Content [1...m] Layer [-]
        ZWT: [mm]
        OF: Water Content for infiltration [-]
        OS: Water Content for evaporation [-]
        Psi_s_H: Soil water potential for first layer of vegetation [MPa]
        Psi_s_L: Soil Water Potential  for Second Layer of Vegetation [MPa]
        gsr_H: Root soil conductance [mmol H20 m^-2 ground s^-1 MPa^-1]
        gsr_L: Root soil conductance [mmol H20 m^-2 ground s^-1 MPa^-1]
        Exwat_H: Maximum extractable water [mm m^2 m^-2 ground s^-1]
        Exwat_L: Maximum extractable water [mm m^2 m^-2 ground s^-1]
        Rd: Water table rise that reaches the surface level (Dunne Runoff) [mm]
        WTR: Water table rise [mm]
        POT: Soil water potential [mm]
        OH: Water content for first layer of vegetation [-]
        OL: Water content for second layer of vegetation [-]
        """

        # Re-define input parameter which is overwritten in this function
        n_local = copy.copy(n)

        # Water Content [-] (n Layer)
        O = numpy.ones(n_local)
        # Water Table Rise [mm]
        WTR = numpy.zeros(n_local)

        for i in range(n_local-1,-1,-1):
            if i == n_local-1:
                O[i] = (V[i] / dz[i]) + Ohy[i]
                WTR[i] = (O[i] - Osat[i]) * dz[i] * int(O[i] > Osat[i])

            else:
                O[i] = (V[i] + WTR[i + 1]) / dz[i] + Ohy[i]
                WTR[i] = (O[i] - Osat[i]) * dz[i] * int(O[i] > Osat[i])

            if O[i] < Ohy[i]:
                O[i] = Ohy[i]
            if O[i] > Osat[i]:
                O[i] = Osat[i]

        i = n_local-1
        while O[i] > Osat[i]-1e-5 and i > -1:
            i = i-1
        ZWT = Zs[i + 1]
        PHead = [(Zs[j] + dz[j]/2)-ZWT for j in range(0,len(dz))]
        for i in range(0,len(PHead)):
            if PHead[i] < 0:
                PHead[i] = 0

        # Compute First Potential
        if SPAR == 1:
            # Hydraulic Head
            Se = [(O[j] - Ohy[j]) / (Osat[j] - Ohy[j]) for j in range(0,len(O))] # [-]
            mVG = [1 - 1 / nVG[j] for j in range(0,len(nVG))] # [mm]
            POT = [PHead[j] + (1 / alpVG[j]) * ((Se[j])**(-1 / mVG[j]) - 1)**(1 / nVG[j]) for j in range(0,len(nVG))]
        else:
            # Hydraulic Head
            Ptem = numpy.zeros(n_local)
            for jk in range(0,n_local):
                CondSuc = self.Conductivity_Suction(2,Ks_Zs[jk],Osat[jk],Ohy[jk],L[jk],Pe[jk],O33[jk],alpVG[jk],nVG[jk],O[jk])
                Ptem[jk] = CondSuc[1] # [mm]
            POT = [PHead[j]-Ptem[j] for j in range(0,len(nVG))] # [mm]

        if i == n_local-1:
            ZWT = Zs[i + 1]
        else:
            if i == 0:
                ZWT = 0
            else:
                # Find Water Table Depth
                CZ = [Zs[j]+dz[j]/2 for j in range(0,len(dz))] # [mm]
                ZWT = CZ[i + 1] + (CZ[i] - CZ[i + 1]) * (-POT[i + 1]) / (POT[i] - POT[i + 1]) # [mm]
                # Recompute Potential
                PHead = [Zs[j]+dz[j]/2-ZWT for j in range(0,len(dz))]
                for i in range(0, len(PHead)):
                    if PHead[i] < 0:
                        PHead[i] = 0

        # Recompute Potential
        if SPAR == 1:
            # Hydraulic Head
            Se = [(O[j]-Ohy[j])/(Osat[j]-Ohy[j]) for j in range(0,len(O))] # [-]
            mVG = [1-1/nVG[j] for j in range(0,len(nVG))] # [mm]
            POT = [PHead[j] + (1/alpVG[j])*((Se[j])**(-1/mVG[j])-1)**(1/nVG[j]) for j in range(0,len(nVG))] # [mm]
        else:
            # Hydraulic Head
            Ptem = numpy.zeros(n_local)
            for jk in range(0,n_local):
                CondSuc = self.Conductivity_Suction(2,Ks_Zs[jk],Osat[jk],Ohy[jk],L[jk],Pe[jk],O33[jk],alpVG[jk],nVG[jk],O[jk])
                Ptem[jk] = CondSuc[1]
            POT = [PHead[j] - Ptem[j] for j in range(0, len(nVG))] # [mm]

        # Dunne Runoff [mm]
        Rd = WTR[0]
        # Evaporation Layer Water Content
        # Evaporation Bare Soil WC [-]
        OS = sum([EvL_Zs[j]*O[j] for j in range(0, len(O))])
        # First layer
        # Infiltration Water Content [-]
        OF = sum([Inf_Zs[j]*O[j] for j in range(0, len(O))])

        # VEGETATION
        cc = 1
        n_local = len(RfH_Zs)

        for i in range(0,cc):
            OH = sum([RfH_Zs[j]*O[j] for j in range(0,len(O))])
            OL = sum([RfL_Zs[j]*O[j] for j in range(0,len(O))])
            CondSuc = self.Conductivity_Suction(SPAR,sum([RfH_Zs[j]*Ks_Zs[j] for j in range(0,len(Ks_Zs))]),
                                               sum([RfH_Zs[j]*Osat[j] for j in range(0,len(Osat))]),
                                               sum([RfH_Zs[j]*Ohy[j] for j in range(0,len(Ohy))]),
                                               sum([RfH_Zs[j]*L[j] for j in range(0,len(L))]),
                                               sum([RfH_Zs[j]*Pe[j] for j in range(0,len(Pe))]),
                                               sum([RfH_Zs[j]*O33[j] for j in range(0,len(O33))]),
                                               sum([RfH_Zs[j]*alpVG[j] for j in range(0,len(alpVG))]),
                                               sum([RfH_Zs[j]*nVG[j] for j in range(0,len(nVG))]),
                                               OH)
            Psi_s_H = CondSuc[1] # [mm]
            CondSuc = self.Conductivity_Suction(SPAR, sum([RfL_Zs[j] * Ks_Zs[j] for j in range(0, len(Ks_Zs))]),
                                               sum([RfL_Zs[j] * Osat[j] for j in range(0, len(Osat))]),
                                               sum([RfL_Zs[j] * Ohy[j] for j in range(0, len(Ohy))]),
                                               sum([RfL_Zs[j] * L[j] for j in range(0, len(L))]),
                                               sum([RfL_Zs[j] * Pe[j] for j in range(0, len(Pe))]),
                                               sum([RfL_Zs[j] * O33[j] for j in range(0, len(O33))]),
                                               sum([RfL_Zs[j] * alpVG[j] for j in range(0, len(alpVG))]),
                                               sum([RfL_Zs[j] * nVG[j] for j in range(0, len(nVG))]),
                                               OL)
            Psi_s_L = CondSuc[1] # [mm]

        Psi_s_H = -(Psi_s_H/1000)*1000*9.81/1e+6  # [MPa]
        Psi_s_L = -(Psi_s_L/1000)*1000*9.81/1e+6  # [MPa]

        # Water density [mmol H20 m^-3]
        rho2 = 55555555
        # Radius cylinder of soil to which root has access to [m]
        rcyl= 2.0*1e-3
        # Root radius [m]
        rroot = 0.5 * 1e-3
        Psi_s = numpy.zeros(n_local)
        gsr_L = numpy.zeros(n_local)
        gsr_H = numpy.zeros(n_local)
        for jk in range(0,n_local):
            CondSuc = self.Conductivity_Suction(SPAR,Ks_Zs[jk],Osat[jk],Ohy[jk],L[jk],Pe[jk],O33[jk],alpVG[jk],nVG[jk],O[jk])
            Ko = CondSuc[0] # [mm s^-1]
            Psi_s[jk] = CondSuc[1] # [mm]
            gsr_L[jk] = self.Root_Soil_Conductance(Ko,RfL_Zs[jk]*Rrootl_L,rcyl,rroot,dz[jk]) # [mmol H20 m^-2 ground s^-1 MPa^-1]

            gsr_H[jk] = self.Root_Soil_Conductance(Ko,RfH_Zs[jk]*Rrootl_H,rcyl,rroot,dz[jk]) # [mmol H20 m^-2 ground s^-1 MPa^-1]

        Psi_s = [-(Psi_s[j] / 1000) * 1000 * 9.81 / 1e+6 for j in range(0,len(Psi_s))] # [MPa]

        Psi_minH = min(PsiX50_H, PsiL50_H)
        Psi_minL = min(PsiX50_L, PsiL50_L)
        # Max extractable water [mm m^2 m^-2 ground s^-1]
        Exwat_L = gsr_L / rho2 * 1000 * (-numpy.tile(Psi_minL,(n_local)) + numpy.tile(Psi_s,(cc)))
        Exwat_H = gsr_H / rho2 * 1000 * (-numpy.tile(Psi_minH,(n_local)) + numpy.tile(Psi_s,(cc)))
        Exwat_L[Exwat_L < 0] = 0
        Exwat_H[Exwat_H < 0] = 0
        gsr_L = numpy.sum(gsr_L,axis=0) # [mmol H20 m^-2 ground s^-1 MPa^-1]
        gsr_H = numpy.sum(gsr_H,axis=0) # [mmol H20 m^-2 ground s^-1 MPa^-1]

        self.SoilWaterMultLay = SoilWaterMultiLayerDef()
        self.SoilWaterMultLay.O = O             # []
        self.SoilWaterMultLay.ZWT = ZWT         # [mm]
        self.SoilWaterMultLay.OF = OF           # []
        self.SoilWaterMultLay.OS = OS           # []
        self.SoilWaterMultLay.Psi_s_H = Psi_s_H # [MPa]
        self.SoilWaterMultLay.Psi_s_L = Psi_s_L # [MPa]
        self.SoilWaterMultLay.gsr_H = gsr_H     # [mmol H20 m^-2 ground s^-1 MPa^-1]
        self.SoilWaterMultLay.gsr_L = gsr_L     # [mmol H20 m^-2 ground s^-1 MPa^-1]
        self.SoilWaterMultLay.Exwat_H = Exwat_H # [mm m^2 m^-2 ground s^-1]
        self.SoilWaterMultLay.Exwat_L = Exwat_L # [mm m^2 m^-2 ground s^-1]
        self.SoilWaterMultLay.Rd = Rd           # [mm]
        self.SoilWaterMultLay.WTR = WTR         # [mm]
        self.SoilWaterMultLay.POT = POT         # [mm]
        self.SoilWaterMultLay.OH = OH           # []
        self.SoilWaterMultLay.OL = OL           # []

        return O,ZWT,OF,OS,Psi_s_H,Psi_s_L,gsr_H,gsr_L,Exwat_H,Exwat_L,Rd,WTR,POT,OH,OL

    def Soil_Parameters_Total(self,Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs):
        """
        ------
        INPUT:
        Pcla: Fraction of clay in the soil [-]
        Psan: Fraction of sand in the soil [-]
        Porg: Fraction of organic material in the soil [-]
        Kfc: Conductivity at field capacity [mm h^-1]
        Phy: Suction at the residual/hygroscopic water content [kPa]
        SPAR: Soil parameter type (1:VanGenuchten, 2:Saxton-Rawls)
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
        -------
        OUTPUT:
        Zs: Soil layer discretization [mm]
        dz: Thickness of soil layers [mm]
        ms: Number of soil layers [-]
        Osat: Water content at saturation, saturation moisture 0 kPa  [-]
        Ohy:  MoistHygroscopicure Evaporation cessation 10000 kPa - 10 MPa - 1000 m [-]
        nVG: n parameter Van-Genuchten soil water retention curve [mm^-1]
        alpVG: Alpha parameter Van-Genuchten soil water retention curve [mm^-1]
        Ks_Zs: Hydraulic conductivity at saturation for each soil layer [mm s^-1]
        L: Slope of logarithmic tension-moisture curve
        Pe: Tension at air entry (bubbling pressure) [kPa]
        O33: 33 kPa Moisture [-]
        SPAR: Soil parameter type (1:VanGenuchten, 2:Saxton-Rawls)
        EvL_Zs: Fraction of evaporation depth in a specific soil layer [-]
        Inf_Zs: Fraction of infiltration depth in a specific soil layer [-]
        RfH_Zs: Root Fraction for High Vegetation [1...m] [%]
        RfL_Zs: Root Fraction for Low Vegetation [1...m] [%]
        Zinf: Depth of infiltration layer (=first layer) [mm]
        Kbot: Conductivity at the bedrock layer [mm s^-1]
        Slo_pot: [fraction dy/dx]
        Dz: Delta Depth Between First Middle Layer and soil surface [mm]
        aR: Anisotropy ratio
        aTop: Ratio between Area and Contour-Lenght [mm]
        rsd: Normal density dry soil [kg m^-3]
        lan_dry: Thermal conductivity dry soil [W m^-1 K^-1]
        lan_s: Thermal conductivity solid soil [W m^-1 K^-1]
        cv_s: Volumetric heat capacity solid soil [J m^-3 K^-1]
        """

        ms = len(Zs) - 1
        # Thickness of the Layers [mm]
        dz = numpy.diff(Zs)
        Dz = numpy.zeros(ms)

        for i in range(0,ms):
            if i > 0:
                # Delta Depth Between Middle Layer  [mm]
                Dz[i] = (dz[i]+dz[i-1])/2
            else:
                # Delta Depth Between First Middle Layer and soil surface [mm]
                Dz[i] = dz[0]/2

        # Depth of desorption, Depth of evaporation layer (=first layer) [mm]
        Zdes = Zs[1] - Zs[0]
        # Depth of infiltration layer (=first layer) [mm]
        Zinf = Zs[1] - Zs[0]

        # Fraction of evaporation depth in a specific soil layer [-]
        EvL_Zs = self.Evaporation_Layers(Zs,Zdes)
        # Fraction of infiltration depth in a specific soil layer [-]
        Inf_Zs = self.Evaporation_Layers(Zs,Zinf)

        # [fraction dy/dx]
        Slo_pot = numpy.zeros(ms)
        # Anisotropy ratio
        aR = 1
        cellsize = 1
        # Ratio between Area and Contour-Lenght [mm]
        aTop = 1000*cellsize**2/cellsize

        #---------------------------------------------------------------------------------------
        # Calculation of total thermal capacity out of soil composition for force restore method
        #---------------------------------------------------------------------------------------
        # Soil Parameters: characterization of soil parameters out of soil composition
        SoilParam = self.Soil_Parameters(Psan,Pcla,Porg)
        Osat = SoilParam[0]   # Saturation moisture 0 kPa  []
        L = SoilParam[1]      # Slope of logaritimc tension-moisture curve
        Pe = SoilParam[2]     # Tension at air entry (bubbling pressure) [kPa]
        Ks = SoilParam[3]     # saturation conductivity [mm s^-1]
        O33 = SoilParam[4]    # 33 kPa Moisture []
        rsd = SoilParam[5]    # density density dry soil [kg m^-3]
        lan_dry = SoilParam[6]# Thermal conductivity dry soil [W m^-1 K^-1]
        lan_s = SoilParam[7]  # Thermal conductivity soil solid [W m^-1 K^-1]
        cv_s = SoilParam[8]   # Volumetric heat capacity soil solid [J m^-3 K^-1]

        # Normal density dry soil [kg m^-3]
        rsd = rsd*numpy.ones(ms)
        # Thermal conductivity dry soil [W m^-1 K^-1]
        lan_dry	= lan_dry*numpy.ones(ms)
        # Thermal conductivity solid soil [W m^-1 K^-1]
        lan_s =	lan_s*numpy.ones(ms)
        # Volumetric heat capacity solid soil [J m^-3 K^-1]
        cv_s = cv_s*numpy.ones(ms)

        # Alpha parameter Van-Genuchten soil water retention curve [mm^-1]
        p = 3+2/L
        m = 2/(p-1)
        nVG = 1/(1-m)
        alpVG = (((-101.9368*Pe)*(2*p*(p-1))/(p+3))*((55.6+7.4*p+p**2)/(147.8+8.1*p+0.092*p**2)))**(-1) # [mm^-1]

        # Initializing variables for each soil layer
        # Water content at saturation, saturation moisture 0 kPa  [-]
        Osat = Osat*numpy.ones(ms)
        # Slope of logarithmic tension-moisture curve [-]
        L = L*numpy.ones(ms)
        # Tension at air antry (bubbling pressure) [kPa]
        Pe = Pe*numpy.ones(ms)
        # Hydraulic conductivity at saturation for each soil layer [mm s^-1]
        Ks_Zs = Ks*numpy.ones(ms)
        # Soil water content at -33 [kPa] of water potential
        O33 = O33*numpy.ones(ms)
        # n parameter Van-Genuchten soil water retention curve [mm^-1]
        nVG = nVG*numpy.ones(ms)
        # Alpha parameter Van-Genuchten soil water retention curve [mm^-1]
        alpVG = alpVG*numpy.ones(ms)

        SoilParamII = self.Soil_ParametersII(ms,Osat,L,Pe,Ks_Zs,O33,nVG,alpVG,Kfc,numpy.nan,numpy.nan,Phy,numpy.nan,2)
        Ohy = SoilParamII[3]
        # Root Distribution ZR_H -ZR_L ---  Root Fraction in a given soil layer [1...m]
        if CASE_ROOT_H == CASE_ROOT_L:
            CASE_ROOT = copy.copy(CASE_ROOT_H)
        else:
            CASE_ROOT = copy.copy(CASE_ROOT_H)
            print('CASE_ROOT_H and CASE_ROOT_L are not the same. CASE_ROOT_H is taken for the calculation')

        RootFracGeneral = self.Root_Fraction_General(Zs,CASE_ROOT,ZR95_H,ZR50_H,ZR95_L,ZR50_L,ZRmax_H,ZRmax_L)
        RfH_Zs = RootFracGeneral[0] # Root Fraction for High Vegetation [1...m] [%]
        RfL_Zs = RootFracGeneral[1] # Root Fraction for Low Vegetation [1...m] [%]

        Kbot_s = Kbot/3600 # [mm s^-1]

        self.SoilParamTotal = SoilParametersTotal()
        self.SoilParamTotal.Zs = Zs          # [mm]
        self.SoilParamTotal.dz = dz          # [mm]
        self.SoilParamTotal.ms = ms          # []
        self.SoilParamTotal.Osat = Osat      # []
        self.SoilParamTotal.Ohy = Ohy        # []
        self.SoilParamTotal.nVG = nVG        # [mm^-1]
        self.SoilParamTotal.alpVG = alpVG    # [mm^-1]
        self.SoilParamTotal.Ks_Zs = Ks_Zs    # [mm s^-1]
        self.SoilParamTotal.L = L            # []
        self.SoilParamTotal.Pe = Pe          # [kPa]
        self.SoilParamTotal.O33 = O33        # []
        self.SoilParamTotal.SPAR = SPAR
        self.SoilParamTotal.EvL_Zs = EvL_Zs # []
        self.SoilParamTotal.Inf_Zs = Inf_Zs # []
        self.SoilParamTotal.RfH_Zs = RfH_Zs # []
        self.SoilParamTotal.RfL_Zs = RfL_Zs # []
        self.SoilParamTotal.Zinf = Zinf     # [mm]
        self.SoilParamTotal.Kbot = Kbot_s   # [mm S^-1]
        self.SoilParamTotal.Slo_pot = Slo_pot
        self.SoilParamTotal.Dz = Dz         # [mm]
        self.SoilParamTotal.aR = aR
        self.SoilParamTotal.aTop = aTop     # [mm]
        self.SoilParamTotal.rsd = rsd       # [kg m^-3]
        self.SoilParamTotal.lan_dry = lan_dry # [W m^-1 K^-1]
        self.SoilParamTotal.lan_s = lan_s   # [W m^-1 K^-1]
        self.SoilParamTotal.cv_s = cv_s     # [J m^-3 K^-1]


        return Zs,dz,ms,Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,SPAR,EvL_Zs,Inf_Zs,RfH_Zs,RfL_Zs,Zinf,Kbot_s,Slo_pot,Dz,aR,aTop,\
               rsd,lan_dry,lan_s,cv_s

    def Evaporation_Layers(self,Zs,Zdes):

        """
        ------
        INPUT:
        Zs: Soil layer discretization [mm]
        Zdes: Depth of desorption, Depth of evaporation layer (=first layer) [mm]
        -------
        OUTPUT:
        EvL_Zs : Evaporation Layer fraction [1...m] [%]
        """

        n = len(Zs) - 1
        EvL_Zs = numpy.zeros(n)
        if Zdes < Zs[0]:
            print('ERROR FIRST LAYER TOO DEPTH')

        if Zdes > Zs[n]:
            print('ERROR LAST LAYER TOO SHALLOW')

        # Number of interested Layer
        NiL = sum(i < Zdes for i in Zs)
        # Thickness of the Layers [mm]
        dz = numpy.diff(Zs)
        # Thickness of last layer [mm]
        dz[NiL] = Zdes - Zs[NiL]
        # Interested Layer
        dz = dz[0:NiL]
        EvL_Zs[0: NiL] = dz / Zdes
        return EvL_Zs

    def Volume_Correction(self,V,EvL_Zs,RfH_Zs,RfL_Zs,EG,T_H,T_L,Lk):

        """
        ------
        INPUT:
        V: Water content in each soil layer [mm]
        EvL_Zs: Fraction of evaporation depth in a specific soil layer [-]
        RfH_Zs: Root Fraction for High Vegetation [1...m] [%]
        RfL_Zs: Root Fraction for Low Vegetation [1...m] [%]
        EG: Water flux from soil under vegetation [mm]
        T_H: Root water uptake from different soil layers of high vegetation [mm]
        T_L: Root water uptake from different soil layers of low vegetation [mm]
        Lk: Leakage at bedrock [mm]
        -------
        OUTPUT:
        V_local: Updated water content in each soil layer [mm]
        T_H_local: Updated root water uptake from different soil layers of high vegetation [mm]
        T_L_local: Updated root water uptake from different soil layers of low vegetation [mm]
        EG_local: Updated water flux from soil under vegetation [mm]
        Lk_local: Updated leakage at bedrock [mm]
        """

        # Re-define input parameters which are overwritten in this function
        V_local = copy.copy(V)
        EG_local = copy.copy(EG)
        T_H_local = copy.copy(T_H)
        T_L_local = copy.copy(T_L)
        Lk_local = copy.copy(Lk)

        # Number of active layer for transpiration H
        lay_H = sum(i > 0 for i in RfH_Zs)
        # Number of active layer for transpiration L
        lay_L = sum(i > 0 for i in RfL_Zs)
        # Number of active layer for Evaporation
        lay_G = sum(i > 0 for i in EvL_Zs)

        #----------------------------------------
        # Compensatory Mechanism on deeper layers
        #----------------------------------------
        for i in range(0,lay_G):
            if (V_local[i]< 0) and (i < lay_G-1):
                V_local[i+1] = V_local[i+1]+V_local[i]
                V_local[i] = 0

        for i in range(0,lay_L):
            if (V_local[i] < 0) and (i < lay_L-1):
                V_local[i+1] = V_local[i+1]+V_local[i]
                V_local[i] = 0

        for i in range(0,lay_H):
            if (V_local[i] < 0) and (i < lay_H-1):
                V_local[i+1] = V_local[i+1]+V_local[i]
                V_local[i] = 0

        # Compensatory Mechanism on shallow layers
        if sum(i < 0 for i in V_local) > 0:
            for i in range(max([idx for idx, val in enumerate(V_local) if val < 0]),1,-1):
                if (V_local[i] < 0) and (V_local[i-1] >= 0):
                    V_local[i-1] = V_local[i-1]+V_local[i]
                    V_local[i] = 0

        # -----------------
        # Brutal correction
        # -----------------
        if sum(i < 0 for i in V_local) > 0:
            for i in range(0,lay_G):
                if V_local[i] < 0:
                    if abs(V_local[i]) > EG_local:
                        V_local[i] = V_local[i]+EG_local
                        EG_local = 0
                    else:
                        EG_local = EG_local+V_local[i]
                        V_local[i] = 0

            for i in range(0,lay_L):
                if V_local[i] < 0:
                    if abs(V_local[i]) > sum(T_L_local):
                        V_local[i] = V_local[i]+sum(T_L_local)
                        T_L_local = [0 for j in range(0,len(T_L_local))]
                    else:
                        T_L_local = [T_L_local[j]+V_local[i]*T_L_local[j]/sum(T_L_local) for j in range(0,len(T_L_local))]
                        V_local[i] = 0

            for i in range(0,lay_H):
                if V_local[i] < 0:
                    if abs(V_local[i]) > sum(T_H_local):
                        V_local[i] = V_local[i] + sum(T_H_local)
                        T_H_local = [0 for j in range(0, len(T_H_local))]
                    else:
                        T_H_local = [T_H_local[j] + V_local[i] * T_H_local[j] / sum(T_H_local) for j in range(0, len(T_H_local))]
                        V_local[i] = 0

            if V_local[-1] < 0:
                i = 1
                Lk_local = Lk_local+V_local[i]
                V_local[i] = 0

        if sum(i < 0 for i in V_local) > 0:
            EG_local = EG_local + sum(i for i in V_local if i < 0)
            for i in range(0,len(V_local)):
                if V_local[i] < 0:
                    V_local[i] = 0

        self.VolCorr = VolumeCorrectionDef()
        self.VolCorr.V = V_local
        self.VolCorr.T_H = T_H_local
        self.VolCorr.T_L = T_L_local
        self.VolCorr.EG = EG_local
        self.VolCorr.Lk = Lk_local

        return V_local,T_H_local,T_L_local,EG_local,Lk_local

    def Soil_Moistures_Rich_Comp(self,V,t,Lk,f,EG,T_H,T_L,Qi_in,Slo_pot,IS,SPAR,Osat,Ohy,O33,dz,Ks_Zs,Dz,numn,L,Pe,aR,aT,
                                 alpVG,nVG,Zs,cosalp,sinalp,SN):

        """
        ------
        INPUT:
        V: Volume of water in each soil layer [mm]
        t: Time [s]
        Lk: Leakage at bedrock [mm s^-1]
        f: Infiltration rate into first soil layer [mm s^-1]
        EG: Distribution of evaporation from bare soil [mm s^-1]
        T_H: Distribution of root water uptake from different soil layers of high vegetation [mm s^-1]
        T_L: Distribution of root water uptake from different soil layers of low vegetation [mm s^-1]
        Qi_in: Lateral water flux [mm s^-1]
        Slo_pot:
        IS: ones vector
        SPAR: oil parameter type
        Osat: Water content at saturation, saturation moisture 0 kPa  [-]
        Ohy: e Evaporation cessation 10000 kPa - 10 MPa - 1000 m [-]
        O33: 33 kPa Moisture [-]
        dz: Thickness of soil layers [mm]
        Ks_Zs: Hydraulic conductivity at saturation for each soil layer [mm s^-1]
        Dz: Delta Depth Between First Middle Layer and soil surface [mm]
        numn: Number of soil layers [-]
        L: Slope of logarithmic tension-moisture curve
        Pe: Tension at air entry (bubbling pressure) [kPa]
        aR: Anisotropy ratio
        aT: Ratio between Area/ContourLenght [mm]
        alpVG: Alpha parameter Van-Genuchten soil water retention curve [mm^-1]
        nVG: n parameter Van-Genuchten soil water retention curve [mm^-1]
        Zs: Soil layer discretization [mm]
        cosalp:
        sinalp:
        SN:
        -------
        OUTPUT:
        dV: Change in volume of water in each soil layer [mm s^-1]
        """

        O = [V[i]/dz[i]+Ohy[i] for i in range(0,len(V))]
        for i in range(0,len(O)):
            if O[i] >= Osat[i]-1e-5:
                O[i] = Osat[i]
            if O[i] <= Ohy[i]+1e-5:
                O[i] = Ohy[i]+1e-5

        numnm1 = numn - 1
        vap = Osat[0:numnm1]

        # Soil water
        dV = numpy.zeros(len(V))
        K = numpy.zeros(len(Osat))
        P = numpy.zeros(len(Osat))
        Khalf = numpy.zeros(len(vap))
        q = numpy.zeros(len(vap))
        qf = numpy.zeros(len(vap))

        if SPAR == 1:
            mVG = [1 - 1 / nVG[i] for i in range(len(nVG))] # [mm]
            Se = [(O[i] - Ohy[i]) / (Osat[i] - Ohy[i]) for i in range(len(O))] # [-]
            P = [(1 / alpVG[i]) * ((Se[i])**(-1 / mVG[i]) - 1)** (1 / nVG[i]) for i in range(0,len(nVG))] # [mm]
            K = [Ks_Zs[i] * ((Se[i])**(0.5)) * (1 - (1 - (Se[i])**(1 / mVG[i]))** mVG[i])** 2 for i in range(0,len(Se))] # [mm s^-1]
        else:
            B = [1/L[i] for i in range(0,len(L))]
            # Coefficient of moisture tension
            A = [numpy.exp(math.log(33) + B[i] * math.log(O33[i])) for i in range(0,len(B))]
            for i in range(0,numn):
                K[i] = Ks_Zs[i] * (O[i] / Osat[i])**(3 + (2 / L[i])) # [mm s^-1]
                if O[i] < O33[i]:
                    P[i] = A[i]*(O[i]**(-B[i])) # [kPa]
                else:
                    P[i] = 33 - ((O[i] - O33[i]) * (33 - Pe[i]) / (Osat[i] - O33[i]))

            P = [-101.9368*P[i] for i in range(0,len(P))] # [mm]

        # Trasmissivity  unsaturated/saturated [mm^2 s^-1]
        To = [K[i] * aR * dz[i] for i in range(0,len(dz))]

        # Slo_pot: Slope Total Hydraulic Head [mm s^-1]
        if SN == 1:
            Qi_out = [IS[i] * (To[i] / aT) * sinalp for i in range(0,len(To))]
        else:
            Qi_out = [(To[i] / aT) * sinalp for i in range(0,len(To))]

        #---------------
        # Richards Model
        #---------------
        Khalf = [0.5*(K[i]+K[i+1]) for i in range(0,len(K)-1)] # [mm s^-1]
        # Flux positive downward from layer i+1 (above) to i (below) [mm s^-1]
        q = [Khalf[i]*(1*cosalp - (P[i+1]-P[i])/Dz[i+1]) for i in range(0,len(Khalf))]

        qf = copy.copy(q)
        for i in range(1,numn):
            if O[i] >= Osat[i]-1e-5:
                qf[i-1] = 0

        for i in range(0,len(q)):
            if q[i] > 0:
                q[i] = qf[i]

        # Soils water balance without other terms
        dV[0] = f - q[0] - T_H[0] - T_L[0] - EG[0] + Qi_in[0] - Qi_out[0]
        Test = numpy.zeros(numn)
        Test[0] = f - q[0] - T_H[0] - T_L[0] - EG[0] + Qi_in[0] - Qi_out[0] - dV[0]
        for i in range(1,numn-1):
            dV[i] = q[i-1] - q[i] - T_H[i] - T_L[i] - EG[i] + Qi_in[i] - Qi_out[i]
            Test[i] = q[i-1] - q[i] - T_H[i] - T_L[i] - EG[i] + Qi_in[i] - Qi_out[i] - dV[i]

        dV[numn-1] = q[numnm1-1] - Lk - T_H[numn-1] - T_L[numn-1] - EG[numn-1] + Qi_in[numn-1] - Qi_out[numn-1]
        Test[numn-1] = q[numnm1-1] - Lk - T_H[numn-1] - T_L[numn-1] - EG[numn-1] + Qi_in[numn-1] - Qi_out[numn-1] - dV[numn-1]

        return dV

    def Soil_Moistures_Rich_Comp_Lat2(self,Vlat,t,dz,SPAR,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,C1,C2,f1,f2,Wcan):
        """
        ------
        INPUT:
        Vlat: Water volume in each layer [mm]
        t: time [s]
        dz: Thickness of the layers [mm]
        SPAR: Soil parameter type
        Ks: Hydraulic conductivity at saturation for each soil layer [mm s^-1]
        Osat: Water content at saturation, saturation moisture 0 kPa [-]
        Ohy: MoistHygroscopicure Evaporation cessation 10000 kPa - 10 MPa - 1000 m [-]
        L: Slope of logarithmic tension-moisture curve
        Pe: Tension at air entry (bubbling pressure) [kPa]
        O33: 33 kPa Moisture [-]
        alpVG: Alpha parameter Van-Genuchten soil water retention curve [mm^-1]
        nVG: n parameter Van-Genuchten soil water retention curve [mm^-1]
        C1: Boolean operator for presence and absence of bare ground
        C2: Boolean operator for presence and absence of vegetated ground
        f1: Fraction of ground coverd by bare soil
        f2: Fraction of ground coverd by vegetation
        Wcan: Width of the canyon [m]
        -------
        OUTPUT:
        dVlat: Soil water balance without other terms [mm s^-1]
        """

        Olat = [Vlat[i] / dz + Ohy for i in range(0, len(Vlat))]
        for i in range(0, len(Olat)):
            if Olat[i] >= Osat - 1e-5:
                Olat[i] = Osat - 1e-5
            if Olat[i] <= Ohy + 1e-5:
                Olat[i] = Ohy + 1e-5

        # Hydraulic conductivity and soil water potential
        self.Conductivity_Suction(SPAR,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,Olat)
        Ko = self.CondSuc.Ko # [mm s^-1]
        Po = self.CondSuc.Po # [mm]

        # Lateral water re-distribution
        # Assumption: horizontal and vertical unsaturated conductivity is the same
        a = 15
        dxsoil = 1000 # [mm] = 1 [m]

        # Calculate lateral water flow [mm s^-1]
        # The higher the soil water potential the drier the soil. Hence, I put a minus to change flux direction.
        Qlat_1to2 = -a * (Ko[0] + Ko[1]) / 2 * (Po[0] - Po[1]) / dxsoil
        Qlat_2to1 = -a * (Ko[0] + Ko[1]) / 2 * (Po[1] - Po[0]) / dxsoil

        # Transmissivity [mm^2 s^-1]
        T_1to2 = Qlat_1to2 * dz
        T_2to1 = Qlat_2to1 * dz

        T_totflux = sum(numpy.nansum([[T_1to2],[T_2to1]],axis=0))
        if T_totflux != 0:
            print('The lateral transmissivities do not add up to 0. Please check Soil_Moistures_Rich_Comp_Lat.m')

        # Re-scale the horizontal water flux over given layer depth [mm s^-1]
        Qin_1to2 = T_1to2 / (f2 * 1000 * Wcan) * C1 * C2
        Qin_2to1 = T_2to1 / (f1 * 1000 * Wcan) * C2 * C1

        # Total flux incoming to one soil column from the two other soil columns [mm s^-1]
        Qin_1 = numpy.nansum(Qin_2to1,axis=0)
        Qin_2 = numpy.nansum(Qin_1to2,axis=0)

        # Soil water balance without other terms [mm s^-1]
        dVlat = [Qin_1,Qin_2]


        return dVlat

    def Soil_Moistures_Rich_Comp_Lat3(self,Vlat,t,dz,SPAR,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,Cimp,Cbare,Cveg,fimp,fbare,fveg,Wcan):
        """
        ------
        INPUT:
        Vlat: Water volume in each layer [mm]
        t: time [s]
        dz: Thickness of the layers [mm]
        SPAR: Soil parameter type
        Ks: Hydraulic conductivity at saturation for each soil layer [mm s^-1]
        Osat: Water content at saturation, saturation moisture 0 kPa [-]
        Ohy: MoistHygroscopicure Evaporation cessation 10000 kPa - 10 MPa - 1000 m [-]
        L: Slope of logarithmic tension-moisture curve
        Pe: Tension at air entry (bubbling pressure) [kPa]
        O33: 33 kPa Moisture [-]
        alpVG: Alpha parameter Van-Genuchten soil water retention curve [mm^-1]
        nVG: n parameter Van-Genuchten soil water retention curve [mm^-1]
        Cimp: Boolean operator for presence and absence of impervious surface
        Cbare: Boolean operator for presence and absence of bare soil
        Cveg: Boolean operator for presence and absence of vegetated ground
        fimp: Fraction of ground covered by impervious surface [-]
        fbare: Fraction of ground covered by bare soil [-]
        fveg: Fraction of ground covered by vegetated ground [-]
        Wcan: Width of the canyon [m]
        -------
        OUTPUT:
        dVlat: Soil water balance without other terms [mm s^-1]
        """


        Olat = [Vlat[i] / dz + Ohy for i in range(0, len(Vlat))]
        for i in range(0, len(Olat)):
            if Olat[i] >= Osat - 1e-5:
                Olat[i] = Osat - 1e-5
            if Olat[i] <= Ohy + 1e-5:
                Olat[i] = Ohy + 1e-5

        # Hydraulic conductivity and soil water potential
        CondSuc = self.Conductivity_Suction(SPAR,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,Olat)
        Ko = CondSuc[0]  # [mm s^-1]
        Po = CondSuc[1]  # [mm]

        # Lateral water re-distribution
        # Assumption: horizontal and vertical unsaturated conductivity is the same
        a = 15
        dxsoil = 1000  # [mm] = 1 [m]

        Kf_gimp = Ko[0]  # [mm s^-1]
        Kf_gbare = Ko[1] # [mm s^-1]
        Kf_gveg = Ko[2]  # [mm s^-1]

        Psi_soil_gimp = Po[0]  # [mm]
        Psi_soil_gbare = Po[1] # [mm]
        Psi_soil_gveg = Po[2]  # [mm]

        # Calculate lateral water flow [mm s^-1]
        # The higher the soil water potential the drier the soil. Hence, I put a minus to change flux direction.
        Qlat_bare2imp = -a * (Kf_gbare + Kf_gimp) / 2 * (Psi_soil_gbare - Psi_soil_gimp) / dxsoil
        Qlat_veg2imp = -a * (Kf_gveg + Kf_gimp) / 2 * (Psi_soil_gveg - Psi_soil_gimp) / dxsoil
        Qlat_veg2bare = -a * (Kf_gbare + Kf_gveg) / 2 * (Psi_soil_gveg - Psi_soil_gbare) / dxsoil
        Qlat_imp2bare = -a * (Kf_gbare + Kf_gimp) / 2 * (Psi_soil_gimp - Psi_soil_gbare) / dxsoil
        Qlat_bare2veg = -a * (Kf_gbare + Kf_gveg) / 2 * (Psi_soil_gbare - Psi_soil_gveg) / dxsoil
        Qlat_imp2veg = -a * (Kf_gveg + Kf_gimp) / 2 * (Psi_soil_gimp - Psi_soil_gveg) / dxsoil

        # Transmissivity [mm^2 s^-1]
        Tveg2imp = Qlat_veg2imp * dz
        Tbare2imp = Qlat_bare2imp * dz
        Tveg2bare = Qlat_veg2bare * dz
        Timp2bare = Qlat_imp2bare * dz
        Tbare2veg = Qlat_bare2veg * dz
        Timp2veg = Qlat_imp2veg * dz

        T_totflux = sum(numpy.nansum([[Tbare2veg],[Tveg2bare],[Tbare2imp],[Timp2bare],[Tveg2imp],[Timp2veg]],axis=0))
        if T_totflux != 0:
            print('The lateral transmissivities do not add up to 0. Please check Soil_Moistures_Rich_Comp_Lat.m')

        # Re-scale the horizontal water flux over given layer depth [mm s^-1]
        # incoming flux 2imp, 2bare, 2veg are positive, outgoing negative
        Qin_veg2imp = Tveg2imp / (fimp * 1000 * Wcan) * Cveg * Cimp
        Qin_bare2imp = Tbare2imp / (fimp * 1000 * Wcan) * Cbare * Cimp
        Qin_veg2bare = Tveg2bare / (fbare * 1000 * Wcan) * Cveg * Cbare
        Qin_imp2bare = Timp2bare / (fbare * 1000 * Wcan) * Cimp * Cbare
        Qin_bare2veg = Tbare2veg / (fveg * 1000 * Wcan) * Cbare * Cveg
        Qin_imp2veg = Timp2veg / (fveg * 1000 * Wcan) * Cimp * Cveg

        # Total flux incoming to one soil column from the two other soil columns  [mm s^-1]
        Qin_imp = numpy.nansum([Qin_veg2imp,Qin_bare2imp], axis=0)
        Qin_bare = numpy.nansum([Qin_veg2bare,Qin_imp2bare], axis=0)
        Qin_veg = numpy.nansum([Qin_bare2veg,Qin_imp2veg], axis=0)

        # Soil water balance without other terms [mm s^-1]
        dVlat = [Qin_imp,Qin_bare,Qin_veg]

        return dVlat

    def Soil_Moisture_Conductivity_Update(self,V,Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,
                                          ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs,Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L):
        """
        ------
        INPUT:
        V: Volume of water in each soil layer per unit area [mm]
        Pcla: Fraction of clay in the soil [-]
        Psan: Fraction of sand in the soil [-]
        Porg: Fraction of organic material in the soil [-]
        Kfc: Conductivity at field capacity [mm h^-1]
        Phy: Suction at the residual/hygroscopic water content [kPa]
        SPAR: Soil parameter type (1:VanGenuchten, 2:Saxton-Rawls)
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
        Rrootl_H: Root length index of high vegetation [m root m^-2 PFT]
        Rrootl_L: Root length index of low vegetation [m root m^-2 PFT]
        PsiL50_H: Water Potential at 50% loss conductivity of high vegetation [MPa]
        PsiL50_L: Water Potential at 50% loss conductivity of low vegetation [MPa]
        PsiX50_H: Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil of high vegetation [MPa]
        PsiX50_L: Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil of low vegetation [MPa]
        -------
        OUTPUT:
        V: Volume of water in each soil layer per unit area [mm]
        O: Water Content [1...m] Layer [-]
        OS: Water Content for evaporation [-]
        Psi_soil: Tension at O [mm]
        Psi_s_H: Soil water potential for first layer of vegetation [MPa]
        Psi_s_L: Soil Water Potential  for Second Layer of Vegetation [MPa]
        Exwat_H: Maximum extractable water [mm m^2 m^-2 ground s^-1]
        Exwat_L: Maximum extractable water [mm m^2 m^-2 ground s^-1]
        Ko: Hydraulic conductivity at O [mm s^-1]
        """

        # Re-define input parameters which are overwritten in this function
        Zs_local = copy.copy(Zs)

        # Calculate soil parameters depending on soil composition
        self.Soil_Parameters_Total(Pcla,Psan,Porg,Kfc,Phy,SPAR,Kbot,CASE_ROOT_H,CASE_ROOT_L,ZR95_H,ZR95_L,ZR50_H,ZR50_L,ZRmax_H,ZRmax_L,Zs_local)
        Zs_local = self.SoilParamTotal.Zs
        dz = self.SoilParamTotal.dz
        ms = self.SoilParamTotal.ms
        Osat = self.SoilParamTotal.Osat
        Ohy = self.SoilParamTotal.Ohy
        nVG = self.SoilParamTotal.nVG
        alpVG = self.SoilParamTotal.alpVG
        Ks_Zs = self.SoilParamTotal.Ks_Zs
        L = self.SoilParamTotal.L
        Pe = self.SoilParamTotal.Pe
        O33 = self.SoilParamTotal.O33
        EvL_Zs = self.SoilParamTotal.EvL_Zs
        Inf_Zs = self.SoilParamTotal.Inf_Zs
        RfH_Zs = self.SoilParamTotal.RfH_Zs
        RfL_Zs = self.SoilParamTotal.RfL_Zs

        # Update soil moisture content in different soil layers.
        self.Soil_Water_MultiLayer(V,Zs_local,dz,ms,Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,SPAR,EvL_Zs,Inf_Zs,RfH_Zs,RfL_Zs,Rrootl_H,
                                   Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L)
        O = self.SoilWaterMultLay.O
        OS = self.SoilWaterMultLay.OS
        Psi_s_H = self.SoilWaterMultLay.Psi_s_H
        Psi_s_L = self.SoilWaterMultLay.Psi_s_L
        Exwat_H = self.SoilWaterMultLay.Exwat_H
        Exwat_L = self.SoilWaterMultLay.Exwat_L

        # Calculate hydraulic conductivity and soil water potential
        Ko = numpy.zeros(len(O))
        Psi_soil = numpy.zeros(len(O))
        for i in range(0,len(O)):
            CondSuc = self.Conductivity_Suction(SPAR,Ks_Zs[i],Osat[i],Ohy[i],L[i],Pe[i],O33[i],alpVG[i],nVG[i],O[i])
            Ko[i] = CondSuc[0]
            Psi_soil[i] = CondSuc[1]



        self.SoilMoistCond = SoilMoistureConductivityUpdate()
        self.SoilMoistCond.V = V
        self.SoilMoistCond.O = O
        self.SoilMoistCond.OS = OS
        self.SoilMoistCond.Psi_soil = Psi_soil  # [mm]
        self.SoilMoistCond.Psi_s_H = Psi_s_H    # [MPa]
        self.SoilMoistCond.Psi_s_L = Psi_s_L    # [MPa]
        self.SoilMoistCond.Exwat_H = Exwat_H    # [mm m^2 m^-2 ground s^-1]
        self.SoilMoistCond.Exwat_L = Exwat_L    # [mm m^2 m^-2 ground s^-1]
        self.SoilMoistCond.Ko = Ko              # [mm s^-1]

        return V,O,OS,Psi_soil,Psi_s_H,Psi_s_L,Exwat_H,Exwat_L,Ko
















