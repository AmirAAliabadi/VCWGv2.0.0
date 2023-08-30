import os
import numpy
import math
from pprint import pprint
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.gridspec as gridspec
import copy
import _pickle as cPickle
from EB_Roof import EnergyBalanceRoof_Def
from EB_Canyon import EnergyBalanceCanyon_Def
from EB_Rural import EnergyBalanceRural_Def
from SurfaceTemperature import Tsurf_Def
from UrbanModel import UCM_Def
from WB_Roof import WaterBalanceRoof_Def
from WB_Canyon import WaterBalanceCanyon_Def
from weather import Weather
from forcing import Forcing
from Simparam import SimParam
from Write_Output import Write_Forcing,Write_EB,Write_Tsurf,Write_WB,Write_TdeepProfiles,Write_1Dprofiles,\
    Write_Ruralprofiles,Write_BEM
from Radiation_Functions import RadiationFunctions
from RSM import RSMDef
from Read_Input import read_VCWG_param,ForcingData,Data_Site
from ReadDOE import readDOE
from Material import Material
from psychrometrics import HumFromRHumTemp
from EPWGenerator import write_epw

"""
Main VCWG script 
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: May 2020
"""

class VCWG_Hydro(object):

    def __init__(self,epwFileName,TopForcingFileName,VCWGParamFileName,ViewFactorFileName,case):
        self.epwFileName = epwFileName
        self.VCWGParamFileName = os.path.join(os.path.join('resources','Parameters'),VCWGParamFileName)
        self.ViewFactorFileName = os.path.join(os.path.join('resources','Parameters'),ViewFactorFileName)
        self.case = case
        self.TopForcingFileName = TopForcingFileName

    def read_input(self):

        # Read the site parameters
        self.Geometry_m, self.ParTree, self.geometry, self.FractionsRoof, self.FractionsGround, self.WallLayers, self.ParSoilRoof, \
        self.ParSoilGround, self.ParInterceptionTree, self.PropOpticalRoof, self.PropOpticalGround, self.PropOpticalWall, \
        self.PropOpticalTree, self.ParThermalRoof, self.ParThermalGround, self.ParThermalWall,self.ParVegRoof,\
        self.ParVegGround,self.ParVegTree,self.Person,self.ColParam,self.RSMParam,self.TimeParam,ViewFactorCal_Param,self.bld,self.zone,\
        self.charLength,self.BEMParam = Data_Site(self.VCWGParamFileName)

        # Calculate view factors
        RadFun = RadiationFunctions()
        self.ViewFactor, ViewFactorPoint = RadFun.VFUrbanCanyon(ViewFactorCal_Param, self.Geometry_m, self.geometry,
                                                                self.Person, self.ParTree,self.ViewFactorFileName)

    def read_epw(self):

        self.simTime = SimParam(self.TimeParam.dts,self.TimeParam.dtWeather,self.TimeParam.Month,self.TimeParam.Day,self.TimeParam.nDay)

        # Build a new epw file using TopForcing dataset, if there is no information from rural site
        if self.epwFileName == None:
            epw_precision = 1
            try:
                write_epw(self.TopForcingFileName, epw_precision, r'rawEPW.epw', self.simTime.timeInitial)
            except:
                raise Exception("Failed to read TopForcing file!")
            self.epwFileName = 'TopForcing.epw'

        weather = Weather(self.epwFileName, self.simTime.timeInitial, self.simTime.timeFinal)
        forcIP = Forcing(weather.staTemp, weather)

        # Deep soil temperature from epw file
        sdd = list([abs(weather.depth_soil[0] - max(self.ParSoilGround.Zs)/1000), abs(weather.depth_soil[1] - max(self.ParSoilGround.Zs)/1000),
                    abs(weather.depth_soil[2] - max(self.ParSoilGround.Zs)/1000)])
        sdd_index = sdd.index(min(sdd))
        Tdeepsoil_z = weather.Tsoil[sdd_index]

        n = len(forcIP.temp)
        Meteo_time_step = int(self.TimeParam.dtWeather)
        Sim_time_step = int(self.TimeParam.dts)
        time_Meteo = [i for i in range(0, n*Meteo_time_step, Meteo_time_step)]
        self.time_n = [i for i in range(0, max(time_Meteo), Sim_time_step)]

        f_LWR_in = interp1d(time_Meteo, forcIP.infra)
        f_Dif_in = interp1d(time_Meteo, forcIP.dif)
        f_Dir_in = interp1d(time_Meteo, forcIP.dir)
        f_T_atm = interp1d(time_Meteo, forcIP.temp)
        f_windspeed_u = interp1d(time_Meteo, forcIP.wind)
        f_pressure_atm = interp1d(time_Meteo, forcIP.pres)
        f_rain = interp1d(time_Meteo, forcIP.prec)
        f_rel_humidity = interp1d(time_Meteo, forcIP.rHum)
        f_Spc_humidity = interp1d(time_Meteo, forcIP.hum)
        f_uDir = interp1d(time_Meteo, forcIP.uDir)

        class MeteoDataRaw_intp_Def():
            pass
        self.MeteoDataRaw_intp = MeteoDataRaw_intp_Def()
        self.MeteoDataRaw_intp.LWR_in = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.Dif_in = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.Dir_in = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.T_atm = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.windspeed_u = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.uDir = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.pressure_atm = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.rain = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.rel_humidity = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.Spc_humidity = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.Tdeepsoil = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.Year = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.Month = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.Day = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.Hour = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.Min = numpy.zeros(len(self.time_n))
        self.MeteoDataRaw_intp.Sec = numpy.zeros(len(self.time_n))
        count = 0
        for i in self.time_n:
            self.MeteoDataRaw_intp.LWR_in[count] = f_LWR_in(i)
            self.MeteoDataRaw_intp.Dif_in[count] = f_Dif_in(i)
            self.MeteoDataRaw_intp.Dir_in[count] = f_Dir_in(i)
            self.MeteoDataRaw_intp.T_atm[count] = f_T_atm(i)
            self.MeteoDataRaw_intp.windspeed_u[count] = f_windspeed_u(i)
            self.MeteoDataRaw_intp.uDir[count] = f_uDir(i)
            self.MeteoDataRaw_intp.pressure_atm[count] = f_pressure_atm(i)
            self.MeteoDataRaw_intp.rain[count] = f_rain(i)
            self.MeteoDataRaw_intp.rel_humidity[count] = f_rel_humidity(i) / 100
            self.MeteoDataRaw_intp.Spc_humidity[count] = f_Spc_humidity(i)
            self.MeteoDataRaw_intp.Year[count] = forcIP.Year[0]
            self.MeteoDataRaw_intp.Month[count] = forcIP.Month[0]
            self.MeteoDataRaw_intp.Day[count] = forcIP.Day[i // Meteo_time_step]
            self.MeteoDataRaw_intp.Hour[count] = forcIP.Hour[i // Meteo_time_step]
            self.MeteoDataRaw_intp.Min[count] = int((i % Meteo_time_step) / 60)
            self.MeteoDataRaw_intp.Sec[count] = 0
            self.MeteoDataRaw_intp.lat = weather.lat
            self.MeteoDataRaw_intp.lon = weather.lon
            self.MeteoDataRaw_intp.GMT = weather.GMT
            count = count + 1
        self.MeteoDataRaw_intp.Tdeepsoil = [Tdeepsoil_z[int(self.MeteoDataRaw_intp.Month[i]-1)] for i in range(len(self.time_n))]

    def is_near_zero(self,num,eps=1e-10):
        return abs(float(num)) < eps

    def instantiate_input(self):

        # -------------------------
        # Initialize energy balance
        # -------------------------
        CiCO2Roof_Sun_init = 400.0
        CiCO2Roof_Shade_init = 400.0
        self.EBRoof = EnergyBalanceRoof_Def(CiCO2Roof_Sun_init,CiCO2Roof_Shade_init)

        CiCO2Ground_Sun_init = 400.0
        CiCO2Ground_Shade_init = 400.0
        CiCO2Tree_Sun_init = 400.0
        CiCO2Tree_Shade_init = 400.0
        Tdepth_gimp_init = 300.0
        Tdepth_gbare_init = 300.0
        Tdepth_gveg_init = 300.0
        Ttree_init = 300.0
        self.EBCanyon = EnergyBalanceCanyon_Def(CiCO2Tree_Sun_init,CiCO2Tree_Shade_init,CiCO2Ground_Sun_init,CiCO2Ground_Shade_init,
                                                Tdepth_gimp_init,Tdepth_gbare_init,Tdepth_gveg_init,Ttree_init)

        Tdepth_rural_init = 300.0
        CiCO2Ground_Sun_Rural_init = 400.0
        CiCO2Ground_Shade_Rural_init = 400.0
        self.EBRural = EnergyBalanceRural_Def(Tdepth_rural_init,CiCO2Ground_Sun_Rural_init,CiCO2Ground_Shade_Rural_init)

        # --------------------------
        # Initialize building energy
        # --------------------------
        # Define BEM for each DOE type (read the fraction)
        # Open pickle file in binary form
        readDOE_file = open('readDOE.PKL', 'rb')
        refDOE = cPickle.load(readDOE_file)
        refBEM = cPickle.load(readDOE_file)
        refSchedule = cPickle.load(readDOE_file)
        readDOE_file.close()

        k = 0
        # Glazing ratio for total building stock
        r_glaze = 0
        # SHGC addition for total building stock
        SHGC = 0
        # total building floor area
        total_urban_bld_area = math.pow(self.charLength, 2)*self.Geometry_m.lambdap*self.Geometry_m.Height_canyon/self.BEMParam.h_floor

        self.BEM = []
        self.Sch = []
        # Loop over multiple building types and built eras and compute weighted average building parameters
        for i in range(len(self.bld)):
            for j in range(3):
                if self.bld[i][j] > 0.:
                    # Add to BEM list
                    self.BEM.append(refBEM[i][j][self.zone])
                    self.BEM[k].frac = self.bld[i][j]
                    self.BEM[k].fl_area = self.bld[i][j] * total_urban_bld_area

                    # Overwrite with optional parameters if provided
                    if self.BEMParam.glzR:
                        self.BEM[k].building.glazingRatio = self.BEMParam.glzR

                    # Keep track of total urban r_glaze, SHGC, and alb_wall for UCM model
                    r_glaze = r_glaze + self.BEM[k].frac * self.BEM[k].building.glazingRatio
                    SHGC = SHGC + self.BEM[k].frac * self.BEM[k].building.shgc

                    # Add to schedule list
                    self.Sch.append(refSchedule[i][j][self.zone])
                    k += 1

        for i in range(len(self.BEM)):
            if self.is_near_zero(self.BEMParam.autosize)==False:
                self.BEM[i].building.coolCap = 9999.
                self.BEM[i].building.heatCap = 9999.

        # ------------------------------
        # Initialize surface temperature
        # ------------------------------
        # number of soil layers
        ng = len(self.ParSoilGround.Zs) - 1
        nrur = len(self.RSMParam.Zs) - 1
        # Define layer thickness for ground
        layerThickness_GroundImp = [numpy.diff(self.ParSoilGround.Zs)[i]/1000 for i in range(ng)]
        layerThickness_GroundBare = [numpy.diff(self.ParSoilGround.Zs)[i]/1000 for i in range(ng)]
        layerThickness_GroundVeg = [numpy.diff(self.ParSoilGround.Zs)[i]/1000 for i in range(ng)]
        layerThickness_Rural = [numpy.diff(self.RSMParam.Zs)[i]/1000 for i in range(nrur)]

        Tlayers_Rural_init = self.MeteoDataRaw_intp.T_atm[0]
        Material_Rural = Material(self.RSMParam.lan_rural,self.RSMParam.cv_s_rural,'RuralMat')
        mat = [Material_Rural for i in range(len(layerThickness_Rural))]
        self.Rural = Tsurf_Def(layerThickness_Rural,mat,Tlayers_Rural_init,"Rural")

        Tlayers_GroundImp_init = 300.0
        Material_GImp = Material(self.ParThermalGround.lan_dry_imp, self.ParThermalGround.cv_s_imp, 'GImpMat')
        mat = [Material_GImp for i in range(len(layerThickness_GroundImp))]
        self.GroundImp = Tsurf_Def(layerThickness_GroundImp, mat, Tlayers_GroundImp_init, "GroundImp")

        Tlayers_GroundVeg_init = 300.0
        Material_GVeg = Material(self.ParThermalGround.lan_dry_veg, self.ParThermalGround.cv_s_veg, 'GVegMat')
        mat = [Material_GVeg for i in range(len(layerThickness_GroundVeg))]
        self.GroundVeg = Tsurf_Def(layerThickness_GroundVeg, mat, Tlayers_GroundVeg_init, "GroundVeg")

        Tlayers_GroundBare_init = 300.0
        Material_GBare = Material(self.ParThermalGround.lan_dry_bare, self.ParThermalGround.cv_s_bare, 'GBareMat')
        mat = [Material_GBare for i in range(len(layerThickness_GroundBare))]
        self.GroundBare = Tsurf_Def(layerThickness_GroundBare, mat, Tlayers_GroundBare_init, "GroundBare")

        # -----------------------------
        # Initialize MOST in rural site
        # -----------------------------
        Tprof_rural_init = self.MeteoDataRaw_intp.T_atm[0]
        Pprof_rural_init = self.MeteoDataRaw_intp.pressure_atm[0]
        S_rural_init = self.MeteoDataRaw_intp.windspeed_u[0]
        self.RSM = RSMDef(self.RSMParam,self.Geometry_m.z,self.Geometry_m.nz,self.Geometry_m.dz,Tprof_rural_init,Pprof_rural_init,S_rural_init)

        # ---------------------------------
        # Initialize 1D model in urban site
        # ---------------------------------
        Vx_urban_init = 0.1
        Vy_urban_init = 0.1
        TKE_urban_init = 0.15
        T_urban_init = 300.0
        Qn_urban_init = 0.01
        Pr_urban_init = self.MeteoDataRaw_intp.pressure_atm[0]
        rho_urban_init = 1.225
        Tref_urban = 300
        self.UCM = UCM_Def(Vx_urban_init,Vy_urban_init,TKE_urban_init,T_urban_init,Qn_urban_init,Pr_urban_init,rho_urban_init,Tref_urban,self.Geometry_m.nz)

        # ------------------------
        # Initialize water balance
        # ------------------------
        IntRoof_init = 0.0
        TERoof_init = 0.0
        ERoof_init = 0.0
        ExWaterRoof_init = 0.0
        if self.FractionsRoof.fimp == 1:
            OwaterRoof_init = 0
            VRoof_init = [0*self.ParSoilRoof.dz[i] for i in range(len(self.ParSoilRoof.dz))]
        else:
            OwaterRoof_init = self.ParSoilRoof.O33
            VRoof_init = [self.ParSoilRoof.O33*self.ParSoilRoof.dz[i] for i in range(len(self.ParSoilRoof.dz))]
        SoilPotWRoof_init = 0.0
        self.WBRoof = WaterBalanceRoof_Def(IntRoof_init, TERoof_init, OwaterRoof_init, ERoof_init, self.ParSoilRoof,
                                           ExWaterRoof_init,VRoof_init,SoilPotWRoof_init)

        IntGround_init = numpy.float64(0.0)
        RunonGround_init = 0.0
        if self.FractionsGround.fimp == 1:
            OwaterGround_init = 0
            VCanyon_init = [0* self.ParSoilGround.dz[i] for i in range(len(self.ParSoilGround.dz))]
        else:
            OwaterGround_init = self.ParSoilGround.O33
            VCanyon_init = [self.ParSoilGround.O33 * self.ParSoilGround.dz[i] for i in range(len(self.ParSoilGround.dz))]
        ExWaterCanyon_init = 0.0
        self.WBCanyon = WaterBalanceCanyon_Def(OwaterGround_init,IntGround_init,RunonGround_init,self.ParSoilGround,VCanyon_init,ExWaterCanyon_init)

    def CheckInputs(self):
        # Check validity of input parameters
        if self.RSMParam.u_star_min_MOST < 0.1:
            print('Error : Minimum friction velocity in the rural area is less than 0.1 [m s^-1]. Please check "u_star_min_MOST" in the input file.')
            quit()
        if self.RSMParam.zToverz0_MOST < 0.1 or self.RSMParam.zToverz0_MOST > 1:
            print('Error : Thermodynamic roughness length over aerodynamic roughness length is unrealistic. Please check "zToverz0_MOST" in the input file.')
            quit()
        if self.RSMParam.dispoverh_MOST < 0.2 or self.RSMParam.dispoverh_MOST > 1:
            print('Error : Displacement height over obstacle height is unrealistic. Please check "dispoverh_MOST" in the input file.')
            quit()
        if self.Geometry_m.nz*self.Geometry_m.dz > 4*self.Geometry_m.Height_canyon or self.Geometry_m.nz*self.Geometry_m.dz < 2*self.Geometry_m.Height_canyon:
            print('Error : Domain height can not be higher than four times of building height or less than two times of building height. Please check "Height_canyon", "nz", and "dz" in the input file. ')
            quit()
        if self.Geometry_m.nz*self.Geometry_m.dz > 150:
            print('Error : Rural model may not be valid for this domain height. Please check "Height_canyon", "nz", and "dz" in the input file. ')
            quit()
        if self.TimeParam.dts < 60 or self.TimeParam.dts > 600:
            print('Error : Simulation time step is too fine or coarse. Please check "dts" in the input file.')
            quit()
        if self.Geometry_m.Width_roof/self.Geometry_m.Width_canyon > 3 or self.Geometry_m.Width_roof/self.Geometry_m.Width_canyon < 0.3:
            print('Error: building width to street width ratio may be out of range. Please check "Width_roof" and "Width_canyon" in the input file.')
            quit()
        if self.Geometry_m.theta_canyon > 90 or self.Geometry_m.theta_canyon < -90:
            print('Error: Canyon orientation must be between -90 and 90 degrees. Please check "theta_canyon" in the input file.')
            quit()
        if self.Geometry_m.Height_canyon > 0.5*self.Geometry_m.nz*self.Geometry_m.dz or self.Geometry_m.Height_canyon < 3:
            print('Error: Building height is out of range. Please check "Height_canyon" in the input file.')
            quit()
        if self.RSMParam.z0overh_MOST > 0.5 or self.RSMParam.z0overh_MOST < 0.05:
            print('Error: Aerodynamic roughness length over obstacle height is out of range. Please check "z0overh_MOST" in the input file.')
            quit()
        if self.RSMParam.WindMin_MOST > 0.7 or self.RSMParam.WindMin_MOST < 0.05:
            print('Error: Minimum wind speed for MOST is out of range. Please check "WindMin_MOST" in the input file.')
            quit()
        if self.RSMParam.h_obs > 10 or self.RSMParam.h_obs < 0.1:
            print('Error: Rural average obstacle height is out of range. Please check "h_obs" in the input file.')
            quit()
        if self.RSMParam.Bowen > 10 or self.RSMParam.Bowen < -10:
            print('Error: Bowen ratio in the rural area is out of range. Please check "Bowen" in the input file.')
            quit()
        if self.RSMParam.MinWind_rural > 0.7 or self.RSMParam.MinWind_rural < 0.05:
            print('Error: Minimum wind speed for rural energy balance is out of range. Please check "MinWind_rural" in the input file.')
            quit()
        if self.Geometry_m.Radius_tree > self.Geometry_m.Distance_tree:
            print('Error: Radius of the tree is greater than the distance of tree from the wall. Please check tree parameters in the input file.')
            quit()
        if self.ColParam.h_LAD[-1] > self.Geometry_m.Height_canyon:
            print('Error: Tree is higher than the building. Please check tree parameters in the input file.')
            quit()
        if 2*2*self.Geometry_m.Radius_tree > self.Geometry_m.Width_canyon:
            print('Error: Tree does not fit in the canyon. Please check tree parameters in the input file.')
            quit()
        if self.ColParam.WindMin_Urban > 0.7 or self.ColParam.WindMin_Urban < 0.01:
            print('Error: Minimum wind speed for urban site is out of range. Please check "WindMin_Urban" in the input file.')
            quit()
        if self.Geometry_m.dz > 5 or self.Geometry_m.dz < 1:
            print('Error: Vertical discretization should be between 1 and 5 m. Please check "dz" in the input file.')
            quit()
        if self.ParSoilRoof.Zs[-1] < 50 or self.ParSoilRoof.Zs[-1] > 500:
            print('Error: Soil depth at the roof is too shallow for accomodatig roots or too deep. Please check "Zs_R" in the input file.')
            quit()
        if self.ParSoilGround.Zs[-1] < 50 or self.ParSoilGround.Zs[-1] > 2500:
            print('Error: Soil depth at the ground is too shallow for accomodatig roots or too deep. Please check "Zs_G" in the input file.')
            quit()

    def Simulate(self):

        # total number of hours in simulation
        self.N = int(self.simTime.days * 24)
        # weather time step counter
        n = 0
        # Define output data
        self.EBRoofData = [None for x in range(self.N)]
        self.EBCanyonData = [None for x in range(self.N)]
        self.EBRuralData = [None for x in range(self.N)]
        self.RuralData = [None for x in range(self.N)]
        self.RoofImpData = [None for x in range(self.N)]
        self.RoofVegData = [None for x in range(self.N)]
        self.GroundImpData = [None for x in range(self.N)]
        self.GroundVegData = [None for x in range(self.N)]
        self.GroundBareData = [None for x in range(self.N)]
        self.WallSunData = [None for x in range(self.N)]
        self.WallShadeData = [None for x in range(self.N)]
        self.RuralGroundImpData = [None for x in range(self.N)]
        self.RuralGroundBareData = [None for x in range(self.N)]
        self.RuralGroundVegData = [None for x in range(self.N)]
        self.RSMData = [None for x in range(self.N)]
        self.UBLData = [None for x in range(self.N)]
        self.UCMData = [None for x in range(self.N)]
        self.WBRoofData = [None for x in range(self.N)]
        self.WBCanyonData = [None for x in range(self.N)]
        self.WBRuralData = [None for x in range(self.N)]
        self.ForcingData = [None for x in range(self.N)]
        self.BEMData = [None for x in range(self.N)]
        self.time = [None for x in range(self.N)]

        # Define surface temperature of two steps back
        TwallSun_ext = 0
        TwallShade_ext = 0
        TwallSun_int = 0
        TwallShade_int = 0
        # Compute weighted average wall temperatures for various building types
        for i in range(len(self.BEM)):
            TwallSun_ext = TwallSun_ext + self.BEM[i].wallSun.Text * self.BEM[i].frac
            TwallShade_ext = TwallShade_ext + self.BEM[i].wallShade.Text * self.BEM[i].frac
            TwallSun_int = TwallSun_int + self.BEM[i].wallSun.Tint * self.BEM[i].frac
            TwallShade_int = TwallShade_int + self.BEM[i].wallShade.Tint * self.BEM[i].frac
        class Tsurf_2back_Def():
            pass
        Tsurf_2back = Tsurf_2back_Def()
        Tsurf_2back.TGroundImp = self.GroundImp.Text
        Tsurf_2back.TGroundBare = self.GroundBare.Text
        Tsurf_2back.TGroundVeg = self.GroundVeg.Text
        Tsurf_2back.TTree = self.EBCanyon.Ttree
        Tsurf_2back.TWallSun = TwallSun_ext
        Tsurf_2back.TWallShade = TwallShade_ext
        Tsurf_2back.TWallIntSun = TwallSun_int
        Tsurf_2back.TWallIntShade = TwallShade_int
        TGroundImp_2last = numpy.zeros(2)
        TGroundBare_2last = numpy.zeros(2)
        TGroundVeg_2last = numpy.zeros(2)
        TTree_2last = numpy.zeros(2)
        TWallSun_2last = numpy.zeros(2)
        TWallShade_2last = numpy.zeros(2)
        TWallIntSun_2last = numpy.zeros(2)
        TWallIntShade_2last = numpy.zeros(2)

        # Start simulation
        for it in range(0,self.simTime.nt-1,1):
            print(r'Progress [%]', numpy.round(100 * it / self.simTime.nt, 2))

            # Simulation time increment raised to weather time step
            SunPosition,MeteoData,Anthropogenic,location,ParCalculation = \
                ForcingData(self.MeteoDataRaw_intp,it, self.WBCanyon.SoilPotW, self.VCWGParamFileName,self.simTime)
            self.simTime.UpdateDate()

            #----------------------
            # Update energy balance
            # ---------------------
            # Roof
            TroofImp_ext = 0
            TroofVeg_ext = 0
            # Compute weighted average roof temperatures for various building types
            for i in range(len(self.BEM)):
                TroofImp_ext = TroofImp_ext + self.BEM[i].roofImp.Text*self.BEM[i].frac
                TroofVeg_ext = TroofVeg_ext + self.BEM[i].roofVeg.Text*self.BEM[i].frac
            TemperatureR = [TroofImp_ext,TroofVeg_ext]
            self.EBRoof.EBSolver_Roof(TemperatureR,MeteoData,self.WBRoof.Int,self.WBRoof.ExWaterRoof_L,self.WBRoof.VRoofSoil,
                                      self.WBRoof.Owater_OwRoofSoilVeg,self.WBRoof.SoilPotWRoof_L,self.Geometry_m,self.FractionsRoof,
                                      self.ParSoilRoof,self.PropOpticalRoof,self.ParVegRoof,ParCalculation,self.UCM.VerticalProfUrban)

            # Canyon
            # Wall temperatures are weighted-averaged according to building type
            TwallSun_ext = 0
            TwallShade_ext = 0
            TwallSun_int = 0
            TwallShade_int = 0
            for i in range(len(self.BEM)):
                TwallSun_ext = TwallSun_ext + self.BEM[i].wallSun.Text*self.BEM[i].frac
                TwallShade_ext = TwallShade_ext + self.BEM[i].wallShade.Text*self.BEM[i].frac
                TwallSun_int = TwallSun_int + self.BEM[i].wallSun.Tint*self.BEM[i].frac
                TwallShade_int = TwallShade_int + self.BEM[i].wallShade.Tint*self.BEM[i].frac

            TemperatureC = [self.GroundImp.Text,self.GroundBare.Text,self.GroundVeg.Text,TwallSun_ext,TwallShade_ext,
                            self.EBCanyon.Ttree[0]]

            # Update surface temperature from two steps back if the force-restore method is used for ground
            if it > 1:
                if it % 2 == 0:
                    Tsurf_2back.TGroundImp = TGroundImp_2last[0]
                    Tsurf_2back.TGroundBare = TGroundBare_2last[0]
                    Tsurf_2back.TGroundVeg = TGroundVeg_2last[0]
                    Tsurf_2back.TTree = TTree_2last[0]
                    Tsurf_2back.TWallSun = TWallSun_2last[0]
                    Tsurf_2back.TWallShade = TWallShade_2last[0]
                    Tsurf_2back.TWallIntSun = TWallIntSun_2last[0]
                    Tsurf_2back.TWallIntShade = TWallIntShade_2last[0]
                else:
                    Tsurf_2back.TGroundImp = TGroundImp_2last[1]
                    Tsurf_2back.TGroundBare = TGroundBare_2last[1]
                    Tsurf_2back.TGroundVeg = TGroundVeg_2last[1]
                    Tsurf_2back.TTree = TTree_2last[1]
                    Tsurf_2back.TWallSun = TWallSun_2last[1]
                    Tsurf_2back.TWallShade = TWallShade_2last[1]
                    Tsurf_2back.TWallIntSun = TWallIntSun_2last[1]
                    Tsurf_2back.TWallIntShade = TWallIntShade_2last[1]
            self.EBCanyon.EBSolver_Canyon(TemperatureC,Tsurf_2back,MeteoData,self.WBCanyon.Int,self.WBCanyon.ExWater,
                                          self.WBCanyon.Vwater,self.WBCanyon.Owater,self.WBCanyon.SoilPotW,self.ViewFactor,
                                          self.Geometry_m,self.ParTree,self.geometry,self.FractionsGround,
                                          self.ParSoilGround,self.ParInterceptionTree,self.PropOpticalGround,self.PropOpticalWall,
                                          self.PropOpticalTree,self.ParThermalGround,self.ParVegGround,
                                          self.ParVegTree, SunPosition, Anthropogenic, ParCalculation,
                                          self.UCM.VerticalProfUrban,self.ColParam)

            if it % 2 == 0:
                TGroundImp_2last[0] = self.GroundImp.Text
                TGroundBare_2last[0] = self.GroundBare.Text
                TGroundVeg_2last[0] = self.GroundVeg.Text
                TTree_2last[0] = self.EBCanyon.Ttree
                TWallSun_2last[0] = TwallSun_ext
                TWallShade_2last[0] = TwallShade_ext
                TWallIntSun_2last[0] = TwallSun_int
                TWallIntShade_2last[0] = TwallShade_int
            else:
                TGroundImp_2last[1] = self.GroundImp.Text
                TGroundBare_2last[1] = self.GroundBare.Text
                TGroundVeg_2last[1] = self.GroundVeg.Text
                TTree_2last[1] = self.EBCanyon.Ttree
                TWallSun_2last[1] = TwallSun_ext
                TWallShade_2last[1] = TwallShade_ext
                TWallIntSun_2last[1] = TwallSun_int
                TWallIntShade_2last[1] = TwallShade_int

            # Rural
            self.EBRural.EBSolver_Rural(MeteoData,self.RSMParam,self.Rural.Text,SunPosition,self.simTime,ParCalculation,self.RSM,self.ParThermalGround)
            # Update rural surface temperature [K]
            self.Rural.Element(self.EBRural.EnergyFlux.SWRabsRural,self.EBRural.EnergyFlux.LWRabsRural,self.EBRural.EnergyFlux.LEfluxRural,
                               self.EBRural.EnergyFlux.HfluxRural,self.TimeParam.dts,self.EBRural.Td,2)
            Text_Rrual = self.Rural.Text

            # -----------------------------------------------------
            # Update rural model and boundary conditions at the top
            # -----------------------------------------------------
            if self.RSMParam.Rural_Model_name == 'MOST':
                # Option 1 (Rural_Model_name): Monin-Obukhov Similarity Theory (MOST) is used
                self.RSM.MOST(self.EBRural,Text_Rrual,MeteoData)
                # Update boundary conditions at the top of domain
                TOPBC_Urban_th = self.RSM.T_rural[-1]
                TOPBC_Urban_q = self.RSM.q_rural[-1]
                TOPBC_Urban_WindSpeed = MeteoData.Uatm
                TOPBC_Urban_WindDir = MeteoData.uDir
                Ustar = copy.copy(self.RSM.u_star)

            elif self.RSMParam.Rural_Model_name == 'Forcing_extFile':
                # Option 2 (Rural_Model_name): Input forcing variables are used at the top of the domain and no rural model is used.
                # Update boundary conditions at the top of domain
                TOPBC_Urban_th = MeteoData.Tatm
                TOPBC_Urban_q = MeteoData.q_atm
                TOPBC_Urban_WindSpeed = MeteoData.Uatm
                TOPBC_Urban_WindDir = MeteoData.uDir
                Ustar = 0

            # ------------------
            # Update urban model
            # ------------------
            self.UCM.UCMCal(TemperatureR,TemperatureC,self.FractionsGround,self.FractionsRoof,TOPBC_Urban_th,TOPBC_Urban_q,
                            TOPBC_Urban_WindSpeed,TOPBC_Urban_WindDir,self.Geometry_m,ParCalculation,self.ParVegGround,
                            self.ParVegRoof,self.EBRoof.Src.thb,self.EBCanyon.Src.thb,self.EBCanyon.Src.tvb,self.EBRoof.Src.qhb,
                            self.EBCanyon.Src.qhb,self.EBCanyon.Hflux,self.EBCanyon.LEflux,self.EBRoof.Hflux,self.EBRoof.LEflux,
                            self.geometry,self.ParTree,self.ColParam,self.BEM,Ustar,self.RSMParam.Rural_Model_name)

            # ----------------------------
            # Update building energy model
            # ----------------------------
            # Update building & traffic schedule
            # Assign day type (1 = weekday, 2 = sat, 3 = sun/other)
            if self.is_near_zero(self.simTime.julian % 7):
                self.dayType = 3  # Sunday
            elif self.is_near_zero(self.simTime.julian % 7 - 6.):
                self.dayType = 2  # Saturday
            else:
                self.dayType = 1  # Weekday
            # Update the energy components for building types
            for i in range(len(self.BEM)):
                # Set point temperature [K]
                # Add from temperature schedule for cooling
                self.BEM[i].building.coolSetpointDay = self.Sch[i].Cool[self.dayType - 1][self.simTime.hourDay] + 273.15
                self.BEM[i].building.coolSetpointNight = self.BEM[i].building.coolSetpointDay
                # Add from temperature schedule for heating
                self.BEM[i].building.heatSetpointDay = self.Sch[i].Heat[self.dayType - 1][self.simTime.hourDay] + 273.15
                self.BEM[i].building.heatSetpointNight = self.BEM[i].building.heatSetpointDay

                # Internal Heat Load Schedule per unit floor area [W m^-2]
                # Electricity consumption per unit floor area [W m^-2] = max for electrical plug process * electricity fraction for the day
                self.BEM[i].Elec = self.Sch[i].Qelec * self.Sch[i].Elec[self.dayType - 1][self.simTime.hourDay]
                # Lighting per unit floor area [W m^-2] = max for light * light fraction for the day
                self.BEM[i].Light = self.Sch[i].Qlight * self.Sch[i].Light[self.dayType - 1][self.simTime.hourDay]
                # Number of occupants x occ fraction for day
                self.BEM[i].Nocc = self.Sch[i].Nocc * self.Sch[i].Occ[self.dayType - 1][self.simTime.hourDay]
                # Sensible Q occupant * fraction occupant sensible Q * number of occupants
                self.BEM[i].Qocc = self.BEMParam.sensOcc*(1-self.BEMParam.LatFOcc)*self.BEM[i].Nocc

                # SWH and ventilation schedule
                # Solar water heating per unit floor area [W m^-2] = Peak Service Hot Water per unit floor [kg hr^-1 m^-2] * SWH fraction for the day
                self.BEM[i].SWH = self.Sch[i].Vswh * self.Sch[i].SWH[self.dayType - 1][self.simTime.hourDay]
                # Ventilation rate per unit floor area [m^3 s^-1 m^-2]
                self.BEM[i].building.vent = self.Sch[i].Vent
                # Gas consumption per unit floor area [W m^-2] = max for gas * Gas fraction for the day
                self.BEM[i].Gas = self.Sch[i].Qgas * self.Sch[i].Gas[self.dayType - 1][self.simTime.hourDay]

                # This is quite messy, should update
                # Update internal heat and corresponding fractional loads per unit floor area [W m^-2]
                intHeat = self.BEM[i].Light + self.BEM[i].Elec + self.BEM[i].Qocc
                self.BEM[i].building.intHeatDay = intHeat
                self.BEM[i].building.intHeatNight = intHeat
                # Fraction of radiant heat from light and equipment of whole internal heat per unit floor area [W m^-2]
                self.BEM[i].building.intHeatFRad = (self.BEMParam.RadFLight*self.BEM[i].Light+self.BEMParam.RadFEquip*self.BEM[i].Elec)/intHeat
                # fraction of latent heat (from occupants) of whole internal heat per unit floor area [W m^-2]
                self.BEM[i].building.intHeatFLat = self.BEMParam.LatFOcc*self.BEMParam.sensOcc*self.BEM[i].Nocc/intHeat

                # Update envelope temperature layers [K]
                # Wall temperature exposed to outdoor environment [K]
                self.BEM[i].T_wallex = (self.BEM[i].wallSun.Text + self.BEM[i].wallShade.Text)/2
                # Wall temperature exposed to indoor environment [K]
                self.BEM[i].T_wallin = (self.BEM[i].wallSun.Tint + self.BEM[i].wallShade.Tint)/2
                # Roof temperature exposed to outdoor environment [K]
                self.BEM[i].T_roofex = self.FractionsRoof.fimp*self.BEM[i].roofImp.Text+self.FractionsRoof.fveg*self.BEM[i].roofVeg.Text
                # Roof temperature exposed to indoor environment [K]
                self.BEM[i].T_roofin = self.FractionsRoof.fimp*self.BEM[i].roofImp.Tint+self.FractionsRoof.fveg*self.BEM[i].roofVeg.Tint

                # Calculate one-point temperature and humidity in the canyon: Using 1-D profiles in the canyon
                canTemp = numpy.mean(self.UCM.VerticalProfUrban.th[0:self.Geometry_m.nz_u])
                canHum = numpy.mean(self.UCM.VerticalProfUrban.qn[0:self.Geometry_m.nz_u])
                self.BEM[i].building.BEMCalc(canTemp,canHum,self.BEM[i],MeteoData,ParCalculation,self.simTime,self.Geometry_m,
                                             self.FractionsRoof,self.EBCanyon.SWR)

                # Electricity consumption of urban area [W]
                self.BEM[i].ElecTotal = self.BEM[i].building.ElecTotal * self.BEM[i].fl_area

                # Update surface temperature of building surfaces
                # Mass
                self.BEM[i].mass.Element(0,0,0,0,self.TimeParam.dts,0.,1,self.BEM[i].building.fluxMass,self.BEM[i].building.fluxMass)
                # Roof
                if self.FractionsRoof.fimp > 0:
                    self.BEM[i].roofImp.Element(self.EBRoof.SWR.SWRabsRoofImp,self.EBRoof.LWR.LWRabsRoofImp,self.EBRoof.LEflux.LEfluxRoofImp,
                                                self.EBRoof.Hflux.HfluxRoofImp,self.TimeParam.dts,0.,1,None,self.BEM[i].building.fluxRoof)
                if self.FractionsRoof.fveg > 0:
                    self.BEM[i].roofVeg.Element(self.EBRoof.SWR.SWRabsRoofVeg,self.EBRoof.LWR.LWRabsRoofVeg,self.EBRoof.LEflux.LEfluxRoofVeg,
                                                self.EBRoof.Hflux.HfluxRoofVeg,self.TimeParam.dts,0.,1,None,self.BEM[i].building.fluxRoof)
                # Walls
                self.BEM[i].wallSun.Element(self.EBCanyon.SWR.SWRabs.SWRabsWallSun,self.EBCanyon.LWR.LWRabs.LWRabsWallSun,
                                            self.EBCanyon.LEflux.LEfluxWallSun,self.EBCanyon.Hflux.HfluxWallSun,self.TimeParam.dts,
                                            0.,1,None,self.BEM[i].building.fluxWall)
                self.BEM[i].wallShade.Element(self.EBCanyon.SWR.SWRabs.SWRabsWallShade,self.EBCanyon.LWR.LWRabs.LWRabsWallShade,
                                              self.EBCanyon.LEflux.LEfluxWallShade,self.EBCanyon.Hflux.HfluxWallShade,self.TimeParam.dts,
                                              0.,1,None,self.BEM[i].building.fluxWall)

            # -----------------------------------
            # Update outdoor surface temperatures
            # -----------------------------------
            if self.FractionsGround.fimp > 0:
                self.GroundImp.Element(self.EBCanyon.SWR.SWRabs.SWRabsGroundImp,self.EBCanyon.LWR.LWRabs.LWRabsGroundImp,
                                       self.EBCanyon.LEflux.LEfluxGroundImp,self.EBCanyon.Hflux.HfluxGroundImp,self.TimeParam.dts,
                                       self.EBCanyon.Tdepth.TDampGroundImp,2)
            if self.FractionsGround.fveg > 0:
                self.GroundVeg.Element(self.EBCanyon.SWR.SWRabs.SWRabsGroundVeg,self.EBCanyon.LWR.LWRabs.LWRabsGroundVeg,
                                       self.EBCanyon.LEflux.LEfluxGroundVeg,self.EBCanyon.Hflux.HfluxGroundVeg,self.TimeParam.dts,
                                       self.EBCanyon.Tdepth.TDampGroundVeg,2)
            if self.FractionsGround.fbare > 0:
                self.GroundBare.Element(self.EBCanyon.SWR.SWRabs.SWRabsGroundBare,self.EBCanyon.LWR.LWRabs.LWRabsGroundBare,
                                        self.EBCanyon.LEflux.LEfluxGroundBare,self.EBCanyon.Hflux.HfluxGroundBare,self.TimeParam.dts,
                                        self.EBCanyon.Tdepth.TDampGroundBare,2)

            #-----------------------
            # Update hydrology model
            #-----------------------
            # Roof
            self.WBRoof.WBSolver_Roof(self.EBRoof.Eflux.EfluxRoofImp,self.EBRoof.Eflux.EfluxRoofVegInt,self.EBRoof.Eflux.EfluxRoofVegPond,
                                      self.EBRoof.Eflux.EfluxRoofVegSoil,self.EBRoof.Eflux.TEfluxRoofVeg,MeteoData,self.FractionsRoof,
                                      self.ParSoilRoof,ParCalculation,self.ParVegRoof,Anthropogenic)
            # Canyon
            self.WBCanyon.WBSolver_Canyon(MeteoData,self.EBCanyon.Eflux.EfluxTreeInt,self.EBCanyon.Eflux.EfluxGroundVegInt,
                                          self.EBCanyon.Eflux.EfluxGroundImp,self.EBCanyon.Eflux.EfluxGroundBarePond,
                                          self.EBCanyon.Eflux.EfluxGroundVegPond, self.EBCanyon.Eflux.EfluxGroundBareSoil,
                                          self.EBCanyon.Eflux.EfluxGroundVegSoil,self.EBCanyon.Eflux.TEfluxGroundVeg,
                                          self.EBCanyon.Eflux.TEfluxTree, self.ParSoilGround, self.ParInterceptionTree,
                                          ParCalculation,self.ParVegGround,self.ParVegTree,self.FractionsGround,self.geometry,
                                          self.ParTree,self.Geometry_m,Anthropogenic)

            #-----------------------
            # Update total variables
            #-----------------------
            # Calculate total shortwave radiations [W m^-2]
            self.EBCanyon.SWR.SWRabsTotalUrban = self.geometry.wroof_norm*self.EBRoof.SWR.SWRabsTotalRoof + \
                                                 self.geometry.wcanyon_norm*self.EBCanyon.SWR.SWRabs.SWRabsTotalCanyon
            self.EBCanyon.SWR.SWRinTotalUrban = self.geometry.wroof_norm*self.EBRoof.SWR.SWRinTotalRoof + \
                                                self.geometry.wcanyon_norm*self.EBCanyon.SWR.SWRin.SWRinTotalCanyon
            self.EBCanyon.SWR.SWRoutTotalUrban = self.geometry.wroof_norm*self.EBRoof.SWR.SWRoutTotalRoof + \
                                                 self.geometry.wcanyon_norm*self.EBCanyon.SWR.SWRout.SWRoutTotalCanyon
            self.EBCanyon.SWR.SWREBTotalUrban = self.geometry.wroof_norm*self.EBRoof.SWR.SWREBTotalRoof + \
                                                self.geometry.wcanyon_norm*self.EBCanyon.SWR.SWREB.SWREBTotalCanyon
            # Calculate total longwave radiations [W m^-2]
            self.EBCanyon.LWR.LWRabsTotalUrban = self.geometry.wroof_norm*self.EBRoof.LWR.LWRabsTotalRoof + \
                                                 self.geometry.wcanyon_norm*self.EBCanyon.LWR.LWRabs.LWRabsTotalCanyon
            self.EBCanyon.LWR.LWRinTotalUrban = self.geometry.wroof_norm*self.EBRoof.LWR.LWRinTotalRoof + \
                                                self.geometry.wcanyon_norm*self.EBCanyon.LWR.LWRin.LWRinTotalCanyon
            self.EBCanyon.LWR.LWRoutTotalUrban = self.geometry.wroof_norm*self.EBRoof.LWR.LWRoutTotalRoof + \
                                                 self.geometry.wcanyon_norm*self.EBCanyon.LWR.LWRout.LWRoutTotalCanyon
            self.EBCanyon.LWR.LWREBTotalUrban = self.geometry.wroof_norm*self.EBRoof.LWR.LWREBTotalRoof + \
                                                self.geometry.wcanyon_norm*self.EBCanyon.LWR.LWREB.LWREBTotalCanyon
            # Calculate total runon in the urban [mm s^-1]
            self.WBCanyon.Runon.RunonUrban = self.geometry.wroof_norm*self.WBRoof.RunonRoofTot + self.geometry.wcanyon_norm*self.WBCanyon.Runon.RunonGroundTot
            # Calculate total runoff in the urban [mm s^-1]
            self.WBCanyon.Runoff.RunoffUrban = self.geometry.wroof_norm*self.WBRoof.RunoffRoofTot + self.geometry.wcanyon_norm*self.WBCanyon.Runoff.RunoffGroundTot

            # Save simulation results in an hourly basis
            if self.is_near_zero(self.simTime.secDay % self.simTime.timePrint) and n < self.N:

                self.EBRoofData[n] = copy.deepcopy(self.EBRoof)
                self.EBCanyonData[n] = copy.deepcopy(self.EBCanyon)
                self.EBRuralData[n] = copy.deepcopy(self.EBRural)
                self.RoofImpData[n] = copy.deepcopy(self.BEM[0].roofImp)
                self.RoofVegData[n] = copy.deepcopy(self.BEM[0].roofVeg)
                self.GroundImpData[n] = copy.deepcopy(self.GroundImp)
                self.GroundVegData[n] = copy.deepcopy(self.GroundVeg)
                self.GroundBareData[n] = copy.deepcopy(self.GroundBare)
                self.WallSunData[n] = copy.deepcopy(self.BEM[0].wallSun)
                self.WallShadeData[n] = copy.deepcopy(self.BEM[0].wallShade)
                self.RSMData[n] = copy.deepcopy(self.RSM)
                self.UCMData[n] = copy.deepcopy(self.UCM)
                self.WBRoofData[n] = copy.deepcopy(self.WBRoof)
                self.WBCanyonData[n] = copy.deepcopy(self.WBCanyon)
                self.ForcingData[n] = copy.deepcopy(MeteoData)
                self.BEMData[n] = copy.deepcopy(self.BEM)
                self.time[n] = self.time_n[it]
                self.RuralData[n] = copy.deepcopy(self.Rural)

                n += 1

    def write_output(self):

        Output_dir =  "Results"

        Write_Forcing(self.case,self.ForcingData,self.time,Output_dir)

        Write_EB(self.case,self.FractionsRoof, self.FractionsGround, self.ParTree, self.RSMParam, self.EBRoofData,self.EBCanyonData,self.EBRuralData,self.UCMData,self.time,Output_dir)

        Write_Tsurf(self.case, self.FractionsRoof, self.FractionsGround, self.ParTree,self.RoofImpData,self.RoofVegData,self.GroundImpData,self.GroundBareData,self.GroundVegData,self.WallSunData,
                    self.WallShadeData,self.EBCanyonData,self.RuralData,self.RuralGroundImpData,self.RuralGroundBareData,self.RuralGroundVegData,
                    self.RSMParam,self.time,Output_dir)

        Write_WB(self.case,self.FractionsRoof, self.FractionsGround, self.ParTree,self.WBRoofData,self.WBCanyonData,self.time,Output_dir)

        Write_BEM(self.BEMData,self.time,self.case,Output_dir)

        # Generate output text file for deep ground temperature profiles in the urban area
        # SurfType = 1 (GroundImp), 2 (GroundVeg), 3 (GroundBare)
        Write_TdeepProfiles("Tdeep_imp",self.FractionsGround,1,self.GroundImpData,self.GroundImpData[0].z_depth[1:],self.time,self.case,Output_dir)
        # Generate output text file for deep vegetated ground temperature profile
        Write_TdeepProfiles("Tdeep_veg",self.FractionsGround,2,self.GroundVegData,self.GroundVegData[0].z_depth[1:],self.time,self.case,Output_dir)
        # Generate output text file for deep bare ground temperature profile
        Write_TdeepProfiles("Tdeep_bare",self.FractionsGround,3,self.GroundBareData,self.GroundBareData[0].z_depth[1:],self.time,self.case,Output_dir)

        # Generate output text file for Qn profile
        Write_1Dprofiles("Qn",self.UCMData,'VerticalProfUrban','qn',self.Geometry_m.z,self.time,self.case,Output_dir)
        # Generate output text file for tke profile
        Write_1Dprofiles("TKE",self.UCMData,'VerticalProfUrban','tke',self.Geometry_m.z,self.time,self.case,Output_dir)
        # Generate output text file for th profile
        Write_1Dprofiles("th",self.UCMData,'VerticalProfUrban','th',self.Geometry_m.z,self.time,self.case,Output_dir)
        # Generate output text file for vx profile
        Write_1Dprofiles("vx",self.UCMData,'VerticalProfUrban','vx',self.Geometry_m.z,self.time,self.case,Output_dir)
        # Generate output text file for vy profile
        Write_1Dprofiles("vy",self.UCMData,'VerticalProfUrban','vy',self.Geometry_m.z,self.time,self.case,Output_dir)
        # Generate output text file for S profile
        Write_1Dprofiles("S",self.UCMData,'VerticalProfUrban','s',self.Geometry_m.z,self.time,self.case,Output_dir)
        # Generate output text file for LE profile
        Write_1Dprofiles("LEfluxProf", self.UCMData, 'VerticalProfUrban', 'LEflux', self.Geometry_m.z, self.time, self.case,Output_dir)
        # Generate output text file for H profile
        Write_1Dprofiles("HfluxProf", self.UCMData, 'VerticalProfUrban', 'Hflux', self.Geometry_m.z, self.time, self.case,Output_dir)

        # Generate output text file for th profile in rural area
        Write_Ruralprofiles("th_rural",self.RSMParam.Rural_Model_name,self.RSMData,'T_rural',self.RSMData[0].z,self.time,self.case,Output_dir)
        # Generate output text file for q profile in rural area
        Write_Ruralprofiles("Qn_rural",self.RSMParam.Rural_Model_name,self.RSMData,'q_rural',self.RSMData[0].z,self.time,self.case,Output_dir)
        # Generate output text file for P profile in rural area
        Write_Ruralprofiles("P_rural",self.RSMParam.Rural_Model_name,self.RSMData,'presProf',self.RSMData[0].z,self.time,self.case,Output_dir)

    def run(self):

        self.read_input()
        self.read_epw()
        self.instantiate_input()
        self.CheckInputs()
        self.Simulate()
        self.write_output()

