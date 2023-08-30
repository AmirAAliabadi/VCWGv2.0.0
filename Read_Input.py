import os
import numpy
import math
from pprint import pprint
from scipy.interpolate import interp1d
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import matplotlib.gridspec as gridspec
import pandas as pd
from Soil_Functions import Soil_Calculations
from datetime import datetime
from datetime import datetime

"""
Read input variables
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: February 2021
"""

def read_VCWG_param(VCWG_param_file_path):
    # Open .uwg file and feed csv data to initializeDataFile

    with open(VCWG_param_file_path) as f:
        lines = f.readlines()
    VCWG_param = []
    for i in range(len(lines)):
        VCWG_param.append(list(lines[i].split(",")))
    # The initialize.uwg is read with a dictionary so that users changing
    # line endings or line numbers doesn't make reading input incorrect
    _init_param_dict = {}
    count = 0
    while count < len(VCWG_param):
        row = VCWG_param[count]
        row = [row[i].replace(" ", "") for i in range(len(row))]  # strip white spaces
        if len(row) == 1 and row == ['\n']:
            count += 1
            continue
        # Optional parameters might be empty so handle separately
        if row == [] or "#" in row[0]:
            count += 1
            continue
        elif row[0] == "SchTraffic":
            # SchTraffic: 3 x 24 matrix
            trafficrows = VCWG_param[count + 1:count + 4]
            trafficrows_float = []
            for i in range(len(trafficrows)):
                trafficrows_float.append(list(map(float, trafficrows[i][:24])))
            _init_param_dict[row[0]] = trafficrows_float
            count += 4
        elif row[0] == "bld":
            # bld: 17 x 3 matrix
            bldrows = VCWG_param[count + 1:count + 17]
            bldrows_float = []
            for i in range(len(bldrows)):
                bldrows_float.append(list(map(float, bldrows[i][:3])))
            _init_param_dict[row[0]] = bldrows_float
            count += 17
        elif row[0] == "LAD":
            # LAD profile
            LADrows = VCWG_param[count + 1:count + 3]
            LADrows_float = []
            for i in range(len(LADrows)):
                LADrows_float.append(list(map(float, LADrows[i][:-1])))
            _init_param_dict[row[0]] = LADrows_float
            count += 3
        elif row[0] == "Zs_R":
            # Zs_R profile
            Zs_Rrows = VCWG_param[count + 1:count + 2]
            Zs_Rrows_float = []
            for i in range(len(Zs_Rrows)):
                Zs_Rrows_float.append(list(map(float, Zs_Rrows[i][:-1])))
            _init_param_dict[row[0]] = Zs_Rrows_float
            count += 2
        elif row[0] == "Zs_G":
            # Zs_G profile
            Zs_Grows = VCWG_param[count + 1:count + 2]
            Zs_Grows_float = []
            for i in range(len(Zs_Grows)):
                Zs_Grows_float.append(list(map(float, Zs_Grows[i][:-1])))
            _init_param_dict[row[0]] = Zs_Grows_float
            count += 2
        elif row[0] == "Zs_W":
            # Zs_G profile
            Zs_Wrows = VCWG_param[count + 1:count + 2]
            Zs_Wrows_float = []
            for i in range(len(Zs_Wrows)):
                Zs_Wrows_float.append(list(map(float, Zs_Wrows[i][:-1])))
            _init_param_dict[row[0]] = Zs_Wrows_float
            count += 2
        elif row[0] == "Zs_rural":
            # Zs_G profile
            Zs_ruralrows = VCWG_param[count + 1:count + 2]
            Zs_ruralrows_float = []
            for i in range(len(Zs_ruralrows)):
                Zs_ruralrows_float.append(list(map(float, Zs_ruralrows[i][:-1])))
            _init_param_dict[row[0]] = Zs_ruralrows_float
            count += 2
        elif row[0] == "lan_dry_imp_W_layers":
            # Zs_G profile
            lan_Wrows = VCWG_param[count + 1:count + 2]
            lan_Wrows_float = []
            for i in range(len(lan_Wrows)):
                lan_Wrows_float.append(list(map(float, lan_Wrows[i][:-1])))
            _init_param_dict[row[0]] = lan_Wrows_float
            count += 2
        elif row[0] == "cv_s_imp_W_layers":
            # Zs_G profile
            cv_Wrows = VCWG_param[count + 1:count + 2]
            cv_Wrows_float = []
            for i in range(len(cv_Wrows)):
                cv_Wrows_float.append(list(map(float, cv_Wrows[i][:-1])))
            _init_param_dict[row[0]] = cv_Wrows_float
            count += 2
        else:
            if row[1] == 'NaN':
                _init_param_dict[row[0]] = numpy.NaN
            elif row[1] == 'inf':
                _init_param_dict[row[0]] = numpy.inf
            else:
                _init_param_dict[row[0]] = float(row[1])
            count += 1

    ipd_vcwg = _init_param_dict

    return ipd_vcwg

def ForcingData(MeteoDataRaw, itt, varargin,VCWG_param_file_path,SimTime):
    """
    ------
    INPUT:
    MeteoDataRaw: Forced variables at all time
    itt: Current time step
    varargin: Soil water potential [MPa]
    VCWG_param_file_path: VCWG parameters file path
    -------
    OUTPUT:
    SunPosition: Sun angles
    MeteoData: Forcing variables
    Anthropogenic: Anthropogenic parameters for water and heat
    Location: Location summary
    ParCalculation: General calculation parameters
    """

    def SetSunVariables(Datam, DeltaGMT, Lon, Lat, t_bef, t_aft):
        # Determine the julian day of the current time
        days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        nowYR = int(Datam[0])
        nowMO = int(Datam[1])
        nowDA = int(Datam[2])
        nowHR = Datam[3] + Datam[4] / 60 + Datam[5] / 3600

        if nowMO == 1:
            jDay = nowDA
        elif nowMO == 2:
            jDay = days[0] + nowDA
        else:
            jDay = sum(days[0:(nowMO - 1)]) + nowDA
            if nowYR % 4 == 0:
                if nowYR % 400 == 0:
                    jDay = jDay + 1
                elif nowYR % 100 != 0:
                    jDay = jDay + 1

        # Compute solar declination [rad]
        delta_S = 23.45 * numpy.pi / 180 * math.cos(2 * numpy.pi / 365 * (172 - jDay))
        # Compute time difference between standard and local meridian
        if Lon < 0:
            Delta_TSL = -1 / 15. * (15 * abs(DeltaGMT) - abs(Lon))
        else:
            Delta_TSL = 1 / 15. * (15 * abs(DeltaGMT) - abs(Lon))

        t = numpy.arange(nowHR - t_bef, nowHR + t_aft, 0.0166666)
        tau_S = numpy.zeros(len(t))
        for i in range(0, len(t)):
            # Compute hour angle of the sun [rad]
            if (t[i] < (12 + Delta_TSL)):
                tau_S[i] = 15 * numpy.pi / 180 * (t[i] + 12 - Delta_TSL)
            else:
                tau_S[i] = 15 * numpy.pi / 180 * (t[i] - 12 - Delta_TSL)

        # Compute solar altitude [rad]
        Lat_rad = Lat * numpy.pi / 180
        sinh_S = [
            math.sin(Lat_rad) * math.sin(delta_S) + math.cos(Lat_rad) * math.cos(delta_S) * math.cos(tau_S[i])
            for i in range(0, len(tau_S))]
        h_S = [math.asin(sinh_S[i]) for i in range(0, len(tau_S))]
        h_S = numpy.mean(h_S)

        # Compute Sun's azimuth [rad]
        zeta_S = [math.atan(-math.sin(tau_S[i]) / (
                math.tan(delta_S) * math.cos(Lat_rad) - math.sin(Lat_rad) * math.cos(tau_S[i]))) for i in
                  range(0, len(tau_S))]

        for i in range(0, len(t)):
            if (tau_S[i] > 0 and tau_S[i] <= numpy.pi):
                if (zeta_S[i] > 0.):
                    zeta_S[i] = zeta_S[i] + numpy.pi
                else:
                    zeta_S[i] = zeta_S[i] + (2. * numpy.pi)
            elif tau_S[i] >= numpy.pi and tau_S[i] <= 2 * numpy.pi:
                if zeta_S[i] < 0.:
                    zeta_S[i] = zeta_S[i] + numpy.pi

        zeta_S = numpy.mean(zeta_S)

        # Compute sunrise time, sunset time, and total day length
        T_sunrise = 180 / (15 * numpy.pi) * (
                2 * numpy.pi - math.acos(-math.tan(delta_S) * math.tan(Lat_rad))) - 12
        T_sunset = 180 / (15 * numpy.pi) * math.acos(-math.tan(delta_S) * math.tan(Lat_rad)) + 12
        L_day = 360 / (15 * numpy.pi) * math.acos(-math.tan(delta_S) * math.tan(Lat_rad))
        T_sunrise = numpy.real(T_sunrise)
        T_sunset = numpy.real(T_sunset)
        L_day = numpy.real(L_day)

        return h_S, delta_S, zeta_S, T_sunrise, T_sunset, L_day, jDay

    def is_near_zero(self,num,eps=1e-10):
        return abs(float(num)) < eps

    ipd = read_VCWG_param(VCWG_param_file_path)

    # Input weather variables
    LWR_in = MeteoDataRaw.LWR_in             # [W m^-2]
    DirectRad = MeteoDataRaw.Dir_in          # [W m^-2]
    DiffusiveRad = MeteoDataRaw.Dif_in       # [W m^-2]
    T_atm = MeteoDataRaw.T_atm               # [K]
    windspeed_u = MeteoDataRaw.windspeed_u   # [m s^-1]
    pressure_atm = MeteoDataRaw.pressure_atm # [Pa]
    rain = MeteoDataRaw.rain                 # [mm h^-1]
    rel_humidity = MeteoDataRaw.rel_humidity # [%]
    Year = MeteoDataRaw.Year
    Month = MeteoDataRaw.Month
    Day = MeteoDataRaw.Day
    Hour = MeteoDataRaw.Hour
    Min = MeteoDataRaw.Min
    Sec = MeteoDataRaw.Sec
    uDir = MeteoDataRaw.uDir                  # wind direction [deg]

    # Location variables
    # latitude positive north [deg]
    phi = MeteoDataRaw.lat
    # longitude positive east [deg]
    Lambda = MeteoDataRaw.lon
    # canyon orientation [rad]
    theta_canyon = ipd['theta_canyon'] * numpy.pi / 180
    # difference with Greenwich Meridian Time [h]
    DeltaGMT = MeteoDataRaw.GMT

    class Location_Def():
        pass
    Location = Location_Def()
    Location.phi = phi
    Location.Lambda = Lambda
    Location.theta_canyon = theta_canyon
    Location.DeltaGMT = DeltaGMT

    Datam = [Year[itt], Month[itt], Day[itt], Hour[itt], Min[itt], Sec[itt]]

    t_bef = 0.5
    t_aft = 0.5

    h_S, _a_, zeta_S, _b_, _c_, _d_, _e_ = SetSunVariables(Datam, DeltaGMT, Lambda, phi, t_bef, t_aft)

    # Solar zenith angle [rad]
    theta_Z = numpy.pi / 2 - h_S
    if theta_Z <= -numpy.pi / 2 or theta_Z >= numpy.pi / 2:
        theta_Z = numpy.pi / 2

    # difference between solar azimuth angle and canyon orientation [rad]
    theta_n = zeta_S - theta_canyon

    class SunPosition_Def():
        pass
    SunPosition = SunPosition_Def()
    SunPosition.Datam = Datam
    SunPosition.t_bef = t_bef
    SunPosition.t_aft = t_aft
    SunPosition.theta_Z = theta_Z
    SunPosition.theta_n = theta_n
    SunPosition.zeta_S = zeta_S

    # Radiation at the time of interest
    # Direct incoming shortwave radiation [W m^-2]
    SW_dir = DirectRad[itt]
    # Diffuse incoming shortwave radiation [W m^-2]
    SW_diff = DiffusiveRad[itt]
    # Atmospheric longwave radiation [W m^-2]
    LWR = LWR_in[itt]

    if abs(math.cos(theta_Z)) < 0.1:
        SW_dir = 0

    # Meteorological data
    # Atmospheric reference height [m]
    Zatm = ipd['dz']*ipd['nz']
    # Air Temperature at atmospheric reference level [K]
    Tatm = T_atm[itt]
    # Wind speed at atmospheric reference level [m s^-1]
    Uatm = windspeed_u[itt]
    if Uatm == 0:
        Uatm = ipd['WindMin_Urban']
    # Saturation vapor pressure at Tatm [Pa]
    esat_Tatm = 611 * numpy.exp(17.27 * (Tatm - 273.16) / (237.3 + (Tatm - 273.16)))
    # Relative humidity
    rel_hum = rel_humidity[itt]
    # vapor pressure [Pa]
    ea = esat_Tatm * rel_hum
    # air pressure [Pa]
    Pre = pressure_atm[itt]
    # Specifc humidity of air at reference height [kg kg^-1]
    q_atm = round(0.622 * ea / (Pre - 0.378 * ea),4)
    # Atmospheric CO2 mixing ratio 2017 [ppm]-[umolCO2 mol^-1]
    Catm_CO2 = ipd['Catm_CO2']
    # Intercellular Partial Pressure Oxygen [ppm] - [umolO2 mol^-1]
    Catm_O2 = ipd['Catm_O2']
    # Precipiation [mm s^-1]
    Rain = rain[itt]/3600
    Time = Datam

    class MeteoData_Def():
        pass
    MeteoData = MeteoData_Def()
    MeteoData.SW_dir = SW_dir
    MeteoData.SW_diff = SW_diff
    MeteoData.LWR = LWR
    MeteoData.Zatm = Zatm
    MeteoData.Tatm = Tatm
    MeteoData.Uatm = Uatm
    MeteoData.uDir = uDir[itt]
    MeteoData.esat_Tatm = esat_Tatm
    MeteoData.rel_hum = rel_hum
    MeteoData.ea = ea
    MeteoData.Pre = Pre
    MeteoData.q_atm = q_atm
    MeteoData.Catm_CO2 = Catm_CO2
    MeteoData.Catm_O2 = Catm_O2
    MeteoData.Rain = Rain
    MeteoData.Time = Time
    MeteoData.TdeepSoil = MeteoDataRaw.Tdeepsoil[itt]
    MeteoData.waterTemp = MeteoData.TdeepSoil

    # ANTHROPOGENIC FACTORS
    # Anthropogenic heat input into the canyon air [W m^-2]
    if is_near_zero(SimTime.julian%7,1e-10):
        dayType = 3 # Sunday
    elif is_near_zero(SimTime.julian%7-6,1e-10):
        dayType = 2 # Saturday
    else:
        dayType = 1 # Weekday

    Qf_canyon = ipd['SchTraffic'][dayType-1][SimTime.hourDay]

    # Anthropogenic water
    # This irrigates so no water stress of plants occur.
    if itt == 0:
        # applied on the vegetated ground surface area [mm h^-1]
        Waterf_canyonVeg = 0
    else:
        if varargin.SoilPotWGroundTot_H <= -0.1 or varargin.SoilPotWGroundVeg_L <= -0.3:
            # applied on the vegetated ground surface area [mm h^-1]
            Waterf_canyonVeg = 0
        else:
            # applied on the vegetated ground surface area [mm h^-1]
            Waterf_canyonVeg = 0

    # applied on the vegetated ground surface area [mm h^-1]
    Waterf_canyonBare = 0
    Waterf_roof = 0

    class Anthropogenic_Def():
        pass
    Anthropogenic = Anthropogenic_Def()
    Anthropogenic.Qf_canyon = Qf_canyon
    Anthropogenic.Waterf_canyonVeg = Waterf_canyonVeg
    Anthropogenic.Waterf_canyonBare = Waterf_canyonBare
    Anthropogenic.Waterf_roof = Waterf_roof

    # GENERAL PARAMETERS FOR CALCULATION
    # Time steps other than dth=1 and tds=3600 are not extensively tested in this version of the code.
    # time step of calculation [h]
    dth = ipd['dtSim']/3600
    # time step of calculation [s]
    dts = ipd['dtSim']
    # density of water [kg m^-3]
    row = 1000
    # specific heat air  [J kg^-1 K^-1]
    cp_atm = 1005 + (((Tatm - 273.15) + 23.15) ** 2) / 3364
    # dry air density at atmosphere [kg m^-3]
    rho_atm = (Pre / (287.04 * Tatm)) * (1 - (ea / Pre) * (1 - 0.622))
    # Latent heat vaporization/condensation [J kg^-1]
    L_heat = 1000 * (2501.3 - 2.361 * (Tatm - 273.15))

    class ParCalculation_Def():
        pass
    ParCalculation = ParCalculation_Def()
    ParCalculation.dts = dts
    ParCalculation.dth = dth
    ParCalculation.rhow = row
    ParCalculation.cp_atm = cp_atm
    ParCalculation.rho_atm = rho_atm
    ParCalculation.Lv = L_heat
    ParCalculation.nightStart = ipd['nightStart']
    ParCalculation.nightEnd = ipd['nightEnd']
    ParCalculation.Elevation = ipd['Elevation']

    return SunPosition, MeteoData, Anthropogenic, Location, ParCalculation

def Data_Site(InputFile):

    SoilCal = Soil_Calculations()

    ipd = read_VCWG_param(InputFile)

    # Rural model parameters
    class RSMParam_Def():
        pass
    RSMParam = RSMParam_Def()
    RSMParam.h_obs = ipd['h_obs']
    RSMParam.lv = 2.26e6
    RSMParam.cp = 1004.
    RSMParam.vk = 0.4
    RSMParam.WindMin_MOST = ipd['WindMin_MOST']
    RSMParam.L_Pos_min = ipd['L_Pos_min']
    RSMParam.L_Pos_max = ipd['L_Pos_max']
    RSMParam.L_Neg_max = ipd['L_Neg_max']
    RSMParam.L_Neg_min = ipd['L_Neg_min']
    RSMParam.ZL_Pos_cutoff = ipd['ZL_Pos_cutoff']
    RSMParam.ZL_Neg_cutoff = ipd['ZL_Neg_cutoff']
    RSMParam.u_star_min_MOST = ipd['u_star_min_MOST']
    RSMParam.z0overh_MOST = ipd['z0overh_MOST']
    RSMParam.zToverz0_MOST = ipd['zToverz0_MOST']
    RSMParam.dispoverh_MOST = ipd['dispoverh_MOST']
    RSMParam.h_wind = ipd['h_wind']
    RSMParam.h_temp = ipd['h_temp']
    RSMParam.Bowen = ipd['BowenRatio_rural']
    RSMParam.e_rural = ipd['e_rural']
    RSMParam.a_rural = ipd['a_rural']
    RSMParam.a_veg = ipd['aveg_rural']
    RSMParam.vegStart = ipd['vegStart']
    RSMParam.vegEnd = ipd['vegEnd']
    RSMParam.rurVegCover = ipd['rurVegCover']
    RSMParam.aveg_rural = ipd['aveg_rural']
    RSMParam.cdmin = ipd['cdmin_rural']
    RSMParam.prandtl = ipd['prandtl_rural']
    RSMParam.schmidt = ipd['schmidt_rural']
    RSMParam.Zs = ipd['Zs_rural'][0]
    RSMParam.lan_rural = ipd['lan_rural']
    RSMParam.cv_s_rural = ipd['cv_s_rural']
    RSMParam.r = 287
    RSMParam.z0_Louis = ipd['z0_Louis']
    RSMParam.MinWind_rural = ipd['MinWind_rural']
    if ipd['Rural_Model_name'] == 1:
        RSMParam.Rural_Model_name = 'MOST'
    elif ipd['Rural_Model_name'] == 2:
        RSMParam.Rural_Model_name = 'Forcing_extFile'
    if ipd['EB_RuralModel_name'] == 1:
        RSMParam.EB_RuralModel_name = 'Louis'
    elif ipd['EB_RuralModel_name'] == 2:
        RSMParam.EB_RuralModel_name = 'Penman_Monteith'

    # Column model parameters
    class ColParam_Def():
        pass
    ColParam = ColParam_Def()
    ColParam.prandtl = ipd['prandtl']
    ColParam.schmidt = ipd['schmidt']
    ColParam.h_LAD = ipd['LAD'][0]
    ColParam.LAD = ipd['LAD'][1]
    ColParam.HVAC_atm_frac = ipd['HVAC_atm_frac']
    ColParam.HVAC_street_frac = ipd['HVAC_street_frac']
    ColParam.Ri_b_cr = ipd['Ri_b_cr']
    ColParam.WindMin_Urban = ipd['WindMin_Urban']
    ColParam.cdmin = ipd['cdmin']
    ColParam.omega = ipd['omega']
    ColParam.omega_drag = ipd['omega_drag']

    # GEOMETRY OF URBAN AREA
    # Canyon (building) height [m]
    Height_canyon = ipd['Height_canyon']
    # Canyon (road) width [m]
    Width_canyon = ipd['Width_canyon']
    # Roof width [m], calculated from the land cover fraction and the street width.
    Width_roof = ipd['Width_roof']
    # Tree-crown radius [m], Calculated out of rescaled tree fraction and street width (assuming two uniform strips of tree rows)
    Radius_tree = ipd['Radius_tree']
    # Tree height [m] from ground to the middle of the crown
    Height_tree = ipd['LAD'][0][-1]-Radius_tree-0.1
    # Tree-to-wall distance [m]
    Distance_tree = ipd['distance_tree']

    # asy switch to include (=1) and exclude (=0) trees in the urban canyon
    trees = ipd['trees']
    # DO NOT CHANGE: Tree fraction along canyon axis
    ftree = ipd['ftree']

    if numpy.isnan(Radius_tree):
        Radius_tree = 0

    # normalized canyon height[-]
    hcanyon = Height_canyon / Width_canyon
    # normalized canyon width [-]
    wcanyon = Width_canyon / Width_canyon
    # normalized roof width [-]
    wroof = Width_roof / Width_canyon
    # normalized tree height [-]
    htree = Height_tree / Width_canyon
    # normalized tree radius [-]
    radius_tree = Radius_tree / Width_canyon
    # normalized tree-to-wall distance [-]
    distance_tree = Distance_tree / Width_canyon
    # height to width ratio [-]
    ratio = hcanyon / wcanyon

    # normalized canyon width overall [-]
    wcanyon_norm = wcanyon / (wcanyon + wroof)
    # normalized roof width overall [-]
    wroof_norm = wroof / (wcanyon + wroof)

    class Gemeotry_m_Def():
        pass
    Gemeotry_m = Gemeotry_m_Def()
    Gemeotry_m.Height_canyon = Height_canyon
    Gemeotry_m.Width_canyon = Width_canyon
    Gemeotry_m.Width_roof = Width_roof
    Gemeotry_m.Height_tree = Height_tree
    Gemeotry_m.Radius_tree = Radius_tree
    Gemeotry_m.Distance_tree = Distance_tree
    Gemeotry_m.nz = int(ipd['nz'])
    Gemeotry_m.nz_u = int(ipd['nz_u'])
    Gemeotry_m.dz = ipd['dz']
    Gemeotry_m.theta_canyon = ipd['theta_canyon']
    Gemeotry_m.z = numpy.linspace(Gemeotry_m.dz/2, Gemeotry_m.nz*Gemeotry_m.dz+Gemeotry_m.dz/2, Gemeotry_m.nz+1)
    Gemeotry_m.lambdap = Width_roof/(Width_roof+Width_canyon)
    Gemeotry_m.lambdaf = Height_canyon/(Width_roof+Width_canyon)
    Gemeotry_m.z_vdm = numpy.arange(RSMParam.h_temp, Gemeotry_m.nz*Gemeotry_m.dz+Gemeotry_m.dz/2, Gemeotry_m.dz)
    Gemeotry_m.nz_vdm = len(Gemeotry_m.z_vdm)-1

    class ParTree_Def():
        pass
    ParTree = ParTree_Def()
    ParTree.trees = trees
    ParTree.ftree = ftree

    class geometry_Def():
        pass
    geometry = geometry_Def()
    geometry.hcanyon = hcanyon
    geometry.wcanyon = wcanyon
    geometry.wroof = wroof
    geometry.htree = htree
    geometry.radius_tree = radius_tree
    geometry.distance_tree = distance_tree
    geometry.ratio = ratio
    geometry.wcanyon_norm = wcanyon_norm
    geometry.wroof_norm = wroof_norm

    # SURFACE FRACTIONS
    # Roof
    # Vegetated roof fraction
    fveg_R = ipd['fveg_R']
    # Impervious roof fraction
    fimp_R = ipd['fimp_R']
    # Percentage of excess water that leaves the system as runoff, needs to be between 0-1 [-]
    Per_runoff_R = ipd['Per_runoff_R']
    class FractionsRoof_Def():
        pass
    FractionsRoof = FractionsRoof_Def()
    FractionsRoof.fveg = fveg_R
    FractionsRoof.fimp = fimp_R
    FractionsRoof.Per_runoff = Per_runoff_R

    # Ground
    # Vegetated ground fraction
    fveg_G = ipd['fveg_G']
    # Bare ground fraction
    fbare_G = ipd['fbare_G']
    # Impervious ground fraction
    fimp_G = ipd['fimp_G']
    # Percentage of excess water that leaves the system as runoff, needs to be between 0-1 [-]
    Per_runoff_G = ipd['Per_runoff_G']
    class FractionsGround_Def():
        pass
    FractionsGround = FractionsGround_Def()
    FractionsGround.fveg = fveg_G
    FractionsGround.fbare = fbare_G
    FractionsGround.fimp = fimp_G
    FractionsGround.Per_runoff = Per_runoff_G

    # OPTICAL PROPERTIES
    # Roof
    # Roof vegetation surface albedo [-]
    aveg_R = ipd['aveg_R']
    # Roof impervious albedo [-]
    aimp_R = ipd['aimp_R']
    # equivalent roof surface albedo [-]
    albedo_R = fveg_R * aveg_R + fimp_R * aimp_R
    # Roof vegetation surface emissivity [-]
    eveg_R = ipd['eveg_R']
    # Roof impervious emissivity [-]
    eimp_R = ipd['eimp_R']
    # equivalent roof surface emissivity [-]
    emissivity_R = fveg_R * eveg_R + fimp_R * eimp_R
    class PropOpticalRoof_Def():
        pass
    PropOpticalRoof = PropOpticalRoof_Def()
    PropOpticalRoof.aveg = aveg_R
    PropOpticalRoof.aimp = aimp_R
    PropOpticalRoof.albedo = albedo_R
    PropOpticalRoof.eveg = eveg_R
    PropOpticalRoof.eimp = eimp_R
    PropOpticalRoof.emissivity = emissivity_R

    # Ground
    # Ground vegetation surface albedo [-]
    aveg_G = ipd['aveg_G']
    # Ground vegetation surface albedo [-]
    abare_G = ipd['abare_G']
    # Ground impervious albedo [-]
    aimp_G = ipd['aimp_G']
    # equivalent Ground surface albedo [-]
    albedo_G = fveg_G*aveg_G + fbare_G*abare_G + fimp_G*aimp_G
    # Ground vegetation surface emissivity [-]
    eveg_G = ipd['eveg_G']
    # Ground vegetation surface emissivity [-]
    ebare_G = ipd['ebare_G']
    # Ground impervious emissivity [-]
    eimp_G = ipd['eimp_G']
    # equivalent Ground surface emissivity [-]
    emissivity_G = fveg_G*eveg_G + fbare_G*ebare_G + fimp_G*eimp_G
    class PropOpticalGround_Def():
        pass
    PropOpticalGround = PropOpticalGround_Def()
    PropOpticalGround.aveg = aveg_G
    PropOpticalGround.abare = abare_G
    PropOpticalGround.aimp = aimp_G
    PropOpticalGround.albedo = albedo_G
    PropOpticalGround.eveg = eveg_G
    PropOpticalGround.ebare = ebare_G
    PropOpticalGround.eimp = eimp_G
    PropOpticalGround.emissivity = emissivity_G

    # Wall
    # Wall surface albedo [-]
    albedo_W = ipd['albedo_W']
    # Wall emissivity [-]
    emissivity_W = ipd['emissivity_W']
    class PropOpticalWall_Def():
        pass
    PropOpticalWall = PropOpticalWall_Def()
    PropOpticalWall.albedo = albedo_W
    PropOpticalWall.emissivity = emissivity_W

    # Tree
    # Tree albedo [-]
    albedo_T = ipd['albedo_T']
    # Tree emissivity [-]
    emissivity_T = ipd['emissivity_T']
    class PropOpticalTree_Def():
        pass
    PropOpticalTree = PropOpticalTree_Def()
    PropOpticalTree.albedo = albedo_T
    PropOpticalTree.emissivity = emissivity_T

    # THERMAL PROPERTIES
    # Roof
    # Thermal conductivity dry solid [W m^-1 K^-1]
    lan_dry_imp_R = ipd['lan_dry_imp_R']
    # Volumetric heat capacity solid [J m^-3 K^-1]
    cv_s_imp_R = ipd['cv_s_imp_R']
    class ParThermalRoof_Def():
        pass
    ParThermalRoof = ParThermalRoof_Def()
    ParThermalRoof.lan_dry_imp = lan_dry_imp_R
    ParThermalRoof.cv_s_imp = cv_s_imp_R

    # Ground
    # Thermal conductivity dry solid [W m^-1 K^-1]
    lan_dry_imp_G = ipd['lan_dry_imp_G']
    # Volumetric heat capacity solid [J m^-3 K^-1]
    cv_s_imp_G = ipd['cv_s_imp_G']
    class ParThermalGround_Def():
        pass
    ParThermalGround = ParThermalGround_Def()
    ParThermalGround.lan_dry_imp = lan_dry_imp_G
    ParThermalGround.cv_s_imp = cv_s_imp_G
    ParThermalGround.lan_dry_veg = ipd['lan_dry_veg_G']
    ParThermalGround.cv_s_veg = ipd['cv_s_veg_G']
    ParThermalGround.lan_dry_bare = ipd['lan_dry_bare_G']
    ParThermalGround.cv_s_bare = ipd['cv_s_bare_G']
    if ipd['Tdeep_ctrl'] == 1:
        ParThermalGround.Tdeep_ctrl = 'Force_Restore'
    elif ipd['Tdeep_ctrl'] == 2:
        ParThermalGround.Tdeep_ctrl = 'Climate_Data'

    # Wall
    # Thermal conductivity dry solid [W m^-1 K^-1]
    lan_dry_imp_W = ipd['lan_dry_imp_W']
    lan_dry_imp_W_layers = ipd['lan_dry_imp_W_layers'][0]
    # Volumetric heat capacity solid [J m^-3 K^-1]
    cv_s_imp_W = ipd['cv_s_imp_W']
    cv_s_imp_W_layers = ipd['cv_s_imp_W_layers'][0]
    class ParThermalWall_Def():
        pass
    ParThermalWall = ParThermalWall_Def()
    ParThermalWall.lan_dry = lan_dry_imp_W
    ParThermalWall.cv_s = cv_s_imp_W
    ParThermalWall.lan_dry_layers = lan_dry_imp_W_layers
    ParThermalWall.cv_s_imp_W_layers = cv_s_imp_W_layers

    # VEGETATION
    # ROOF VEGETATION
    # General
    # Leaf area index for the roof vegetation [-]
    LAI_R = ipd['LAI_R']
    # Stem area index for the roof vegetation [-]
    SAI_R = ipd['SAI_R']
    # canopy height roof vegetation	[m]
    hc_R = ipd['hc_R']
    # Zero plane displacement height of roof vegetation [m]
    h_disp_R = 2 / 3 * hc_R
    # Leaf dimension of roof vegetation [cm]
    d_leaf_R = ipd['d_leaf_R']
    # Roof water uptake
    # Type of Root Profile
    CASE_ROOT_R = ipd['CASE_ROOT_R']
    # Root depth 95 percentile [mm]
    ZR95_R = ipd['ZR95_R']
    # Root depth 50 percentile [mm]
    ZR50_R = ipd['ZR50_R']
    # Maximum Root depth [mm]
    ZRmax_R = ipd['ZRmax_R']
    # Root length index [m root m^-2 PFT]
    Rrootl_R = ipd['Rrootl_R']
    # Water Potential at 50% loss conductivity [MPa]
    PsiL50_R = ipd['PsiL50_R']
    # Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil [MPa]
    PsiX50_R = ipd['PsiX50_R']
    # Photosynthesis and Transpiration
    # Intrinsic quantum Efficiency [umolCO2 umolPhotons^-1]
    FI_R = ipd['FI_R']
    # Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis [Pa]
    Do_R = ipd['Do_R']
    # Empirical parameter connecting stomatal aperture and net assimilation []
    a1_R = ipd['a1_R']
    # minimal Stomatal Conductance [mol s^-1 m^-2]
    go_R = ipd['go_R']
    # --> 'CT' == 3  'CT' ==  4 Photosyntesis Typology for Plants, Photosynthetic pathway C3 or C4
    CT_R = ipd['CT_R']
    # Activation Energy - Plant Dependent, Activation Energy in Photosynthesis for Rubisco Capacity [kJ mol^-1]
    DSE_R = ipd['DSE_R']
    # entropy factor - Plant Dependent, Activation energy. [kJ mol^-1 K^-1]
    Ha_R = ipd['Ha_R']
    # Mesophyll conductance, not used [mol CO2 s^-1 m^-2]
    gmes_R = ipd['gmes_R']
    # Scaling factor between Jmax and Vmax
    rjv_R = ipd['rjv_R']
    # Light extinction parameter [-]
    Kopt_R = ipd['Kopt_R']
    # Canopy nitrogen decay coefficient [-]
    Knit_R = ipd['Knit_R']
    # Maximum Rubisco capacity at 25°C leaf level
    Vmax_R = ipd['Vmax_R']
    #
    mSl_R = ipd['mSl_R']
    # Relative Efficiency of the photosynthesis apparatus due to Age/Day-length []
    e_rel_R = ipd['e_rel_R']
    # Relative efficiency of the photosynthesis apparatus due to N limitations []
    e_relN_R = ipd['e_relN_R']
    # Water Potential at PLCs loss conductivity [MPa]
    Psi_sto_00_R = ipd['Psi_sto_00_R']
    # Water Potential at 50% loss conductivity [MPa]
    Psi_sto_50_R = ipd['Psi_sto_50_R']
    # specific leaf area of  biomass [m^2 gC^-1]
    Sl_R = ipd['Sl_R']
    class ParVegRoof_Def():
        pass
    ParVegRoof = ParVegRoof_Def()
    ParVegRoof.LAI = LAI_R
    ParVegRoof.SAI = SAI_R
    ParVegRoof.hc = hc_R
    ParVegRoof.h_disp = h_disp_R
    ParVegRoof.d_leaf = d_leaf_R
    ParVegRoof.CASE_ROOT = CASE_ROOT_R
    ParVegRoof.ZR95 = ZR95_R
    ParVegRoof.ZR50 = ZR50_R
    ParVegRoof.ZRmax = ZRmax_R
    ParVegRoof.Rrootl = Rrootl_R
    ParVegRoof.PsiL50 = PsiL50_R
    ParVegRoof.PsiX50 = PsiX50_R
    ParVegRoof.FI = FI_R
    ParVegRoof.Do = Do_R
    ParVegRoof.a1 = a1_R
    ParVegRoof.go = go_R
    ParVegRoof.CT = CT_R
    ParVegRoof.DSE = DSE_R
    ParVegRoof.Ha = Ha_R
    ParVegRoof.gmes = gmes_R
    ParVegRoof.rjv = rjv_R
    ParVegRoof.Kopt = Kopt_R
    ParVegRoof.Knit = Knit_R
    ParVegRoof.Vmax = Vmax_R
    ParVegRoof.mSl = mSl_R
    ParVegRoof.e_rel = e_rel_R
    ParVegRoof.e_relN = e_relN_R
    ParVegRoof.Psi_sto_00 = Psi_sto_00_R
    ParVegRoof.Psi_sto_50 = Psi_sto_50_R
    ParVegRoof.Sl = Sl_R

    # GROUND VEGETATION
    # Grass
    # General
    # Leaf area index for the ground vegetation [-]
    LAI_G = ipd['LAI_G']
    # Stem area index for the ground vegetation [-]
    SAI_G = ipd['SAI_G']
    # canopy height ground vegetation [m]
    hc_G = ipd['hc_G']
    # Zero plane displacement height of ground vegetation [m]
    h_disp_G = 2 / 3 * hc_G
    # Leaf dimension of ground vegetation [cm]
    d_leaf_G = ipd['d_leaf_G']
    # ground water uptake
    # Type of Root Profile
    CASE_ROOT_G = ipd['CASE_ROOT_G']
    # Root depth 95 percentile [mm]
    ZR95_G = ipd['ZR95_G']
    # Root depth 50 percentile [mm]
    ZR50_G = ipd['ZR50_G']
    # Maximum Root depth [mm]
    ZRmax_G = ipd['ZRmax_G']
    # Root length index [m root m^-2 PFT]
    Rrootl_G = ipd['Rrootl_G']
    # Water Potential at 50% loss conductivity [MPa]
    PsiL50_G = ipd['PsiL50_G']
    # Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil [MPa]
    PsiX50_G = ipd['PsiX50_G']
    # Photosynthesis and Transpiration
    # Intrinsic quantum Efficiency [umolCO2 umolPhotons^-1]
    FI_G = ipd['FI_G']
    # Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis [Pa]
    Do_G = ipd['Do_G']
    # Empirical parameter connecting stomatal aperture and net assimilation []
    a1_G = ipd['a1_G']
    # minimal Stomatal Conductance [mol s^-1 m^-2]
    go_G = ipd['go_G']
    # --> 'CT' == 3  'CT' ==  4 Photosyntesis Typology for Plants, Photosynthetic pathway C3 or C4
    CT_G = ipd['CT_G']
    # Activation Energy - Plant Dependent, Activation Energy in Photosynthesis for Rubisco Capacity [kJ mol^-1]
    DSE_G = ipd['DSE_G']
    # entropy factor - Plant Dependent, Activation energy. [kJ mol^-1 K^-1]
    Ha_G = ipd['Ha_G']
    # Mesophyll conductance, not used [mol CO2 s^-1 m^-2]
    gmes_G = ipd['gmes_G']
    # Scaling factor between Jmax and Vmax
    rjv_G = ipd['rjv_G']
    # light extinction parameter [-]
    Kopt_G = ipd['Kopt_G']
    # Canopy nitrogen decay coefficient [-]
    Knit_G = ipd['Knit_G']
    # Maximum Rubisco capacity at 25°C leaf level
    Vmax_G = ipd['Vmax_G']
    #
    mSl_G = ipd['mSl_G']
    # Relative Efficiency of the photosynthesis apparatus due to Age/Day-length []
    e_rel_G = ipd['e_rel_G']
    # Relative efficiency of the photosynthesis apparatus due to N limitations []
    e_relN_G = ipd['e_relN_G']
    # Water Potential at PLCs loss conductivity [MPa]
    Psi_sto_00_G = ipd['Psi_sto_00_G']
    # Water Potential at 50% loss conductivity [MPa]
    Psi_sto_50_G = ipd['Psi_sto_50_G']
    # specific leaf area of  biomass [m^2 gC^-1]
    Sl_G = ipd['Sl_G']
    class ParVegGround_Def():
        pass
    ParVegGround = ParVegGround_Def()
    ParVegGround.LAI = LAI_G
    ParVegGround.SAI = SAI_G
    ParVegGround.hc = hc_G
    ParVegGround.h_disp = h_disp_G
    ParVegGround.d_leaf = d_leaf_G
    ParVegGround.CASE_ROOT = CASE_ROOT_G
    ParVegGround.ZR95 = ZR95_G
    ParVegGround.ZR50 = ZR50_G
    ParVegGround.ZRmax = ZRmax_G
    ParVegGround.Rrootl = Rrootl_G
    ParVegGround.PsiL50 = PsiL50_G
    ParVegGround.PsiX50 = PsiX50_G
    ParVegGround.FI = FI_G
    ParVegGround.Do = Do_G
    ParVegGround.a1 = a1_G
    ParVegGround.go = go_G
    ParVegGround.CT = CT_G
    ParVegGround.DSE = DSE_G
    ParVegGround.Ha = Ha_G
    ParVegGround.gmes = gmes_G
    ParVegGround.rjv = rjv_G
    ParVegGround.Kopt = Kopt_G
    ParVegGround.Knit = Knit_G
    ParVegGround.Vmax = Vmax_G
    ParVegGround.mSl = mSl_G
    ParVegGround.e_rel = e_rel_G
    ParVegGround.e_relN = e_relN_G
    ParVegGround.Psi_sto_00 = Psi_sto_00_G
    ParVegGround.Psi_sto_50 = Psi_sto_50_G
    ParVegGround.Sl = Sl_G

    # TREE VEGETATION
    # General
    # Leaf area index for the Tree vegetation [-]
    LAI_T = ipd['LAI_T']
    # Stem area index for the Tree vegetation [-]
    SAI_T = ipd['SAI_T']
    # Leaf dimension of Tree vegetation [cm]
    d_leaf_T = ipd['d_leaf_T']
    # Tree water uptake
    # Type of Root Profile
    CASE_ROOT_T = ipd['CASE_ROOT_T']
    # Root depth 95 percentile [mm]
    ZR95_T = ipd['ZR95_T']
    # Root depth 50 percentile [mm]
    ZR50_T = ipd['ZR50_T']
    # Maximum Root depth [mm]
    ZRmax_T = ipd['ZRmax_T']
    # Root length index [m root m^-2 PFT]
    Rrootl_T = ipd['Rrootl_T']
    # Water Potential at 50% loss conductivity [MPa]
    PsiL50_T = ipd['PsiL50_T']
    # Water potential at 50 of xylem hydraulic conductivity and limit for water extraction from soil [MPa]
    PsiX50_T = ipd['PsiX50_T']
    # Photosynthesis and Transpiration
    # Intrinsic quantum Efficiency [umolCO2 umolPhotons^-1]
    FI_T = ipd['FI_T']
    # Empirical coefficient for the role of vapor pressure in the biochemical model of photosynthesis [Pa]
    Do_T = ipd['Do_T']
    # Empirical parameter connecting stomatal aperture and net assimilation []
    a1_T = ipd['a1_T']
    # minimal Stomatal Conductance [mol s^-1 m^-2]
    go_T = ipd['go_T']
    # --> 'CT' == 3  'CT' ==  4 Photosyntesis Typology for Plants, Photosynthetic pathway C3 or C4
    CT_T = ipd['CT_T']
    # Activation Energy - Plant Dependent, Activation Energy in Photosynthesis for Rubisco Capacity [kJ mol^-1]
    DSE_T = ipd['DSE_T']
    # entropy factor - Plant Dependent, Activation energy. [kJ mol^-1 K^-1]
    Ha_T = ipd['Ha_T']
    # Mesophyll conductance, not used [mol CO2 s^-1 m^-2]
    gmes_T = ipd['gmes_T']
    # Scaling factor between Jmax and Vmax
    rjv_T = ipd['rjv_T']
    # Light extinction parameter [-]
    Kopt_T = ipd['Kopt_T']
    # Canopy nitrogen decay coefficient [-]
    Knit_T = ipd['Knit_T']
    # Maximum Rubisco capacity at 25°C leaf level
    Vmax_T = ipd['Vmax_T']
    #
    mSl_T = ipd['mSl_T']
    # Relative Efficiency of the photosynthesis apparatus due to Age/Day-length []
    e_rel_T = ipd['e_rel_T']
    # Relative efficiency of the photosynthesis apparatus due to N limitations []
    e_relN_T = ipd['e_relN_T']
    # Water Potential at PLCs loss conductivity [MPa]
    Psi_sto_00_T = ipd['Psi_sto_00_T']
    # Water Potential at 50% loss conductivity [MPa]
    Psi_sto_50_T = ipd['Psi_sto_50_T']
    # specific leaf area of  biomass [m^2 gC^-1]
    Sl_T = ipd['Sl_T']
    # Tree root distribution: 1 = Tree roots can access all water in the soil (imp, bare, veg) equally 2 =  If the tree
    # crown is smaller than the combined vegetated and bare fraction, then the trees only transpire from these fractions.
    # Otherwise, they also transpire from the impervious ground fraction.
    SPARTREE = ipd['SPARTREE']
    class ParVegTree_Def():
        pass
    ParVegTree = ParVegTree_Def()
    ParVegTree.LAI = LAI_T
    ParVegTree.SAI = SAI_T
    ParVegTree.d_leaf = d_leaf_T
    ParVegTree.CASE_ROOT = CASE_ROOT_T
    ParVegTree.ZR95 = ZR95_T
    ParVegTree.ZR50 = ZR50_T
    ParVegTree.ZRmax = ZRmax_T
    ParVegTree.Rrootl = Rrootl_T
    ParVegTree.PsiL50 = PsiL50_T
    ParVegTree.PsiX50 = PsiX50_T
    ParVegTree.FI = FI_T
    ParVegTree.Do = Do_T
    ParVegTree.a1 = a1_T
    ParVegTree.go = go_T
    ParVegTree.CT = CT_T
    ParVegTree.DSE = DSE_T
    ParVegTree.Ha = Ha_T
    ParVegTree.gmes = gmes_T
    ParVegTree.rjv = rjv_T
    ParVegTree.Kopt = Kopt_T
    ParVegTree.Knit = Knit_T
    ParVegTree.Vmax = Vmax_T
    ParVegTree.mSl = mSl_T
    ParVegTree.e_rel = e_rel_T
    ParVegTree.e_relN = e_relN_T
    ParVegTree.Psi_sto_00 = Psi_sto_00_T
    ParVegTree.Psi_sto_50 = Psi_sto_50_T
    ParVegTree.Sl = Sl_T
    ParVegTree.SPARTREE = SPARTREE

    # SOLID LAYER DISCRETIZATION
    # Roof
    # Layer thickness
    # Soil layer discritization
    # soil layer discretization [mm]
    Zs_R = ipd['Zs_R'][0]
    # number of soil layers [-]
    ms_R = len(Zs_R) - 1

    # Ground
    # Soil layer discritization
    # soil layer discretization [mm]
    Zs_G = ipd['Zs_G'][0]
    # number of soil layers [-]
    ms_G = len(Zs_G) - 1

    # Wall
    # Wall layer discretization
    Zs_W = ipd['Zs_W'][0]
    class WallLayers_Def():
        pass
    WallLayers = WallLayers_Def()
    WallLayers.Zs_W = Zs_W

    # INTERCEPTION AND SOIL PARAMETERS
    # Roof
    # Soil classification: Sandy Loam
    # Fraction of clay in the soil [-]
    Pcla_R = ipd['Pcla_R']
    # Fraction of sand in the soil [-]
    Psan_R = ipd['Psan_R']
    # Fraction of organic material in the soil [-]
    Porg_R = ipd['Porg_R']
    # Interception and soil parameters
    # Maxiumum interception capacity of roof impervious area [mm]
    In_max_imp_R = ipd['In_max_imp_R']
    # Maximum interception capacity of ground under roof vegetation [mm]
    In_max_ground_R = ipd['In_max_ground_R']
    # specific water retained by a vegetated surface [mm m^2 VEG area m^-2 plant area]
    Sp_In_R = ipd['Sp_In_R']
    # Hydraulic conductivity of impervious area [mm h^-1]
    Kimp_R = ipd['Kimp_R']
    # Conductivity at field capacity [mm h^-1]
    Kfc_R = ipd['Kfc_R']
    # Suction at the residual/hygroscopic water content [kPa]
    Phy_R = ipd['Phy_R']
    # Conductivity at the bedrock layer [mm h^-1]
    Kbot_R = ipd['Kbot_R']
    # SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls. Do not use 1-VanGenuchten at the moment as very high soil water potential when dry
    SPAR_R = ipd['SPAR_R']
    SoilCal.Soil_Parameters(Psan_R, Pcla_R, Porg_R)
    class ParSoilRoof_Def():
        pass
    ParSoilRoof = ParSoilRoof_Def()
    ParSoilRoof.Zs = Zs_R
    ParSoilRoof.ms = ms_R
    ParSoilRoof.In_max_imp = In_max_imp_R
    ParSoilRoof.In_max_ground = In_max_ground_R
    ParSoilRoof.Sp_In = Sp_In_R
    ParSoilRoof.Kimp = Kimp_R
    ParSoilRoof.Kfc = Kfc_R
    ParSoilRoof.Phy = Phy_R
    ParSoilRoof.SPAR = SPAR_R
    ParSoilRoof.Kbot = Kbot_R
    ParSoilRoof.Pcla = Pcla_R
    ParSoilRoof.Psan = Psan_R
    ParSoilRoof.Porg = Porg_R
    ParSoilRoof.dz = numpy.diff(ParSoilRoof.Zs)
    ParSoilRoof.O33 = SoilCal.SoilParam.O33

    # Ground
    # Soil classification: Sandy Loam
    # Fraction of clay in the soil [-]
    Pcla_G = ipd['Pcla_G']
    # Fraction of sand in the soil [-]
    Psan_G = ipd['Psan_G']
    # Fraction of organic material in the soil [-]
    Porg_G = ipd['Porg_G']
    # Interception and soil parameters
    # Maxiumum interception capacity of impervious ground area [mm]
    In_max_imp_G = ipd['In_max_imp_G']
    # Maxiumum interception capacity of vegetated ground area [mm]
    In_max_underveg_G = ipd['In_max_underveg_G']
    # Maxiumum interception capacity of bare ground area [mm]
    In_max_bare_G = ipd['In_max_bare_G']
    # specific water retained by a vegetated surface [mm m^2 VEG area m^-2 plant area]
    Sp_In_G = ipd['Sp_In_G']
    # Hydraulic conductivity of impervious area [mm h^-1]
    Kimp_G = ipd['Kimp_G']
    # Conductivity at field capacity [mm h^-1]
    Kfc_G = ipd['Kfc_G']
    # Suction at the residual/hygroscopic water content [kPa]
    Phy_G = ipd['Phy_G']
    # Conductivity at the bedrock layer [mm h^-1]
    Kbot_G = ipd['Kbot_G']
    # SOIL PARAMETER TYPE 1-VanGenuchten 2-Saxton-Rawls. Do not use 1-VanGenuchten at the moment as very high soil water potential when dry
    SPAR_G = ipd['SPAR_G']
    SoilCal.Soil_Parameters(Psan_G, Pcla_G, Porg_G)
    class ParSoilGround_Def():
        pass
    ParSoilGround = ParSoilGround_Def()
    ParSoilGround.Zs = Zs_G
    ParSoilGround.ms = ms_G
    ParSoilGround.In_max_imp = In_max_imp_G
    ParSoilGround.In_max_underveg = In_max_underveg_G
    ParSoilGround.In_max_bare = In_max_bare_G
    ParSoilGround.Sp_In = Sp_In_G
    ParSoilGround.Kimp = Kimp_G
    ParSoilGround.Kfc = Kfc_G
    ParSoilGround.Phy = Phy_G
    ParSoilGround.SPAR = SPAR_G
    ParSoilGround.Kbot = Kbot_G
    ParSoilGround.Pcla = Pcla_G
    ParSoilGround.Psan = Psan_G
    ParSoilGround.Porg = Porg_G
    ParSoilGround.dz = numpy.diff(ParSoilGround.Zs)
    ParSoilGround.O33 = SoilCal.SoilParam.O33

    # Tree
    # Interception tree
    # specific water retained by the tree [mm m^2 VEG area m^-2 plant area]
    Sp_In_T = ipd['Sp_In_T']
    class ParInterceptionTree_Def():
        pass
    ParInterceptionTree = ParInterceptionTree_Def()
    ParInterceptionTree.Sp_In = Sp_In_T

    # PERSON for MRT calculation
    # position within canyon [m]
    PositionPx = Gemeotry_m.Width_canyon / 2
    # height of centre of person, usually choose 1.1 [m]
    PositionPz = ipd['PositionPz']
    # PersonWidth and PersonHeight are not used at the moment
    # horizontal radius of ellipse describing person (=hip width / 2)
    PersonWidth = ipd['PersonWidth']
    # Vertical radius of ellipse describing person (= height / 2)
    PersonHeight = ipd['PersonHeight']
    # Automatic wind speed calculation at user speficied height
    # height for wind speed to calculate OTC [m]
    HeightWind = ipd['HeightWind']
    class Person_Def():
        pass
    Person = Person_Def()
    Person.PositionPx = PositionPx
    Person.PositionPz = PositionPz
    Person.PersonWidth = PersonWidth
    Person.PersonHeight = PersonHeight
    Person.HeightWind = HeightWind

    class TimeParam_Def():
        pass
    TimeParam = TimeParam_Def()
    TimeParam.dts = ipd['dtSim']
    TimeParam.dth = ipd['dtSim']/3600
    TimeParam.nDay = ipd['nDay']
    TimeParam.Day = ipd['Day']
    TimeParam.Month = ipd['Month']
    TimeParam.dtWeather = ipd['dtWeather']
    TimeParam.nDay_spinup = ipd['nDay_spinup']

    class ViewFactorCal_Param_Def():
        pass
    ViewFactorCal_Param = ViewFactorCal_Param_Def()
    ViewFactorCal_Param.OPTION_RAY = ipd['OPTION_RAY']
    ViewFactorCal_Param.MCSampleSize = ipd['MCSampleSize']
    ViewFactorCal_Param.NRays = ipd['NRays']

    bld = ipd['bld']
    zone = int(ipd['zone'])-1

    charLength = ipd['charLength']

    class BEMParam_Def():
        pass
    BEMParam = BEMParam_Def()
    BEMParam.glzR = ipd['glzR']
    BEMParam.autosize = ipd['autosize']
    BEMParam.sensOcc = ipd['sensOcc']
    BEMParam.LatFOcc = ipd['LatFOcc']
    BEMParam.RadFOcc = ipd['RadFOcc']
    BEMParam.RadFEquip = ipd['RadFEquip']
    BEMParam.RadFLight = ipd['RadFLight']
    BEMParam.hvac = ipd['hvac']
    BEMParam.h_floor = ipd['h_floor']

    RSMParam.fimp = 0
    RSMParam.fveg = 1
    RSMParam.fbare = 0
    RSMParam.LAI = ParVegGround.LAI
    RSMParam.SAI = ParVegGround.SAI
    RSMParam.Kopt = ParVegGround.Kopt
    RSMParam.Psi_sto_50 = ParVegGround.Psi_sto_50
    RSMParam.Psi_sto_00 = ParVegGround.Psi_sto_00
    RSMParam.CT = ParVegGround.CT
    RSMParam.Vmax = ParVegGround.Vmax
    RSMParam.DSE = ParVegGround.DSE
    RSMParam.Ha = ParVegGround.Ha
    RSMParam.FI = ParVegGround.FI
    RSMParam.Do = ParVegGround.Do
    RSMParam.a1 = ParVegGround.a1
    RSMParam.go = ParVegGround.go
    RSMParam.e_rel = ParVegGround.e_rel
    RSMParam.e_relN = ParVegGround.e_relN
    RSMParam.gmes = ParVegGround.gmes
    RSMParam.rjv = ParVegGround.rjv
    RSMParam.mSl = ParVegGround.mSl
    RSMParam.Sl = ParVegGround.Sl
    RSMParam.CASE_ROOT = ParVegGround.CASE_ROOT
    RSMParam.ZR95 = ParVegGround.ZR95
    RSMParam.ZR50 = ParVegGround.ZR50
    RSMParam.ZRmax = ParVegGround.ZRmax
    RSMParam.d_leaf = ParVegGround.d_leaf
    RSMParam.Knit = ParVegGround.Knit
    RSMParam.Pcla = ParSoilGround.Pcla
    RSMParam.Psan = ParSoilGround.Psan
    RSMParam.Porg = ParSoilGround.Porg
    RSMParam.Kfc = ParSoilGround.Kfc
    RSMParam.Phy = ParSoilGround.Phy
    RSMParam.SPAR = ParSoilGround.SPAR
    RSMParam.Kbot = ParSoilGround.Kbot
    RSMParam.Zs = ParSoilGround.Zs
    RSMParam.Sp_In = ParSoilGround.Sp_In
    RSMParam.In_max_imp = ParSoilGround.In_max_imp
    RSMParam.In_max_bare = ParSoilGround.In_max_bare
    RSMParam.In_max_underveg = ParSoilGround.In_max_underveg
    RSMParam.Kimp = ParSoilGround.Kimp
    RSMParam.trees = 0
    RSMParam.SPARTREE = ParVegTree.SPARTREE
    RSMParam.eimp = PropOpticalGround.eimp
    RSMParam.ebare = PropOpticalGround.ebare
    RSMParam.eveg = PropOpticalGround.eveg
    RSMParam.aimp = PropOpticalGround.aimp
    RSMParam.abare = PropOpticalGround.abare
    RSMParam.aveg = PropOpticalGround.aveg
    RSMParam.Rrootl = ParVegGround.Rrootl
    RSMParam.PsiL50 = ParVegGround.PsiL50
    RSMParam.PsiX50 = ParVegGround.PsiX50
    RSMParam.Per_runoff = FractionsGround.Per_runoff


    return Gemeotry_m,ParTree,geometry,FractionsRoof,FractionsGround,WallLayers,ParSoilRoof,ParSoilGround,ParInterceptionTree,\
           PropOpticalRoof,PropOpticalGround,PropOpticalWall,PropOpticalTree,ParThermalRoof,ParThermalGround,ParThermalWall,\
           ParVegRoof,ParVegGround,ParVegTree,Person,ColParam,RSMParam,TimeParam,ViewFactorCal_Param,bld,zone,charLength,BEMParam