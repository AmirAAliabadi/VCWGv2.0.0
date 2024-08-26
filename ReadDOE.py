import sys
import os
import _pickle as cPickle

from BuildingEnergy import Building
from Material import Material
# from Element import Element
from BEMDef import BEMDef
from schdef import SchDef
from SurfaceTemperature import Tsurf_Def
from Utilities import read_csv, str2fl
import Utilities
import pandas
import numpy

"""
Developed by Bruno Bueno
Building Technology, Massachusetts Institute of Technology (MIT), Cambridge, U.S.A.
Last update: 2012
"""


DIR_CURR = os.path.abspath(os.path.dirname(__file__))
DIR_DOE_PATH = os.path.join(DIR_CURR,"resources","DOERefBuildings")

# Define standards: 16 building types, 3 built eras, 17 climate zones

# DOE Building Types
BLDTYPE = [
    'FullServiceRestaurant',    # 1
    'Hospital',                 # 2
    'LargeHotel',               # 3
    'LargeOffice',              # 4
    'MedOffice',                # 5
    'MidRiseApartment',         # 6
    'OutPatient',               # 7
    'PrimarySchool',            # 8
    'QuickServiceRestaurant',   # 9
    'SecondarySchool',          # 10
    'SmallHotel',               # 11
    'SmallOffice',              # 12
    'StandAloneRetail',         # 13
    'StripMall',                # 14
    'SuperMarket',              # 15
    'WareHouse']                # 16

BUILTERA = [
    'Pre80',                    # 1
    'Pst80',                    # 2
    'New'                       # 3
    ]

ZONETYPE = [
    '1A (Miami)',               # 1
    '2A (Houston)',             # 2
    '2B (Phoenix)',             # 3
    '3A (Atlanta)',             # 4
    '3B-CA (Los Angeles)',      # 5
    '3B (Las Vegas)',           # 6
    '3C (San Francisco)',       # 7
    '4A (Baltimore)',           # 8
    '4B (Albuquerque)',         # 9
    '4C (Seattle)',             # 10
    '5A (Chicago)',             # 11
    '5B (Boulder)',             # 12
    '6A (Minneapolis)',         # 13
    '6B (Helena)',              # 14
    '7 (Duluth)',               # 15
    '8 (Fairbanks)',            # 16
    'Custom'                    # 17
    ]

def readDOE(serialize_output=True):

    """
    Read csv files of DOE buildings
    Sheet 1 = BuildingSummary
    Sheet 2 = ZoneSummary
    Sheet 3 = LocationSummary
    Sheet 4 = Schedules
    Note BLD8 & 10 = school


    Then make matrix of ref data as nested nested lists [16, 3, 17]:
    matrix refDOE = Building objs
    matrix Schedule = SchDef objs
    matrix refBEM (16,3,17) = BEMDef
    where:
        [16,3,17] is Type = 1-16, Era = 1-3, climate zone = 1-17
        i.e.
        Type: FullServiceRestaurant, Era: Pre80, Zone: 6A Minneapolis
    Nested tree:
    [TYPE_1:
        ERA_1:
            CLIMATE_ZONE_1
            ...
            CLIMATE_ZONE_17
        ERA_2:
            CLIMATE_ZONE_1
            ...
            CLIMATE_ZONE_17
        ...
        ERA_3:
            CLIMATE_ZONE_1
            ...
            CLIMATE_ZONE_17]

    """

    #Nested, nested lists of Building, SchDef, BEMDef objects
    refDOE = list(map(lambda j_: list(map (lambda k_: [None]*17,[None]*3)), [None]*16))     #refDOE(17,3,16) = Building;
    Schedule = list(map(lambda j_: list(map (lambda k_: [None]*17,[None]*3)), [None]*16))   #Schedule (17,3,16) = SchDef;
    refBEM = list(map(lambda j_: list(map (lambda k_: [None]*17,[None]*3)), [None]*16))     #refBEM (17,3,16) = BEMDef;

    #Purpose: Loop through every DOE reference csv and extract building data
    #Nested loop = 16 types, 3 era, 17 zones 

    for i in range(16):

        #i = 16 types of buildings
        #print "\tType: {} @i={}".format(BLDTYPE[i], i)

        # Read building summary (Sheet 1)
        file_doe_name_bld = os.path.join("{}".format(DIR_DOE_PATH), "BLD{}".format(i+1),"BLD{}_BuildingSummary.csv".format(i+1))
        df_BuildingSummary = pandas.read_csv(file_doe_name_bld,header=1,usecols=[0,3,4,5],index_col=[0])
        #listof(listof 3 era values)
        nFloor      = df_BuildingSummary.loc['nFloor'].values      # Number of Floors, this will be list of floats and str if "basement"
        glazing     = df_BuildingSummary.loc['glazing'].values     # Ratio Total
        hCeiling    = df_BuildingSummary.loc['hCeiling'].values    # [m] Ceiling height
        ver2hor     = df_BuildingSummary.loc['ver2hor'].values     # Wall to Skin Ratio
        AreaRoof    = df_BuildingSummary.loc['areaRoof'].values    # [m^2] Gross Dimensions - Total area

        # Read zone summary (Sheet 2)
        file_doe_name_zone = os.path.join("{}".format(DIR_DOE_PATH), "BLD{}".format(i+1),"BLD{}_ZoneSummary.csv".format(i+1))
        df_ZoneSummary = pandas.read_csv(file_doe_name_zone, header=1, usecols=[0,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20], index_col=[0])
        #listof(listof 3 eras)
        AreaFloor   = df_ZoneSummary.loc[:,'Area (m2)'].values                   # [m^2]
        Volume      = df_ZoneSummary.loc[:,'Volume (m3)'].values                 # [m^3]
        AreaWall    = df_ZoneSummary.loc[:,'Gross Wall Area  (m2)'].values       # [m^2]
        AreaWindow  = df_ZoneSummary.loc[:,'Window Glass Area (m2)'].values      # [m^2]
        Occupant    = df_ZoneSummary.loc[:,'People'].values                      # Number of People
        Light       = df_ZoneSummary.loc[:,'Lights (W/m2)'].values               # [W m^-2]
        Elec        = df_ZoneSummary.loc[:,'Elec Plug and Process (W/m2)'].values# [W m^-2] Electric Plug and Process
        Gas         = df_ZoneSummary.loc[:,'Gas Plug and Process (W/m2)'].values # Gas Plug and Process per unit floor area [W m^-2]
        SHW         = df_ZoneSummary.loc[:,'SWH (L/h)'].values                   # Peak Service Hot Water per unit floor [kg hr^-1 m^-2]
        Vent        = df_ZoneSummary.loc[:,'Ventilation (L/s/m2)'].values        # [L s^-1 m^-2] Ventilation rate per unit floor area
        Infil       = df_ZoneSummary.loc[:,'Infiltration (ACH)'].values          # Infiltration Air Change per Hour (ACH) [hr^-1]

        # Read location summary (Sheet 3)
        file_doe_name_location = os.path.join("{}".format(DIR_DOE_PATH), "BLD{}".format(i+1),"BLD{}_LocationSummary.csv".format(i+1))
        df_LocationSummary = pandas.read_csv(file_doe_name_location,header=1,usecols=[1,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20],index_col=[0],encoding = "ISO-8859-1")
        #(listof (listof 3 eras (listof 16 climate types)))
        TypeWall    = df_LocationSummary.loc['TypeWall'].values                     # Construction type
        RvalWall    = df_LocationSummary.loc['RvalWall'].values.astype(numpy.float) # [m^2 K W^-1] R-value
        TypeRoof    = df_LocationSummary.loc['TypeRoof'].values                     # Construction type
        RvalRoof    = df_LocationSummary.loc['RvalRoof'].values.astype(numpy.float) # [m^2 K W^-1] R-value
        Uwindow     = df_LocationSummary.loc['Uwindow'].values.astype(numpy.float)  # [W m^-2 K^-1] U-factor
        SHGC        = df_LocationSummary.loc['SHGC'].values.astype(numpy.float)     # [-] coefficient
        HVAC        = df_LocationSummary.loc['HVAC'].values.astype(numpy.float)     # [kW] Air Conditioning
        HEAT        = df_LocationSummary.loc['HEAT'].values.astype(numpy.float)     # [kW] Heating
        COP         = df_LocationSummary.loc['COP'].values.astype(numpy.float)      # [-] Air Conditioning COP
        EffHeat     = df_LocationSummary.loc['EffHeat'].values.astype(numpy.float)  # [%] Heating Efficiency
        FanFlow     = df_LocationSummary.loc['Fan'].values.astype(numpy.float)      # [m^3 s^-1] Fan Max Flow Rate

        # Read Schedules (Sheet 4)
        file_doe_name_schedules = os.path.join("{}".format(DIR_DOE_PATH), "BLD{}".format(i+1),"BLD{}_Schedules.csv".format(i+1))
        df_Schedules = pandas.read_csv(file_doe_name_schedules,header=0,usecols=[0,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29],index_col=[0])

        #listof(listof weekday, sat, sun (list of 24 fractions)))
        SchEquip    = df_Schedules.iloc[0:3].values      # Equipment Schedule 24 hrs
        SchLight    = df_Schedules.iloc[3:6].values      # Light Schedule 24 hrs; Wkday=Sat=Sun=Hol
        SchOcc      = df_Schedules.iloc[6:9].values      # Occupancy Schedule 24 hrs
        SetCool     = df_Schedules.iloc[9:12].values     # Cooling Setpoint Schedule 24 hrs
        SetHeat     = df_Schedules.iloc[12:15].values    # Heating Setpoint Schedule 24 hrs; summer design
        SchGas      = df_Schedules.iloc[15:18].values    # Gas Equipment Schedule 24 hrs; wkday=sat
        SchSWH      = df_Schedules.iloc[18:21].values    # Solar Water Heating Schedule 24 hrs; wkday=summerdesign, sat=winterdesgin


        for j in range(3):

            # j = 3 built eras
            #print"\tEra: {} @j={}".format(BUILTERA[j], j)

            for k in range(17):

                # k = 17 climate zones
                #print "\tClimate zone: {} @k={}".format(ZONETYPE[k], k)

                B = Building(
                    hCeiling[j],                        # floorHeight by area
                    1,                                  # intHeatNight
                    1,                                  # intHeatDay
                    0.1,                                # intHeatFRad
                    0.1,                                # intHeatFLat
                    Infil[j],                           # infiltration rate Air Change per Hour (ACH) [hr^-1]
                    Vent[j]/1000.,                      # ventilation rate per unit floor area converted from liters to cubic meter [m^3 s^-1 m^-2]
                    glazing[j],                         # glazing ratio by area
                    Uwindow[j][k],                      # uValue by area, by climate type
                    SHGC[j][k],                         # SHGC, by area, by climate type
                    'AIR',                              # cooling condensation system type: AIR, WATER
                    COP[j][k],                          # cop by area, climate type
                    297,                                # coolSetpointDay = 24 C
                    297,                                # coolSetpointNight
                    293,                                # heatSetpointDay = 20 C
                    293,                                # heatSetpointNight
                    (HVAC[j][k]*1000.0)/AreaFloor[j],   # coolCap converted from kW per entire floor area to Watt per unit floor area [W m^-2]
                    EffHeat[j][k],                      # heatEff by area, climate type
                    293)                                # initialTemp at 20 C

                #Not defined in the constructor
                B.heatCap = (HEAT[j][k]*1000.0)/AreaFloor[j]         # heating Capacity converted to [W m^-2] by area, climate type
                B.Type = BLDTYPE[i]
                B.Era = BUILTERA[j]
                B.Zone = ZONETYPE[k]
                refDOE[i][j][k] = B

                # Define wall, mass(floor), roof
                # Reference from E+ for conductivity, thickness (reference below)

                # Material: (thermalCond, volHeat = specific heat * density)
                # Concrete = Material (1.311, 836.8 * 2240,"Concrete")
                # Insulation = Material (0.049, 836.8 * 265.0, "Insulation")
                # Gypsum = Material (0.16, 830.0 * 784.9, "Gypsum")
                # Wood = Material (0.11, 1210.0 * 544.62, "Wood")
                # Stucco = Material(0.6918,  837.0 * 1858.0, "Stucco")

                # Wall (1 in stucco, concrete, insulation, gypsum)
                # Check TypWall by area, by climate
                if TypeWall[j][k] == "MassWall":
                    #Construct wall based on R value of Wall from refDOE and properties defined above
                    # 1" stucco, 8" concrete, tbd insulation, 1/2" gypsum
                    Concrete = Material(1.311, 836.8 * 2240, "Concrete")
                    Insulation = Material (0.049, 836.8 * 265.0, "Insulation")
                    Gypsum = Material(0.17, 870.0 * 625.9, "Gypsum")
                    Wood = Material(0.11, 1210.0 * 544.62, "Wood")
                    Stucco = Material(0.6918, 837.0 * 1858.0, "Stucco")
                    Brick = Material(0.52, 900 * 899, "Brick")
                    Rbase = 0.271087 # R val based on stucco, concrete, gypsum
                    Rins = RvalWall[j][k] - Rbase #find insulation value
                    D_ins = Rins * Insulation.thermalCond # depth of insulation [m^2 K W^-1] * [W m^-1 K^-1] = m
                    if D_ins > 0.01:
                        thickness = [0.0254,0.0508,0.0508,0.0508,0.0508,D_ins,0.0127]
                        layers = [Brick,Concrete,Concrete,Concrete,Concrete,Insulation,Gypsum]
                    else:
                        #if it's less than 1 cm don't include in layers
                        thickness = [0.0254,0.0508,0.0508,0.0508,0.0508,0.0127]
                        layers = [Brick,Concrete,Concrete,Concrete,Concrete,Gypsum]

                    # wall = Element(0.08,0.92,thickness,layers,0.,293.,0.,"MassWall")
                    # thickness = [0.015,0.06,0.05,0.12,0.015]
                    # lay1 = Material(0.8,1820*863.90,"LimePlaster")
                    # lay2 = Material(0.52, 900 * 899, "brick1")
                    # lay3 = Material (0.049, 836.8 * 265.0, "Insulation")
                    # lay4 = Material(1.1, 900 * 2000, "brick2")
                    # lay5 = Material(0.16, 830.0 * 784.9, "Gypsum")
                    # layers = [lay1,lay2,lay3,lay4,lay5]
                    # wall = Element(thickness, layers, 0., 293., 0., "MassWall")
                    wallSun = Tsurf_Def(thickness, layers,300.,"MassWall")
                    wallShade = Tsurf_Def(thickness, layers,300.,"MassWall")

                    # If mass wall, assume mass floor (4" concrete)
                    # Mass (assume 4" concrete);
                    alb = 0.2
                    emis = 0.9
                    thickness = [0.054,0.054]
                    concrete = Material (1.31, 2240.0*836.8)
                    # mass = Element(thickness, [concrete, concrete], 0, 300, 1, "MassFloor")
                    mass = Tsurf_Def(thickness, [concrete, concrete],300,"MassFloor")

                elif TypeWall[j][k] == "WoodFrame":
                    # 0.01m wood siding, tbd insulation, 1/2" gypsum
                    Concrete = Material(1.311, 836.8 * 2240, "Concrete")
                    Insulation = Material(0.049, 836.8 * 265.0, "Insulation")
                    Gypsum = Material(0.16, 830.0 * 784.9, "Gypsum")
                    Wood = Material(0.11, 1210.0 * 544.62, "Wood")
                    Stucco = Material(0.6918, 837.0 * 1858.0, "Stucco")
                    Rbase = 0.170284091    # based on wood siding, gypsum
                    Rins = RvalWall[j][k] - Rbase
                    D_ins = Rins * Insulation.thermalCond # depth of insulation

                    if D_ins > 0.01:
                        thickness = [0.01,D_ins,0.0127]
                        layers = [Wood,Insulation,Gypsum]
                    else:
                        thickness = [0.01,0.0127]
                        layers = [Wood,Gypsum]
                    # wall = Element(thickness, layers, 0., 293., 0., "WoodFrameWall")
                    wallSun = Tsurf_Def(thickness, layers,300.,"WoodFrameWall")
                    wallShade = Tsurf_Def(thickness, layers,300.,"WoodFrameWall")

                    # If wood frame wall, assume wooden floor
                    alb = 0.2
                    emis = 0.9
                    thickness = [0.05,0.05]
                    wood = Material(1.31, 2240.0*836.8)
                    # mass = Element(thickness, [wood, wood], 0., 300., 1., "WoodFloor")
                    mass = Tsurf_Def(thickness, [wood, wood],300.,"WoodFloor")


                elif TypeWall[j][k] == "SteelFrame":
                    # 1" stucco, 8" concrete, tbd insulation, 1/2" gypsum
                    Concrete = Material(1.311, 836.8 * 2240, "Concrete")
                    Insulation = Material(0.049, 836.8 * 265.0, "Insulation")
                    Gypsum = Material(0.16, 830.0 * 784.9, "Gypsum")
                    Wood = Material(0.11, 1210.0 * 544.62, "Wood")
                    Stucco = Material(0.6918, 837.0 * 1858.0, "Stucco")
                    Rbase = 0.271087 # based on stucco, concrete, gypsum
                    Rins = RvalWall[j][k] - Rbase
                    D_ins = Rins * Insulation.thermalCond
                    if D_ins > 0.01:
                        thickness = [0.0254,0.0508,0.0508,0.0508,0.0508,D_ins,0.0127]
                        layers = [Stucco,Concrete,Concrete,Concrete,Concrete,Insulation,Gypsum]
                    else:    # If insulation is too thin, assume no insulation
                        thickness = [0.0254,0.0508,0.0508,0.0508,0.0508,0.0127]
                        layers = [Stucco,Concrete,Concrete,Concrete,Concrete,Gypsum]
                    # wall = Element(thickness, layers, 0., 293., 0., "SteelFrame")
                    wallSun = Tsurf_Def(thickness, layers,300.,"SteelFrame")
                    wallShade = Tsurf_Def(thickness, layers,300.,"SteelFrame")

                    # If mass wall, assume mass floor
                    # Mass (assume 4" concrete),
                    alb = 0.2
                    emis = 0.93
                    thickness = [0.05,0.05]
                    # mass = Element(thickness, [Concrete, Concrete], 0., 300., 1., "MassFloor")
                    mass = Tsurf_Def(thickness, [Concrete, Concrete],300.,"MassFloor")


                elif TypeWall[j][k] == "MetalWall":
                    # metal siding, insulation, 1/2" gypsum
                    alb = 0.2
                    emis = 0.9
                    Concrete = Material(1.311, 836.8 * 2240, "Concrete")
                    Insulation = Material(0.049, 836.8 * 265.0, "Insulation")
                    Gypsum = Material(0.16, 830.0 * 784.9, "Gypsum")
                    Wood = Material(0.11, 1210.0 * 544.62, "Wood")
                    Stucco = Material(0.6918, 837.0 * 1858.0, "Stucco")
                    D_ins = max((RvalWall[j][k] * Insulation.thermalCond)/2, 0.01) #use derived insulation thickness or 0.01 based on max
                    thickness = [D_ins,D_ins,0.0127]
                    materials = [Insulation,Insulation,Gypsum]
                    # wall = Element(thickness, materials, 0, 293, 0, "MetalWall")
                    wallSun = Tsurf_Def(thickness, materials,300,"MetalWall")
                    wallShade = Tsurf_Def(thickness, materials,300,"MetalWall")

                    # Mass (assume 4" concrete);
                    alb = 0.2
                    emis = 0.9
                    Concrete = Material(1.311, 836.8 * 2240, "Concrete")
                    Insulation = Material(0.049, 836.8 * 265.0, "Insulation")
                    Gypsum = Material(0.16, 830.0 * 784.9, "Gypsum")
                    Wood = Material(0.11, 1210.0 * 544.62, "Wood")
                    Stucco = Material(0.6918, 837.0 * 1858.0, "Stucco")
                    thickness = [0.05, 0.05]
                    concrete = Material(1.31, 2240.0*836.8)
                    # mass = Element(thickness, [concrete, concrete], 0., 300., 1., "MassFloor")

                    mass = Tsurf_Def(thickness, [concrete, concrete],300.,"MassFloor")


                elif TypeWall[j][k] == "BaselWall":

                    mat1 = Material(0.6918,1555146,'Mat1')
                    mat2 = Material(1.311,1874432,'Mat2')
                    mat3 = Material(0.16,651467,'Mat3')
                    materials = [mat1,mat2,mat2,mat2,mat2,mat3]
                    thickness = [0.0254,0.0508,0.0508,0.0508,0.0508,0.0127]
                    wallSun = Tsurf_Def(thickness, materials,300,"BaselWall")
                    wallShade = Tsurf_Def(thickness, materials,300,"BaselWall")

                    Concrete = Material(1.311, 836.8 * 2240, "Concrete")
                    thickness = [0.05, 0.05]
                    # mass = Element(thickness, [Concrete, Concrete], 0., 300., 1., "MassFloor")

                    mass = Tsurf_Def(thickness, [Concrete, Concrete],300.,"MassFloor")


                # Roof
                if TypeRoof[j][k] == "IEAD": #Insulation Entirely Above Deck
                    # IEAD-> membrane, insulation, decking
                    alb = 0.2
                    emis = 0.93
                    Concrete = Material(1.311, 836.8 * 2240, "Concrete")
                    Insulation = Material(0.049, 836.8 * 265.0, "Insulation")
                    Gypsum = Material(0.16, 830.0 * 784.9, "Gypsum")
                    Wood = Material(0.11, 1210.0 * 544.62, "Wood")
                    Stucco = Material(0.6918, 837.0 * 1858.0, "Stucco")
                    D_ins = max(RvalRoof[j][k] * Insulation.thermalCond/2.,0.01)
                    # roof = Element([D_ins, D_ins], [Insulation, Insulation], 0., 293., 0., "IEAD")
                    roofImp = Tsurf_Def([D_ins, D_ins], [Insulation, Insulation],300.,"IEAD")

                    Soil = Material(6.84, 2.2137*10**6,'RoofVegMat')
                    roofVeg = Tsurf_Def([D_ins, D_ins], [Soil, Soil],300.,"IEAD")

                elif TypeRoof[j][k] == "Attic":
                    # IEAD-> membrane, insulation, decking
                    alb = 0.2
                    emis = 0.9
                    # D_ins = max(RvalRoof[j][k] * Insulation.thermalCond/2.,0.01)
                    Concrete = Material(1.311, 836.8 * 2240, "Concrete")
                    Insulation = Material(1.3, 1800000, "Insulation")
                    Gypsum = Material(0.16, 830.0 * 784.9, "Gypsum")
                    Wood = Material(0.11, 1210.0 * 544.62, "Wood")
                    Stucco = Material(0.6918, 837.0 * 1858.0, "Stucco")
                    D_ins = 0.2
                    # roof = Element([D_ins, D_ins], [Insulation, Insulation], 0., 293., 0., "Attic")
                    roofImp = Tsurf_Def([D_ins, D_ins], [Insulation, Insulation],300.,"Attic")

                    Soil = Material(6.84, 2.2137 * 10 ** 6, 'RoofVegMat')
                    roofVeg = Tsurf_Def([D_ins, D_ins], [Soil, Soil], 300., "IEAD")


                elif TypeRoof[j][k] == "MetalRoof":
                    # IEAD-> membrane, insulation, decking
                    alb = 0.2
                    emis = 0.9
                    Concrete = Material(1.311, 836.8 * 2240, "Concrete")
                    Insulation = Material(0.049, 836.8 * 265.0, "Insulation")
                    Gypsum = Material(0.16, 830.0 * 784.9, "Gypsum")
                    Wood = Material(0.11, 1210.0 * 544.62, "Wood")
                    Stucco = Material(0.6918, 837.0 * 1858.0, "Stucco")
                    D_ins = max(RvalRoof[j][k] * Insulation.thermalCond/2.,0.01)
                    # roof = Element([D_ins, D_ins], [Insulation, Insulation], 0., 293., 0., "MetalRoof")
                    roofImp = Tsurf_Def([D_ins, D_ins], [Insulation, Insulation],300.,"MetalRoof")

                    Soil = Material(6.84, 2.2137 * 10 ** 6, 'RoofVegMat')
                    roofVeg = Tsurf_Def([D_ins, D_ins], [Soil, Soil], 300., "IEAD")


                elif TypeRoof[j][k] == "BaselRoof":

                    mat = Material(0.94,1400000,'RoofMat1')
                    thickness = [0.05819,0.05819]
                    # roof = Element(thickness, [mat, mat], 0., 293., 0., "BaselRoof")
                    roofImp = Tsurf_Def(thickness, [mat, mat],300.,"BaselRoof")

                    Soil = Material(6.84, 2.2137 * 10 ** 6, 'RoofVegMat')
                    roofVeg = Tsurf_Def(thickness, [Soil, Soil], 300., "IEAD")



                # Define building energy model, set fraction of the urban floor space of this typology to zero
                refBEM[i][j][k] = BEMDef(B, mass, wallSun, wallShade, roofImp, roofVeg, 0.0)
                #refBEM[i][j][k].building.FanMax = FanFlow[j][k] # max fan flow rate [m^3 s^-1] per DOE

                Schedule[i][j][k] = SchDef()

                Schedule[i][j][k].Elec = SchEquip   # 3x24 matrix of schedule for fraction electricity (WD,Sat,Sun)
                Schedule[i][j][k].Light = SchLight  # 3x24 matrix of schedule for fraction light (WD,Sat,Sun)
                Schedule[i][j][k].Gas = SchGas      # 3x24 matrix of schedule for fraction gas (WD,Sat,Sun)
                Schedule[i][j][k].Occ = SchOcc      # 3x24 matrix of schedule for fraction occupancy (WD,Sat,Sun)
                Schedule[i][j][k].Cool = SetCool    # 3x24 matrix of schedule for fraction cooling temp (WD,Sat,Sun)
                Schedule[i][j][k].Heat = SetHeat    # 3x24 matrix of schedule for fraction heating temp (WD,Sat,Sun)
                Schedule[i][j][k].SWH = SchSWH      # 3x24 matrix of schedule for fraction SWH (WD,Sat,Sun

                Schedule[i][j][k].Qelec = Elec[j]                   # [W m^-2] (max) for electrical plug process
                Schedule[i][j][k].Qlight = Light[j]                 # [W m^-2] (max) for light
                Schedule[i][j][k].Nocc = Occupant[j]/AreaFloor[j]   # [Person m^-2]
                Schedule[i][j][k].Qgas = Gas[j]                     # [W m^-2] (max) for gas
                Schedule[i][j][k].Vent = Vent[j]/1000.0             # [m^3 m^-2] per person
                Schedule[i][j][k].Vswh = SHW[j]/AreaFloor[j]        # litres per hour per m^2 of floor


    # if not test serialize refDOE,refBEM,Schedule and store in resources
    if not serialize_output:

        # create a binary file for serialized obj
        pkl_file_path = os.path.join(DIR_CURR,'readDOE.pkl')
        pickle_readDOE = open(pkl_file_path, 'wb')

        # dump in ../resources
        # Pickle objects, protocol 1 b/c binary file
        cPickle.dump(refDOE, pickle_readDOE,1)
        cPickle.dump(refBEM, pickle_readDOE,1)
        cPickle.dump(Schedule, pickle_readDOE,1)

        pickle_readDOE.close()

    return refDOE, refBEM, Schedule

if __name__ == "ReadDOE": # "__main__":#

    # Set to True only if you want create new .pkls of DOE refs
    # Use --serialize switch to serialize the readDOE data
    print(sys.argv)
    if len(sys.argv)> 1 and sys.argv[1]=="--serialize":
        refDOE, refBEM, Schedule = readDOE(True)
    else:
        refDOE, refBEM, Schedule = readDOE(False)


# Material ref from E+
#     1/2IN Gypsum,            !- Name
#     Smooth,                  !- Roughness
#     0.0127,                  !- Thickness [m]
#     0.1600,                  !- Conductivity [W m^-1 K^-1]
#     784.9000,                !- Density [kg m^-3]
#     830.0000,                !- Specific Heat [J kg^-1 K^-1]
#     0.9000,                  !- Thermal Absorptance
#     0.9200,                  !- Solar Absorptance
#     0.9200;                  !- Visible Absorptance
#
# Material,
#     1IN Stucco,              !- Name
#     Smooth,                  !- Roughness
#     0.0253,                  !- Thickness
#     0.6918,                  !- Conductivity
#     1858.0000,               !- Density
#     837.0000,                !- Specific Heat
#     0.9000,                  !- Thermal Absorptance
#     0.9200,                  !- Solar Absorptance
#     0.9200;                  !- Visible Absorptance
#
# Material,
#     8IN CONCRETE HW,  !- Name
#     Rough,                   !- Roughness
#     0.2032,                  !- Thickness
#     1.3110,                  !- Conductivity
#     2240.0000,               !- Density
#     836.8000,                !- Specific Heat
#     0.9000,                  !- Thermal Absorptance
#     0.7000,                  !- Solar Absorptance
#     0.7000;                  !- Visible Absorptance
#
# Material,
#     Mass NonRes Wall Insulation, !- Name
#     MediumRough,             !- Roughness
#     0.0484268844343858,      !- Thickness
#     0.049,                   !- Conductivity
#     265.0000,                !- Density
#     836.8000,                !- Specific Heat
#     0.9000,                  !- Thermal Absorptance
#     0.7000,                  !- Solar Absorptance
#     0.7000;                  !- Visible Absorptance
#
# Material,
#     Std Wood 6inch,          !- Name
#     MediumSmooth,            !- Roughness
#     0.15,                    !- Thickness
#     0.12,                    !- Conductivity
#     540.0000,                !- Density
#     1210,                    !- Specific Heat
#     0.9000000,               !- Thermal Absorptance
#     0.7000000,               !- Solar Absorptance
#     0.7000000;               !- Visible Absorptance! Common Materials
#
# Material,
#     Wood Siding,             !- Name
#     MediumSmooth,            !- Roughness
#     0.0100,                  !- Thickness
#     0.1100,                  !- Conductivity
#     544.6200,                !- Density
#     1210.0000,               !- Specific Heat
#     0.9000,                  !- Thermal Absorptance
#     0.7800,                  !- Solar Absorptance
#     0.7800;                  !- Visible Absorptance
