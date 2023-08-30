import numpy
import pandas as pd
from psychrometrics import psychrometrics,HumFromRHumTemp
import functools
import os

"""
For top forcing generate an EPW file using the excel user input file
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: May 2021
"""

def write_epw(TopForcingFile,epw_precision,climate_file,timeInitial,climate_file_ERA5=None):
    """
    Generate new EPW file based on TopForcing file
    """
    EPW_path = os.path.join('resources','epw')
    TopForcing_path = os.path.join('resources','TopForcing')
    rawEPW = os.path.join(EPW_path,climate_file)
    TopForcing = os.path.join(TopForcing_path,TopForcingFile)

    # Read lat and lon
    with open(TopForcing) as Forcing_file:
        lines = Forcing_file.readlines()
    line0 = lines[0].split(',')
    lon = float(line0[1])
    lat = float(line0[2])
    GMT = float(line0[3])

    # open the TopForcingFile and clean it if there is any NaN or extra blank space
    df = pd.read_csv(TopForcing,header=1,parse_dates=['Time'])
    df_obj = df.select_dtypes(['object'])
    df[df_obj.columns] = df_obj.apply(lambda x: x.str.strip())
    index_nan_LWR = df[df['IncomingLongwaveRadiation[Wm-2]']=='NaN'].index.values
    index_nan_SWRDir = df[df['IncomingDirShortwaveRadiation[Wm-2]']=='NaN'].index.values
    index_nan_SWRDif = df[df['IncomingDifShortwaveRadiation[Wm-2]']=='NaN'].index.values
    index_nan_T = df[df['AirTemperature[C]']=='NaN'].index.values
    index_nan_Pressure = df[df['BarometricPressure[KPa]']=='NaN'].index.values
    index_nan_RH = df[df['RelativeHumidity[%]']=='NaN'].index.values
    index_nan_WindDir = df[df['WindDirection[Deg]']=='NaN'].index.values
    index_nan_WindSpeed = df[df['WindSpeed[ms-1]']=='NaN'].index.values
    for i_nan_LWR in range(0,len(index_nan_LWR)):
        df.loc[index_nan_LWR[i_nan_LWR],'IncomingLongwaveRadiation[Wm-2]'] = df.loc[index_nan_LWR[i_nan_LWR]-1,'IncomingLongwaveRadiation[Wm-2]']
    for i_nan_SWRDir in range(0,len(index_nan_SWRDir)):
        df.loc[index_nan_SWRDir[i_nan_SWRDir],'IncomingDirShortwaveRadiation[Wm-2]'] = df.loc[index_nan_SWRDir[i_nan_SWRDir]-1,'IncomingDirShortwaveRadiation[Wm-2]']
    for i_nan_SWRDif in range(0,len(index_nan_SWRDif)):
        df.loc[index_nan_SWRDif[i_nan_SWRDif],'IncomingDifShortwaveRadiation[Wm-2]'] = df.loc[index_nan_SWRDif[i_nan_SWRDif]-1,'IncomingDifShortwaveRadiation[Wm-2]']
    for i_nan_T in range(0,len(index_nan_T)):
        df.loc[index_nan_T[i_nan_T],'AirTemperature[C]'] = df.loc[index_nan_T[i_nan_T]-1,'AirTemperature[C]']
    for i_nan_Pre in range(0,len(index_nan_Pressure)):
        df.loc[index_nan_Pressure[i_nan_Pre],'BarometricPressure[KPa]'] = df.loc[index_nan_Pressure[i_nan_Pre]-1,'BarometricPressure[KPa]']
    for i_nan_RH in range(0,len(index_nan_RH)):
        df.loc[index_nan_RH[i_nan_RH],'RelativeHumidity[%]'] = df.loc[index_nan_RH[i_nan_RH]-1,'RelativeHumidity[%]']
    for i_nan_WDir in range(0,len(index_nan_WindDir)):
        df.loc[index_nan_WindDir[i_nan_WDir],'WindDirection[Deg]'] = df.loc[index_nan_WindDir[i_nan_WDir]-1,'WindDirection[Deg]']
    for i_nan_WS in range(0,len(index_nan_WindSpeed)):
        df.loc[index_nan_WindSpeed[i_nan_WS],'WindSpeed[ms-1]'] = df.loc[index_nan_WindSpeed[i_nan_WS]-1,'WindSpeed[ms-1]']

    # Build hourly data
    df_hourly = df.loc[df['Time'].dt.minute == 0]
    df_hourly = df_hourly.fillna(method='bfill')

    # Calculate specific humidity
    SpecHum = []
    for i in range(0,len(df_hourly)):
        SpecHum.append(HumFromRHumTemp(float(df_hourly['RelativeHumidity[%]'].values[i]),float(df_hourly['AirTemperature[C]'].values[i]),
                                       float(df_hourly['BarometricPressure[KPa]'].values[i])*1000))

    # Calculate dew point temperature
    Tdp = []
    for i in range(0, len(df_hourly)):
        Tdp.append(psychrometrics(float(df_hourly['AirTemperature[C]'].values[i])+273.15,SpecHum[i],
                                  float(df_hourly['BarometricPressure[KPa]'].values[i])*1000)[4])
    # Insert Tdp in the dataframe
    df_hourly.insert(loc=8,column='DewPointTemperature[C]',value=Tdp,allow_duplicates=True)


    # Open .epw file
    with open(rawEPW) as f:
        lines = f.readlines()
    climate_data = []
    for i in range(len(lines)):
        climate_data.append(list(lines[i].split(",")))
    _header = climate_data[0:8]
    _header[0][6] = str(lat)
    _header[0][7] = str(lon)
    _header[0][8] = str(GMT)
    epwinput = climate_data[8:]
    epw_prec = epw_precision  # precision of epw file input

    for iJ in range(len(df_hourly)):
        # [iJ+timeInitial[Month]-8] = increments along every weather timestep in epw
        # [6 to 21]                       = column data of epw
        # Year
        epwinput[iJ + timeInitial - 8][0] = "{0:.{1}f}".format(df['Time'].dt.year.values[iJ], 0)
        # dry bulb temperature [C]
        epwinput[iJ + timeInitial - 8][6] = "{0:.{1}f}".format(df_hourly['AirTemperature[C]'].values[iJ], epw_prec)
        # dew point temperature [C]
        epwinput[iJ + timeInitial - 8][7] = "{0:.{1}f}".format(df_hourly['DewPointTemperature[C]'].values[iJ], epw_prec)
        # relative humidity [%]
        epwinput[iJ + timeInitial - 8][8] = "{0:.{1}f}".format(float(df_hourly['RelativeHumidity[%]'].values[iJ]), epw_prec)
        # Pressure [Pa]
        epwinput[iJ + timeInitial - 8][9] = "{0:.{1}f}".format(float(df_hourly['BarometricPressure[KPa]'].values[iJ])*1000.0, epw_prec)
        # LWR [W m^-2]
        epwinput[iJ + timeInitial - 8][12] = "{0:.{1}f}".format(float(df_hourly['IncomingLongwaveRadiation[Wm-2]'].values[iJ]), epw_prec)
        # SWR_Diff [W m^-2]
        epwinput[iJ + timeInitial - 8][15] = "{0:.{1}f}".format(float(df_hourly['IncomingDifShortwaveRadiation[Wm-2]'].values[iJ]), epw_prec)
        # SWR_Dir [W m^-2]
        epwinput[iJ + timeInitial - 8][14] = "{0:.{1}f}".format(float(df_hourly['IncomingDirShortwaveRadiation[Wm-2]'].values[iJ]), epw_prec)
        # Wind Direction [deg]
        epwinput[iJ + timeInitial - 8][20] = "{0:.{1}f}".format(float(df_hourly['WindDirection[Deg]'].values[iJ]), epw_prec)
        # wind speed [m s^-1]
        epwinput[iJ + timeInitial - 8][21] = "{0:.{1}f}".format(float(df_hourly['WindSpeed[ms-1]'].values[iJ]), epw_prec)
        # Precipitation [mm h^-1]
        epwinput[iJ + timeInitial - 8][33] = "{0:.{1}f}".format(float(df_hourly['RainFall[mm]'].values[iJ]), 5)


    # Writing new EPW file
    TopForcing_EPW = EPW_path+r'\TopForcing.epw'
    epw_new_id = open(TopForcing_EPW, "w")

    for i in range(8):
        new_epw_line = '{}'.format(functools.reduce(lambda x, y: x + "," + y, _header[i]))
        epw_new_id.write(new_epw_line)

    for i in range(len(epwinput)):
        printme = ""
        for ei in range(34):
            printme += "{}".format(epwinput[i][ei]) + ','
        printme = printme + "{}".format(epwinput[i][ei])
        new_epw_line = "{0}\n".format(printme)
        epw_new_id.write(new_epw_line)

    epw_new_id.close()


