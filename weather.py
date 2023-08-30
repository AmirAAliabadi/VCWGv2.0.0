from Utilities import read_csv, str2fl
from math import pow, log, exp
from psychrometrics import HumFromRHumTemp
import numpy
import os

"""
Developed by Bruno Bueno
Building Technology, Massachusetts Institute of Technology (MIT), Cambridge, U.S.A.
Last update: 2012
"""

class Weather(object):

    """
    Weather
    Read epw file
    http://bigladdersoftware.com/epx/docs/8-2/auxiliary-programs/epw-csv-format-inout.html
    properties
        location  # location name
        staTemp   % air temperature [C}
        staTdp    % dewpoint temperature [C]
        staRhum   % air relative humidity [%]
        staPres   % air pressure [Pa]
        staInfra  % horizontal Infrared Radiation Intensity [W m^-2]
        staHor    % horizontal radiation
        staDir    % normal solar direct radiation [W m^-2]
        staDif    % horizontal solar diffuse radiation [W m^-2]
        staUdir   % wind direction [deg]
        staUmod   % wind speed [m s^-1]
        staRobs   % Precipitation [mm h^-1]
        staHum    % specific humidty [kg kg^-1]
    """

    def __init__(self,EPW_file,HI,HF):
        #HI: Julian start date
        #HF: Julian final date
        #H1 and HF define the row we want
        EPW_path = os.path.join('resources','epw')
        Climate_file = os.path.join(EPW_path,EPW_file)
        # Open .epw file and feed csv data to self.climate_data
        with open(Climate_file) as f:
            lines = f.readlines()
        self.climate_data = []
        for i in range(len(lines)):
            self.climate_data.append(list(lines[i].split(",")))

        # Read header lines (1 to 8) from EPW and ensure TMY2 format.
        self._header = self.climate_data[0:8]
        # Read Lat, Long (line 1 of EPW)
        self.lat = float(self._header[0][6])
        self.lon = float(self._header[0][7])
        self.GMT = float(self._header[0][8])

        # Read in soil temperature data (assumes this is always there)
        # ref: http://bigladdersoftware.com/epx/docs/8-2/auxiliary-programs/epw-csv-format-inout.html
        soilData = self._header[3]
        self.nSoil = int(soilData[1])  # Number of ground temperature depths
        self.Tsoil = numpy.zeros((self.nSoil, 12))  # nSoil x 12 matrix for soil temperture (K)
        self.depth_soil = numpy.zeros((self.nSoil, 1))  # nSoil x 1 matrix for soil depth (m)

        # Read monthly data for each layer of soil from EPW file
        for i in range(self.nSoil):
            self.depth_soil[i][0] = float(soilData[2 + (i * 16)])  # get soil depth for each nSoil
            # Monthly data
            for j in range(12):
                self.Tsoil[i][j] = float(soilData[6 + (i * 16) + j]) + 273.15  # 12 months of soil T for specific depth

        self.location = self.climate_data[0][1]
        cd = self.climate_data[HI:HF+1]
        self.staYear = numpy.array([cd[i][0] for i in range(len(cd))],dtype=float)
        self.staMonth = numpy.array([cd[i][1] for i in range(len(cd))],dtype=float)
        self.staDay = numpy.array([cd[i][2] for i in range(len(cd))],dtype=float)
        self.staHour = numpy.array([cd[i][3] for i in range(len(cd))],dtype=float)
        self.staMin = numpy.array([cd[i][4] for i in range(len(cd))],dtype=float)
        self.staTemp = numpy.array([cd[i][6] for i in range(len(cd))],dtype=float)           # drybulb [C]
        self.staTdp = numpy.array([cd[i][7] for i in range(len(cd))],dtype=float)            # dewpoint [C]
        self.staRhum = numpy.array([cd[i][8] for i in range(len(cd))],dtype=float)           # air relative humidity [%]
        self.staPres = numpy.array([cd[i][9] for i in range(len(cd))],dtype=float)           # air pressure [Pa]
        self.staInfra = numpy.array([cd[i][12] for i in range(len(cd))],dtype=float)         # horizontal Infrared Radiation Intensity [W m^-2]
        self.staHor = numpy.array([cd[i][13] for i in range(len(cd))],dtype=float)           # horizontal radiation [W m^-2]
        self.staDir = numpy.array([cd[i][14] for i in range(len(cd))],dtype=float)           # normal solar direct radiation [W m^-2]
        self.staDif = numpy.array([cd[i][15] for i in range(len(cd))],dtype=float)           # horizontal solar diffuse radiation [W m^-2]
        self.staUdir = numpy.array([cd[i][20] for i in range(len(cd))],dtype=float)          # wind direction [deg]
        self.staUmod = numpy.array([cd[i][21] for i in range(len(cd))],dtype=float)          # wind speed [m s^-1]
        self.staRobs = numpy.array([cd[i][33] for i in range(len(cd))],dtype=float)          # Precipitation [mm h^-1]
        self.staHum = [0.0] * len(self.staTemp)                              # specific humidity [kg kg^-1]
        for i in range(len(self.staTemp)):
            self.staHum[i] = HumFromRHumTemp(self.staRhum[i], self.staTemp[i], self.staPres[i])
        print(self.staHum)
        self.staTemp = [s+273.15 for s in self.staTemp]                      # air temperature [K]
        print(HI)
        print(HF)
    def __repr__(self):
        return "Weather: {a}, HI Tdb:{b}, HF Tdb:{c}".format(
            a=self.location,
            b=self.staTemp[0]-273.15,
            c=self.staTemp[-1]-273.15
            )
