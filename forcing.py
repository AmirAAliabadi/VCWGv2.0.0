import numpy

"""
Developed by Bruno Bueno
Building Technology, Massachusetts Institute of Technology (MIT), Cambridge, U.S.A.
Last update: 2012
"""

class Forcing (object):
    def __init__(self,staTemp=None,weather=None):
        # weather: Weather obj
        # staTemp: list of hourly air temperature for simulation period
        # Define default values for instance variables, when the type can be mutable:
        if staTemp==None and weather==None:
            self.deepTemp = None    # deep soil temperature [K]
            self.waterTemp = None   # ground water temp, set to temp at 2m
            self.infra = None       # horizontal Infrared Radiation Intensity [W m^-2]
            self.uDir = None        # wind direction
            self.hum  = None        # specific humidity [kg kg^-1]
            self.pres = None        # Pressure [Pa]
            self.temp = None        # air temperature [C]
            self.rHum = None        # Relative humidity [%]
            self.dir = None         # Direct solar radiation [W m^-2]
            self.dif = None         # Diffusive solar radiation [W m^-2]. Amount of solar radiation received from the sky (excluding the solar disk) on a horizontal surface
            self.prec = None        # Precipitation [mm h^-1]
            self.wind = None        # wind speed [m s^-1]
        else:
            self.deepTemp = sum(staTemp)/float(len(staTemp))
            self.waterTemp = sum(staTemp)/float(len(staTemp))
            self.infra = weather.staInfra
            self.uDir = weather.staUdir
            self.hum  = weather.staHum
            self.pres = weather.staPres
            self.temp = weather.staTemp
            self.rHum = weather.staRhum
            self.dir = weather.staDir
            self.dif = weather.staDif
            self.prec = [p for p in weather.staRobs]
            self.wind = weather.staUmod
            self.Year = weather.staYear
            self.Month = weather.staMonth
            self.Day = weather.staDay
            self.Hour = weather.staHour
            self.Min = weather.staMin

            self.infra = numpy.append(self.infra, self.infra[-1])
            self.uDir = numpy.append(self.uDir, self.uDir[-1])
            self.hum = numpy.append(self.hum, self.hum[-1])
            self.pres = numpy.append(self.pres, self.pres[-1])
            self.temp = numpy.append(self.temp, self.temp[-1])
            self.rHum = numpy.append(self.rHum, self.rHum[-1])
            self.dir = numpy.append(self.dir, self.dir[-1])
            self.dif = numpy.append(self.dif, self.dif[-1])
            self.prec = numpy.append(self.prec, self.prec[-1])
            self.wind = numpy.append(self.wind, self.wind[-1])
            self.Year = numpy.append(self.Year, self.Year[-1])
            self.Month = numpy.append(self.Month, self.Month[-1])
            self.Day = numpy.append(self.Day, self.Day[-1])
            self.Hour = numpy.append(self.Hour, 24)
            self.Min = numpy.append(self.Min, self.Min[-1])


    def __repr__(self):
        return "Forcing: deepT={a}, waterT={b}".format(
            a=int(self.deepTemp) if self.deepTemp else None,
            b=int(self.waterTemp) if self.waterTemp else None
            )
