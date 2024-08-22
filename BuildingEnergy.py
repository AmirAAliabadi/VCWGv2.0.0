
from psychrometrics import psychrometrics, moist_air_density
import logging
import numpy
import copy
"""
Calculate building characteristics
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: February 2020
Originally developed by Bruno Bueno
"""
class Building(object):
    """

    properties
        % Building parameters
        floorHeight         % floor height [m]
        intHeat;            % time step internal heat gains per unit floor area [W m^-2] (bld) (sensible only)
        intHeatNight;       % nighttime internal heat gains per unit floor area [W m^-2] (floor)
        intHeatDay;         % daytime internal heat gains per unit floor area [W m^-2] (floor)
        intHeatFRad;        % radiant fraction of internal gains
        intHeatFLat;        % latent fraction of internal gains
        infil;              % Infiltration Air Change per Hour (ACH) [hr^-1]
        vent;               % Ventilation rate per unit floor area [m^3 s^-1 m^-2]
        glazingRatio;       % glazing ratio
        uValue;             % window U-value [W m^-2 K^-1] (including film coeff)
        shgc;               % window Solar Heat Gain Coefficient (SHGC)
        condType;           % cooling condensation system type {'AIR', 'WATER'}
        cop;                % COP of the cooling system (nominal)
        coolSetpointDay;    % daytime indoor cooling set-point [K]
        coolSetpointNight;  % nighttime indoor cooling set-point [K]
        heatSetpointDay;    % daytime indoor heating set-point [K]
        heatSetpointNight;  % nighttime indoor heating set-point [K]
        coolCap;            % rated cooling system capacity [W m^-2]
        heatCap;            % rated heating system capacity [W m^-2]
        heatEff;            % heating system efficiency (-)
        canyon_fraction     # fraction of waste heat released to canyon, default = 1
        mSys;               % HVAC supply mass flowrate [kg s^-1 m^-2]
        indoorTemp;         % indoor air temperature [K]
        indoorHum;          % indoor specific humidity [kgv kga^-1]
        Twb;                % wetbulb temperature [C]
        Tdp;                % dew point [C]
        indoorRhum;         % indoor relative humidity [%]

        area_floor;         % total floor space of the BEM
        FanMax;             % max fan flow rate [m^3 s^-1] per DOE
        nFloor;             % number of floors
        RadFOcc;            % Radiant fraction of occupant
        LatFOcc;            % Latent fraction of occupant
        RadFEquip;          % Radiant fraction of equipment
        RadFLight;          % Radiant fraction of light

        Type;               % DOE reference building type
        Era;                % PRE80, PST80, NEW
        Zone;               % Climate zone number

        % Calculated values
        sensCoolDemand;     % building sensible cooling demand per unit building footprint area [W m^-2]
        sensHeatDemand;     % building sensible heating demand per unit building footprint area [W m^-2]
        copAdj;             % adjusted COP per temperature
        dehumDemand;        % Latent heat demand for dehumidification of air per unit building footprint area [W m^-2]
        coolConsump;        % cooling energy consumption per unit building footprint area OR per unit floor area [W m^-2]
        heatConsump;        % heating energy consumption per unit floor area [W m^-2]
        sensWaste;          % sensible waste heat per unit building footprint area [W m^-2]
        latWaste;           % lat waste heat per unit building footprint area [W m^-2]
        fluxMass;           % mass surface heat flux [W m^-2] (mass to indoor air)
        fluxWall;           % wall surface heat flux [W m^-2] (wall to inside)
        fluxRoof;           % roof surface heat flux [W m^-2] (roof to inside)
        fluxSolar;          % solar heat gain per unit floor area [W m^-2] through window (SHGC)
        fluxWindow;         % heat gain/loss from window per unit floor area [W m^-2] (U-value)
        fluxInterior;       % internal heat gain adjusted for latent/LW heat per unit floor area [W m^-2]
        fluxInfil;          % heat flux from infiltration per unit floor area [W m^-2]
        fluxVent;           % heat flux from ventilation per unit floor area [W m^-2]
        ElecTotal;          % total electricity consumption per unit floor area [W m^-2]
        GasTotal;           % total gas consumption per unit floor area [W m^-2]
        Qhvac;              % total heat removed (sensible + latent) per unit building footprint area [W m^-2] (calculated in cooling system)
        Qheat;              % total heat added (sensible only) per unit building footprint area [W m^-2] (calculated in heating system)
    """

    TEMPERATURE_COEFFICIENT_CONFLICT_MSG = "FATAL ERROR!"

    def __init__(self,floorHeight,intHeatNight,intHeatDay,intHeatFRad,\
            intHeatFLat,infil,vent,glazingRatio,uValue,shgc,\
            condType,cop,coolSetpointDay,coolSetpointNight,\
            heatSetpointDay,heatSetpointNight,coolCap,heatEff,initialTemp):

            self.floorHeight =float(floorHeight)        # floor height
            self.intHeat = intHeatNight                 # timestep internal sensible heat gain per unit floor area [W m^-2]
            self.intHeatNight = intHeatNight            # nighttime internal heat gain per unit floor area [W m^-2]
            self.intHeatDay = intHeatDay                # daytime internal heat gain per unit floor area [W m^-2]
            self.intHeatFRad = intHeatFRad              # internal gain radiant fraction
            self.intHeatFLat = intHeatFLat              # internal gain latent fraction
            self.infil = infil                          # Infiltration Air Change per Hour (ACH) [hr^-1]
            self.vent = vent                            # Ventilation rate per unit floor area [m^3 s^-1 m^-2]
            self.glazingRatio = glazingRatio            # glazing ratio
            self.uValue = uValue                        # window U-value [W m^-2 K^-1] including film coefficient
            self.shgc = shgc                            # window Solar Heat Gain Coefficient (SHGC), fraction of radiation that is admitted through a window
            self.condType = condType                    # cooling condensation system type: AIR, WATER
            self.cop = cop                              # COP of cooling system (nominal)
            self.coolSetpointDay = coolSetpointDay      # daytime indoor cooling setpoint [K]
            self.coolSetpointNight = coolSetpointNight  # nighttime indoor heating setpoint [K]
            self.heatSetpointDay = heatSetpointDay      # daytimge indoor heating setpoint [K]
            self.heatSetpointNight = heatSetpointNight  # nighttime indoor heating setpoint [K]
            self.coolCap = coolCap                      # rated cooling system capacity per floor area [W m^-2]
            self.heatEff = heatEff                      # heating system capacity (-)
            self.mSys = coolCap/1004./(min(coolSetpointDay,coolSetpointNight)-14-273.15) # HVAC supply mass flowrate [kg s^-1 m^-2]
            self.indoorTemp = initialTemp               # Indoor Air Temperature [K]
            self.indoorHum = 0.006                      # Indoor specific humidity [kgv kga^-1] fixed at 40% RH at 20C
            self.heatCap = 999                          # rated heating system capacity per floor area [W m^-2]
            self.copAdj = cop                           # adjusted COP per temperature
            self.canyon_fraction = 1.0                  # Default canyon fraction
            self.sensWaste = 0

            self.Type = "null"                          # DOE reference building type
            self.Era = "null"                           # pre80, pst80, new
            self.Zone = "null"                          # Climate zone number

            # Logger will be disabled by default unless explicitly called in tests
            self.logger = logging.getLogger(__name__)

    def __repr__(self):
        return "BuildingType: {a}, Era: {b}, Zone: {c}".format(
            a=self.Type,
            b=self.Era,
            c=self.Zone
            )

    def is_near_zero(self,val,tol=1e-14):
        return abs(float(val)) < tol

    def BEMCalc(self,canTemp,canHum,BEM,MeteoData,ParCalculation,simTime,Geometry_m,FractionsRoof,SWR):

        """
        ------
        INPUT:
        canTemp: Average canyon temperature [K]
        canHum: Average canyon specific humidity [kg kg^-1]
        BEM: Building energy parameters
        MeteoData: Forcing variables
        ParCalculation: General calculation parameters
        simTime: Simulation time parameters
        Geometry_m: Geometric parameters
        FractionsRoof:
        SWR: Shortwave radiation fluxes [W m^-2]
        -------
        OUTPUT:
        sensCoolDemand: building sensible cooling demand per unit building footprint area [W m^-2]
        sensHeatDemand: building sensible heating demand per unit building footprint area [W m^-2]
        copAdj: adjusted COP per temperature
        dehumDemand: Latent heat demand for dehumidification of air per unit building footprint area [W m^-2]
        coolConsump: cooling energy consumption per unit building footprint area OR per unit floor area [W m^-2]
        heatConsump: heating energy consumption per unit floor area [W m^-2]
        sensWaste: sensible waste heat per unit building footprint area [W m^-2]
        sensWasteCoolHeatDehum: Sensible waste heat per unit building footprint area only including cool, heat, and dehum [W m-2]
        latWaste: lat waste heat per unit building footprint area [W m^-2]
        fluxMass: mass surface heat flux [W m^-2] (mass to indoor air)
        fluxWall: wall surface heat flux [W m^-2] (wall to inside)
        fluxRoof: roof surface heat flux [W m^-2] (roof to inside)
        fluxSolar: solar heat gain per unit floor area [W m^-2] through window (SHGC)
        fluxWindow: heat gain/loss from window per unit floor area [W m^-2] (U-value)
        fluxInterior: internal heat gain adjusted for latent/LW heat per unit floor area [W m^-2]
        fluxInfil: heat flux from infiltration per unit floor area [W m^-2]
        fluxVent: heat flux from ventilation per unit floor area [W m^-2]
        ElecTotal: total electricity consumption per unit floor area [W m^-2]
        GasTotal: total gas consumption per unit floor area [W m^-2]
        Qhvac: total heat removed (sensible + latent) per unit building footprint area [W m^-2] (calculated in cooling system)
        Qheat: total heat added (sensible only) per unit building footprint area [W m^-2] (calculated in heating system)
        nFloor: Number of floors
        indoorTemp: Indoor air temperature [K]
        indoorHum: Indoor specific humidity [kg kg^-1]
        QWindowSolar: Solar Heat Gain on windows per building footprint area [W m^-2]
        QWall: Wall load per unit building footprint area [W m^-2]
        QMass: Other surfaces load per unit building footprint area [W m^-2]
        QWindow: window load due to temperature difference per unit building footprint area [W m^-2]
        QCeil: ceiling load per unit building footprint area [W m^-2]
        QInfil: infiltration load per unit building footprint area [W m^-2]
        QVen: ventilation load per unit building footprint area [W m^-2]
        QWater: energy consumption for domestic hot water [W m^-2]
        QGas: energy consumption for gas [W m^-2]
        """

        self.logger.debug("Logging at {} {}".format(__name__, self.__repr__()))

        # Building Energy Model
        self.ElecTotal = 0.0                            # total electricity consumption - (W/m^2) of floor
        self.nFloor = max(Geometry_m.Height_canyon/float(self.floorHeight),1)   # At least one floor
        self.Qheat = 0.0                                # total sensible heat added (or heating demand) per unit building footprint area [W m^-2]
        self.sensCoolDemand = 0.0                       # building sensible cooling demand per unit building footprint area [W m^-2]
        self.sensHeatDemand = 0.0                       # building sensible heating demand per unit building footprint area [W m^-2]
        self.sensWaterHeatDemand = 0.0                  # building sensible water heating demand per unit building footprint area [W m^-2]
        self.coolConsump  = 0.0                         # cooling energy consumption per unit building footprint area OR per unit floor area [W m^-2]
        self.heatConsump  = 0.0                         # heating energy consumption per unit floor area [W m^-2]
        self.sensWaste = 0.0                            # Total Sensible waste heat per unit building footprint area including cool, heat, dehum, water, and gas [W m^-2]
        self.sensWasteCoolHeatDehum = 0.0               # Sensible waste heat per unit building footprint area only including cool, heat, and dehum [W m-2]
        self.dehumDemand  = 0.0                         # Latent heat demand for dehumidification of air per unit building footprint area [W m^-2]
        self.Qhvac = 0.0                                # Total heat removed (sensible + latent)
        self.elecDomesticDemand = 0.0                   # Electricity demand for appliances and lighting (not for energy) per building footprint area [W m^-2]

        Qdehum = 0.0
        # Moist air density given dry bulb temperature, humidity ratio, and pressure [kgv m^-3]
        dens =  moist_air_density(MeteoData.Pre,self.indoorTemp,self.indoorHum)
        # evaporation efficiency in the condenser for evaporative cooling devices
        evapEff = 1.
        # total ventilation volumetric flow rate per building footprint area [m^3 s^-1 m^-2]
        volVent = self.vent * self.nFloor
        # total infiltration volumetric flow rate per building footprint area [m^3 s^-1 m^-2]
        volInfil = self.infil * Geometry_m.Height_canyon / 3600.
        # Interior wall temperature [K]
        T_wall = (BEM.wallSun.Tint+BEM.wallShade.Tint)/2
        # Solar water heating per building footprint area per hour [kg s^-1 m^-2] (Change of units [hr^-1] to [s^-1]
        massFlowRateSWH = BEM.SWH * self.nFloor/3600.
        # Interior roof temperature [K]
        T_ceil = FractionsRoof.fimp*BEM.roofImp.Tint+FractionsRoof.fveg*BEM.roofVeg.Tint
        T_mass = BEM.mass.Text              # Outer layer [K]
        T_indoor = self.indoorTemp          # Indoor temp (initial) [K]
        T_can = canTemp                     # Canyon temperature [K]

        # Normalize areas to building foot print [m^2/m^2(bld)]
        # Facade (exterior) area per unit building footprint area [m^2 m^-2]
        facArea = 2*Geometry_m.Height_canyon/numpy.sqrt(Geometry_m.Width_roof*Geometry_m.Width_roof)
        wallArea = facArea*(1.-self.glazingRatio)       # Wall area per unit building footprint area [m^2 m^-2]
        winArea = facArea*self.glazingRatio             # Window area per unit building footprint area [m^2 m^-2]
        massArea = 2*self.nFloor-1                      # ceiling and floor (top & bottom) per unit building footprint area [m^2 m^-2]
        ceilingArea = 1                                 # ceiling area per unit building footprint area [m^2 m^-2]; must be equal to 1

        # Set temperature set points according to night/day set points in building schedule & simTime; need the time in [hr]
        isEqualNightStart = self.is_near_zero((simTime.secDay/3600.) - ParCalculation.nightStart)
        if simTime.secDay/3600. < ParCalculation.nightEnd or (simTime.secDay/3600. > ParCalculation.nightStart or isEqualNightStart):
            self.logger.debug("{} Night set points @{}".format(__name__,simTime.secDay/3600.))

            # Set point temperatures in [K]
            T_cool = self.coolSetpointNight
            T_heat = self.heatSetpointNight

            # Internal heat per unit building footprint area [W m^-2]
            self.intHeat = self.intHeatNight * self.nFloor
        else:
            self.logger.debug("{} Day set points @{}".format(__name__,simTime.secDay/3600.))

            # Set point temperatures in [K]
            T_cool = self.coolSetpointDay
            T_heat = self.heatSetpointDay

            # Internal heat per unit building footprint area [W m^-2]
            self.intHeat = self.intHeatDay*self.nFloor

        # Indoor convection heat transfer coefficients
        # wall convective heat transfer coefficient [W m^-2 K^-1]
        zac_in_wall = 3.076
        # other surfaces convective heat transfer coefficient [W m^-2 K^-1]
        zac_in_mass = 0.948

        # Option 1 (zac_in_ceil): assume the same convective heat transfer coefficient regardless of temperature difference
        zac_in_ceil = 0.948

        # Option 2 (zac_in_ceil): make convective heat transfer coefficient dependent on ceiling-indoor temperature difference
        '''
        # If ceiling temperature is greater than indoor temperature use a different convective heat transfer coefficient
        if T_ceil > T_indoor:
            zac_in_ceil  = 0.948
        # If ceiling temperature is less than indoor temperature use a different convective heat transfer coefficient
        elif (T_ceil < T_indoor) or self.is_near_zero(T_ceil-T_indoor):
            zac_in_ceil  = 4.040
        else:
            print T_ceil, T_indoor
            raise Exception(self.TEMPERATURE_COEFFICIENT_CONFLICT_MSG)
            return
        '''

        # -------------------------------------------------------------
        # Heat fluxes [W m^-2]
        # -------------------------------------------------------------
        # Solar Heat Gain on windows per building footprint area [W m^-2]:
        # = radiation intensity [W m^-2] * Solar Heat Gain Coefficient (SHGC) * window area per unit building foot print area [m^2 m^-2]
        SWRinWall = (SWR.SWRin.SWRinWallSun + SWR.SWRin.SWRinWallShade)/2
        self.QWindowSolar = (SWRinWall * self.shgc * winArea)

        # QL: Latent heat per unit floor area [W m^-2] from infiltration & ventilation from
        # volInfil and volVent: volumetric rate of infiltration or ventilation per unit area [m^3 s^-1 m^-2]
        # ParCalculation.Lv: latent heat of evaporation [J kgv^-1]
        # dens: density [kga m^-3]
        # canHum: canyon specific humidity [kgv kga^-1]
        # indoorHum: indoor specific humidity [kgv kga^-1]
        # Note: at the moment the infiltration and system specific humidity are considered to be the same
        # This is a serious limitation.
        # Future versions of the UWG must calculate the system specific humidity based on HVAC system parameters

        # Latent heat per building footprint area [W m^-2]
        QLinfil = volInfil * dens * ParCalculation.Lv * (canHum - self.indoorHum)
        QLvent = volVent * dens * ParCalculation.Lv * (canHum - self.indoorHum)

        # Latent heat load per unit building footprint area [W m^-2]
        QLintload = self.intHeat * self.intHeatFLat

        # Note: at the moment the infiltration and system air temperatures are considered to be the same
        # This is a serious limitation.
        # Future versions of UWG must calculate the system temperature based on HVAC system parameters
        # wall load per unit building footprint area [W m^-2]
        self.QWall = wallArea * zac_in_wall * (T_wall - T_cool)
        # other surfaces load per unit building footprint area [W m^-2]
        self.QMass = massArea * zac_in_mass * (T_mass - T_cool)
        # window load due to temperature difference per unit building footprint area [W m^-2]
        self.QWindow = winArea * self.uValue * (T_can - T_cool)
        # ceiling load per unit building footprint area [W m^-2]
        self.QCeil = ceilingArea * zac_in_ceil * (T_ceil - T_cool)
        # infiltration load per unit building footprint area [W m^-2]
        self.QInfil = volInfil * dens * ParCalculation.cp_atm * (T_can - T_cool)
        # ventilation load per unit building footprint area [W m^-2]
        self.QVen = volVent * dens * ParCalculation.cp_atm * (T_can - T_cool)
        # Heat/Cooling load per unit building footprint area [W m^-2], if any
        self.sensCoolDemand = max(self.QWall+self.QMass+self.QWindow+self.QCeil+self.intHeat+self.QInfil+self.QVen+self.QWindowSolar,0.)

        self.sensHeatDemand = max(
            -(wallArea*zac_in_wall*(T_wall-T_heat) +               # wall load per unit building footprint area [W m^-2]
            massArea*zac_in_mass*(T_mass-T_heat) +                 # other surfaces load per unit building footprint area [W m^-2]
            winArea*self.uValue*(T_can-T_heat) +                   # window load due to temperature difference per unit building footprint area [W m^-2]
            ceilingArea*zac_in_ceil*(T_ceil-T_heat) +              # ceiling load per unit building footprint area [W m^-2]
            self.intHeat +                                         # internal load per unit building footprint area [W m^-2]
            volInfil*dens*ParCalculation.cp_atm*(T_can-T_heat) +   # infiltration load per unit building footprint area [W m^-2]
            volVent*dens*ParCalculation.cp_atm*(T_can-T_heat) +    # ventilation load per unit building footprint area [W m^-2]
            self.QWindowSolar),                                    # solar load through window per unit building footprint area [W m^-2]
            0.)

        # -------------------------------------------------------------
        # HVAC system (cooling demand = [W m^-2] bld footprint)
        # -------------------------------------------------------------
        # If the canyon air temperature is greater than 288 K building energy system is under cooling mode
        if self.sensCoolDemand > 0. and canTemp > 288.:
            # Option 1 (dehumDemand): consider a dehumidifiction scenario
            '''
            # Energy is used to dehumidify volumetric flow rate of air per unit building footprint area that flows through dehumidifier,
            # Calculate an arbitrary Volumetric Flow rate based on sensible Cooling demand
            # Assume air is cooled to 10C and fraction f = 0.02 is dehumidified
            # equal to f * sensCoolDemand / (dens * Cp * (T_indoor - (273.15+10)))
            # Fraction of volumetric flow rate of air to cool per unit building footprint area [m^3 s^-1 m^-2]
            # VolDehum = 0.02 * self.sensCoolDemand / (dens*ParCalculation.cp_atm*(T_indoor - (273.15+10)))
            # Energy is used to dehumidify this volumetric flow rate of air per unit building footprint area,
            # Assume air is cooled to 10C and conditioned to 70% RH (5.5 [gv kg^-1])
            # This energy is equal to VolDehum * dens * (self.indoorHum - 0.0055)*ParCalculation.Lv
            # Latent heat demand for dehumidification of air per unit building footprint area [W m^-2]
            # self.dehumDemand = max(VolDehum * dens * (self.indoorHum - 0.0055)*ParCalculation.Lv, 0.)
            '''
            # Option 2 (dehumDemand): set dehumidification demand as a fraction of sensible cooling demenad
            self.dehumDemand = 0.1 * self.sensCoolDemand

            # Save dehumidification demand for later use because it may be modified
            Qdehum = copy.copy(self.dehumDemand)

            # Calculate total cooling demand in per unit building footprint area [W m^-2]
            # if cooling energy demand is greater then HVAC cooling capacity
            if (self.dehumDemand + self.sensCoolDemand) > (self.coolCap * self.nFloor):
                self.Qhvac = self.coolCap * self.nFloor
                # Part load ratio
                PLR = (self.coolCap * self.nFloor) / (self.dehumDemand + self.sensCoolDemand)
                # For option 1 above, we need to discount VolDehum
                # VolDehum = VolDehum * PLR
                self.sensCoolDemand = self.sensCoolDemand * PLR
                self.dehumDemand = self.dehumDemand * PLR
            else:
                self.Qhvac = self.dehumDemand + self.sensCoolDemand

            # Calculate input work required by the refrigeration cycle per unit building footprint area [W m^-2]
            # COP = QL/Win or Win = QL/COP
            self.coolConsump = (max(self.sensCoolDemand+self.dehumDemand,0.0))/self.copAdj

            # Calculate waste heat from HVAC system per unit building footprint area [W m^-2]
            # Using 1st law of thermodynamics QH = Win + QL
            if (self.condType == 'AIR'):
                self.sensWasteCoolHeatDehum = max(self.sensCoolDemand+self.dehumDemand,0)+self.coolConsump
                self.latWaste = 0.0
            # We have not tested this option; it must be investigated further
            elif (self.condType == 'WAT'):
                self.sensWasteCoolHeatDehum = max(self.sensCoolDemand+self.dehumDemand,0)+self.coolConsump*(1.-evapEff)
                self.latWaste = max(self.sensCoolDemand+self.dehumDemand,0)+self.coolConsump*evapEff

            self.sensHeatDemand = 0.

        # -------------------------------------------------------------
        # HVAC system (heating demand = [W m^-2] bld footprint)
        # -------------------------------------------------------------
        # If the canyon air temperature is less than 288 K building energy system is under heating mode
        # Under heating mode, there is no dehumidification
        elif self.sensHeatDemand > 0. and canTemp < 288.:
            # Calculate total heating demand in per unit building footprint area [W m^-2]
            # Heating demand must be less than or equal to heating capacity
            self.Qheat = min(self.sensHeatDemand, self.heatCap*self.nFloor)
            # Calculate the energy consumption of the heating system per unit building footprint area [W m^-2] from heating demand divided by efficiency
            self.heatConsump  = self.Qheat / self.heatEff
            # Calculate waste heat from HVAC system per unit building footprint area [W m^-2]
            # Using 1st law of thermodynamics QL = Win - QH
            self.sensWasteCoolHeatDehum = self.heatConsump - self.Qheat
            # The heating system model assumes that the indoor air humidity is not controlled
            Qdehum = 0.0
            self.sensCoolDemand = 0.0


        # -------------------------------------------------------------
        # Evolution of the internal temperature and humidity
        # -------------------------------------------------------------
        # Solve sensible heat balance equation for indoor air, considering effect of heat fluxes from wall, mass, roof,
        # window, solar heat gain on windows, internal heat, infiltration, ventilation and HVAC (cooling or heating)
        # per unit building footprint area [W m^-2]
        # Explicit terms in eq. 2 which either do not contain Tin or contain Tin from previous iteration (Bueno et al., 2012)
        Q = self.intHeat + self.QWindowSolar + self.Qheat - self.sensCoolDemand

        H1 = (T_wall*wallArea*zac_in_wall +
            T_mass*massArea*zac_in_mass +
            T_ceil*zac_in_ceil +
            T_can*winArea*self.uValue +
            T_can*volInfil * dens * ParCalculation.cp_atm +
            T_can*volVent * dens * ParCalculation.cp_atm)
        # Implicit terms in eq. 2 which directly contain coefficient for newest Tin to be solved (Bueno et al., 2012)
        H2 = (wallArea*zac_in_wall +
            massArea*zac_in_mass +
            ceilingArea*zac_in_ceil +
            winArea*self.uValue +
            volInfil * dens * ParCalculation.cp_atm +
            volVent * dens * ParCalculation.cp_atm)

        # Assumes air temperature of control volume is sum of surface boundary temperatures
        # weighted by area and heat transfer coefficient + generated heat
        # Calculate indoor air temperature [K]
        self.indoorTemp = (H1 + Q)/H2
        # Solve the latent heat balance equation for indoor air, considering effect of internal, infiltration and
        # ventilation latent heat and latent heat demand for dehumidification per unit building footprint area [W m^-2];
        # eq. 3 (Bueno et al., 2012)
        # Calculate indoor specific humidity [kgv kga^-1]
        self.indoorHum = self.indoorHum + (simTime.dt/(dens * ParCalculation.Lv * Geometry_m.Height_canyon)) * \
            (QLintload + QLinfil + QLvent - Qdehum)

        # Calculate relative humidity ((Pw/Pws)*100) using pressure, indoor temperature, humidity
        _Tdb, _w, _phi, _h, _Tdp, _v = psychrometrics(self.indoorTemp, self.indoorHum, MeteoData.Pre)
        # Indoor relative humidity
        self.indoorRhum = _phi

        # Heat fluxes of elements [W m^-2]
        # (will be used for element calculation)
        # Wall heat flux per unit wall area [W m^-2]
        self.fluxWall = zac_in_wall * (T_indoor - T_wall)
        # Top ceiling heat flux per unit ceiling or building footprint area [W m^-2]
        self.fluxRoof = zac_in_ceil * (T_indoor - T_ceil)
        # Inner horizontal heat flux per unit floor area [W m^-2]
        self.fluxMass = zac_in_mass * (T_indoor - T_mass) + self.intHeat * self.intHeatFRad/massArea


        # Calculate heat fluxes per unit floor area [W m^-2] (These are for record keeping only)
        self.fluxSolar = self.QWindowSolar/self.nFloor
        self.fluxWindow = winArea * self.uValue *(T_can - T_indoor)/self.nFloor
        self.fluxInterior = self.intHeat * self.intHeatFRad *(1.-self.intHeatFLat)/self.nFloor
        self.fluxInfil= volInfil * dens * ParCalculation.cp_atm *(T_can - T_indoor)/self.nFloor
        self.fluxVent = volVent * dens * ParCalculation.cp_atm *(T_can - T_indoor)/self.nFloor

        # Total Electricity consumption per unit floor area [W m^-2] which is equal to
        # cooling consumption + electricity consumption + lighting
        self.ElecTotal = self.coolConsump/self.nFloor + BEM.Elec + BEM.Light
        # electricity demand other than cooling consumption per building footprint area [W m^-2]
        self.elecDomesticDemand = self.nFloor * (BEM.Elec + BEM.Light)
        # Sensible hot water heating demand
        CpH20 = 4200.           # heat capacity of water [J Kg^-1 K^-1]
        T_hot = 49 + 273.15     # Service water temp (assume no storage) [K]
        self.sensWaterHeatDemand = massFlowRateSWH * CpH20 * (T_hot - MeteoData.waterTemp)

        # Calculate total sensible waste heat to canyon per unit building footprint area [W m^-2]
        # which can be determined from sensible waste to canyon, energy consumption for domestic hot water and gas consumption
        # Sensible hot water heating demand
        self.sensWaterHeatDemand = massFlowRateSWH * CpH20 * (T_hot - MeteoData.waterTemp)
        # Waste heat of water heating
        self.QWater = (1 / self.heatEff - 1.) * self.sensWaterHeatDemand
        self.QGas = BEM.Gas * (1 - self.heatEff) * self.nFloor
        self.sensWaste = self.sensWasteCoolHeatDehum + self.QWater + self.QGas
        # Calculate total gas consumption per unit floor area [W m^-2] which is equal to gas consumption per unit floor area +
        # energy consumption for domestic hot water per unit floor area + energy consumption of the heating system per unit floor area
        self.GasTotal = BEM.Gas + (massFlowRateSWH*CpH20*(T_hot - MeteoData.waterTemp)/self.nFloor)/self.heatEff + self.heatConsump/self.nFloor
