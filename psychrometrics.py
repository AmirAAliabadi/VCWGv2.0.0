from math import log, pow, exp

"""
Calculate dew point temperature, saturation pressure, specific humidity
Developed by Bruno Bueno
Building Technology Lab; Massachusetts Institute of Technology, Cambridge, USA
Last update: March 2012
"""

# Modified version of Psychometrics by Tea Zakula
# MIT Building Technology Lab
# Input: Tdb_in, w_in, P
# Output: Tdb, w, phi, h, Tdp, v
#
# where:
# Tdb_in = dry bulb temperature [K]
# w_in = Humidity Ratio also known as specific humidity (q) [kgv kgda^-1]
# P = Atmospheric Station Pressure [Pa]
#
# Tdb: dry bulb temperature [C]
# w: Humidity Ratio also known as specific humidity (q) [kgv kgda^-1]
# phi: relative humidity [(Pw/Pws)*100]
# Tdp: dew point temperature [C]
# h: enthalpy [J kga^-1]
# v: specific volume also equal to inverse of density [m^3 kga^-1]

def psychrometrics (Tdb_in, w_in, P):

    """
    ------
    INPUT:
    Tdb_in: Dry bulb temperature [K]
    w_in: Specific humidity [kg kg^-1]
    P: Pressure [Pa]
    -------
    OUTPUT:
    Tdb: Dry bulb temperature [C]
    w: Specific humidity [kg kg^-1]
    phi: Relative humidity [%]
    h: Enthalpy [J kga^-1]
    Tdp: Dew point temperature [C]
    v: Specific volume calculation [m^3 kga^-1]
    """

    # Change units
    c_air = 1006.   # air heat capacity, value from ASHRAE Fundamentals [J kg^-1]
    hlg = 2501000.  # latent heat, value from ASHRAE Fundamentals [J kg^-1]
    cw  = 1860.     # water vapor heat capacity, value from ASHRAE Fundamentals [J kg^-1]
    P = P/1000.     # convert from [Pa] to [kPa]

    # Dry bulb temperature [C]
    Tdb = Tdb_in - 273.15
    w = w_in


    # phi (RH) calculation from Tdb and w
    Pw = (w*P)/(0.621945 + w)                             # partial pressure of water vapor [kPa]
    Pws = saturation_pressure(Tdb)                        # Get saturation pressure for given Tdb [kPa]

    phi = (Pw/Pws)*100.0

    # enthalpy calculation from Tdb and w [J kga^-1]
    h = c_air*Tdb + w*(hlg+cw*Tdb)

    # specific volume calculation from Tdb and w [m^3 kga^-1]
    v = 0.287042 * (Tdb+273.15)*(1+1.607858*w)/P

    # dew point calculation from w
    # water vapor partial pressure in [kPa]
    _pw = (w*P)/(0.621945 + w)
    alpha = log(_pw)

    # Dew point temperature valid for Tdp between 0 C and 93 C [C]
    Tdp = 6.54 + 14.526*alpha + pow(alpha,2)*0.7389 + pow(alpha,3)*0.09486 + pow(_pw,0.1984)*0.4569

    return  Tdb, w, phi, h, Tdp, v

def saturation_pressure(Tdb_):

    """
    ------
    INPUT:
    Tdb_: Dry bulb temperature [C]
    -------
    OUTPUT:
    _Pws: Saturation pressure [kPa]
    """

    # Temperature [K]
    T = Tdb_ + 273.15

    # N.B In Matlab, negative values are converted to complex values.
    # log(-x) = log(x) + log(-1) = log(x) + i*pi
    # Python will throw an exception. Negative value occurs here if
    # simulation timestep (dtSim) is large, i.e 3600s.
    # Saturation pressure [Pa]
    _Pws = exp(-1*(5.8002206e3) / T+1.3914993 + (4.8640239e-2)*T*(-1.) + (4.1764768e-5)*pow(T,2) - (1.4452093e-8)*pow(T,3) + 6.5459673*log(T))
    # Saturation pressure [kPa]
    _Pws = _Pws/1000.
    return _Pws

def moist_air_density(P,Tdb,H):
    # Moist air density [kgv m^-3] given dry bulb temperature, humidity ratio, and pressure.
    # ASHRAE Fundamentals (2005) ch. 6 eqn. 28
    # ASHRAE Fundamentals (2009) ch. 1 eqn. 28
    # from: https://github.com/psychrometrics/Libraries/blob/master/Psychrometrics_SI.cpp
    moist_air_density = P/(1000*0.287042*Tdb*(1.+1.607858*H))
    return moist_air_density

def HumFromRHumTemp(RH,T,P):

    """
    ------
    INPUT:
    RH: Relative humidity [%]
    T: Temperature [C]
    P: Pressure [Pa]
    -------
    OUTPUT:
    W: Specific humidity [kg kg^-1]
    """

    # Derive Specific Humidity [kgh20 kgn202^-1] from RH, T and Pa
    # Saturation vapour pressure from ASHRAE
    C8 = -5.8002206e3
    C9 = 1.3914993
    C10 = -4.8640239e-2
    C11 = 4.1764768e-5
    C12 = -1.4452093e-8
    C13 = 6.5459673

    # Convert temperature from Celcius [C] to Kelvin [K]
    T += 273.15

    PWS = exp(C8/T + C9 + C10*T + C11 * pow(T,2) + C12 * pow(T,3) + C13 * log(T))
    PW = RH*PWS/100.0        # Vapour pressure
    W = 0.62198*PW/(P-PW)    # Specific humidity
    return W

"""
function psat = psat(temp,parameter)
    gamw  = (parameter.cl - parameter.cpv) / parameter.rv;
    betaw = (parameter.lvtt/parameter.rv) + (gamw * parameter.tt);
    alpw = log(parameter.estt) + (betaw /parameter.tt) + (gamw *log(parameter.tt));
    psat = zeros(size(temp));
    for jj=1:size(temp)
        psat = exp(alpw - betaw/temp - gamw*log(temp));
    end
end

% Not used for this release but saved for possible future use
function Twb = wet_bulb(Tdb,Tdp,pres)

    % Copyright (c) 2015, Rolf Henry Goodwin
    % All rights reserved.
    %
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are
    % met:
    %
    %     * Redistributions of source code must retain the above copyright
    %       notice, this list of conditions and the following disclaimer.
    %     * Redistributions in binary form must reproduce the above copyright
    %       notice, this list of conditions and the following disclaimer in
    %       the documentation and/or other materials provided with the distribution
    %
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    % POSSIBILITY OF SUCH DAMAGE.

    % Code modified to merge into a single file - Joseph Yang, 2016


    % Tdb, Tdp, Twb in K
    % p in Pa (obtained function uses hPa, so /100 needed)
    global T;
    global T_d;
    global p;
    T = Tdb;
    T_d = Tdp;
    p = pres/100;

    Twb = root_finder(@Delta_q,T_d,T);
end

function dQTw = Delta_q(T_w)
    %Delta_q finds the value of function dq(Tw)
    %INPUT wet bulb temperature T_w
    %OUTPUT dq(Tw)
    global T;
    global T_d;
    global p;

    Cp = 1005; % Heat capacity of water vapor in J/(kg*K)
    L = 2.501e6; % Latent heat of water vapor at 0 degC in J/kg
    w1 = mixing_ratio(T_d,p); % Mixing ratio corresponding to T_d and p
    w2 = mixing_ratio(T_w,p); % Mixing ratio corresponding to T_w and p

    dQTw = (L*(w2-w1))/(1+w2)-Cp*(T-T_w)*(1+0.8*w2); % Finds deltaq(Tw)

end
"""
