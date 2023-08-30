import numpy
import math
from BuildingColumnModel import BuildingCol
from DragTurb import TurbCoeff
from shear import ShearProd
from Buoyancy import BuoProd
from NumericalSolver import Diff
from Invert import Invert
import copy

"""
Column Model for momentum, turbulent kinetic energy, temperature, and specific humidity in the urban environment
Developed by Mohsen Moradi and Amir A. Aliabadi
Atmospheric Innovations Research (AIR) Laboratory, University of Guelph, Guelph, Canada
Last update: May 2021
Originally developed by Alberto Martilli, Scott Krayenhoff, and Negin Nazarian
"""

def ColumnModelCal(z0_road,z0_roof,Ceps,Cdrag,Ck,thb,qhb,tvb,FractionsGround,FractionsRoof,TemperatureC,TemperatureR,
                   ForcingVariable,VerticalProfUrban,Geometry_m,geometry,ColParam,f_LAD,dlk,dls,pb,ss,sf,vol,Cp,Ustar,
                   SensHt_HVAC,dts,Rural_Model_name):
    """
    ------
    INPUT:
    z0_road: Road aerodynamic roughness length [m]
    z0_roof: Roof aerodynamic roughness length [m]
    Ceps: Coefficient for the destruction of turbulent dissipation rate [-]
    Cdrag: Drag coefficient due to buildings (sectional drag coefficient) [-]
    Ck: Coefficient used in the equation of diffusion coefficient [-]
    thb: Sink/source terms in temperature equation in 1-D model caused by roof and ground
    qhb: Sink/source terms in specific humidity equation in 1-D model
    tvb: Sink/source terms in temperature equation in 1-D model caused by walls
    FractionsGround: Fractions of ground covered by vegetation, impervious,and bare surfaces
    FractionsRoof: Fractions of roof covered by vegetation and impervious surface
    TemperatureC: Temperature of the exterior surfaces of urban canyon [K]
    TemperatureR: Temperature of the exterior surfaces of roof [K]
    ForcingVariable: Forcing variable at the top of the domain
    VerticalProfUrban: Vertical profile of temperature, humidity, wind speed, and turbulent kinetic energy from previous time step
    Geometry_m: Geometric parameters
    geometry: Normalized geometric parameters
    ColParam: Column model parameters
    f_LAD: Interpolate LAD profile [m^2 m^-3]
    dlk: Turbulent length scale [m]
    dls: Dissipation length scale [m]
    pb: Probability that a building has a height greater or equal to z [-]
    ss: Probability that a building has a height equal to z [-]
    sf: Fraction of air at the interface between cells [-]
    vol: Volume fraction of air in each urban unit cell [-]
    Cp: Air Specific heat [J kg^-1 K^-1]
    Ustar: Friction velocity in the rural site [m s^-1]
    SensHt_HVAC: Total waste heat from buildings [W m^-2]
    dts: time step [s]
    Rural_Model_name: Rural model name
    -------
    OUTPUT:
    vx: Wind speed profile in x direction [m s^-1]
    vy: Wind speed profile in y direction [m s^-1]
    tke: Turbulent kinetic energy profile [m^2 s^-2]
    th: Potential temperature profile [K]
    qn: Specific humidity profile [kg kg^-1]
    presProf: Pressure profile [Pa]
    th_eq: Terms in temperature equation
    vx_eq: Terms in wind speed in x direction equation
    vy_eq: Terms in wind speed in y direction equation
    tke_eq: Terms in turbulent kinetic energy equation
    qn_eq: Terms in specific humidity equation
    ColumnModel_OtherParam:
    """

    TGround = FractionsGround.fimp*TemperatureC[0] + FractionsGround.fbare*TemperatureC[1] + FractionsGround.fveg*TemperatureC[2]
    TWallSun = TemperatureC[3]
    TWallShade = TemperatureC[4]
    TRoof = FractionsRoof.fimp*TemperatureR[0] + FractionsRoof.fveg*TemperatureR[1]

    vx = copy.copy(VerticalProfUrban.vx)             # x component of horizontal wind speed [m s^-1]
    vy = copy.copy(VerticalProfUrban.vy)             # y component of horizontal wind speed [m s^-1]
    tke = copy.copy(VerticalProfUrban.tke)           # Turbulent kinetic energy [m^2 s^-2]
    th = copy.copy(VerticalProfUrban.th)             # Potential temperature [K]
    qn = copy.copy(VerticalProfUrban.qn)             # Specific humidity [kgv kga^-1]
    presProf = copy.copy(VerticalProfUrban.presProf) # Pressure profile [Pa]
    rho = copy.copy(VerticalProfUrban.rho)
    th_ref = copy.copy(VerticalProfUrban.th_ref)     # Reference potential temperature [K]

    # gravitational acceleration [m s^-2]
    g = 9.81

    # Define explicit and implicit parts of source and sink terms
    srex_vx = numpy.zeros(Geometry_m.nz)       # Explicit part of x component of horizontal wind speed [m s^-2]
    srim_vx = numpy.zeros(Geometry_m.nz)       # Implicit part of x component of horizontal wind speed [s^-1]
    srex_vy = numpy.zeros(Geometry_m.nz)       # Explicit part of y component of horizontal wind speed [m s^-2]
    srim_vy = numpy.zeros(Geometry_m.nz)       # Implicit part of y component of horizontal wind speed [s^-1]
    srex_tke = numpy.zeros(Geometry_m.nz)      # Explicit part of turbulent kinetic energy [m^2 s^-3]
    srim_tke = numpy.zeros(Geometry_m.nz)      # Implicit part of turbulent kinetic energy [s^-1]
    srex_th = numpy.zeros(Geometry_m.nz)       # Explicit part of potential temperature [K s^-1]
    srim_th = numpy.zeros(Geometry_m.nz)       # Implicit part of potential temperature [s^-1]
    srex_qn = numpy.zeros(Geometry_m.nz)       # Explicit part of specific humidity [kg kg^-1 s^-1]
    srim_qn = numpy.zeros(Geometry_m.nz)       # Implicit part of specific humidity [s^-1]
    dissip_tke_veg = numpy.zeros(Geometry_m.nz)# Dissipation of tke caused by vegetation [s^-1]
    srim_vxy_veg = numpy.zeros(Geometry_m.nz)  # Implicit term in x/y momentum equation due to vegetation [s^-1]
    srex_tke_veg = numpy.zeros(Geometry_m.nz)  # Implicit term in turbulent kinetic energy equation due to vegetation [m^2 s^-3]
    srex_th_veg = numpy.zeros(Geometry_m.nz)   # Explicit part of potential temperature caused by vegetation [K s^-1]
    srex_qn_veg = numpy.zeros(Geometry_m.nz)   # Explicit part of specific humidity caused by vegetation [kg kg^-1 s^-1]

    #---------------------------
    # Define Boundary Conditions
    #---------------------------
    # Define boundary conditions for wind speed (1:Neumann boundary condition (Flux), 2:Dirichlet boundary condition (Constant value))
    if Rural_Model_name == 'Forcing_extFile':
        Wind_bc_bottom = 1
        Wind_bc_top = 2
    elif Rural_Model_name == 'MOST':
        Wind_bc_bottom = 1
        Wind_bc_top = 1

    kappa = 0.4

    if Wind_bc_top == 2 and Wind_bc_bottom == 1:
        # Use externally forced wind component on the top of the urban column
        vx[Geometry_m.nz - 1] = ForcingVariable[1] * math.cos(math.radians(ForcingVariable[2]))
        vy[Geometry_m.nz - 1] = ForcingVariable[1] * math.sin(math.radians(ForcingVariable[2]))

        dpdx = 0
        dpdy = 0

    elif Wind_bc_top == 1 and Wind_bc_bottom == 1:
        # Calculate pressure gradient [kg m^-2 s^-2]
        dpdx = rho[Geometry_m.nz-1]*(Ustar**2)*(math.cos(math.radians(ForcingVariable[2])))/(Geometry_m.dz*Geometry_m.nz)
        dpdy = rho[Geometry_m.nz-1]*(Ustar**2)*(math.sin(math.radians(ForcingVariable[2])))/(Geometry_m.dz*Geometry_m.nz)


    # Define boundary conditions for temperature (1:Neumann boundary condition (Flux), 2:Dirichlet boundary condition (Constant value))
    T_bc_bottom = 1
    T_bc_top = 2
    th[Geometry_m.nz - 1] = ForcingVariable[0]

    # Define boundary conditions for humidity (1:Neumann boundary condition (Flux), 2:Dirichlet boundary condition (Constant value))
    q_bc_bottom = 1
    q_bc_top = 2
    qn[Geometry_m.nz - 1] = ForcingVariable[3]

    # Define boundary conditions for tke (1:Neumann boundary condition (Flux), 2:Dirichlet boundary condition (Constant value))
    tke_bc_bottom = 1
    tke_bc_top = 1

    # -----------------------------------------
    # Calculate Turbulent Diffusion Coefficient
    # -----------------------------------------
    # Calculate turbulent diffusion coefficient (Km) [m^2 s^-1]
    Km = TurbCoeff(Geometry_m.nz, Ck, tke, dlk)

    # Calculate shear production [m^2 s^-3] in TKE equation. (Term II of equation 5.2, Krayenhoff 2014, PhD thesis)
    sh = ShearProd(ColParam.cdmin,Geometry_m.nz, Geometry_m.dz, vx, vy, Km)

    # Calculate buoyant production [m^2 s^-3] in TKE equation. (Term IX of equation 5.2, Krayenhoff 2014, PhD thesis)
    bu = BuoProd(ColParam.cdmin,Geometry_m.nz, Geometry_m.dz, th, Km, th_ref, ColParam.prandtl)

    # Calculate dissipation (td) [s^-1] in TKE equation. (Term VI of equation 5.2, Krayenhoff 2014, PhD thesis)
    # parameterization of dissipation is based on Nazarian's code. (https://github.com/nenazarian/MLUCM/blob/master/Column_Model/column_lkPro.f90)
    td = numpy.zeros(Geometry_m.nz)
    for i in range(0, Geometry_m.nz):
        if dls[i] != 0:
            td[i] = -Ceps*(math.sqrt(tke[i]))/dls[i]
        else:
            td[i] = 0
        sh[i] = sh[i]*sf[i]
        bu[i] = bu[i]*sf[i]

    # Calculate sink and source terms in momentum, temperature and turbulent kinetic energy (TKE) equations which are caused by building
    BuildingCoef = BuildingCol(Geometry_m.nz, Geometry_m.dz, dts, vol,Geometry_m.lambdap, Geometry_m.lambdaf,
                               Geometry_m.Height_canyon, Ck, Cp, th_ref, vx, vy, th, Cdrag, rho,Geometry_m.nz_u, pb, ss,
                               z0_road,z0_roof,SensHt_HVAC,ColParam.HVAC_street_frac,ColParam.HVAC_atm_frac,ColParam.WindMin_Urban)
    BuildingCoef.BuildingDrag_UTC(thb,qhb,tvb,FractionsGround,FractionsRoof,TWallSun,TWallShade,TGround,TRoof)


    # Drag coefficient for vegetation foliage
    cdv = 0.2

    # Calculate source and sink terms caused by trees and then calculate total source and sink terms
    for i in range(0, Geometry_m.nz):

        # source/sink terms of specific humidity
        wind = numpy.sqrt(vx[i] ** 2 + vy[i] ** 2)

        # Calculate terms in transport equations caused by trees. "wt" is term in temperature equation and "wt_drag"
        # is term in TKE and momentum equations. It is assumed there is no vegetation above average building height
        if Geometry_m.dz*i+Geometry_m.dz/2 > max(ColParam.h_LAD):
            wt = 0         # [m^2 m^-3]
            wt_drag = 0    # [m^2 m^-3]
        else:
            wt = f_LAD(Geometry_m.dz*i+Geometry_m.dz/2) * ColParam.omega * (1 - Geometry_m.lambdap) / vol[i]             # [m^2 m^-3]
            wt_drag = f_LAD(Geometry_m.dz*i+Geometry_m.dz/2) * ColParam.omega_drag * (1. - Geometry_m.lambdap) / vol[i]  # [m^2 m^-3]

        # Calculate explicit terms in temperature and humidity equations caused by trees.
        srex_th_veg[i] = wt*thb.tree*2*2*geometry.radius_tree # [K s^-1]
        srex_qn_veg[i] = wt*qhb.tree*2*2*geometry.radius_tree # [kg kg^-1 s^-1]

        # Calculate total explicit terms
        # Explicit term in x momentum equation [m s^-2] = fluxes from horizontal surfaces + pressure gradient
        srex_vx[i] = BuildingCoef.srex_vx_h[i]+dpdx

        # Explicit term in y momentum equation [m s^-2] = fluxes from horizontal surfaces + pressure gradient
        srex_vy[i] = BuildingCoef.srex_vy_h[i]+dpdy

        # Explicit term in TKE equation [m^2 s^-3] = terms from urban horizontal surfaces +
        # terms from walls [m^2 s^-3] + shear production [m^2 s^-3] + buoyant production [m^2 s^-3] +
        # term caused by vegetation [m^2 s^-3]
        srex_tke_veg[i] = cdv*wind**3.*wt_drag
        srex_tke[i] = BuildingCoef.srex_tke_h[i] + BuildingCoef.srex_tke_v[i] + sh[i] + bu[i] + srex_tke_veg[i]

        # Explicit term in temperature equation [K s^-1] = term from urban horizontal surfaces [K s^-1] +
        # term from walls [K s^-1] + term caused by vegetation [K s^-1]
        srex_th[i] = BuildingCoef.srex_th_h[i] + BuildingCoef.srex_th_v[i] + srex_th_veg[i]

        # Explicit term in humidity equation [K s^-1] = term caused by latent heat from ground and vegetation [K s^-1]
        srex_qn[i] = BuildingCoef.srex_qn_h[i] + srex_qn[i] + srex_qn_veg[i]

        # Calculate total Implicit terms
        # Implicit term in x momentum equation [s^-1] = term from walls [s^-1] - term caused by vegetation [s^-1]
        srim_vxy_veg[i] = cdv*wind*wt_drag
        srim_vx[i] = BuildingCoef.srim_vx_v[i]-srim_vxy_veg[i]

        # Implicit term in y momentum equation [s^-1] = term from walls [s^-1] - term caused by vegetation [s^-1]
        srim_vy[i] = BuildingCoef.srim_vy_v[i]-srim_vxy_veg[i]

        # Implicit term in TKE equation [s^-1] = dissipation [s^-1] - term caused by vegetation [s^-1]
        dissip_tke_veg[i] = 6.5*cdv*wind*wt_drag
        srim_tke[i] = td[i]-dissip_tke_veg[i]

        # Implicit term in temperature equation [s^-1] = term from wall [s^-1] - term caused by vegetation [s^-1]
        srim_th[i] = BuildingCoef.srim_th_v[i]

        # Implicit term in humidity equation [s^-1] = term caused by latent heat from vegetation [s^-1]
        srim_qn[i] = srim_qn[i]


    # Solve transport equations
    Sol = Diff(Geometry_m.nz, dts, sf, vol, Geometry_m.dz, rho)
    # Solve x component of momentum equation
    vx_new,uw,duwdz = Sol.Solver(Geometry_m.nz,Geometry_m.nz,Wind_bc_bottom,Wind_bc_top,dts,rho,vx,Km,srim_vx,srex_vx,sf,vol,Geometry_m.dz)
    # Solve y component of momentum equation
    vy_new,vw,dvwdz = Sol.Solver(Geometry_m.nz,Geometry_m.nz,Wind_bc_bottom,Wind_bc_top,dts,rho,vy,Km,srim_vy,srex_vy,sf,vol,Geometry_m.dz)
    # Solve TKE equation
    tke_new,wtke,dwtkedz = Sol.Solver(Geometry_m.nz,Geometry_m.nz,tke_bc_bottom,tke_bc_top,dts,rho,tke,Km,srim_tke,srex_tke,sf,vol,Geometry_m.dz)
    # Solve temperature equation
    th_new,wth,dwthdz = Sol.Solver(Geometry_m.nz,Geometry_m.nz,T_bc_bottom,T_bc_top,dts,rho,th,Km/ColParam.prandtl,srim_th,srex_th,sf,vol,Geometry_m.dz)
    # Solve specific humidity equation
    qn_new,wqn,dwqndz = Sol.Solver(Geometry_m.nz,Geometry_m.nz,q_bc_bottom,q_bc_top,dts,rho,qn,Km/ColParam.schmidt,srim_qn,srex_qn,sf,vol,Geometry_m.dz)

    vx = copy.copy(vx_new)
    vy = copy.copy(vy_new)
    th = copy.copy(th_new)
    tke = copy.copy(tke_new)
    qn = copy.copy(qn_new)
    # Set a minimum value for kinetic energy which avoid trapping of heat at street level
    for i in range(0, Geometry_m.nz):
        if tke[i] < 1e-3:
           tke[i] = 1e-3

    # Calculate pressure profile [Pa]
    for iz in reversed(list(range(Geometry_m.nz))[1:]):
        presProf[iz-1] = (math.pow(presProf[iz], 287 / Cp) +
                                 9.81 / Cp * (math.pow(presProf[0], 287 / Cp)) *
                                 (1. / th[iz] + 1. / th[iz - 1]) * 0.5 * Geometry_m.dz) ** (1. / (287 / Cp))

    class th_eq_Def():
        pass
    th_eq = th_eq_Def()
    th_eq.srex = srex_th
    th_eq.srex_h = BuildingCoef.srex_th_h
    th_eq.srex_v = BuildingCoef.srex_th_v
    th_eq.srex_veg = srex_th_veg
    th_eq.srim = srim_th
    th_eq.srim_v = BuildingCoef.srim_th_v
    th_eq.diffusion = Km/ColParam.prandtl

    class vx_eq_Def():
        pass
    vx_eq = vx_eq_Def()
    vx_eq.srex = srex_vx
    vx_eq.srex_h = BuildingCoef.srex_vx_h
    vx_eq.srex_dpdx = [dpdx for i in range(Geometry_m.nz)]
    vx_eq.srim = srim_vx
    vx_eq.srim_v = BuildingCoef.srim_vx_v
    vx_eq.srim_veg = srim_vxy_veg
    vx_eq.diffusion = Km

    class vy_eq_Def():
        pass
    vy_eq = vy_eq_Def()
    vy_eq.srex = srex_vy
    vy_eq.srex_h = BuildingCoef.srex_vy_h
    vy_eq.srex_dpdx = [dpdy for i in range(Geometry_m.nz)]
    vy_eq.srim = srim_vy
    vy_eq.srim_v = BuildingCoef.srim_vy_v
    vy_eq.srim_veg = srim_vxy_veg
    vy_eq.diffusion = Km

    class tke_eq_Def():
        pass
    tke_eq = tke_eq_Def()
    tke_eq.srex = srex_tke
    tke_eq.srex_h = BuildingCoef.srex_tke_h
    tke_eq.srex_v = BuildingCoef.srex_tke_v
    tke_eq.srex_sh = sh
    tke_eq.srex_bu = bu
    tke_eq.srex_veg = srex_tke_veg
    tke_eq.srim = srim_tke
    tke_eq.srim_td = td
    tke_eq.srim_veg = dissip_tke_veg
    tke_eq.diffusion = Km

    class qn_eq_Def():
        pass
    qn_eq = qn_eq_Def()
    qn_eq.srex = srex_qn
    qn_eq.srex_h = BuildingCoef.srex_qn_h
    qn_eq.srex_veg = srex_qn_veg
    qn_eq.srim = srim_qn
    qn_eq.diffusion = Km/ColParam.schmidt

    class ColumnModel_OtherParam_Def():
        pass
    ColumnModel_OtherParam = ColumnModel_OtherParam_Def()
    ColumnModel_OtherParam.Km = Km
    ColumnModel_OtherParam.dlk = dlk
    ColumnModel_OtherParam.dls = dls
    ColumnModel_OtherParam.ExpImp_th = [srex_th[i]+srim_th[i]*th[i] for i in range(0,Geometry_m.nz)]
    ColumnModel_OtherParam.wth = wth
    ColumnModel_OtherParam.uw = uw
    ColumnModel_OtherParam.vw = vw
    ColumnModel_OtherParam.wtke = wtke
    ColumnModel_OtherParam.wqn = wqn
    ColumnModel_OtherParam.duwdz = duwdz
    ColumnModel_OtherParam.dvwdz = dvwdz
    ColumnModel_OtherParam.dwtkedz = dwtkedz
    ColumnModel_OtherParam.dwthdz = dwthdz
    ColumnModel_OtherParam.dwqndz = dwqndz

    return vx,vy,tke,th,qn,presProf, th_eq, vx_eq, vy_eq, tke_eq, qn_eq, ColumnModel_OtherParam

