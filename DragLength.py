import numpy
import math
import copy

"""
Calculate mixing length scale for the 1-D column model in the urban area
Developed by Negin Nazarian(1) and Scott Krayenhoff(2)
1.University of New South Wales, Sydney, Australia
2.University of Guelph, Guelph, Canada
Last update: February 2020
"""


# Calculate sectional drag coefficient (Cdrag) for buildings w/o trees(eq. 4.17, Krayenhoff, PhD thesis)
def Drag_Coef(nz,lambdaf,pb):
    """
    -----
    INPUT:
    nz: Number of grid points in the urban area
    lambdaf: Frontal area density
    pb(z): Probability that a building has a height greater or equal to z
    -------
    OUTPUT:
    Cdrag: Drag coefficient due to buildings [-]
    """
    Cdrag = numpy.zeros(nz)
    for i in range(0,nz):
        if (lambdaf*pb[i+1]) <= 0.33:
            Cdrag[i] = 7.3*((lambdaf*pb[i+1])**(0.62))
        else:
            Cdrag[i] = 3.67
    return Cdrag

# Calculate turbulent and dissipation length scales
def Length_Scale(nz,z,lambdap,bldHeight,Ceps,Cmu,Ck):
    """
    ------
    INPUT:
    nz: Number of grid points in the urban area
    z: Vertical distribution of grids [m]
    lambdap: Plane area density
    bldHeight: Building height [m]
    Ceps: Coefficient for the destruction of turbulent dissipation rate
    Cmu: Model constant
    Ck: Model constant
    -------
    OUTPUT:
    dls: Dissipation length scale [m]
    dlk: Turbulent length scale [m]
    """

    # Option 1:
    # Eqs. 4.15 and 4.18, Krayenhoff, PhD thesis
    #a1 = 1.95
    #a2 = 1.07
    # Option 2:
    # Eq. 12 in "N. Nazarian et al., 2020"
    a1 = 4
    a2 = min(5,max(2,1.3*lambdap**(-0.45)))
    # Calculate displacement height (eq. 4.19, Krayenhoff 2014, PhD thesis)
    disp = bldHeight * (lambdap ** (0.15))
    # Dissipation length scale [m]
    dls = numpy.zeros(nz)
    # Turbulent length scale [m]
    dlk = numpy.zeros(nz)
    for i in range(0,nz):
        zc = (z[i]+z[i+1])/2

        if bldHeight == 0:
            dls[i] = Ceps*a2*zc[i]
        elif (zc/bldHeight) <= 1:
            dls[i] = Ceps*a1*(bldHeight-disp)
        elif (zc/bldHeight) > 1 and (zc/bldHeight) <= 1.5:
            dls[i] = Ceps*a1*(zc-disp)
        elif (zc/bldHeight) > 1.5:
            d2 = (1 - a1 / a2) * 1.5 * bldHeight + (a1 / a2) * disp
            dls[i] = Ceps*a2*(zc-d2)
        dlk[i] = Cmu[i]*dls[i]/(Ceps*Ck)

    return dls,dlk

# Calculate turbulent and dissipation length scales
def Length_Scale_StabilityCorrection(nz,z,bldHeight,Ceps,Ck,hfx,vx_eq,vy_eq,VerticalProfUrban,lambdaf,disp,dz,vl):

    """
    ------
    INPUT:
    nz: Number of grid points in the urban area
    z: Vertical distribution of grids [m]
    bldHeight: Building height [m]
    Ceps: Coefficient for the destruction of turbulent dissipation rate
    Ck: Model constant
    hfx: Total urban heat flux per unit flat area of the earth [W m^-2]
    vx_eq: Terms in x-momentum equation
    vy_eq: Terms in y-momentum equation
    VerticalProfUrban: Vertical profile of temperature, humidity, wind speed, and turbulent kinetic energy from previous time step
    lambdaf: Frontal area density
    disp: Displacement height [m]
    dz: Vertical resolution [m]
    vl: Volume fraction of air in each urban unit cell [-]
    -------
    OUTPUT:
    dls: Dissipation length scale [m]
    dlk: Turbulent length scale [m]
    """

    srim_vx_sav = copy.copy(vx_eq.srim)
    srim_vy_sav = copy.copy(vy_eq.srim)
    srex_vx_sav = copy.copy(vx_eq.srex)
    srex_vy_sav = copy.copy(vy_eq.srex)
    vx_sav = copy.copy(VerticalProfUrban.vx)
    vy_sav = copy.copy(VerticalProfUrban.vy)
    th0 = copy.copy(VerticalProfUrban.th_ref[0])

    g = 9.81

    dls = numpy.zeros(nz)
    dlk = numpy.zeros(nz)
    dls_neut = numpy.zeros(nz)
    if hfx > 0:

        sum_drag = 0
        for iz in range(0,nz):
            sum_drag = sum_drag+dz*vl[iz]*(abs(srim_vx_sav[iz]*vx_sav[iz]) + abs(srim_vy_sav[iz]*vy_sav[iz]) +
                                           abs(srex_vx_sav[iz]) + abs(srex_vy_sav[iz]))
        u_tau = numpy.sqrt(abs(sum_drag))
        H_L = bldHeight/(u_tau**3.)/(g/th0)/abs(hfx)

        a1 = 5.00
        if lambdaf == 0:
            a2 = 5.
        else:
            a2 = min(5.,max(2.,1.23*lambdaf**(-0.442)))

        cmu_can = max(0.06,-1.3678*lambdaf**2.+0.658*lambdaf+0.0303)
        cmu_above = 0.037

        for iz in range(0,nz):
            zc = (z[iz]+z[iz+1])/2
            if bldHeight == 0:
                dls[iz] = Ceps*a2*zc
                cmu = cmu_can
                dls_neut[iz] = dls[iz]
                dls[iz] = dls[iz]*(1.+0.5*min(H_L,3.))
            else:
                if zc/bldHeight <= 1:
                    dls[iz] = Ceps*a1*(bldHeight-disp)
                    cmu = cmu_can

                    if zc/bldHeight > 0.9:
                        cmu = (cmu_can*(1.1-zc/bldHeight)+cmu_above*(zc/bldHeight-0.9))/0.2
                    dls_neut[iz] = dls[iz]
                    dls[iz] = dls[iz]*(1+(0.2+0.3*zc/bldHeight)*min(H_L,3.))
                elif zc/bldHeight > 1 and zc/bldHeight <= 1.5:
                    dls[iz] = Ceps*a1*(zc-disp)
                    cmu = cmu_above

                    if zc/bldHeight < 1.25:
                        cmu = (cmu_can*(1.25-zc/bldHeight)+cmu_above*(zc/bldHeight-0.75))/0.5
                    dls_neut[iz] = dls[iz]
                    dls[iz] = dls[iz]*(1.+0.5*min(H_L,3.))

                elif zc/bldHeight > 1.5:
                    d2 = (1.-a1/a2)*1.5*bldHeight+a1/a2*disp
                    dls[iz] = Ceps*a2*(zc-d2)
                    cmu = cmu_above
                    dls_neut[iz] = dls[iz]
                    dls[iz] = dls[iz]*(1.+0.5*min(H_L,3.))

            dlk[iz] = cmu/(Ceps*Ck)*dls[iz]
    else:

        a1 = 5.00
        if lambdaf == 0:
            a2 = 5.
        else:
            a2 = min(5., max(2., 1.23 * lambdaf ** (-0.442)))

        cmu_can = max(0.06, -1.3678 * lambdaf ** 2. + 0.658 * lambdaf + 0.0303)
        cmu_above = 0.037

        for iz in range(0, nz):
            zc = (z[iz] + z[iz + 1]) / 2
            if bldHeight == 0:
                dls[iz] = Ceps * a2 * zc
                cmu = cmu_can
            else:
                if zc / bldHeight <= 1:
                    dls[iz] = Ceps * a1 * (bldHeight - disp)
                    cmu = cmu_can
                    if zc / bldHeight > 0.9:
                        cmu = (cmu_can*(1.1-zc/bldHeight) + cmu_above*(zc/bldHeight-0.9))/0.2
                elif zc / bldHeight > 1 and zc / bldHeight <= 1.5:
                    dls[iz] = Ceps * a1 * (zc - disp)
                    cmu = cmu_above
                    if zc / bldHeight < 1.1:
                        cmu = (cmu_can*(1.1-zc/bldHeight) + cmu_above*(zc/bldHeight-0.9))/0.2
                elif zc / bldHeight > 1.5:
                    d2 = (1. - a1 / a2) * 1.5 * bldHeight + a1 / a2 * disp
                    dls[iz] = Ceps * a2 * (zc - d2)
                    cmu = cmu_above

            dlk[iz] = cmu / (Ceps * Ck) * dls[iz]

    return dls, dlk






